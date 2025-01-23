#pragma once

// template<typename T>
struct MatPoly{
  using T = CuC;

  std::vector<std::vector<const LinOp*>> vec_mats;
  std::vector<T> coeffs;

  cublasHandle_t handle;
  cudaStream_t stream;
  const bool is_external;

  MatPoly()
    : is_external(false)
  {
    handle = NULL;
    stream = 0;
    CUBLAS_CHECK(cublasCreate(&handle));
    // CUDA_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    CUBLAS_CHECK(cublasSetStream(handle, stream));
  }

  MatPoly( cublasHandle_t handle_, cudaStream_t stream_ )
    : handle(handle_)
    , stream(stream_)
    , is_external(true)
  {}


  MatPoly( const std::initializer_list<int> struc )
    : vec_mats(struc.size())
    , coeffs(struc.size())
    , is_external(false)
  {
    auto itr = struc.begin();
    for(int i=0; i<struc.size(); i++) {
      vec_mats[i].resize( *itr );
      ++itr;
    }

    handle = NULL;
    stream = 0;
    CUBLAS_CHECK(cublasCreate(&handle));
    CUBLAS_CHECK(cublasSetStream(handle, stream));
  }

  ~MatPoly(){
    if(!is_external) CUBLAS_CHECK(cublasDestroy(handle));
  }

  int size() const { return vec_mats.size(); }
  std::vector<const LinOp*> operator[](const int i) const { return vec_mats[i]; }

  void set_coeff( const int i, const T coeff ){ coeffs[i] = coeff; }

  void set_matrix( const int i, const int j,
		   const LinOp* mat_ptr ){
    vec_mats[i][j] = mat_ptr;
  }

  void push_back( const T coeff,
		  const std::initializer_list<const LinOp*> struc={} ){
    coeffs.push_back(coeff);
    vec_mats.push_back(std::vector<const LinOp*>(0));
    for( auto itr=struc.begin(); itr!=struc.end(); ++itr ) vec_mats.back().push_back( *itr );
  }


  // template<Idx N>
  // void on_cpu( std::vector<Complex>& v, const std::vector<Complex>& v0 ) const {
  //   for(Idx j=0; j<N; j++) v[j] = 0.0;
  //   std::vector<Complex> tmp(N), Mv0(N);

  //   for(int i=0; i<vec_mats.size(); i++){
  //     tmp = v0;
  //     Mv0 = v0;
  //     for(int j=0; j<vec_mats[i].size(); j++){
  //       vec_mats[i][j]->act_cpu( Mv0, tmp );
  //       tmp = Mv0;
  //     }
  //     for(Idx j=0; j<N; j++) v[j] += coeffs[i] * Mv0[j];
  //   }
  // }


  template<Idx N>
  void from_cpu( std::vector<Complex>& v, const std::vector<Complex>& v0 ) const {
    CuC *d_v, *d_v0;
    CUDA_CHECK(cudaMalloc(&d_v, N*CD));
    CUDA_CHECK(cudaMalloc(&d_v0, N*CD));

    CUDA_CHECK(cudaMemcpy(d_v0, reinterpret_cast<const CuC*>(v0.data()), N*CD, H2D));

    on_gpu<N>( d_v, d_v0 );

    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(v.data()), d_v, N*CD, D2H));

    CUDA_CHECK(cudaFree(d_v));
    CUDA_CHECK(cudaFree(d_v0));
  }

  template<Idx N> __host__
  void on_gpu(CuC* d_v, const CuC* d_v0) const {
    CuC *d_tmp, *d_Mv0, *d_coeffs;
    CUDA_CHECK(cudaMalloc(&d_tmp, N*CD));
    CUDA_CHECK(cudaMalloc(&d_Mv0, N*CD));
    CUDA_CHECK(cudaMalloc(&d_coeffs, coeffs.size()*CD));

    CUDA_CHECK(cudaMemset(d_v, 0, N*CD));
    CUDA_CHECK(cudaMemcpy(d_coeffs, coeffs.data(), coeffs.size()*CD, H2D));

    for(int i=0; i<vec_mats.size(); i++){
      CUDA_CHECK(cudaMemcpy(d_tmp, d_v0, N*CD, D2D));
      CUDA_CHECK(cudaMemcpy(d_Mv0, d_v0, N*CD, D2D));
      for(int j=0; j<vec_mats[i].size(); j++){
	vec_mats[i][j]->operator()(d_Mv0, d_tmp);
	CUDA_CHECK(cudaMemcpy(d_tmp, d_Mv0, N*CD, D2D));
      }
      Taxpy<CuC, N><<<NBlocks, NThreadsPerBlock>>>(d_v,
        					   d_coeffs+i,
        					   d_Mv0,
        					   d_v);
      // zaxpy<N>( d_coeffs+i, d_Mv0, d_v );
    }

    CUDA_CHECK(cudaFree(d_tmp));
    CUDA_CHECK(cudaFree(d_Mv0));
    CUDA_CHECK(cudaFree(d_coeffs));
  }

  // template<Idx N> __host__ // use taxpy
  // void zaxpy( const CuC* a,
  //             const CuC* x, CuC* y) const {
  //   const int incx=1, incy=1;
  //   CUBLAS_CHECK( cublasZaxpy(handle, N,
  //       		      a,
  //       		      x, incx,
  //       		      y, incy) );
  // }

  template<Idx N> __host__
  void Zdscal( const double alpha,
               CuC* x) const {
    const int incx=1, incy=1;
    CUBLAS_CHECK( cublasZdscal(handle, N,
                               &alpha,
                               x, incx) );
  }

  template<Idx N> __host__
  inline void dot( CuC* result, const CuC* x, const CuC* y) const {
    CUBLAS_CHECK( cublasZdotc(handle, N,
			      x, 1,
			      y, 1,
                              result) );
  }

  template<Idx N> __host__
  void dot2self( double* result, // CuC* d_scalar,
                 const CuC* x, const double TOL=1.0e-12) const {
    CuC dummy;
    cublasZdotc(handle, N,
		x, 1,
		x, 1,
		&dummy );

    double crit = abs( imag(dummy)/real(dummy) );
    if( isnan(crit) || isinf(crit) ){
      crit = abs( imag(dummy) );
    }
    //#ifdef IsVerbose
    if( crit >= TOL*std::sqrt(N) || isnan(crit) || isinf(crit)  ){
      std::clog << crit << std::endl;
    }
    // #endif
    assert( crit < TOL*std::sqrt(N) );
    *result = real(dummy);
  }


  template<Idx N> __host__
  void solve(std::vector<Complex>& x, const std::vector<Complex>& b,
	     const double tol=1.0e-13, const int maxiter=1e8) const {
    // CG
    CuC *d_x, *d_r; // , *d_tmp, *d_tmp2;
    CUDA_CHECK(cudaMalloc(&d_x, N*CD));
    CUDA_CHECK(cudaMalloc(&d_r, N*CD));
    CUDA_CHECK(cudaMemcpy(d_r, reinterpret_cast<const CuC*>(b.data()), N*CD, H2D));

    solve<N>(d_x, d_r, tol, maxiter);

    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(x.data()), d_x, N*CD, D2H));
    CUDA_CHECK(cudaFree(d_x));
    CUDA_CHECK(cudaFree(d_r));
  }



  template<Idx N> __host__
  void solve(CuC* d_x, const CuC* d_b,
	     const double tol=1.0e-13, const int maxiter=1e8) const {
    // CG
    CuC *d_p, *d_q, *d_r; // , *d_tmp, *d_tmp2;
    CUDA_CHECK(cudaMalloc(&d_r, N*CD));
    CUDA_CHECK(cudaMalloc(&d_p, N*CD));
    CUDA_CHECK(cudaMalloc(&d_q, N*CD));
    //
    CUDA_CHECK(cudaMemset(d_x, 0, N*CD));
    CUDA_CHECK(cudaMemset(d_q, 0, N*CD));

    CUDA_CHECK(cudaMemcpy(d_r, d_b, N*CD, D2D));
    CUDA_CHECK(cudaMemcpy(d_p, d_r, N*CD, D2D));

    double mu;
    dot2self<N>(&mu, d_r);
    assert(mu>=0.0);
    double mu_old = mu;

    double b_norm_sq;
    dot2self<N>(&b_norm_sq, d_r);
    assert(b_norm_sq>=0.0);
    double mu_crit = tol*tol*b_norm_sq;

    if(mu<mu_crit) {
#ifdef IsVerbose
      std::clog << "NO SOLVE" << std::endl;
#endif
    }
    else{
      int k=0;
      CuC gam;

      for(; k<maxiter; ++k){
	this->on_gpu<N>(d_q, d_p);

	dot<N>(&gam, d_p, d_q);

	CuC al = mu/gam;
        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock>>>(d_x, al, d_p, d_x);

	// al = -al;
        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock>>>(d_r, -al, d_q, d_r);

	dot2self<N>(&mu, d_r);
	assert(mu>=0.0);

	if(mu<mu_crit || std::isnan(mu)) break;
	CuC bet = cplx(mu/mu_old);
	mu_old = mu;

        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock>>>(d_p, bet, d_p, d_r);

	if(k%100==0) {
#ifdef IsVerbose
	  std::clog << "SOLVER:       #iterations: " << k << ", mu =         " << mu << std::endl;
#endif
	}
      }
#ifdef IsVerbose
      std::clog << "SOLVER:       #iterations: " << k << std::endl;
      std::clog << "SOLVER:       mu =         " << mu << std::endl;
#endif
    }

    CUDA_CHECK(cudaFree(d_r));
    CUDA_CHECK(cudaFree(d_p));
    CUDA_CHECK(cudaFree(d_q));
  }



};


