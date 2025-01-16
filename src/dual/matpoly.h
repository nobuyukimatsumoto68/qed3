#pragma once

// template<typename T>
struct MatPoly{
  using T = CuC;

  std::vector<std::vector<const SparseMatrix*>> vec_mats;
  std::vector<T> coeffs;

  cublasHandle_t handle;
  cudaStream_t stream;

  MatPoly()
  {
    handle = NULL;
    stream = NULL;

    CUBLAS_CHECK(cublasCreate(&handle));
    CUDA_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    CUBLAS_CHECK(cublasSetStream(handle, stream));
  }

  MatPoly( const std::initializer_list<int> struc )
    : vec_mats(struc.size())
    , coeffs(struc.size())
  {
    auto itr = struc.begin();
    for(int i=0; i<struc.size(); i++) {
      vec_mats[i].resize( *itr );
      ++itr;
    }

    handle = NULL;
    stream = NULL;

    CUBLAS_CHECK(cublasCreate(&handle));
    CUDA_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    CUBLAS_CHECK(cublasSetStream(handle, stream));
  }

  ~MatPoly(){
    CUBLAS_CHECK(cublasDestroy(handle));
    CUDA_CHECK(cudaStreamDestroy(stream));
  }

  int size() const { return vec_mats.size(); }
  std::vector<const SparseMatrix*> operator[](const int i) const { return vec_mats[i]; }

  void set_coeff( const int i, const T coeff ){ coeffs[i] = coeff; }

  void set_matrix( const int i, const int j,
		   const SparseMatrix* mat_ptr ){
    vec_mats[i][j] = mat_ptr;
  }
  
  void push_back( const T coeff,
		  const std::initializer_list<const SparseMatrix*> struc={} ){
    coeffs.push_back(coeff);
    vec_mats.push_back(std::vector<const SparseMatrix*>(0));
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
    CUDA_CHECK(cudaMemset(d_v, 0, N*CD));

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
    CUDA_CHECK(cudaMemset(d_tmp, 0, N*CD));
    CUDA_CHECK(cudaMemset(d_Mv0, 0, N*CD));
    CUDA_CHECK(cudaMemset(d_coeffs, 0, coeffs.size()*CD));
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
      // zaxpy( coeffs[i], d_Mv0, d_v );
    }
    
    CUDA_CHECK(cudaFree(d_tmp));
    CUDA_CHECK(cudaFree(d_Mv0));
    CUDA_CHECK(cudaFree(d_coeffs));
  }

  // __host__ // use taxpy
  // void zaxpy( const CuC a,
  // 	      const CuC* x, CuC* y) const {
  //   // std::cout << "debug 1." << std::endl;
  //   const int incx=1, incy=1;
  //   CUBLAS_CHECK( cublasZaxpy(handle, CompilationConst::N,
  // 			      &a,
  // 			      x, incx,
  // 			      y, incy) );
  //   // std::cout << "debug 2." << std::endl;
  // }


  template<Idx N> __host__
  void dot( CuC& result, CuC* d_scalar,
	    const CuC* x, const CuC* y) const {
    CUBLAS_CHECK( cublasZdotc(handle, N,
			      x, 1,
			      y, 1,
			      d_scalar) );
    CUDA_CHECK(cudaMemcpy(&result, d_scalar, CD, D2H));
  }

  template<Idx N> __host__
  void dot2self( double& result, CuC* d_scalar, const CuC* x, const double TOL=1.0e-12) const {
    cublasZdotc(handle, N,
		x, 1,
		x, 1,
		d_scalar);
    CuC dummy;
    CUDA_CHECK(cudaMemcpy(&dummy, d_scalar, CD, D2H));
#ifdef IsVerbose
    if(abs( imag(dummy)/real(dummy) )>=TOL*std::sqrt(N)){
      std::clog << abs( imag(dummy)/real(dummy) ) << std::endl;
    }
#endif
    assert( abs( imag(dummy)/real(dummy) )<TOL*std::sqrt(N) );
    result = real(dummy);
  }


  template<Idx N> __host__
  void solve(std::vector<Complex>& x, const std::vector<Complex>& b,
	     const double tol=1.0e-13, const int maxiter=1e8) const {
    // CG
    CuC *d_x, *d_r, *d_p, *d_q; // , *d_tmp, *d_tmp2;
    CUDA_CHECK(cudaMalloc(&d_x, N*CD));
    CUDA_CHECK(cudaMalloc(&d_r, N*CD));
    CUDA_CHECK(cudaMalloc(&d_p, N*CD));
    CUDA_CHECK(cudaMalloc(&d_q, N*CD));
    //
    CUDA_CHECK(cudaMemset(d_x, 0, N*CD)); // added @@
    CUDA_CHECK(cudaMemset(d_q, 0, N*CD)); // added @@

    CuC *d_scalar;
    CUDA_CHECK(cudaMalloc(&d_scalar, CD));
    CUDA_CHECK(cudaMemset(d_scalar, 0, CD)); // added @@

    CUDA_CHECK(cudaMemcpy(d_r, reinterpret_cast<const CuC*>(b.data()), N*CD, H2D));
    // CUDA_CHECK(cudaMemcpy(d_r, reinterpret_cast<const CuC*>(b), N*CD, H2D));
    CUDA_CHECK(cudaMemcpy(d_p, d_r, N*CD, D2D));

    double mu;
    // dot2self_normalized_wrapper<N>(mu, d_scalar, d_r);
    dot2self<N>(mu, d_scalar, d_r);
    assert(mu>=0.0);
    double mu_old = mu;

    double b_norm_sq;
    // dot2self_normalized_wrapper<N>(b_norm_sq, d_scalar, d_r);
    dot2self<N>(b_norm_sq, d_scalar, d_r);
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

	// dot_normalized_wrapper<N>(gam, d_scalar, d_p, d_q);
	dot<N>(gam, d_scalar, d_p, d_q);

	CuC al = mu/gam;
	CUDA_CHECK(cudaMemcpy(d_scalar, &al, CD, H2D));
	Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock>>>(d_x, d_scalar, d_p, d_x);
	// Zaxpy<N>( d_scalar, d_p, d_x );

	al = -al;
	CUDA_CHECK(cudaMemcpy(d_scalar, &al, CD, H2D));
	Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock>>>(d_r, d_scalar, d_q, d_r);
	// Zaxpy<N>( d_scalar, d_q, d_r );

	// dot2self_normalized_wrapper<N>(mu, d_scalar, d_r);
	dot2self<N>(mu, d_scalar, d_r);
	assert(mu>=0.0);

	if(mu<mu_crit || std::isnan(mu)) break;
	CuC bet = cplx(mu/mu_old);
	mu_old = mu;

	CUDA_CHECK(cudaMemcpy(d_scalar, &bet, CD, H2D));
	Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock>>>(d_p, d_scalar, d_p, d_r);
	// Zaxpy<N>( d_scalar, d_p, d_r );
	// CUBLAS_CHECK(cublasZswap(handle, N,
	// 			 d_p, 1,
	// 			 d_r, 1));

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

    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(x.data()), d_x, N*CD, D2H));
    // CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(x), d_x, N*CD, D2H));

    CUDA_CHECK(cudaFree(d_x));
    CUDA_CHECK(cudaFree(d_r));
    CUDA_CHECK(cudaFree(d_p));
    CUDA_CHECK(cudaFree(d_q));
    //
    CUDA_CHECK(cudaFree(d_scalar));
  }



};


