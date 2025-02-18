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
    stream = NULL;
    CUBLAS_CHECK(cublasCreate(&handle));
    // CUBLAS_CHECK(cublasSetStream(handle, stream));
  }

  MatPoly( cublasHandle_t handle_, cudaStream_t stream_ )
    : is_external(true)
  {
    handle = handle_;
    stream = stream_; // pointer
  }

  ~MatPoly(){
    if(!is_external) CUBLAS_CHECK(cublasDestroy(handle));
  }

  int size() const { return vec_mats.size(); }
  std::vector<const LinOp*> operator[](const int i) const { return vec_mats[i]; }

  void push_back( const T coeff,
		  const std::initializer_list<const LinOp*> struc={} ){
    coeffs.push_back(coeff);
    vec_mats.push_back(std::vector<const LinOp*>(0));
    for( auto itr=struc.begin(); itr!=struc.end(); ++itr ) vec_mats.back().push_back( *itr );
  }


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
    CuC *d_tmp, *d_Mv0;
    CUDA_CHECK(cudaMalloc(&d_tmp, N*CD));
    CUDA_CHECK(cudaMalloc(&d_Mv0, N*CD));
    CUDA_CHECK(cudaMemset(d_v, 0, N*CD));

    for(int i=0; i<vec_mats.size(); i++){
      CUDA_CHECK(cudaMemcpy(d_tmp, d_v0, N*CD, D2D));
      CUDA_CHECK(cudaMemcpy(d_Mv0, d_v0, N*CD, D2D));
      for(int j=0; j<vec_mats[i].size(); j++){
        vec_mats[i][j]->operator()(d_Mv0, d_tmp);
        CUDA_CHECK(cudaMemcpy(d_tmp, d_Mv0, N*CD, D2D));
      }
      Taxpy<CuC, N><<<NBlocks, NThreadsPerBlock>>>(d_v,
                                                   coeffs[i],
        					   d_Mv0,
        					   d_v);
    }

    CUDA_CHECK(cudaFree(d_tmp));
    CUDA_CHECK(cudaFree(d_Mv0));
  }


  template<Idx N> __host__
  void on_gpuAsync(CuC* d_v, const CuC* d_v0) const {
    CuC *d_tmp, *d_Mv0;
    CUDA_CHECK(cudaMallocAsync(&d_tmp, N*CD, stream));
    CUDA_CHECK(cudaMallocAsync(&d_Mv0, N*CD, stream));
    CUDA_CHECK(cudaMemsetAsync(d_v, 0, N*CD, stream));

    for(int i=0; i<vec_mats.size(); i++){
      CUDA_CHECK(cudaMemcpyAsync(d_tmp, d_v0, N*CD, D2D, stream));
      CUDA_CHECK(cudaMemcpyAsync(d_Mv0, d_v0, N*CD, D2D, stream));
      for(int j=0; j<vec_mats[i].size(); j++){
        vec_mats[i][j]->Async(d_Mv0, d_tmp, stream);
        CUDA_CHECK(cudaMemcpyAsync(d_tmp, d_Mv0, N*CD, D2D, stream));
      }
      Taxpy<CuC, N><<<NBlocks, NThreadsPerBlock, 0, stream>>>(d_v,
                                                              coeffs[i],
                                                              d_Mv0,
                                                              d_v);
    }

    CUDA_CHECK(cudaFreeAsync(d_tmp, stream));
    CUDA_CHECK(cudaFreeAsync(d_Mv0, stream));
    CUDA_CHECK( cudaStreamSynchronize(stream) );
  }

  template<Idx N> __host__
  inline void Zdscal( const double alpha,
                      CuC* x) const {
    CUBLAS_CHECK( cublasZdscal(handle, N,
                               &alpha,
                               x, 1) );
  }

  template<Idx N> __host__
  inline void ZdscalAsync( const double alpha,
                      CuC* x) const {
    CUBLAS_CHECK( cublasZdscal(handle, N,
                               &alpha,
                               x, 1) );
    CUDA_CHECK(cudaStreamSynchronize(stream));
  }

  template<Idx N> __host__
  inline void dot( CuC* result, const CuC* x, const CuC* y) const {
    CUBLAS_CHECK( cublasZdotc(handle, N,
        		      x, 1,
        		      y, 1,
                              result) );
  }

  template<Idx N> __host__
  inline void dotAsync( CuC* result, const CuC* x, const CuC* y) const {
    CUBLAS_CHECK( cublasZdotc(handle, N,
			      x, 1,
			      y, 1,
                              result) );
    CUDA_CHECK(cudaStreamSynchronize(stream));
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
    if( crit >= TOL*std::sqrt(N) || isnan(crit) || isinf(crit)  ) std::clog << "CRIT = " << crit << std::endl;
    assert( crit < TOL*std::sqrt(N) );
    *result = real(dummy);
  }

  template<Idx N> __host__
  void dot2selfAsync( double* result, // CuC* d_scalar,
                      const CuC* x, const double TOL=1.0e-12) const {
    CuC dummy;
    cublasZdotc(handle, N,
		x, 1,
		x, 1,
		&dummy );
    CUDA_CHECK( cudaStreamSynchronize(stream) );

    double crit = abs( imag(dummy)/real(dummy) );
    if( isnan(crit) || isinf(crit) ){
      crit = abs( imag(dummy) );
    }
    if( crit >= TOL*std::sqrt(N) || isnan(crit) || isinf(crit)  ) std::clog << "CRIT = " << crit << std::endl;
    assert( crit < TOL*std::sqrt(N) );
    *result = real(dummy);

    CUDA_CHECK(cudaStreamSynchronize(stream));
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


  // necessary for outer loop
  template<Idx N> __host__
  void solve(CuC* d_x, const CuC* d_b,
	     const double tol=1.0e-13, const int maxiter=1e8) const {
    // CG
    CuC *d_p, *d_q, *d_r;
    CUDA_CHECK(cudaMalloc(&d_p, N*CD));
    CUDA_CHECK(cudaMalloc(&d_q, N*CD));
    CUDA_CHECK(cudaMalloc(&d_r, N*CD));
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
#ifdef IsVerbose2
      std::clog << "NO SOLVE" << std::endl;
#endif
    }
    else{
      int k=0;
      for(; k<maxiter; ++k){
	this->on_gpu<N>(d_q, d_p);

        CuC gam; dot<N>(&gam, d_p, d_q);
	const CuC al = mu/gam;
        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock>>>(d_x, al, d_p, d_x);
        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock>>>(d_r, -al, d_q, d_r);

	dot2self<N>(&mu, d_r);
        assert(mu>=0.0);
	if(mu<mu_crit || std::isnan(mu)) break;
	const CuC bet = cplx(mu/mu_old);
	mu_old = mu;

        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock>>>(d_p, bet, d_p, d_r);

	if(k%100==0) {
#ifdef IsVerbose2
	  std::clog << "SOLVER:       #iterations: " << k << ", mu =         " << mu << std::endl;
#endif
	}
      }
#ifdef IsVerbose2
      std::clog << "SOLVER:       #iterations: " << k << std::endl;
      std::clog << "SOLVER:       mu =         " << mu << std::endl;
#endif
    }

    CUDA_CHECK(cudaFree(d_r));
    CUDA_CHECK(cudaFree(d_p));
    CUDA_CHECK(cudaFree(d_q));
  }



  template<Idx N> __host__
  void solveAsync(CuC* d_x, const CuC* d_b,
                  const double tol=1.0e-13, const int maxiter=1e8) const {
    // CG
    CuC *d_p, *d_q, *d_r;
    CUDA_CHECK(cudaMallocAsync(&d_p, N*CD, stream));
    CUDA_CHECK(cudaMallocAsync(&d_q, N*CD, stream));
    CUDA_CHECK(cudaMallocAsync(&d_r, N*CD, stream));
    //
    CUDA_CHECK(cudaMemsetAsync(d_x, 0, N*CD, stream));
    CUDA_CHECK(cudaMemsetAsync(d_q, 0, N*CD, stream));

    CUDA_CHECK(cudaMemcpyAsync(d_r, d_b, N*CD, D2D, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_p, d_r, N*CD, D2D, stream));

    double mu;
    this->dot2selfAsync<N>(&mu, d_r);
    assert(mu>=0.0);
    double mu_old = mu;

    double b_norm_sq;
    this->dot2selfAsync<N>(&b_norm_sq, d_r);
    assert(b_norm_sq>=0.0);
    double mu_crit = tol*tol*b_norm_sq;

    if(mu<mu_crit) {
#ifdef IsVerbose
      std::clog << "NO SOLVE" << std::endl;
#endif
    }
    else{
      int k=0;
      for(; k<maxiter; ++k){
	this->on_gpuAsync<N>(d_q, d_p);

        CuC gam;
        this->dotAsync<N>(&gam, d_p, d_q);
	const CuC al = mu/gam;
        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock, 0, stream>>>(d_x, al, d_p, d_x);
        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock, 0, stream>>>(d_r, -al, d_q, d_r);
        CUDA_CHECK( cudaStreamSynchronize(stream) );

	this->dot2selfAsync<N>(&mu, d_r);
        assert(mu>=0.0);
	if(mu<mu_crit || std::isnan(mu)) break;
	const CuC bet = cplx(mu/mu_old);
	mu_old = mu;

        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock, 0, stream>>>(d_p, bet, d_p, d_r);
        CUDA_CHECK( cudaStreamSynchronize(stream) );

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

    CUDA_CHECK(cudaFreeAsync(d_r, stream));
    CUDA_CHECK(cudaFreeAsync(d_p, stream));
    CUDA_CHECK(cudaFreeAsync(d_q, stream));
    CUDA_CHECK( cudaStreamSynchronize(stream) );
  }


  template<Idx N>
  void bicgstab( std::vector<Complex>& v, const std::vector<Complex>& v0,
                 const std::vector<Complex>& rc,
                 const double tol=1.0e-6, const int maxiter=1e8,
                 const double eps=1.0e-8 ) const {
    CuC *d_v, *d_v0, *d_rc;
    CUDA_CHECK(cudaMalloc(&d_v, N*CD));
    CUDA_CHECK(cudaMalloc(&d_v0, N*CD));
    CUDA_CHECK(cudaMalloc(&d_rc, N*CD));

    CUDA_CHECK(cudaMemcpy(d_v0, reinterpret_cast<const CuC*>(v0.data()), N*CD, H2D));
    CUDA_CHECK(cudaMemcpy(d_rc, reinterpret_cast<const CuC*>(rc.data()), N*CD, H2D));

    bicgstab<N>( d_v, d_v0, d_rc, tol, maxiter, eps );

    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(v.data()), d_v, N*CD, D2H));

    CUDA_CHECK(cudaFree(d_v));
    CUDA_CHECK(cudaFree(d_v0));
    CUDA_CHECK(cudaFree(d_rc));
  }



  template<Idx N> __host__
  void bicgstab(CuC* d_x, const CuC* d_b, CuC* d_rc,
                // double& mu=1.0,
                const double tol=1.0e-4, const int maxiter=1e8,
                const double eps=1.0e-8 ) const {
    CuC *d_p, *d_t, *d_Ap, *d_At, *d_r;
    CUDA_CHECK(cudaMalloc(&d_p, N*CD));
    CUDA_CHECK(cudaMalloc(&d_t, N*CD));
    CUDA_CHECK(cudaMalloc(&d_Ap, N*CD));
    CUDA_CHECK(cudaMalloc(&d_At, N*CD));
    CUDA_CHECK(cudaMalloc(&d_r, N*CD));
    // CUDA_CHECK(cudaMalloc(&d_rc, N*CD));

    CUDA_CHECK(cudaMemset(d_x, 0, N*CD));
    CUDA_CHECK(cudaMemset(d_p, 0, N*CD));
    CUDA_CHECK(cudaMemset(d_t, 0, N*CD));
    CUDA_CHECK(cudaMemset(d_Ap, 0, N*CD));
    CUDA_CHECK(cudaMemset(d_At, 0, N*CD));
    CUDA_CHECK(cudaMemset(d_r, 0, N*CD));
    // CUDA_CHECK(cudaMemset(d_rc, 0, N*CD));

    CUDA_CHECK(cudaMemcpy(d_x, d_b, N*CD, D2D));
    this->on_gpu<N>(d_Ap, d_x);
    Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock, 0>>>(d_r, -1.0, d_Ap, d_b);
    CUDA_CHECK(cudaMemcpy(d_rc, d_b, N*CD, D2D));
    Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock, 0>>>(d_rc, -1.0, d_Ap, d_rc);

    CuC gam;
    this->dot<N>(&gam, d_rc, d_r);
    assert( abs(gam)>tol );

    CUDA_CHECK(cudaMemcpy(d_p, d_r, N*CD, D2D));

    double mu;
    this->dot2self<N>(&mu, d_r);
    assert(mu>=0.0);
    double mu_old = mu;

    double nm;

    double b_norm_sq;
    this->dot2self<N>(&b_norm_sq, d_r);
    assert(b_norm_sq>=0.0);
    double mu_crit = tol*tol*b_norm_sq;

    if(mu<mu_crit) {
#ifdef IsVerbose
      std::clog << "NO SOLVE" << std::endl;
#endif
    }
    else{
      int n = 0;
      for (; n < maxiter; n++) {
        this->on_gpu<N>(d_Ap, d_p);

        CuC gam1, gam2;
        this->dot<N>(&gam1, d_rc, d_r);
        this->dot<N>(&gam2, d_rc, d_Ap);
        const CuC al = gam1/gam2;

        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock, 0>>>(d_t, -al, d_Ap, d_r);
        this->dot2self<N>(&nm, d_t);
        if(nm<eps) break;
        this->on_gpu<N>(d_At, d_t);

        this->dot<N>(&gam1, d_t, d_At);
        this->dot<N>(&gam2, d_At, d_At);
        const CuC om = gam1/gam2;

        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock, 0>>>(d_x, al, d_p, d_x);
        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock, 0>>>(d_x, om, d_t, d_x);

        this->dot<N>(&gam2, d_rc, d_r);
        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock, 0>>>(d_r, -om, d_At, d_t);
        this->dot<N>(&gam1, d_rc, d_r);
        const CuC bet = al/om * gam1/gam2;

        this->dot2self<N>(&mu, d_r);
        assert(mu>=0.0);
        if(mu<mu_crit || std::isnan(mu)) break;
        mu_old = mu;

        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock, 0>>>(d_p, -om, d_Ap, d_p);
        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock, 0>>>(d_p, bet, d_p, d_r);

        // restart
        this->dot<N>(&gam, d_rc, d_r);
        if( abs(gam)<eps || 2.0*mu_old < mu ){
          CUDA_CHECK(cudaMemcpy(d_rc, d_r, N*CD, D2D));
          CUDA_CHECK(cudaMemcpy(d_p, d_r, N*CD, D2D));
#ifdef IsVerbose2
	  std::clog << "SOLVER:      restart " << std::endl;
#endif
        }

        // std::cout << "debug. mu = " << mu << std::endl;

        if(n%10==0) {
#ifdef IsVerbose2
	  std::clog << "SOLVER:       #iterations: " << n << ", mu =         " << mu << std::endl;
#endif
	}

      }

      if (n == maxiter) {
        std::cerr << "BiCGStab iteration did not converge." << std::endl;
        // std::cerr << "v = " << v << std::endl;
      }
    }

    CUDA_CHECK(cudaFree(d_p));
    CUDA_CHECK(cudaFree(d_Ap));
    CUDA_CHECK(cudaFree(d_At));
    CUDA_CHECK(cudaFree(d_r));
    // CUDA_CHECK(cudaFree(d_rc));
    CUDA_CHECK(cudaFree(d_t));
    // CUDA_CH/ECK( cudaStreamSynchronize(stream) );

    // return x;
  }




};


