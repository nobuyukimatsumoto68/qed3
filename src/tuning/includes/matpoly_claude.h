#pragma once

#include "multishift_kernels_claude.h"        // per-pole block kernels for solve_multishift
#include "multishift_block_kernels_claude.h"  // mrhs (nstack) block kernels, used by BlockedMat (blocked_mat_claude.h)

// _claude: global CG-iteration counter for the Osborn-style HMC cost function
// (Cost = C_S/(tau^2 <P_acc>), C_S = <N_CG>/traj). Every CG solve (solve / solve_multishift /
// solveAsync) adds its Krylov-iteration count k on exit. Reset before a trajectory, read after.
// Osborn et al., PoS LATTICE2024 052 "Automated tuning for HMC mass ratios".
inline unsigned long long g_matpoly_cg_iters = 0;
inline void reset_cg_iters(){ g_matpoly_cg_iters = 0; }
inline unsigned long long get_cg_iters(){ return g_matpoly_cg_iters; }

// device memory set
struct DeviceMemorySetN{
  // d_x, ...;
  CuC *d_tmp, *d_Mv0;
  CuC *d_x, *d_r0; // , *d_tmp, *d_tmp2;
  CuC *d_p, *d_q, *d_r;

  void allocate(){
    CUDA_CHECK(cudaMalloc(&d_tmp, Comp::N*CD));
    CUDA_CHECK(cudaMalloc(&d_Mv0, Comp::N*CD));
    CUDA_CHECK(cudaMemset(d_tmp, 0, Comp::N*CD));
    CUDA_CHECK(cudaMemset(d_Mv0, 0, Comp::N*CD));

    CUDA_CHECK(cudaMalloc(&d_x, Comp::N*CD));
    CUDA_CHECK(cudaMalloc(&d_r0, Comp::N*CD));
    CUDA_CHECK(cudaMemset(d_x, 0, Comp::N*CD));
    CUDA_CHECK(cudaMemset(d_r0, 0, Comp::N*CD));

    CUDA_CHECK(cudaMalloc(&d_p, Comp::N*CD));
    CUDA_CHECK(cudaMalloc(&d_q, Comp::N*CD));
    CUDA_CHECK(cudaMalloc(&d_r, Comp::N*CD));
    CUDA_CHECK(cudaMemset(d_p, 0, Comp::N*CD));
    CUDA_CHECK(cudaMemset(d_q, 0, Comp::N*CD));
    CUDA_CHECK(cudaMemset(d_r, 0, Comp::N*CD));
  }
  void deallocate(){
    CUDA_CHECK(cudaFree(d_tmp));
    CUDA_CHECK(cudaFree(d_Mv0));

    CUDA_CHECK(cudaFree(d_x));
    CUDA_CHECK(cudaFree(d_r0));

    CUDA_CHECK(cudaFree(d_r));
    CUDA_CHECK(cudaFree(d_p));
    CUDA_CHECK(cudaFree(d_q));
  }
};


DeviceMemorySetN d_MemorySets[Comp::NSTREAMS];


// template<typename T>
struct MatPoly{
  using T = CuC;

  std::vector<std::vector<const LinOp*>> vec_mats;
  std::vector<T> coeffs;

  cublasHandle_t handle;
  cudaStream_t stream;
  const bool is_external;
  const int istream;

  MatPoly()
    : is_external(false)
    , istream(-1)
  {
    handle = NULL;
    stream = NULL;
    CUBLAS_CHECK(cublasCreate(&handle));
  }

  MatPoly( cublasHandle_t handle_, cudaStream_t stream_, const int istream_)
    : is_external(true)
    , istream(istream_)
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
  void from_cpu( Complex* v, const Complex* v0 ) const {
    CuC *d_v, *d_v0;
    CUDA_CHECK(cudaMalloc(&d_v, N*CD));
    CUDA_CHECK(cudaMalloc(&d_v0, N*CD));
    CUDA_CHECK(cudaMemcpy(d_v0, reinterpret_cast<const CuC*>(v0), N*CD, H2D));

    on_gpu<N>( d_v, d_v0 );

    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(v), d_v, N*CD, D2H));
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


  // PREALLOC variant of on_gpu: the caller supplies scratch d_tmp/d_Mv0 (length N each) so that
  // on_gpu's per-call cudaMalloc/cudaFree -- DEVICE-WIDE syncs, ~10-100 us each -- are hoisted out
  // of the CG loop. Used by solve_multishift's hot, occupancy-starved seed matvec. Body is
  // otherwise byte-identical to on_gpu. Profile-gate A/B: swap the seed matvec call between on_gpu
  // and on_gpu_pre to size the cudaMalloc/cudaFree overhead (see solve_multishift seed matvec).
  template<Idx N> __host__
  void on_gpu_pre(CuC* d_v, const CuC* d_v0, CuC* d_tmp, CuC* d_Mv0) const {
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
  }


  template<Idx N> __host__
  void on_gpuAsync(CuC* d_v, const CuC* d_v0) const {
    CuC *d_tmp = d_MemorySets[istream].d_tmp;
    CuC *d_Mv0 = d_MemorySets[istream].d_Mv0;

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
  void solve(Complex* x, const Complex* b,
             const double tol=1.0e-13, const int maxiter=1e8) const {
    // CG
    CuC *d_x, *d_r; // , *d_tmp, *d_tmp2;
    CUDA_CHECK(cudaMalloc(&d_x, N*CD));
    CUDA_CHECK(cudaMalloc(&d_r, N*CD));
    CUDA_CHECK(cudaMemcpy(d_r, reinterpret_cast<const CuC*>(b), N*CD, H2D));

    solve<N>(d_x, d_r, tol, maxiter);

    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(x), d_x, N*CD, D2H));
    CUDA_CHECK(cudaFree(d_x));
    CUDA_CHECK(cudaFree(d_r));
  }


  // template<Idx N> __host__
  // void solveAsync(std::vector<Complex>& x, const std::vector<Complex>& b,
  //            const double tol=1.0e-13, const int maxiter=1e8) const {
  //   // CG
  //   assert(istream>=0);
  //   CuC *d_x = d_MemorySets[istream].d_x;
  //   CuC *d_r0 = d_MemorySets[istream].d_r0;

  //   CUDA_CHECK(cudaMemcpyAsync(d_r0, reinterpret_cast<const CuC*>(b.data()), N*CD, H2D, stream));

  //   solveAsync<N>(d_x, d_r0, tol, maxiter);

  //   CUDA_CHECK(cudaMemcpyAsync(reinterpret_cast<CuC*>(x.data()), d_x, N*CD, D2H, stream));

  //   CUDA_CHECK( cudaStreamSynchronize(stream) );
  // }


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

    // zero RHS: exact solution is the memset-zero d_x. The relative test
    // mu < mu_crit = tol^2*||b||^2 degenerates to 0 < 0 when ||b||=0, so guard
    // it explicitly to avoid al = mu/gam = 0/0 = NaN in the CG loop.
    if(abs(b_norm_sq) < 1.0e-14 || mu<mu_crit) {
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
      g_matpoly_cg_iters += (unsigned long long)k;   // _claude: cost accounting
#ifdef IsVerbose2
      std::clog << "SOLVER:       #iterations: " << k << std::endl;
      std::clog << "SOLVER:       mu =         " << mu << std::endl;
#endif
    }

    CUDA_CHECK(cudaFree(d_r));
    CUDA_CHECK(cudaFree(d_p));
    CUDA_CHECK(cudaFree(d_q));
  }


  // ======================================================================
  // Multi-shift CG (Jegerlehner, hep-lat/9612014). Solves
  //     (A + sigma_m) X_m = b   for all m = 0..npole-1   in ONE Krylov pass,
  // where A = this MatPoly (the SEED matrix, built WITHOUT any shift, e.g. the
  // (1/lambda_max^2) D_W^dag D_W term) and sigma_m > 0 are scalar shifts. The
  // expensive matvec A p is done ONCE per iteration and shared by every pole;
  // the per-pole solutions ride cheap scalar recurrences. Seed = smallest shift
  // (slowest system) => when the seed meets tol all poles have converged.
  //   d_X   : N*npole COLUMN-MAJOR block [m*N + i] (output; zeroed here).
  //   d_b   : single RHS (length N), shared by all poles.
  //   sigma : HOST array of npole shifts.
  // Host reference (validated, C1): test_multishift_claude.cpp.
  //
  // ASYNC NOTE: this path is intentionally fully SYNCHRONOUS on the default
  // stream. The per-pole coeff host buffers (alm/zeta_new/betm) are overwritten
  // EVERY iteration, so their H2D copies MUST be synchronous cudaMemcpy -- an
  // async copy could race the next iteration's overwrite. The blocking
  // host-pointer dots (dot/dot2self) serialize the kernels on the default
  // stream. Do NOT convert these to *Async without pinning the coeff buffers
  // AND adding a per-copy stream sync before the next overwrite.
  // ======================================================================
  template<Idx N> __host__
  void solve_multishift(CuC* d_X, const CuC* d_b,
                        const double* sigma, const int npole,
                        const double tol=1.0e-13, const int maxiter=1e8) const {
    // seed = smallest shift; relative shifts hat_m = sigma_m - sigma_seed >= 0
    int seed=0; for(int m=1;m<npole;m++) if(sigma[m]<sigma[seed]) seed=m;
    const double sig0 = sigma[seed];
    std::vector<double> hat(npole);
    for(int m=0;m<npole;m++) hat[m] = sigma[m]-sig0;

    // seed work vectors (length N); per-pole directions d_pm (N*npole block);
    // per-pole coeff arrays d_alm/d_zeta/d_betm (length npole).
    CuC *d_p0,*d_q,*d_r,*d_pm;
    CUDA_CHECK(cudaMalloc(&d_p0, N*CD));
    CUDA_CHECK(cudaMalloc(&d_q,  N*CD));
    CUDA_CHECK(cudaMalloc(&d_r,  N*CD));
    CUDA_CHECK(cudaMalloc(&d_pm, (size_t)N*npole*CD));
    double *d_alm,*d_zeta,*d_betm;
    CUDA_CHECK(cudaMalloc(&d_alm,  npole*sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_zeta, npole*sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_betm, npole*sizeof(double)));

    // PREALLOC scratch for the seed matvec (on_gpu_pre): hoists on_gpu's per-call cudaMalloc/Free
    // out of the CG loop. Allocated once per solve. Profile-gate A/B (see the seed matvec below).
    CuC *d_tmp_mv, *d_Mv0_mv;
    CUDA_CHECK(cudaMalloc(&d_tmp_mv, N*CD));
    CUDA_CHECK(cudaMalloc(&d_Mv0_mv, N*CD));

    // grid sized for the N*npole block work (NOT the N-sized NBlocks macro)
    const Idx total = (Idx)N*npole;
    const int nb_blk = (total + NThreadsPerBlock)/NThreadsPerBlock;

    // init: X=0; r=p0=b; p_m=b for every column
    CUDA_CHECK(cudaMemset(d_X, 0, (size_t)N*npole*CD));
    CUDA_CHECK(cudaMemcpy(d_r,  d_b, N*CD, D2D));
    CUDA_CHECK(cudaMemcpy(d_p0, d_b, N*CD, D2D));
    for(int m=0;m<npole;m++) CUDA_CHECK(cudaMemcpy(d_pm + (size_t)m*N, d_b, N*CD, D2D));

    double mu; dot2self<N>(&mu, d_r); assert(mu>=0.0);
    const double b_norm_sq = mu;
    const double mu_crit = tol*tol*b_norm_sq;

    // host per-pole state: zeta_m^j (zeta), zeta_m^{j-1} (zeta_old); seed al/bet history
    std::vector<double> zeta(npole,1.0), zeta_old(npole,1.0);
    std::vector<double> alm(npole), zeta_new(npole), betm(npole);
    double al_old=1.0, bet_old=0.0;   // al_{-1}=1, bet_{-1}=0
    // freeze-converged: once a shift's residual zeta_m^2 * mu < mu_crit it is CONVERGED.
    // Freeze it (stop updating zeta_m/x_m/p_m) so the recurrence never drives zeta_m -> 0,
    // which underflows to 0/0 = NaN for the fast large-shift poles (the bug we diagnosed).
    std::vector<char> frozen(npole, 0);

    // zero RHS guard (cf. solve()): X stays 0 when ||b||=0.
    if(abs(b_norm_sq) < 1.0e-14 || mu<mu_crit){
#ifdef IsVerbose2
      std::clog << "NO SOLVE (multishift)" << std::endl;
#endif
    }
    else{
      int k=0;
      for(; k<maxiter; ++k){
        // seed matvec: q = A p0 + sig0 p0   (one matvec shared by all poles)
        this->on_gpu<N>(d_q, d_p0);
        // PROFILE-GATE A/B (OFF -- free-field wash: only ~370 cudaMalloc removed, dominated by a
        // one-time startup alloc): preallocated-scratch seed matvec hoists the per-matvec
        // cudaMalloc/cudaFree out of the CG loop. May help on long thermalized solves; swap to the
        // line below to re-enable.
        // this->on_gpu_pre<N>(d_q, d_p0, d_tmp_mv, d_Mv0_mv);
        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock>>>(d_q, cplx(sig0), d_p0, d_q);
        CuC gam; dot<N>(&gam, d_p0, d_q);     // host-pointer dot => blocks (sync point)
        // BREAKDOWN GUARD: the seed curvature <p0|A p0> is > 0 for SPD A+sigma_s; roundoff
        // drives it to <=0 only when the orthogonal search direction collapses (seed at its
        // residual floor). Per request this is a HARD FAILURE (core dump), not a silent
        // return -- dump the per-shift state at the breakdown moment, then assert(false).
        const double gam_re = real(gam);
        if(!(gam_re > 0.0) || !std::isfinite(gam_re) || !std::isfinite(mu)){
          std::clog << "MULTISHIFT BREAKDOWN at iter " << k
                    << ": gam = " << gam_re << ", mu = " << mu << std::endl;
#ifdef IsVerbose2
          { const double bnorm = std::sqrt(b_norm_sq);
            for(int m=0; m<npole; m++){
              this->on_gpu<N>(d_q, d_X + (size_t)m*N);                                              // A Z_m
              Taxpy<CuC,N><<<NBlocks,NThreadsPerBlock>>>(d_q, cplx(sigma[m]), d_X+(size_t)m*N, d_q); // + sigma_m Z_m
              Taxpy<CuC,N><<<NBlocks,NThreadsPerBlock>>>(d_q, cplx(-1.0), d_b, d_q);                 // - b
              double rn=0.0, zn=0.0;
              CUBLAS_CHECK( cublasDznrm2(handle, N, d_q, 1, &rn) );
              CUBLAS_CHECK( cublasDznrm2(handle, N, d_X + (size_t)m*N, 1, &zn) );
              std::clog << "MULTISHIFT-CHK(breakdown): m=" << m << " sigma=" << sigma[m]
                        << " true_resid=" << (bnorm>0.0 ? rn/bnorm : rn)
                        << " ||Z_m||=" << zn << std::endl;
            } }
#endif
          assert(false);   // core dump on breakdown (asserts are ON in this project)
        }
        const double al = mu/gam_re;

        // host zeta recurrence + al_m (uses al and the previous al_old/bet_old).
        // Freeze a shift once its residual zeta_m^2 * mu falls below mu_crit: its x_m is
        // already converged, and continuing would drive zeta_m -> 0 -> 0/0 = NaN.
        for(int m=0;m<npole;m++){
          if(!frozen[m] && zeta[m]*zeta[m]*mu < mu_crit) frozen[m]=1;
          if(frozen[m]){ alm[m]=0.0; zeta_new[m]=0.0; continue; }  // no x update; zeta_m -> 0
          const double denom = al*bet_old*(zeta_old[m]-zeta[m])
                             + zeta_old[m]*al_old*(1.0 + hat[m]*al);
          zeta_new[m] = (zeta[m]*zeta_old[m]*al_old)/denom;
          alm[m]      = al * zeta_new[m]/zeta[m];
        }
        CUDA_CHECK(cudaMemcpy(d_alm, alm.data(), npole*sizeof(double), H2D)); // sync (ASYNC NOTE)
        multishift_x_update<N><<<nb_blk, NThreadsPerBlock>>>(d_X, d_alm, d_pm, npole);

        // seed residual + convergence (seed = slowest pole drives all)
        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock>>>(d_r, cplx(-al), d_q, d_r);
        double mu_new; dot2self<N>(&mu_new, d_r); assert(mu_new>=0.0);
        if(mu_new<mu_crit || std::isnan(mu_new)) break;   // X_m already updated this iter => final
        const double bet = mu_new/mu;

        // host bet_m; p_m = zeta_new r + bet_m p_m; seed p0 = r + bet p0.
        // Frozen shifts: zeta_new[m]=0 (above) and bet_m=0 => p_m zeroed (unused), no recurrence.
        for(int m=0;m<npole;m++){
          if(frozen[m]){ betm[m]=0.0; continue; }
          const double ratio=zeta_new[m]/zeta[m]; betm[m]=bet*ratio*ratio;
        }
        CUDA_CHECK(cudaMemcpy(d_zeta, zeta_new.data(), npole*sizeof(double), H2D)); // sync
        CUDA_CHECK(cudaMemcpy(d_betm, betm.data(),     npole*sizeof(double), H2D)); // sync
        multishift_p_update<N><<<nb_blk, NThreadsPerBlock>>>(d_pm, d_zeta, d_r, d_betm, npole);
        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock>>>(d_p0, cplx(bet), d_p0, d_r);

        // roll history
        zeta_old=zeta; zeta=zeta_new; al_old=al; bet_old=bet; mu=mu_new;

#ifdef IsVerbose2
        if(k%100==0) std::clog << "MULTISHIFT:   #iter: " << k << ", mu = " << mu << std::endl;
#endif
      }
      g_matpoly_cg_iters += (unsigned long long)k;   // _claude: cost accounting (one Krylov pass, all poles)
#ifdef IsVerbose2
      std::clog << "MULTISHIFT:   #iter (final): " << k << ", mu = " << mu << std::endl;
#endif
    }

#ifdef IsVerbose2
    {
      // DEBUG: per-shift TRUE residual ||(A+sigma_m) Z_m - b||/||b|| and ||Z_m||, to detect
      // shifted-solution DRIFT/overflow despite seed convergence (the multi-shift weakness:
      // fast poles get over-corrected by the zeta recurrence while the slow seed grinds on).
      // Costs npole extra matvecs -- verbose-only. d_q is reused as scratch (freed below).
      const double bnorm = std::sqrt(b_norm_sq);
      for(int m=0; m<npole; m++){
        this->on_gpu<N>(d_q, d_X + (size_t)m*N);                                              // A Z_m
        Taxpy<CuC,N><<<NBlocks,NThreadsPerBlock>>>(d_q, cplx(sigma[m]), d_X+(size_t)m*N, d_q); // + sigma_m Z_m
        Taxpy<CuC,N><<<NBlocks,NThreadsPerBlock>>>(d_q, cplx(-1.0), d_b, d_q);                 // - b
        double rn=0.0, zn=0.0;
        CUBLAS_CHECK( cublasDznrm2(handle, N, d_q, 1, &rn) );
        CUBLAS_CHECK( cublasDznrm2(handle, N, d_X + (size_t)m*N, 1, &zn) );
        std::clog << "MULTISHIFT-CHK: m=" << m << " sigma=" << sigma[m]
                  << " true_resid=" << (bnorm>0.0 ? rn/bnorm : rn)
                  << " ||Z_m||=" << zn << std::endl;
      }
    }
#endif

    CUDA_CHECK(cudaFree(d_p0)); CUDA_CHECK(cudaFree(d_q));
    CUDA_CHECK(cudaFree(d_r));  CUDA_CHECK(cudaFree(d_pm));
    CUDA_CHECK(cudaFree(d_alm)); CUDA_CHECK(cudaFree(d_zeta)); CUDA_CHECK(cudaFree(d_betm));
    CUDA_CHECK(cudaFree(d_tmp_mv)); CUDA_CHECK(cudaFree(d_Mv0_mv));
  }


  template<Idx N> __host__
  void solveAsync(CuC* d_x, const CuC* d_b,
                  const double tol=1.0e-13, const int maxiter=1e8) const {
    // CG
    assert(istream>=0);
    CuC *d_p = d_MemorySets[istream].d_p;
    CuC *d_q = d_MemorySets[istream].d_q;
    CuC *d_r = d_MemorySets[istream].d_r;

    CUDA_CHECK(cudaMemsetAsync(d_x, 0, N*CD, stream));
    CUDA_CHECK(cudaMemsetAsync(d_q, 0, N*CD, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_r, d_b, N*CD, D2D, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_p, d_r, N*CD, D2D, stream));

    double mu, b_norm_sq;
    this->dot2selfAsync<N>(&mu, d_r);
    this->dot2selfAsync<N>(&b_norm_sq, d_r);
    CUDA_CHECK( cudaStreamSynchronize(stream) );

    assert(mu>=0.0);
    assert(b_norm_sq>=0.0);
    double mu_old = mu;
    double mu_crit = tol*tol*b_norm_sq;

    // zero RHS guard (see solve()): avoid 0/0 = NaN when ||b|| = 0.
    if(abs(b_norm_sq) < 1.0e-14 || mu<mu_crit) {
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
        CUDA_CHECK( cudaStreamSynchronize(stream) );

        const CuC al = mu/gam;
        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock, 0, stream>>>(d_x, al, d_p, d_x);
        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock, 0, stream>>>(d_r, -al, d_q, d_r);

	this->dot2selfAsync<N>(&mu, d_r);
        CUDA_CHECK( cudaStreamSynchronize(stream) );

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
      g_matpoly_cg_iters += (unsigned long long)k;   // _claude: cost accounting
#ifdef IsVerbose
      std::clog << "SOLVER:       #iterations: " << k << std::endl;
      std::clog << "SOLVER:       mu =         " << mu << std::endl;
#endif
    }
  }


//   template<Idx N>
//   void bicgstab( std::vector<Complex>& v, const std::vector<Complex>& v0,
//                  const std::vector<Complex>& rc,
//                  const double tol=1.0e-6, const int maxiter=1e8,
//                  const double eps=1.0e-8 ) const {
//     CuC *d_v, *d_v0, *d_rc;
//     CUDA_CHECK(cudaMalloc(&d_v, N*CD));
//     CUDA_CHECK(cudaMalloc(&d_v0, N*CD));
//     CUDA_CHECK(cudaMalloc(&d_rc, N*CD));

//     CUDA_CHECK(cudaMemcpy(d_v0, reinterpret_cast<const CuC*>(v0.data()), N*CD, H2D));
//     CUDA_CHECK(cudaMemcpy(d_rc, reinterpret_cast<const CuC*>(rc.data()), N*CD, H2D));

//     bicgstab<N>( d_v, d_v0, d_rc, tol, maxiter, eps );

//     CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(v.data()), d_v, N*CD, D2H));

//     CUDA_CHECK(cudaFree(d_v));
//     CUDA_CHECK(cudaFree(d_v0));
//     CUDA_CHECK(cudaFree(d_rc));
//   }



//   template<Idx N> __host__
//   void bicgstab(CuC* d_x, const CuC* d_b, CuC* d_rc,
//                 // double& mu=1.0,
//                 const double tol=1.0e-4, const int maxiter=1e8,
//                 const double eps=1.0e-8 ) const {
//     CuC *d_p, *d_t, *d_Ap, *d_At, *d_r;
//     CUDA_CHECK(cudaMalloc(&d_p, N*CD));
//     CUDA_CHECK(cudaMalloc(&d_t, N*CD));
//     CUDA_CHECK(cudaMalloc(&d_Ap, N*CD));
//     CUDA_CHECK(cudaMalloc(&d_At, N*CD));
//     CUDA_CHECK(cudaMalloc(&d_r, N*CD));
//     // CUDA_CHECK(cudaMalloc(&d_rc, N*CD));

//     CUDA_CHECK(cudaMemset(d_x, 0, N*CD));
//     CUDA_CHECK(cudaMemset(d_p, 0, N*CD));
//     CUDA_CHECK(cudaMemset(d_t, 0, N*CD));
//     CUDA_CHECK(cudaMemset(d_Ap, 0, N*CD));
//     CUDA_CHECK(cudaMemset(d_At, 0, N*CD));
//     CUDA_CHECK(cudaMemset(d_r, 0, N*CD));
//     // CUDA_CHECK(cudaMemset(d_rc, 0, N*CD));

//     CUDA_CHECK(cudaMemcpy(d_x, d_b, N*CD, D2D));
//     this->on_gpu<N>(d_Ap, d_x);
//     Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock, 0>>>(d_r, -1.0, d_Ap, d_b);
//     CUDA_CHECK(cudaMemcpy(d_rc, d_b, N*CD, D2D));
//     Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock, 0>>>(d_rc, -1.0, d_Ap, d_rc);

//     CuC gam;
//     this->dot<N>(&gam, d_rc, d_r);
//     assert( abs(gam)>tol );

//     CUDA_CHECK(cudaMemcpy(d_p, d_r, N*CD, D2D));

//     double mu;
//     this->dot2self<N>(&mu, d_r);
//     assert(mu>=0.0);
//     double mu_old = mu;

//     double nm;

//     double b_norm_sq;
//     this->dot2self<N>(&b_norm_sq, d_r);
//     assert(b_norm_sq>=0.0);
//     double mu_crit = tol*tol*b_norm_sq;

//     if(mu<mu_crit) {
// #ifdef IsVerbose
//       std::clog << "NO SOLVE" << std::endl;
// #endif
//     }
//     else{
//       int n = 0;
//       for (; n < maxiter; n++) {
//         this->on_gpu<N>(d_Ap, d_p);

//         CuC gam1, gam2;
//         this->dot<N>(&gam1, d_rc, d_r);
//         this->dot<N>(&gam2, d_rc, d_Ap);
//         const CuC al = gam1/gam2;

//         Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock, 0>>>(d_t, -al, d_Ap, d_r);
//         this->dot2self<N>(&nm, d_t);
//         if(nm<eps) break;
//         this->on_gpu<N>(d_At, d_t);

//         this->dot<N>(&gam1, d_t, d_At);
//         this->dot<N>(&gam2, d_At, d_At);
//         const CuC om = gam1/gam2;

//         Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock, 0>>>(d_x, al, d_p, d_x);
//         Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock, 0>>>(d_x, om, d_t, d_x);

//         this->dot<N>(&gam2, d_rc, d_r);
//         Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock, 0>>>(d_r, -om, d_At, d_t);
//         this->dot<N>(&gam1, d_rc, d_r);
//         const CuC bet = al/om * gam1/gam2;

//         this->dot2self<N>(&mu, d_r);
//         assert(mu>=0.0);
//         if(mu<mu_crit || std::isnan(mu)) break;
//         mu_old = mu;

//         Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock, 0>>>(d_p, -om, d_Ap, d_p);
//         Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock, 0>>>(d_p, bet, d_p, d_r);

//         // restart
//         this->dot<N>(&gam, d_rc, d_r);
//         if( abs(gam)<eps || 2.0*mu_old < mu ){
//           CUDA_CHECK(cudaMemcpy(d_rc, d_r, N*CD, D2D));
//           CUDA_CHECK(cudaMemcpy(d_p, d_r, N*CD, D2D));
// #ifdef IsVerbose2
// 	  std::clog << "SOLVER:      restart " << std::endl;
// #endif
//         }

//         // std::cout << "debug. mu = " << mu << std::endl;

//         if(n%10==0) {
// #ifdef IsVerbose2
// 	  std::clog << "SOLVER:       #iterations: " << n << ", mu =         " << mu << std::endl;
// #endif
// 	}

//       }

//       if (n == maxiter) {
//         std::cerr << "BiCGStab iteration did not converge." << std::endl;
//         // std::cerr << "v = " << v << std::endl;
//       }
//     }

//     CUDA_CHECK(cudaFree(d_p));
//     CUDA_CHECK(cudaFree(d_Ap));
//     CUDA_CHECK(cudaFree(d_At));
//     CUDA_CHECK(cudaFree(d_r));
//     // CUDA_CHECK(cudaFree(d_rc));
//     CUDA_CHECK(cudaFree(d_t));
//     // CUDA_CH/ECK( cudaStreamSynchronize(stream) );

//     // return x;
//   }



};


