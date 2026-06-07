#pragma once

// blocked_mat_claude.h -- BlockMemPool<N,NSTACK> + BlockedMat<Idx N, int NSTACK, class Op> (C7).
// ALL mrhs (multi-RHS block) functionality. BlockMemPool OWNS the device scratch (RAII; ctor allocates,
// with_pole_blocks gates the three big N*NSTACK*npole arrays); BlockedMat holds a const ref to the
// overlap operator Op (for M_DW/M_DWH/A/C/mass/lambda_max/k/cp/size) + ALIAS pointers into a pool
// (shared via the 2-arg ctor, or an owned full pool via the 1-arg ctor) and provides:
//   - solve_multishift_block : inner Zolotarev sign-function solve (seed (1/lmax^2) M_DWH M_DW)
//   - solve_shift_block      : single-shift block CG over A+sigma (Term B inner pole solve, C6f-b0);
//                              needs ONLY the N*NSTACK scratch => can run on a lean (no-pole-block) pool
//   - mult / adj / DDH       : the overlap block applies (D_ov+m, its adjoint, (D+m)^dag(D+m)) on an
//                              N*NSTACK block (were OverlapWMass::*_deviceAsyncLaunch_ms_block)
//   - solve_sq               : outer per-column block CG over DDH (op_Dmsq), bit-identical per column
// NSTACK is compile-time => scratch exactly sized; npole = Op.size-1 (runtime, from the ctor's Op).
// Created ONCE (the ctor alloc is heavy); the applies run many times per solve so a per-call object
// would re-malloc the scratch -- hence a persistent object, not a local.
//
// Reference (algorithm): B. Jegerlehner, "Krylov space solvers for shifted linear systems",
// hep-lat/9612014 (inner multishift) + standard CG (outer). Block math validated bit-identical
// (C6a-d). Method names (solve_multishift_block/mult/adj/DDH/solve) avoid the free block-kernel names
// (mult_block / block_reduce_poles / axpy_col / axpby_col / axpy_uniform / multishift_x/p_update_block
// / axpy_col_c) in multishift_block_kernels_claude.h.
//
// Include AFTER gpu_header.h, sparse_matrix.h, multishift_block_kernels_claude.h, and the overlap
// header (for the Op type). Op must expose: CSR<N> M_DW, M_DWH; double lambda_max, k, C; int size;
// std::vector<double> cp, A; Complex mass.

// BlockMemPool<Idx N, int NSTACK> -- OWNS the device block scratch so multiple BlockedMat clients can
// SHARE one preallocated pool (pass it by reference at the BlockedMat ctor; its lifetime must outlive
// the BlockedMat). with_pole_blocks gates the THREE big N*NSTACK*npole arrays (d_pm/d_Zblk/d_Yblk):
//   - solve_multishift_block + mult/adj/DDH/solve_sq NEED them (with_pole_blocks=true, the default);
//   - solve_shift_block (single-shift block CG; ConservedCurrent::apply_k_block_t at NSTACK=Nt) does
//     NOT, so its pool is built with with_pole_blocks=false to skip ~3*N*Nt*npole*16 B (~1.9 GB at L1,
//     Nt=128). Such a lean pool may ONLY drive solve_shift_block (not the pole-block methods).
// npole is runtime (from the Op) -- pass D.size-1; a BlockedMat asserts its npole matches the pool's.
template<Idx N, int NSTACK>
struct BlockMemPool {
  const int npole;
  const bool has_pole_blocks;
  // ---- in-use lease (single-thread re-entrancy guard; see PoolGuard) ----
  // catches the overwrite hazard when TWO BlockedMat clients share one pool and interleave. Tracks the
  // OWNER (not just a bool) so legitimate self-nesting (solve_sq -> DDH -> mult -> solve_multishift_block,
  // all the same BlockedMat) does NOT false-positive. mutable: the const block methods take the lease.
  mutable const void* lessee = nullptr;
  mutable int lease_depth = 0;
  CuC *d_p0,*d_q,*d_r,*d_seedtmp;            // N*NSTACK : inner-solve seed/outer vectors
  CuC *d_pm,*d_Zblk,*d_Yblk;                 // N*NSTACK*npole : per-(rhs,pole) blocks (only if pole-blocks)
  CuC *d_mp,*d_mm;                           // N*NSTACK : DDH temporaries
  CuC *d_ocg_p,*d_ocg_q,*d_ocg_r;            // N*NSTACK : outer/shift block-CG p/q/r
  CuC *d_blk_in,*d_blk_out;                  // N*NSTACK : host-block solve I/O staging
  CuC *d_alc;                                // NSTACK : per-column complex al
  double *d_alm,*d_zeta,*d_betm;             // NSTACK*npole : inner per-(c,m) scalars
  double *d_sc,*d_sc2;                       // NSTACK : per-column real scalars

  BlockMemPool(int npole_, bool with_pole_blocks=true)
    : npole(npole_), has_pole_blocks(with_pole_blocks),
      d_pm(nullptr), d_Zblk(nullptr), d_Yblk(nullptr)
  {
    const size_t Ns  = (size_t)N*NSTACK;
    const size_t Nsp = (size_t)N*NSTACK*npole;
    CUDA_CHECK(cudaMalloc(&d_p0, Ns*CD));      CUDA_CHECK(cudaMalloc(&d_q, Ns*CD));
    CUDA_CHECK(cudaMalloc(&d_r,  Ns*CD));      CUDA_CHECK(cudaMalloc(&d_seedtmp, Ns*CD));
    CUDA_CHECK(cudaMalloc(&d_mp, Ns*CD));      CUDA_CHECK(cudaMalloc(&d_mm, Ns*CD));
    CUDA_CHECK(cudaMalloc(&d_ocg_p, Ns*CD));   CUDA_CHECK(cudaMalloc(&d_ocg_q, Ns*CD));
    CUDA_CHECK(cudaMalloc(&d_ocg_r, Ns*CD));   CUDA_CHECK(cudaMalloc(&d_alc, (size_t)NSTACK*CD));
    CUDA_CHECK(cudaMalloc(&d_blk_in, Ns*CD));  CUDA_CHECK(cudaMalloc(&d_blk_out, Ns*CD));
    CUDA_CHECK(cudaMalloc(&d_alm,  (size_t)NSTACK*npole*DB));
    CUDA_CHECK(cudaMalloc(&d_zeta, (size_t)NSTACK*npole*DB));
    CUDA_CHECK(cudaMalloc(&d_betm, (size_t)NSTACK*npole*DB));
    CUDA_CHECK(cudaMalloc(&d_sc,   (size_t)NSTACK*DB));
    CUDA_CHECK(cudaMalloc(&d_sc2,  (size_t)NSTACK*DB));
    if(with_pole_blocks){
      CUDA_CHECK(cudaMalloc(&d_pm, Nsp*CD));   CUDA_CHECK(cudaMalloc(&d_Zblk, Nsp*CD));
      CUDA_CHECK(cudaMalloc(&d_Yblk, Nsp*CD));
    }
  }
  ~BlockMemPool(){
    CUDA_CHECK(cudaFree(d_p0)); CUDA_CHECK(cudaFree(d_q)); CUDA_CHECK(cudaFree(d_r));
    CUDA_CHECK(cudaFree(d_seedtmp)); CUDA_CHECK(cudaFree(d_mp)); CUDA_CHECK(cudaFree(d_mm));
    CUDA_CHECK(cudaFree(d_ocg_p)); CUDA_CHECK(cudaFree(d_ocg_q)); CUDA_CHECK(cudaFree(d_ocg_r));
    CUDA_CHECK(cudaFree(d_blk_in)); CUDA_CHECK(cudaFree(d_blk_out));
    CUDA_CHECK(cudaFree(d_alc)); CUDA_CHECK(cudaFree(d_alm)); CUDA_CHECK(cudaFree(d_zeta));
    CUDA_CHECK(cudaFree(d_betm)); CUDA_CHECK(cudaFree(d_sc)); CUDA_CHECK(cudaFree(d_sc2));
    if(has_pole_blocks){ CUDA_CHECK(cudaFree(d_pm)); CUDA_CHECK(cudaFree(d_Zblk)); CUDA_CHECK(cudaFree(d_Yblk)); }
  }
  BlockMemPool(const BlockMemPool&)=delete;
  BlockMemPool& operator=(const BlockMemPool&)=delete;
};

// RAII lease on a BlockMemPool: assert-fails if a DIFFERENT owner already holds the pool's scratch
// (the shared-pool overwrite hazard). Same-owner re-entry (nested block-method calls) just counts the
// depth. assert() compiles out under NDEBUG => zero release-build cost; never touches buffers, so
// bit-identity is unaffected. Constructed at the top of every public BlockedMat method (who = this).
template<Idx N, int NSTACK>
struct PoolGuard {
  const BlockMemPool<N,NSTACK>& p;
  PoolGuard(const BlockMemPool<N,NSTACK>& p_, const void* who) : p(p_) {
    assert((p.lease_depth==0 || p.lessee==who) &&
           "BlockMemPool already in use by another client -- shared scratch would be overwritten");
    p.lessee = who; ++p.lease_depth;
  }
  ~PoolGuard(){ if(--p.lease_depth==0) p.lessee=nullptr; }
  PoolGuard(const PoolGuard&)=delete;
  PoolGuard& operator=(const PoolGuard&)=delete;
};

template<Idx N, int NSTACK, class Op>
struct BlockedMat {
  const Op& D;
  cublasHandle_t handle;
  const int npole;                          // = D.size - 1
  BlockMemPool<N,NSTACK>* owned_pool;       // non-null only when WE allocated the pool (1-arg ctor)
  BlockMemPool<N,NSTACK>* mempool;          // the active pool (external or owned); for the lease guard
  // ---- block scratch: ALIAS pointers into the shared/owned BlockMemPool (NOT owned here) ----
  CuC *d_p0,*d_q,*d_r,*d_seedtmp;
  CuC *d_pm,*d_Zblk,*d_Yblk;
  CuC *d_mp,*d_mm;
  CuC *d_ocg_p,*d_ocg_q,*d_ocg_r;
  CuC *d_blk_in,*d_blk_out;
  CuC *d_alc;
  double *d_alm,*d_zeta,*d_betm;
  double *d_sc,*d_sc2;

  void bind(BlockMemPool<N,NSTACK>& p){
    d_p0=p.d_p0; d_q=p.d_q; d_r=p.d_r; d_seedtmp=p.d_seedtmp;
    d_pm=p.d_pm; d_Zblk=p.d_Zblk; d_Yblk=p.d_Yblk;
    d_mp=p.d_mp; d_mm=p.d_mm;
    d_ocg_p=p.d_ocg_p; d_ocg_q=p.d_ocg_q; d_ocg_r=p.d_ocg_r;
    d_blk_in=p.d_blk_in; d_blk_out=p.d_blk_out; d_alc=p.d_alc;
    d_alm=p.d_alm; d_zeta=p.d_zeta; d_betm=p.d_betm; d_sc=p.d_sc; d_sc2=p.d_sc2;
  }

  // share an EXTERNAL preallocated pool (lifetime must outlive this BlockedMat)
  BlockedMat(const Op& D_, BlockMemPool<N,NSTACK>& pool_)
    : D(D_), npole(D_.size-1), owned_pool(nullptr)
  {
    assert(pool_.npole==npole);
    CUBLAS_CHECK(cublasCreate(&handle));
    mempool = &pool_;
    bind(pool_);
  }
  // convenience: allocate an OWNED full pool (with_pole_blocks=true). Backward compatible with the
  // original single-arg ctor -- existing call sites (C6 tests, jj_corr_mrhs) are unchanged.
  explicit BlockedMat(const Op& D_)
    : D(D_), npole(D_.size-1), owned_pool(new BlockMemPool<N,NSTACK>(D_.size-1, true))
  {
    CUBLAS_CHECK(cublasCreate(&handle));
    mempool = owned_pool;
    bind(*owned_pool);
  }
  ~BlockedMat(){
    CUBLAS_CHECK(cublasDestroy(handle));
    delete owned_pool;
  }
  BlockedMat(const BlockedMat&) = delete;
  BlockedMat& operator=(const BlockedMat&) = delete;

  static int ng(size_t tot){ return (int)((tot + NThreadsPerBlock - 1)/NThreadsPerBlock); }

  // ===== inner Zolotarev sign-function solve (was MatPoly::solve_multishift_block) =====
  // (A + sigma_m) X_{c,m} = b_c for all RHS c, poles m; seed A = (1/lmax^2) M_DWH M_DW (apply M_DW
  // then M_DWH), shifts sigma_m = -k^2/cp[m] from D. Per-column + per-pole freeze-converged (avoids
  // the gam/zeta 0/0); shared worst-column stop. d_X: N*NSTACK*npole; d_B: N*NSTACK.
  void solve_multishift_block(CuC* d_X, const CuC* d_B, const double tol=1.0e-13, const int maxiter=1e8) const {
    PoolGuard<N,NSTACK> g(*mempool, this);
    const double coeff_seed = 1.0/(D.lambda_max*D.lambda_max);
    std::vector<double> sigma(npole);
    for(int m=1;m<D.size;m++) sigma[m-1] = -D.k*D.k/D.cp[m];
    int seed=0; for(int m=1;m<npole;m++) if(sigma[m]<sigma[seed]) seed=m;
    const double sig0 = sigma[seed];
    std::vector<double> hat(npole); for(int m=0;m<npole;m++) hat[m]=sigma[m]-sig0;

    const size_t Ns=(size_t)N*NSTACK, Nsp=(size_t)N*NSTACK*npole;
    const int ncolp=NSTACK*npole;
    const int g_seed=ng(Ns), g_blk=ng(Nsp);

    CUDA_CHECK(cudaMemset(d_X, 0, Nsp*CD));
    CUDA_CHECK(cudaMemcpy(d_r,  d_B, Ns*CD, D2D));
    CUDA_CHECK(cudaMemcpy(d_p0, d_B, Ns*CD, D2D));
    for(int c=0;c<NSTACK;c++) for(int m=0;m<npole;m++)
      CUDA_CHECK(cudaMemcpy(d_pm + (size_t)(c*npole+m)*N, d_B + (size_t)c*N, N*CD, D2D));

    auto dotself=[&](const CuC* x,int c){ CuC d; CUBLAS_CHECK(cublasZdotc(handle,N,x+(size_t)c*N,1,x+(size_t)c*N,1,&d)); return real(d); };
    auto dotcol =[&](const CuC* x,const CuC* y,int c){ CuC d; CUBLAS_CHECK(cublasZdotc(handle,N,x+(size_t)c*N,1,y+(size_t)c*N,1,&d)); return real(d); };

    std::vector<double> mu(NSTACK), b_norm_sq(NSTACK), mu_crit(NSTACK), mu_new(NSTACK);
    std::vector<double> al(NSTACK,0.0), bet(NSTACK,0.0), al_old(NSTACK,1.0), bet_old(NSTACK,0.0);
    std::vector<char> colfrozen(NSTACK,0);
    std::vector<double> sc(NSTACK), sc2(NSTACK);
    for(int c=0;c<NSTACK;c++){ mu[c]=dotself(d_r,c); b_norm_sq[c]=mu[c]; mu_crit[c]=tol*tol*b_norm_sq[c];
      if(std::abs(b_norm_sq[c])<1.0e-14 || mu[c]<mu_crit[c]) colfrozen[c]=1; }
    std::vector<double> zeta(ncolp,1.0), zeta_old(ncolp,1.0), zeta_new(ncolp), alm(ncolp), betm(ncolp);
    std::vector<char> frozen(ncolp,0);

    int k=0;
    for(; k<maxiter; ++k){
      mult_block<CuC,N,NSTACK><<<g_seed,NThreadsPerBlock>>>(d_seedtmp, d_p0, D.M_DW.val, D.M_DW.cols, D.M_DW.rows);
      mult_block<CuC,N,NSTACK><<<g_seed,NThreadsPerBlock>>>(d_q, d_seedtmp, D.M_DWH.val, D.M_DWH.cols, D.M_DWH.rows);
      for(int c=0;c<NSTACK;c++){ sc[c]=coeff_seed; sc2[c]=sig0; }
      CUDA_CHECK(cudaMemcpy(d_sc, sc.data(),  NSTACK*DB, H2D));
      CUDA_CHECK(cudaMemcpy(d_sc2,sc2.data(), NSTACK*DB, H2D));
      axpby_col<N><<<g_seed,NThreadsPerBlock>>>(d_q, d_sc, d_q, d_sc2, d_p0, NSTACK);

      for(int c=0;c<NSTACK;c++){
        if(colfrozen[c]){ al[c]=0.0; continue; }
        const double gam_re = dotcol(d_p0, d_q, c);
        if(!(gam_re>0.0) || !std::isfinite(gam_re) || !std::isfinite(mu[c])){
          std::clog << "MULTISHIFT-BLOCK BREAKDOWN col "<<c<<" iter "<<k<<": gam="<<gam_re<<" mu="<<mu[c]<<std::endl;
          assert(false);
        }
        al[c]=mu[c]/gam_re;
      }
      for(int c=0;c<NSTACK;c++) for(int m=0;m<npole;m++){
        const int col=c*npole+m;
        if(colfrozen[c]){ alm[col]=0.0; zeta_new[col]=zeta[col]; continue; }
        if(!frozen[col] && zeta[col]*zeta[col]*mu[c] < mu_crit[c]) frozen[col]=1;
        if(frozen[col]){ alm[col]=0.0; zeta_new[col]=0.0; continue; }
        const double denom = al[c]*bet_old[c]*(zeta_old[col]-zeta[col]) + zeta_old[col]*al_old[c]*(1.0 + hat[m]*al[c]);
        zeta_new[col]=(zeta[col]*zeta_old[col]*al_old[c])/denom;
        alm[col]=al[c]*zeta_new[col]/zeta[col];
      }
      CUDA_CHECK(cudaMemcpy(d_alm, alm.data(), ncolp*DB, H2D));
      multishift_x_update_block<N,NSTACK><<<g_blk,NThreadsPerBlock>>>(d_X, d_alm, d_pm, npole);

      for(int c=0;c<NSTACK;c++) sc[c] = -al[c];
      CUDA_CHECK(cudaMemcpy(d_sc, sc.data(), NSTACK*DB, H2D));
      axpy_col<N><<<g_seed,NThreadsPerBlock>>>(d_r, d_sc, d_q, d_r, NSTACK);

      double worst=0.0; bool nan_seen=false;
      for(int c=0;c<NSTACK;c++){
        if(colfrozen[c]){ mu_new[c]=mu[c]; continue; }
        mu_new[c]=dotself(d_r,c); if(std::isnan(mu_new[c])) nan_seen=true;
        worst=std::max(worst, b_norm_sq[c]>0 ? mu_new[c]/b_norm_sq[c] : mu_new[c]);
      }
      if(worst<tol*tol || nan_seen) break;

      for(int c=0;c<NSTACK;c++) bet[c] = colfrozen[c] ? 0.0 : mu_new[c]/mu[c];
      for(int c=0;c<NSTACK;c++) for(int m=0;m<npole;m++){
        const int col=c*npole+m;
        if(colfrozen[c] || frozen[col]){ betm[col]=0.0; continue; }
        const double r2=zeta_new[col]/zeta[col]; betm[col]=bet[c]*r2*r2;
      }
      CUDA_CHECK(cudaMemcpy(d_zeta, zeta_new.data(), ncolp*DB, H2D));
      CUDA_CHECK(cudaMemcpy(d_betm, betm.data(),     ncolp*DB, H2D));
      multishift_p_update_block<N,NSTACK><<<g_blk,NThreadsPerBlock>>>(d_pm, d_zeta, d_r, d_betm, npole);

      CUDA_CHECK(cudaMemcpy(d_sc, bet.data(), NSTACK*DB, H2D));
      axpy_col<N><<<g_seed,NThreadsPerBlock>>>(d_p0, d_sc, d_p0, d_r, NSTACK);

      zeta_old=zeta; zeta=zeta_new;
      for(int c=0;c<NSTACK;c++){ if(colfrozen[c]) continue; al_old[c]=al[c]; bet_old[c]=bet[c]; mu[c]=mu_new[c];
        if(mu_new[c]<mu_crit[c]) colfrozen[c]=1; }
    }
  }

  // ===== single-shift block CG: (A + sigma) X_c = b_c for all RHS c (C6f -- Term B inner pole solve) =====
  // A = (1/lmax^2) M_DWH M_DW (same seed matvec as solve_multishift_block: apply M_DW then M_DWH);
  // sigma = -k^2/cp[m] (the single pole shift, passed explicitly). A+sigma is Hermitian PD => REAL CG
  // coefficients (al=mu/gam, bet=mu_new/mu), mirroring solve_sq's proven-bit-identical flow but with
  // the seed matvec+real shift instead of DDH. Per-column freeze-converged + shared worst-column stop.
  // nstack=1 reproduces the single MatPoly::solve (the per-pole Term B solve at conserved_current
  // _claude.h:354) BIT-FOR-BIT. Used by ConservedCurrent::apply_k_block_t to batch Term B's npole
  // independent inner solves over the Nt-wide t-block. d_X, d_B are N*NSTACK blocks.
  // Reference (algorithm): standard CG; block batching for occupancy, cf. solve_sq.
  void solve_shift_block(CuC* d_X, const CuC* d_B, const double sigma, const double tol=1.0e-13, const int maxiter=1e8) const {
    PoolGuard<N,NSTACK> gd(*mempool, this);
    const double coeff_seed = 1.0/(D.lambda_max*D.lambda_max);
    const size_t Ns=(size_t)N*NSTACK; const int g=ng(Ns);
    CuC *d_p=d_ocg_p, *d_qq=d_ocg_q, *d_rr=d_ocg_r;

    CUDA_CHECK(cudaMemset(d_X, 0, Ns*CD));
    CUDA_CHECK(cudaMemcpy(d_rr, d_B, Ns*CD, D2D));
    CUDA_CHECK(cudaMemcpy(d_p,  d_B, Ns*CD, D2D));

    auto dotself=[&](const CuC* x,int c){ CuC d; CUBLAS_CHECK(cublasZdotc(handle,N,x+(size_t)c*N,1,x+(size_t)c*N,1,&d)); return real(d); };
    auto dotcol =[&](const CuC* x,const CuC* y,int c){ CuC d; CUBLAS_CHECK(cublasZdotc(handle,N,x+(size_t)c*N,1,y+(size_t)c*N,1,&d)); return real(d); };

    std::vector<double> mu(NSTACK), b_norm_sq(NSTACK), mu_crit(NSTACK), mu_new(NSTACK);
    std::vector<double> al(NSTACK,0.0), bet(NSTACK,0.0), sc(NSTACK), sc2(NSTACK);
    std::vector<char> colfrozen(NSTACK,0);
    for(int c=0;c<NSTACK;c++){ mu[c]=dotself(d_rr,c); b_norm_sq[c]=mu[c]; mu_crit[c]=tol*tol*b_norm_sq[c];
      if(std::abs(b_norm_sq[c])<1.0e-14 || mu[c]<mu_crit[c]) colfrozen[c]=1; }

    int k=0;
    for(; k<maxiter; ++k){
      // q = (A + sigma) p = coeff_seed * M_DWH(M_DW p) + sigma p
      mult_block<CuC,N,NSTACK><<<g,NThreadsPerBlock>>>(d_seedtmp, d_p, D.M_DW.val, D.M_DW.cols, D.M_DW.rows);
      mult_block<CuC,N,NSTACK><<<g,NThreadsPerBlock>>>(d_qq, d_seedtmp, D.M_DWH.val, D.M_DWH.cols, D.M_DWH.rows);
      for(int c=0;c<NSTACK;c++){ sc[c]=coeff_seed; sc2[c]=sigma; }
      CUDA_CHECK(cudaMemcpy(d_sc, sc.data(),  NSTACK*DB, H2D));
      CUDA_CHECK(cudaMemcpy(d_sc2,sc2.data(), NSTACK*DB, H2D));
      axpby_col<N><<<g,NThreadsPerBlock>>>(d_qq, d_sc, d_qq, d_sc2, d_p, NSTACK);

      for(int c=0;c<NSTACK;c++){
        if(colfrozen[c]){ al[c]=0.0; continue; }
        const double gam = dotcol(d_p, d_qq, c);
        if(!(gam>0.0) || !std::isfinite(gam) || !std::isfinite(mu[c])){
          std::clog << "SHIFT-BLOCK BREAKDOWN col "<<c<<" iter "<<k<<": gam="<<gam<<" mu="<<mu[c]<<std::endl;
          assert(false);
        }
        al[c]=mu[c]/gam;
      }
      CUDA_CHECK(cudaMemcpy(d_sc, al.data(), NSTACK*DB, H2D));
      axpy_col<N><<<g,NThreadsPerBlock>>>(d_X, d_sc, d_p, d_X, NSTACK);    // X += al p
      for(int c=0;c<NSTACK;c++) sc[c] = -al[c];
      CUDA_CHECK(cudaMemcpy(d_sc, sc.data(), NSTACK*DB, H2D));
      axpy_col<N><<<g,NThreadsPerBlock>>>(d_rr, d_sc, d_qq, d_rr, NSTACK); // r -= al q

      double worst=0.0; bool nan_seen=false;
      for(int c=0;c<NSTACK;c++){
        if(colfrozen[c]){ mu_new[c]=mu[c]; continue; }
        mu_new[c]=dotself(d_rr,c); if(std::isnan(mu_new[c])) nan_seen=true;
        worst=std::max(worst, b_norm_sq[c]>0 ? mu_new[c]/b_norm_sq[c] : mu_new[c]);
      }
      if(worst<tol*tol || nan_seen) break;

      for(int c=0;c<NSTACK;c++) bet[c] = colfrozen[c] ? 0.0 : mu_new[c]/mu[c];
      CUDA_CHECK(cudaMemcpy(d_sc, bet.data(), NSTACK*DB, H2D));
      axpy_col<N><<<g,NThreadsPerBlock>>>(d_p, d_sc, d_p, d_rr, NSTACK);   // p = bet p + r

      for(int c=0;c<NSTACK;c++){ if(colfrozen[c]) continue; mu[c]=mu_new[c]; if(mu_new[c]<mu_crit[c]) colfrozen[c]=1; }
    }
  }

  // copy residues A[1..size-1] -> d_alm (free scratch after the inner solve), for block_reduce_poles
  void load_residues() const {
    std::vector<double> Ah(npole); for(int m=1;m<D.size;m++) Ah[m-1]=D.A[m];
    CUDA_CHECK(cudaMemcpy(d_alm, Ah.data(), (size_t)npole*DB, H2D));
  }

  // ===== overlap block applies (were OverlapWMass::*_deviceAsyncLaunch_ms_block) =====
  void mult(CuC* d_res, const CuC* d_xi) const {
    PoolGuard<N,NSTACK> gd(*mempool, this);
    const size_t Ns=(size_t)N*NSTACK; const int g=ng(Ns);
    solve_multishift_block(d_Zblk, d_xi, Comp::TOL_INNER);                       // Z_m
    load_residues();
    block_reduce_poles<N,NSTACK><<<g,NThreadsPerBlock>>>(d_q, d_xi, d_alm, d_Zblk, npole);  // Zs0 = xi + sum A[m] Z_m
    mult_block<CuC,N,NSTACK><<<g,NThreadsPerBlock>>>(d_seedtmp, d_q, D.M_DW.val, D.M_DW.cols, D.M_DW.rows);
    CUDA_CHECK(cudaMemset(d_res, 0, Ns*CD));
    axpy_uniform<N><<<g,NThreadsPerBlock>>>(d_res, cplx(1.0/D.lambda_max), d_seedtmp, d_res, NSTACK);
    axpy_uniform<N><<<g,NThreadsPerBlock>>>(d_res, cplx(D.C),    d_res, d_xi, NSTACK);       // D_ov = C*. + xi
    axpy_uniform<N><<<g,NThreadsPerBlock>>>(d_res, cplx(D.mass), d_xi,  d_res, NSTACK);      // + m xi
    CUDA_CHECK(cudaDeviceSynchronize());
  }
  void adj(CuC* d_res, const CuC* d_xi) const {
    PoolGuard<N,NSTACK> gd(*mempool, this);
    const size_t Ns=(size_t)N*NSTACK; const int g=ng(Ns);
    mult_block<CuC,N,NSTACK><<<g,NThreadsPerBlock>>>(d_seedtmp, d_xi, D.M_DWH.val, D.M_DWH.cols, D.M_DWH.rows);
    CUDA_CHECK(cudaMemset(d_res, 0, Ns*CD));
    axpy_uniform<N><<<g,NThreadsPerBlock>>>(d_res, cplx(1.0/D.lambda_max), d_seedtmp, d_res, NSTACK);  // Ys0 = (1/lmax) M_DWH xi (in d_res)
    solve_multishift_block(d_Yblk, d_res, Comp::TOL_INNER);                       // Y_m (RHS = Ys0 = d_res, const in solve)
    load_residues();
    block_reduce_poles<N,NSTACK><<<g,NThreadsPerBlock>>>(d_res, d_res, d_alm, d_Yblk, npole);  // d_res = Ys0 + sum A[m] Y_m
    axpy_uniform<N><<<g,NThreadsPerBlock>>>(d_res, cplx(D.C),                 d_res, d_xi, NSTACK);   // D^dag_ov
    axpy_uniform<N><<<g,NThreadsPerBlock>>>(d_res, cplx(std::conj(D.mass)),   d_xi,  d_res, NSTACK);  // + m^* xi
    CUDA_CHECK(cudaDeviceSynchronize());
  }
  void DDH(CuC* d_res, const CuC* d_xi) const {
    PoolGuard<N,NSTACK> gd(*mempool, this);
    const size_t Ns=(size_t)N*NSTACK; const int g=ng(Ns);
    this->mult(d_mp, d_xi); CUDA_CHECK(cudaDeviceSynchronize());
    this->adj (d_mm, d_xi); CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemcpy(d_res, d_mp, Ns*CD, D2D));
    axpy_uniform<N><<<g,NThreadsPerBlock>>>(d_res, cplx(std::conj(D.mass)),   d_mp, d_res, NSTACK);
    axpy_uniform<N><<<g,NThreadsPerBlock>>>(d_res, cplx(Complex(1.0)+D.mass), d_mm, d_res, NSTACK);
    const double scalar = 2.0*D.mass.real() + std::norm(D.mass);
    axpy_uniform<N><<<g,NThreadsPerBlock>>>(d_res, cplx(-scalar),             d_xi, d_res, NSTACK);
    CUDA_CHECK(cudaDeviceSynchronize());
  }

  // ===== outer per-column block CG over DDH = op_Dmsq (was MatPoly::solve_block_cg) =====
  // d_X, d_B are N*NSTACK blocks. nstack=1 reproduces the single MatPoly::solve bit-for-bit.
  void solve_sq(CuC* d_X, const CuC* d_B, const double tol=1.0e-13, const int maxiter=1e8) const {
    PoolGuard<N,NSTACK> gd(*mempool, this);
    const size_t Ns=(size_t)N*NSTACK; const int g=ng(Ns);
    CuC *d_p=d_ocg_p, *d_qq=d_ocg_q, *d_rr=d_ocg_r;

    CUDA_CHECK(cudaMemset(d_X, 0, Ns*CD));
    CUDA_CHECK(cudaMemcpy(d_rr, d_B, Ns*CD, D2D));
    CUDA_CHECK(cudaMemcpy(d_p,  d_B, Ns*CD, D2D));

    auto dotself=[&](const CuC* x,int c){ CuC d; CUBLAS_CHECK(cublasZdotc(handle,N,x+(size_t)c*N,1,x+(size_t)c*N,1,&d)); return real(d); };
    auto dotc   =[&](const CuC* x,const CuC* y,int c){ CuC d; CUBLAS_CHECK(cublasZdotc(handle,N,x+(size_t)c*N,1,y+(size_t)c*N,1,&d)); return d; };

    std::vector<double> mu(NSTACK), b_norm_sq(NSTACK), mu_crit(NSTACK), mu_new(NSTACK), bet(NSTACK);
    std::vector<char> colfrozen(NSTACK,0);
    std::vector<CuC> alc(NSTACK);
    for(int c=0;c<NSTACK;c++){ mu[c]=dotself(d_rr,c); b_norm_sq[c]=mu[c]; mu_crit[c]=tol*tol*b_norm_sq[c];
      if(std::abs(b_norm_sq[c])<1.0e-14 || mu[c]<mu_crit[c]) colfrozen[c]=1; }

    int k=0;
    for(; k<maxiter; ++k){
      this->DDH(d_qq, d_p);                                  // q = (op_Dmsq) p  (block)
      for(int c=0;c<NSTACK;c++) alc[c] = colfrozen[c] ? cplx(0.0) : (cplx(mu[c]) / dotc(d_p,d_qq,c));
      CUDA_CHECK(cudaMemcpy(d_alc, alc.data(), NSTACK*CD, H2D));
      axpy_col_c<N><<<g,NThreadsPerBlock>>>(d_X, d_alc, d_p, d_X, NSTACK);
      for(int c=0;c<NSTACK;c++) alc[c] = -alc[c];
      CUDA_CHECK(cudaMemcpy(d_alc, alc.data(), NSTACK*CD, H2D));
      axpy_col_c<N><<<g,NThreadsPerBlock>>>(d_rr, d_alc, d_qq, d_rr, NSTACK);

      double worst=0.0; bool nan_seen=false;
      for(int c=0;c<NSTACK;c++){
        if(colfrozen[c]){ mu_new[c]=mu[c]; continue; }
        mu_new[c]=dotself(d_rr,c); if(std::isnan(mu_new[c])) nan_seen=true;
        worst=std::max(worst, b_norm_sq[c]>0 ? mu_new[c]/b_norm_sq[c] : mu_new[c]);
      }
      if(worst<tol*tol || nan_seen) break;

      for(int c=0;c<NSTACK;c++) bet[c] = colfrozen[c] ? 0.0 : mu_new[c]/mu[c];
      CUDA_CHECK(cudaMemcpy(d_sc, bet.data(), NSTACK*DB, H2D));
      axpy_col<N><<<g,NThreadsPerBlock>>>(d_p, d_sc, d_p, d_rr, NSTACK);

      for(int c=0;c<NSTACK;c++){ if(colfrozen[c]) continue; mu[c]=mu_new[c]; if(mu_new[c]<mu_crit[c]) colfrozen[c]=1; }
    }
  }

  // ===== host-block convenience: x_host = op_Dmsq^{-1} b_host (N*NSTACK host blocks). Hides the
  // device I/O via the owned d_blk_in/d_blk_out staging. In-place safe (x_host may alias b_host). =====
  void solve_sq_from_cpu(Complex* x_host, const Complex* b_host, const double tol=1.0e-13, const int maxiter=1e8) const {
    PoolGuard<N,NSTACK> gd(*mempool, this);
    const size_t Ns=(size_t)N*NSTACK;
    CUDA_CHECK(cudaMemcpy(d_blk_in, reinterpret_cast<const CuC*>(b_host), Ns*CD, H2D));
    this->solve_sq(d_blk_out, d_blk_in, tol, maxiter);
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(x_host), d_blk_out, Ns*CD, D2H));
  }
};
