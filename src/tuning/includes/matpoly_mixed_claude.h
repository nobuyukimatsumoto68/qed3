#pragma once

// matpoly_mixed_claude.h -- Chunk 1 of the mixed-precision inner-solver plan
// (mixedprec_impl_plan_claude.md). Single-shift mixed-precision reliable-update CG for the shifted
// seed system  (scale * D_W^dag D_W + sigma) x = b , mirroring Grid's
// ConjugateGradientMixedPrec.h (MixedPrecisionConjugateGradient, C. Kelly): a RESTARTED DEFECT
// CORRECTION -- the OUTER loop recomputes the true fp64 residual r = b - A x, casts it to fp32, runs
// a WHOLE fp32 CG on the correction system A delta = r, casts delta back and accumulates into the
// fp64 solution x; a final fp64 patch-up CG hits the exact tol. The Krylov iteration lives entirely
// in fp32 (matvec + axpy + dots), so the bulk cost runs at the ~2x fp32 bandwidth measured in the
// Chunk-0 gate, while the converged solution is fp64-accurate to `tol` (reliable update => the force
// tolerance is NOT loosened; reversibility-clean).
//
// Algorithm source: M.A. Clark, R. Babich, K. Barros, R.C. Brower, C. Rebbi, arXiv:0911.3191;
// reliable update = G.L.G. Sleijpen, H.A. van der Vorst, Computing 56 (1996) 141.
//
// A = scale * D_W^dag D_W = scale * M_DWH ( M_DW v ),  scale = 1/lambda_max^2. The fp32 mirrors of
// M_DW / M_DWH are held in CSRf and re-cast from the fp64 CSR via refresh() after every gauge update
// (in the FORCE, every MD step) -- mirroring Grid's per-solve precisionChange(Umu_f, Umu_d).

#include "sparse_matrix_mixed_claude.h"

template<Idx N>
struct MixedSolver {
  const CSR<N>& M_DW;    // fp64 D_W  CSR (borrowed from the OverlapWMass operator)
  const CSR<N>& M_DWH;   // fp64 D_W^dag CSR
  double scale;          // 1/lambda_max^2
  Idx nnz;

  CSRf<N> dw_f, dwh_f;   // owned fp32 mirrors (val cast from the fp64 CSR; cols/rows shared)
  cublasHandle_t hb;     // dots/norms (fp32 Cdotc/Scnrm2 + fp64 Zdotc/Dznrm2 share one handle)

  // fp32 inner-Krylov scratch
  CuCf *d_xf=nullptr, *d_rf=nullptr, *d_pf=nullptr, *d_qf=nullptr, *d_tf=nullptr;
  // fp64 outer/patch-up scratch. d_r = outer residual (also the patch-up RHS); d_rcg = the patch-up
  // CG's OWN residual (kept distinct from d_r so the RHS is never mutated during the patch-up).
  CuC *d_tmp=nullptr, *d_Ax=nullptr, *d_r=nullptr, *d_rcg=nullptr, *d_dx=nullptr, *d_p=nullptr, *d_q=nullptr;
  // nested-mixed per-shift cleanup scratch (residual + correction increment)
  CuC *d_res=nullptr, *d_corr=nullptr;

  // diagnostics (last solve): last_stage = # fp32 stages / reliable updates, last_f = fp32-phase iters,
  // last_d = fp64-phase iters, last_clean_lo/hi = per-shift cleanup fp32/fp64 iters (multishift only).
  int last_stage=0, last_f=0, last_d=0, last_clean_lo=0, last_clean_hi=0;

  MixedSolver(const CSR<N>& dw, const CSR<N>& dwh, const double scale_, const Idx nnz_)
    : M_DW(dw)
    , M_DWH(dwh)
    , scale(scale_)
    , nnz(nnz_)
  {
    dw_f .associate( M_DW.cols,  M_DW.rows,  nnz );
    dwh_f.associate( M_DWH.cols, M_DWH.rows, nnz );
    refresh();

    CUBLAS_CHECK( cublasCreate(&hb) );

    CUDA_CHECK(cudaMalloc(&d_xf, N*CDf));
    CUDA_CHECK(cudaMalloc(&d_rf, N*CDf));
    CUDA_CHECK(cudaMalloc(&d_pf, N*CDf));
    CUDA_CHECK(cudaMalloc(&d_qf, N*CDf));
    CUDA_CHECK(cudaMalloc(&d_tf, N*CDf));

    CUDA_CHECK(cudaMalloc(&d_tmp, N*CD));
    CUDA_CHECK(cudaMalloc(&d_Ax,  N*CD));
    CUDA_CHECK(cudaMalloc(&d_r,   N*CD));
    CUDA_CHECK(cudaMalloc(&d_rcg, N*CD));
    CUDA_CHECK(cudaMalloc(&d_dx,  N*CD));
    CUDA_CHECK(cudaMalloc(&d_p,   N*CD));
    CUDA_CHECK(cudaMalloc(&d_q,   N*CD));
    CUDA_CHECK(cudaMalloc(&d_res, N*CD));
    CUDA_CHECK(cudaMalloc(&d_corr,N*CD));
  }

  ~MixedSolver(){
    CUBLAS_CHECK( cublasDestroy(hb) );
    CUDA_CHECK(cudaFree(d_xf)); CUDA_CHECK(cudaFree(d_rf)); CUDA_CHECK(cudaFree(d_pf));
    CUDA_CHECK(cudaFree(d_qf)); CUDA_CHECK(cudaFree(d_tf));
    CUDA_CHECK(cudaFree(d_tmp)); CUDA_CHECK(cudaFree(d_Ax)); CUDA_CHECK(cudaFree(d_r));
    CUDA_CHECK(cudaFree(d_rcg)); CUDA_CHECK(cudaFree(d_dx)); CUDA_CHECK(cudaFree(d_p)); CUDA_CHECK(cudaFree(d_q));
    CUDA_CHECK(cudaFree(d_res)); CUDA_CHECK(cudaFree(d_corr));
  }

  // re-cast the fp32 mirror values from the current fp64 CSR (call after every gauge update).
  void refresh(){
    dw_f .cast_from( M_DW.val  );
    dwh_f.cast_from( M_DWH.val );
  }

  // ---- fp64 apply: out = (scale * D_W^dag D_W + sigma) in -------------------------------------
  void applyA_fp64(CuC* out, const CuC* in, const double sigma) const {
    M_DW ( d_tmp, in  );                 // d_tmp = D_W in
    M_DWH( out,   d_tmp );               // out   = D_W^dag D_W in
    axpby_shift_d<N><<<NBlocks,NThreadsPerBlock>>>(out, out, in, scale, sigma);
  }

  // ---- fp32 apply: out = (scale * D_W^dag D_W + sigma) in -------------------------------------
  void applyA_fp32(CuCf* out, const CuCf* in, const float sigma) const {
    dw_f .apply( d_tf, in   );           // d_tf = D_W in   (fp32)
    dwh_f.apply( out,  d_tf );           // out  = D_W^dag D_W in
    axpby_shift_f<N><<<NBlocks,NThreadsPerBlock>>>(out, out, in, (float)scale, sigma);
  }

  inline double norm2_fp64(const CuC* v) const {
    double nrm=0.0;
    CUBLAS_CHECK( cublasDznrm2(hb, N, v, 1, &nrm) );
    return nrm*nrm;
  }
  inline double norm2_fp32(const CuCf* v) const {
    float nrm=0.0f;
    CUBLAS_CHECK( cublasScnrm2(hb, N, v, 1, &nrm) );
    return (double)nrm*(double)nrm;
  }

  // ---- fp32 inner CG: solve (A+sigma) x = rhs from x=0, to relative tol tol_in. RHS is in d_rf on
  //      entry (residual = RHS, p = RHS); correction returned in d_xf. Returns iteration count. ----
  int cg_fp32(const float sigma, const double tol_in, const int maxit){
    CUDA_CHECK(cudaMemset(d_xf, 0, N*CDf));
    CUDA_CHECK(cudaMemcpy(d_pf, d_rf, N*CDf, D2D));
    double mu = norm2_fp32(d_rf);
    const double crit = tol_in*tol_in*mu;
    if(mu<=0.0 || mu<crit) return 0;
    int k=0;
    for(; k<maxit; ++k){
      applyA_fp32(d_qf, d_pf, sigma);
      cuComplex gamc;
      CUBLAS_CHECK( cublasCdotc(hb, N, d_pf, 1, d_qf, 1, &gamc) );
      const double gam = (double)gamc.x;
      if(!(gam>0.0)) break;             // fp32 curvature floor -> let the fp64 outer/patch recover
      const double al = mu/gam;
      Taxpy_f<N><<<NBlocks,NThreadsPerBlock>>>(d_xf, cplxf((float) al), d_pf, d_xf);
      Taxpy_f<N><<<NBlocks,NThreadsPerBlock>>>(d_rf, cplxf((float)-al), d_qf, d_rf);
      const double mu_new = norm2_fp32(d_rf);
      if(mu_new<crit || std::isnan(mu_new)) { ++k; break; }
      const double bet = mu_new/mu;
      Taxpy_f<N><<<NBlocks,NThreadsPerBlock>>>(d_pf, cplxf((float)bet), d_pf, d_rf);
      mu = mu_new;
    }
    return k;
  }

  // ---- fp64 CG: solve (A+sigma) x = d_b from x=0, to an ABSOLUTE squared-residual target abs_crit.
  //      Patch-up: d_b is the outer residual r (= b - A x_mix), so ||r - A dx||^2 < abs_crit = stop
  //      makes the TOTAL residual ||b - A(x_mix+dx)|| < tol*||b|| -- reducing r by only ~1 order
  //      (few iters), NOT re-solving the whole system (a relative tol would over-solve by ~1/tol). ---
  int cg_fp64(CuC* d_x, const CuC* d_b, const double sigma, const double abs_crit, const int maxit){
    CUDA_CHECK(cudaMemset(d_x, 0, N*CD));
    CUDA_CHECK(cudaMemcpy(d_rcg, d_b, N*CD, D2D));   // own residual (d_b never mutated -- may alias d_r)
    CUDA_CHECK(cudaMemcpy(d_p, d_rcg, N*CD, D2D));
    double mu = norm2_fp64(d_rcg);
    const double crit = abs_crit;
    if(mu<=0.0 || mu<crit) return 0;
    int k=0;
    for(; k<maxit; ++k){
      applyA_fp64(d_q, d_p, sigma);
      cuDoubleComplex gamc;
      CUBLAS_CHECK( cublasZdotc(hb, N, d_p, 1, d_q, 1, &gamc) );
      const double gam = real(gamc);
      if(!(gam>0.0)) break;
      const double al = mu/gam;
      Taxpy<CuC,N><<<NBlocks,NThreadsPerBlock>>>(d_x,   cplx( al), d_p, d_x);
      Taxpy<CuC,N><<<NBlocks,NThreadsPerBlock>>>(d_rcg, cplx(-al), d_q, d_rcg);
      const double mu_new = norm2_fp64(d_rcg);
      if(mu_new<crit || std::isnan(mu_new)) { ++k; break; }
      const double bet = mu_new/mu;
      Taxpy<CuC,N><<<NBlocks,NThreadsPerBlock>>>(d_p, cplx(bet), d_p, d_rcg);
      mu = mu_new;
    }
    return k;
  }

  // ---- nested-MIXED per-shift cleanup: solve (A+sigma) x = rhs from 0 to the ABSOLUTE target abs_crit,
  //      with n_nested fp32 defect-correction stages (each to relative tol_f_clean) then a fp64 finish.
  //      Used to refine a stalled pole from the fp32 floor (~1e-7) to tol_d: the leftover is only a few
  //      orders, so ~1 fp32 stage suffices and runs at ~2x. iters_lo/hi += the fp32/fp64 iteration counts.
  void cleanup_pole_mixed(CuC* d_x, const CuC* d_rhs, const double sigma, const double abs_crit,
                          const double tol_f_clean, const int n_nested, const int maxit,
                          int& iters_lo, int& iters_hi){
    CUDA_CHECK(cudaMemset(d_x, 0, N*CD));
    CUDA_CHECK(cudaMemcpy(d_res, d_rhs, N*CD, D2D));   // residual = rhs (x=0)
    for(int s=0; s<n_nested; s++){
      const double rr = norm2_fp64(d_res);
      if(rr < abs_crit) return;
      // fp32 stage targets the ABSOLUTE crit_d, not a fixed relative tol: rel = sqrt(abs_crit/||res||^2)
      // (~1e-2 for a pole at ~1e-7). A fixed relative tol_f_clean over-solves by extra orders (bug: 8x
      // more iters). Floor at tol_f_clean so we never ask fp32 below what it can reach in one stage.
      const double rel_tol = std::max( std::sqrt(abs_crit/rr), tol_f_clean );
      cast_z2c_launch(d_rf, d_res, N);                                        // fp32 stage RHS = residual
      iters_lo += cg_fp32((float)sigma, rel_tol, maxit);                      // fp32 correction -> d_xf
      cast_c2z_launch(d_corr, d_xf, N);
      Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_x, 1.0, d_corr, d_x);       // x += corr
      applyA_fp64(d_Ax, d_x, sigma);                                          // (A+sigma) x
      Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_res, -1.0, d_Ax, d_rhs);    // res = rhs - (A+sigma)x
    }
    if(norm2_fp64(d_res) >= abs_crit){
      iters_hi += cg_fp64(d_corr, d_res, sigma, abs_crit, maxit);            // fp64 finish on the remainder
      Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_x, 1.0, d_corr, d_x);
    }
  }

  // ============== mixed defect-correction solve (A+sigma) x = b: TWO explicit tolerances ==========
  // tol_f = fp32 STAGE tolerance (each fp32 CG solves the fp64 residual system A dx = r to RELATIVE
  //         tol_f -- a loose, fp32-ACHIEVABLE target, e.g. 1e-4..1e-5); the fp64 residual is recomputed
  //         and re-cast between stages (restarted defect correction, Clark et al. arXiv:0911.3191).
  // tol_d = final fp64 tolerance (= physics target). A closing fp64 CG reduces the residual to the
  //         ABSOLUTE target stop_d = tol_d^2||b||^2 -> total residual < tol_d*||b||.
  // The two are INDEPENDENT knobs (NM): loose tol_f -> more cheap fp32 stages + tiny fp64 finish;
  // tight tol_f (below the fp32 floor ~1e-6..1e-7) wastes fp32 iters spinning at the floor. Tune tol_f
  // (Chunk 1b) to minimize total cost. Krylov is discarded between stages (restart penalty) -- this
  // single-shift path is the beachhead; the multishift path (Chunk 2) keeps one continuous pass.
  void solve(CuC* d_x, const CuC* d_b, const double sigma,
             const double tol_f=1.0e-5, const double tol_d=1.0e-9,
             const int maxit=100000, const int maxstage=30){
    last_stage=0; last_f=0; last_d=0;
    const double bnorm2 = norm2_fp64(d_b);
    const double stop_d = tol_d*tol_d*bnorm2;

    CUDA_CHECK(cudaMemset(d_x, 0, N*CD));
    if(bnorm2<=0.0){ return; }   // zero RHS

    // fp32 stages: solve the current fp64 residual system to RELATIVE tol_f, accumulate into x.
    for(int stage=0; stage<maxstage; ++stage){
      applyA_fp64(d_Ax, d_x, sigma);                                     // Ax = (A+sigma) x
      Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_r, -1.0, d_Ax, d_b);  // r = b - Ax (TRUE)
      const double nrm2 = norm2_fp64(d_r);
      if(nrm2 < stop_d) break;                                          // reached tol_d in fp32 alone
      cast_z2c_launch(d_rf, d_r, N);                                    // fp32 stage RHS = residual
      last_f += cg_fp32((float)sigma, tol_f, maxit);                    // fp32 correction -> d_xf
      last_stage = stage+1;
      cast_c2z_launch(d_dx, d_xf, N);
      Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_x, 1.0, d_dx, d_x);   // x += dx
    }

    // fp64 finish: reduce the fresh residual to the ABSOLUTE target stop_d -> total < tol_d*||b||.
    applyA_fp64(d_Ax, d_x, sigma);
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_r, -1.0, d_Ax, d_b);
    last_d = cg_fp64(d_dx, d_r, sigma, stop_d, maxit);
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_x, 1.0, d_dx, d_x);
    CUDA_CHECK(cudaDeviceSynchronize());
  }

  // ============ mixed MULTISHIFT: (A + sigma_m) X_m = b for all m, ONE Krylov pass ==============
  // Ports MatPoly::solve_multishift (Jegerlehner hep-lat/9612014). Search directions (d_pm), shifted
  // solutions (d_X), and ALL scalar recurrences stay fp64 -> NO restart, CG superlinearity preserved
  // (Grid ConjugateGradientMultiShiftMixedPrec, C. Kelly 2020). Two phases in ONE continuous recurrence:
  //   fp32 phase  : seed matvec in fp32 (~2x) + frequent (every `freq`) fp64 reliable updates of the
  //                 true seed residual, down to the crossover tol_f. The reliable update MUST be
  //                 frequent (SMALL corrections) so r and the directions stay consistent.
  //   fp64 tail   : below tol_f the fp32 matvec floors (~1e-7: q noise kills gam positivity -> CG
  //                 breakdown), so switch the matvec to fp64 -- a PLAIN precision switch, NO big
  //                 residual reset (a single big reset stalls the POLES at ~1e-7, verified) -- and
  //                 continue the SAME recurrence to tol_d. Seed drives; residual-proportionality
  //                 r_m = zeta_m r_seed (|zeta_m| <= 1) => seed at tol_d => all poles at tol_d.
  // tol_f = fp32->fp64 crossover (~1e-5, ABOVE the ~1e-7 fp32 floor); tol_d = physics target; freq =
  // reliable-update cadence. Diagnostics: last_f/last_d = fp32-/fp64-phase matvec iters, last_stage =
  // # reliable-update fp64 matvecs.
  //   d_X : N*npole COLUMN-MAJOR block [m*N+i] (output, zeroed here). d_b : single shared RHS.
  //   sigma : HOST array of npole shifts (all > 0 for SPD A+sigma).
  // Defaults = the TUNED L4 config (2026-07-15): tol_f=1e-6 crossover, n_nested=1 fp32 cleanup stage,
  // tol_f_clean=1e-4 floor, freq=10 -> net ~1.4-1.45x at L4 (see mixedprec_impl_plan_claude.md).
  void solve_multishift(CuC* d_X, const CuC* d_b, const double* sigma, const int npole,
                        const double tol_d=1.0e-9, const double tol_f=1.0e-6, const int freq=10,
                        const bool cleanup=true, const int n_nested=1, const double tol_f_clean=1.0e-4,
                        const int maxit=100000){
    last_stage=0; last_f=0; last_d=0; last_clean_lo=0; last_clean_hi=0;   // last_f = fp32-phase matvec iters, last_d = fp64-phase matvec iters

    int seed=0;
    for(int m=1;m<npole;m++) if(sigma[m]<sigma[seed]) seed=m;
    const double sig0 = sigma[seed];
    std::vector<double> hat(npole);
    for(int m=0;m<npole;m++) hat[m] = sigma[m]-sig0;

    // per-pole search directions (N*npole) + coeff arrays; seed work vectors reuse members d_p/d_q/d_r.
    CuC* d_pm;
    CUDA_CHECK(cudaMalloc(&d_pm, (size_t)N*npole*CD));
    double *d_alm,*d_zeta,*d_betm;
    CUDA_CHECK(cudaMalloc(&d_alm,  npole*sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_zeta, npole*sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_betm, npole*sizeof(double)));
    const Idx total = (Idx)N*npole;
    const int nb_blk = (total + NThreadsPerBlock)/NThreadsPerBlock;

    CUDA_CHECK(cudaMemset(d_X, 0, (size_t)N*npole*CD));
    CUDA_CHECK(cudaMemcpy(d_r, d_b, N*CD, D2D));
    CUDA_CHECK(cudaMemcpy(d_p, d_b, N*CD, D2D));   // d_p = seed direction p0
    for(int m=0;m<npole;m++) CUDA_CHECK(cudaMemcpy(d_pm+(size_t)m*N, d_b, N*CD, D2D));

    double mu = norm2_fp64(d_r);
    const double bnorm2 = mu;
    const double crit_d = tol_d*tol_d*bnorm2;
    const double crit_f = tol_f*tol_f*bnorm2;   // fp32->fp64 crossover (ABOVE the fp32 matvec floor)

    std::vector<double> zeta(npole,1.0), zeta_old(npole,1.0);
    std::vector<double> alm(npole), zeta_new(npole), betm(npole);
    double al_old=1.0, bet_old=0.0;
    std::vector<char> frozen(npole,0);   // freeze a pole once converged (stops zeta_m -> 0/0 NaN)

    bool hi=false;   // false = fp32 matvec phase, true = fp64 matvec tail (below the fp32 floor)
    if(bnorm2>0.0 && mu>=crit_d){
      int k=0;
      for(; k<maxit; ++k){
        // seed matvec q = (A + sig0) p0 -- fp32 until the crossover tol_f, fp64 for the sub-floor tail.
        // The crossover is a plain precision switch of the matvec: NO big residual reset (that broke
        // r-p consistency and stalled the poles); frequent small reliable updates in the fp32 phase
        // keep r honest, and the fp64 matvec below the floor descends without breakdown.
        if(!hi){
          cast_z2c_launch(d_pf, d_p, N);
          applyA_fp32(d_qf, d_pf, (float)sig0);
          cast_c2z_launch(d_q, d_qf, N);
          last_f++;
        }
        else{
          applyA_fp64(d_q, d_p, sig0);
          last_d++;
        }

        CuC gamc;
        CUBLAS_CHECK( cublasZdotc(hb, N, d_p, 1, d_q, 1, &gamc) );
        const double gam_re = real(gamc);
        if(!(gam_re>0.0) || !std::isfinite(gam_re) || !std::isfinite(mu)) break;   // breakdown
        const double al = mu/gam_re;

        for(int m=0;m<npole;m++){
          if(!frozen[m] && zeta[m]*zeta[m]*mu < crit_d) frozen[m]=1;
          if(frozen[m]){ alm[m]=0.0; zeta_new[m]=0.0; continue; }
          const double denom = al*bet_old*(zeta_old[m]-zeta[m]) + zeta_old[m]*al_old*(1.0 + hat[m]*al);
          zeta_new[m] = (zeta[m]*zeta_old[m]*al_old)/denom;
          alm[m]      = al * zeta_new[m]/zeta[m];
        }
        CUDA_CHECK(cudaMemcpy(d_alm, alm.data(), npole*sizeof(double), H2D));
        multishift_x_update<N><<<nb_blk, NThreadsPerBlock>>>(d_X, d_alm, d_pm, npole);

        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock>>>(d_r, cplx(-al), d_q, d_r);     // r -= al q (recurrence)
        double mu_new = norm2_fp64(d_r);

        // RELIABLE UPDATE (fp64 true seed residual r = b - (A+sig0)X_seed): in the FP32 phase every
        // `freq` iters (bound the fp32 drift; frequent SMALL corrections keep r-p consistent, Grid
        // CGMultiShiftMixedPrec), AND whenever the recurrence CLAIMS convergence -- so the stop test
        // sees a TRUE residual. In the fp64 phase the matvec is exact -> only the convergence check RUs.
        if( (!hi && (k+1)%freq==0) || mu_new<crit_d ){
          applyA_fp64(d_Ax, d_X + (size_t)seed*N, sig0);                              // (A+sig0) X_seed
          Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_r, -1.0, d_Ax, d_b); // r = b - (A+sig0)X_seed
          mu_new = norm2_fp64(d_r);
          last_stage++;   // count reliable-update fp64 matvecs
        }

        // crossover: switch the matvec to fp64 once the (honest) seed residual crosses tol_f -- ABOVE
        // the fp32 breakdown floor, so no big reset is needed and the recurrence continues seamlessly.
        if(!hi && mu_new < crit_f) hi = true;

        if(mu_new<crit_d || std::isnan(mu_new)) break;   // seed converged => all poles converged
        const double bet = mu_new/mu;

        for(int m=0;m<npole;m++){
          if(frozen[m]){ betm[m]=0.0; continue; }
          const double ratio = zeta_new[m]/zeta[m];
          betm[m] = bet*ratio*ratio;
        }
        CUDA_CHECK(cudaMemcpy(d_zeta, zeta_new.data(), npole*sizeof(double), H2D));
        CUDA_CHECK(cudaMemcpy(d_betm, betm.data(),     npole*sizeof(double), H2D));
        multishift_p_update<N><<<nb_blk, NThreadsPerBlock>>>(d_pm, d_zeta, d_r, d_betm, npole);
        Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock>>>(d_p, cplx(bet), d_p, d_r);     // p0 = r + bet p0

        zeta_old=zeta; zeta=zeta_new; al_old=al; bet_old=bet; mu=mu_new;
      }
    }

    // per-shift fp64 CLEANUP (PoS LATTICE2022 338 Sec 3.2): the fp32 matvec breaks the shifted-residual
    // co-linearity, so the poles stall at the fp32 floor (~1e-7) even when the seed reaches tol_d. Bring
    // each unconverged pole to tol_d with its OWN fp64 CG on (A+sigma_m) dX = r_m -- different RHS per
    // pole, so NO shared Krylov (this is the un-shareable cost that mixed multishift pays).
    if(cleanup){
      for(int m=0; m<npole; m++){
        applyA_fp64(d_Ax, d_X + (size_t)m*N, sigma[m]);                                // (A+sigma_m) X_m
        Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_r, -1.0, d_Ax, d_b);   // r_m = b - (A+sigma_m)X_m
        const double rn2 = norm2_fp64(d_r);
        if(rn2 < crit_d) continue;                                                     // pole already at tol_d
        // (A+sigma_m) dX = r_m -> abs crit_d: nested-mixed (n_nested fp32 stages + fp64 finish) if
        // n_nested>0, else pure fp64. Warm from the multishift's X_m via the defect-correction form.
        if(n_nested>0) cleanup_pole_mixed(d_dx, d_r, sigma[m], crit_d, tol_f_clean, n_nested, maxit,
                                          last_clean_lo, last_clean_hi);
        else           last_clean_hi += cg_fp64(d_dx, d_r, sigma[m], crit_d, maxit);
        Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_X + (size_t)m*N, 1.0, d_dx, d_X + (size_t)m*N);
      }
    }

    CUDA_CHECK(cudaFree(d_pm));
    CUDA_CHECK(cudaFree(d_alm)); CUDA_CHECK(cudaFree(d_zeta)); CUDA_CHECK(cudaFree(d_betm));
    CUDA_CHECK(cudaDeviceSynchronize());
  }
};
