#pragma once

#include <vector>
#include <utility>

// conserved_current_block_claude.h -- C6f-b1: t-BLOCKED sink K apply (the ~45%-of-jj lever).
// ConservedCurrentBlockT<OverlapOp,Gauge,NSTACK> wraps a ConservedCurrent<OverlapOp,Gauge> and computes
// K(t,fixed) xi for ALL t=0..NSTACK-1 (NSTACK = Nt) in ONE pass, sharing xi=phi'. For a fixed dual site
// n (temporal) or spatial link lk, the Nt link elements el(t)=(t,fixed) share xi but differ in W^{wz}.
//   - Step 1 (Z_m = R_m xi, xi-ONLY) is done ONCE via single-vector MatPoly::solve_multishift -- the
//     C6f-b1 INVARIANT: Step 1 stays multishift (npole shifts share one Krylov sequence), NOT
//     solve_shift_block, so the npole->1 collapse is preserved.
//   - Term A (C W(t) Z_0) and Term B (-C sum_m A[m] X R_m u_m(t)) are blocked over the Nt el's. Term B's
//     per-el inner pole solve R_m u_m(t), u_m(t)=X^dag W(t) Z_m - W^dag(t) X Z_m, batches over t at
//     FIXED pole m via BlockedMat::solve_shift_block(sigma_m) (single-shift block CG; the occupancy
//     win). Multishift (Step 1) and block (Term B) are ORTHOGONAL -- see mrhs plan .md C6f REFINEMENT.
// Numerically equals per-t ConservedCurrent::apply_k_ms (validate: test_apply_k_block_t_claude.cu).
//
// SEPARATE header (NOT folded into conserved_current_claude.h) so that header's many consumers stay
// unaffected: this references BlockedMat + the block kernels (mult_block/axpy_uniform), so include it
// AFTER blocked_mat_claude.h (-> multishift_block_kernels_claude.h) and conserved_current_claude.h.
// It only reads the (public) ConservedCurrent members (D, d_Zs, d_Zblock, d_tmp1/2/3, build_W).
//
// Reference (algorithm): conserved current K eq.(3.34); inner shifted solves = B. Jegerlehner,
// "Krylov space solvers for shifted linear systems", hep-lat/9612014.

// K-block working scratch: the three N*NSTACK blocks Term B needs, plus the t-independent X Z_m vector.
// (Per-t single-vector temporaries reuse the wrapped ConservedCurrent's d_tmp1/d_tmp2.)
template<Idx N, int NSTACK>
struct KBlockScratch {
  CuC *d_B;     // N*NSTACK : Term B RHS block, column t = u_m(t)
  CuC *d_V;     // N*NSTACK : solve result, column t = R_m u_m(t)
  CuC *d_MV;    // N*NSTACK : M_DW V block (X V = (1/lmax) M_DW V)
  CuC *d_XZ;    // N        : X Z_m (per pole m, t-independent)
  KBlockScratch(){
    CUDA_CHECK(cudaMalloc(&d_B,  (size_t)N*NSTACK*CD));
    CUDA_CHECK(cudaMalloc(&d_V,  (size_t)N*NSTACK*CD));
    CUDA_CHECK(cudaMalloc(&d_MV, (size_t)N*NSTACK*CD));
    CUDA_CHECK(cudaMalloc(&d_XZ, (size_t)N*CD));
  }
  ~KBlockScratch(){
    CUDA_CHECK(cudaFree(d_B)); CUDA_CHECK(cudaFree(d_V));
    CUDA_CHECK(cudaFree(d_MV)); CUDA_CHECK(cudaFree(d_XZ));
  }
  KBlockScratch(const KBlockScratch&)=delete;
  KBlockScratch& operator=(const KBlockScratch&)=delete;
};

template<typename OverlapOp, typename Gauge, int NSTACK>
struct ConservedCurrentBlockT {
  static constexpr Idx N = Comp::N;
  const ConservedCurrent<OverlapOp,Gauge>& kop;
  const OverlapOp& D;
  // FULL pool (with_pole_blocks=true): apply_k_block_t (K) only needs solve_shift_block, but
  // apply_k_dag_block_t (K^dag) Term A^dag is a per-t MULTISHIFT (w0(t)=W^dag(t) xi shared across poles,
  // varying over t) => solve_multishift_block, which needs the pole-blocks. Cost ~3*N*Nt*npole*16 B
  // (~189 MB at L1, Nt=128) -- negligible on a 12 GB card, so one full pool serves both K and K^dag.
  BlockMemPool<N,NSTACK> pool;
  BlockedMat<N,NSTACK,OverlapOp> blk;          // shares the pool
  KBlockScratch<N,NSTACK> scr;

  explicit ConservedCurrentBlockT(const ConservedCurrent<OverlapOp,Gauge>& kop_)
    : kop(kop_), D(kop_.D), pool(kop_.D.size-1, /*with_pole_blocks=*/true), blk(kop_.D, pool)
  {}

  static int ng(size_t tot){ return (int)((tot + NThreadsPerBlock - 1)/NThreadsPerBlock); }

  // build the el=(t,fixed) link element with the right pair type (matches apply_k / apply_k_dag).
  static std::pair<int,Idx>      makeEl(int t, const Idx& n)       { return {t, n}; }
  static std::pair<int,BaseLink> makeEl(int t, const BaseLink& lk) { return {t, lk}; }

  // adjoint single-link COO W^dag (replicates ConservedCurrent::apply_k's inline WH build; see :158-166)
  template<typename LinkEl>
  void build_WH(COO<N>& coo, const Gauge& U, const LinkEl& el) const {
    D.DW.d_coo_format(coo.en, U, el);
    for(auto& e : coo.en){
      const Complex z(cuCreal(e.v), cuCimag(e.v));
      e.v = cplx(Complex(-z.imag()/D.lambda_max, -z.real()/D.lambda_max));
      std::swap(e.i, e.j);
    }
    coo.do_it();
  }

  // K(t,fixed) xi for all t -> d_result_block (N*NSTACK; column t = K(t,fixed) xi).
  // FixedLink = Idx (temporal site) or BaseLink (spatial link); el(t) = {t, fixed}.
  template<typename FixedLink>
  void apply_k_block_t(CuC* d_result_block, const CuC* d_xi, const Gauge& U, const FixedLink& fixed) const {
    const int npole = D.size-1;
    const size_t Ns = (size_t)N*NSTACK;
    const int g = (int)((Ns + NThreadsPerBlock - 1)/NThreadsPerBlock);

    // --- Step 1 (ONCE, xi-only): Z_m = R_m xi via single-vector multishift; Z_0 = xi + sum A[m] Z_m ---
    MatPoly Aseed;
    Aseed.push_back(cplx(1.0/(D.lambda_max*D.lambda_max)), {&D.M_DW, &D.M_DWH});
    std::vector<double> sigma(npole);
    for(int m=1;m<D.size;m++) sigma[m-1] = -D.k*D.k/D.cp[m];
    Aseed.solve_multishift<N>(kop.d_Zblock, d_xi, sigma.data(), npole, Comp::TOL_INNER);
    for(int m=1;m<D.size;m++)
      CUDA_CHECK(cudaMemcpy(kop.d_Zs[m], kop.d_Zblock + (size_t)(m-1)*N, N*CD, D2D));
    CUDA_CHECK(cudaMemcpy(kop.d_Zs[0], d_xi, N*CD, D2D));
    for(int m=1;m<D.size;m++)
      Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(kop.d_Zs[0], D.A[m], kop.d_Zs[m], kop.d_Zs[0]);
    CUDA_CHECK(cudaDeviceSynchronize());

    // --- precompute the NSTACK single-link COOs W(t), W^dag(t) once (reused in Term A + Term B) ---
    // count-ctor value-initializes each COO IN PLACE (is_set=false); the vector is never resized, so no
    // COO copy/move is invoked (COO's freeing dtor + default shallow copy make by-value storage unsafe
    // under reallocation -- avoided here). Each COO's dtor frees its device CSR exactly once.
    std::vector<COO<N>> coo_W(NSTACK), coo_WH(NSTACK);
    for(int t=0;t<NSTACK;t++){
      const auto el = makeEl(t, fixed);
      kop.build_W(coo_W[t], U, el);
      build_WH(coo_WH[t], U, el);
    }

    // --- Term A: result[:,t] = C W(t) Z_0 ---
    CUDA_CHECK(cudaMemset(d_result_block, 0, Ns*CD));
    for(int t=0;t<NSTACK;t++){
      coo_W[t](kop.d_tmp1, kop.d_Zs[0]);
      CUDA_CHECK(cudaDeviceSynchronize());
      Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_result_block + (size_t)t*N, D.C, kop.d_tmp1, d_result_block + (size_t)t*N);
    }
    CUDA_CHECK(cudaDeviceSynchronize());

    // --- Term B: result[:,t] -= C A[m] X R_m u_m(t),  u_m(t) = X^dag W(t) Z_m - W^dag(t) X Z_m ---
    for(int m=1;m<D.size;m++){
      // XZ_m = X Z_m = (1/lmax) M_DW Z_m   (t-independent)
      { MatPoly X; X.push_back(cplx(1.0/D.lambda_max), {&D.M_DW}); X.on_gpu<N>(scr.d_XZ, kop.d_Zs[m]); }
      // assemble the Nt-wide RHS block u_m(t)
      for(int t=0;t<NSTACK;t++){
        coo_W[t](kop.d_tmp1, kop.d_Zs[m]);                    // W(t) Z_m
        CUDA_CHECK(cudaDeviceSynchronize());
        { MatPoly XH; XH.push_back(cplx(1.0/D.lambda_max), {&D.M_DWH}); XH.on_gpu<N>(kop.d_tmp2, kop.d_tmp1); } // X^dag W(t) Z_m
        coo_WH[t](kop.d_tmp1, scr.d_XZ);                      // W^dag(t) X Z_m
        CUDA_CHECK(cudaDeviceSynchronize());
        Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(scr.d_B + (size_t)t*N, -1.0, kop.d_tmp1, kop.d_tmp2); // u_m(t)
      }
      CUDA_CHECK(cudaDeviceSynchronize());
      // V_m(t) = R_m u_m(t): single-shift block CG over A+sigma_m, batched over t
      blk.solve_shift_block(scr.d_V, scr.d_B, sigma[m-1], Comp::TOL_INNER);
      // result -= C A[m] X V_m = result + (-C A[m]/lmax) M_DW V_m
      mult_block<CuC,N,NSTACK><<<g,NThreadsPerBlock>>>(scr.d_MV, scr.d_V, D.M_DW.val, D.M_DW.cols, D.M_DW.rows);
      axpy_uniform<N><<<g,NThreadsPerBlock>>>(d_result_block, cplx(-D.C*D.A[m]/D.lambda_max), scr.d_MV, d_result_block, NSTACK);
      CUDA_CHECK(cudaDeviceSynchronize());
    }
  }

  // K^dag(t,fixed) xi for all t -> d_result_block. Operator-ADJOINT mirror of apply_k_block_t; equals
  // per-t ConservedCurrent::apply_k_dag_ms. With X=D_W/lmax, R_m=(X X^dag - k^2/cp[m])^{-1}, C=D.C,
  // S=I+sum_m A_m R_m:  K^dag = C S W^dag - C sum_m A_m R_m (W^dag X - X^dag W) R_m X^dag.
  //   - Step 1 (ONCE, xi-only): Y_m = R_m X^dag xi (single-vector multishift; the b1 invariant).
  //   - Term A^dag: result(t) = C (w0(t) + sum_m A_m R_m w0(t)), w0(t)=W^dag(t) xi (t-DEPENDENT) =>
  //     batches over t as solve_multishift_block (shared RHS across poles per t) + block_reduce_poles.
  //   - Term B^dag: result(t) -= C A_m R_m (W^dag(t) X - X^dag W(t)) Y_m => Term-B-style block over t via
  //     solve_shift_block(sigma_m).
  template<typename FixedLink>
  void apply_k_dag_block_t(CuC* d_result_block, const CuC* d_xi, const Gauge& U, const FixedLink& fixed) const {
    const size_t Ns = (size_t)N*NSTACK;
    const int g = ng(Ns);
    std::vector<double> sigma(D.size-1);
    for(int m=1;m<D.size;m++) sigma[m-1] = -D.k*D.k/D.cp[m];

    // --- Step 1 (ONCE, xi-only): Y_m = R_m (X^dag xi) via single-vector multishift; Y_m -> kop.d_Zs[m] ---
    { MatPoly XH; XH.push_back(cplx(1.0/D.lambda_max), {&D.M_DWH}); XH.on_gpu<N>(kop.d_tmp1, d_xi); } // X^dag xi
    MatPoly Aseed;
    Aseed.push_back(cplx(1.0/(D.lambda_max*D.lambda_max)), {&D.M_DW, &D.M_DWH});
    Aseed.solve_multishift<N>(kop.d_Zblock, kop.d_tmp1, sigma.data(), D.size-1, Comp::TOL_INNER);
    for(int m=1;m<D.size;m++)
      CUDA_CHECK(cudaMemcpy(kop.d_Zs[m], kop.d_Zblock + (size_t)(m-1)*N, N*CD, D2D)); // Y_m
    CUDA_CHECK(cudaDeviceSynchronize());

    // --- precompute the NSTACK single-link COOs W(t), W^dag(t) once (count-ctor; no COO copy/move) ---
    std::vector<COO<N>> coo_W(NSTACK), coo_WH(NSTACK);
    for(int t=0;t<NSTACK;t++){
      const auto el = makeEl(t, fixed);
      kop.build_W(coo_W[t], U, el);
      build_WH(coo_WH[t], U, el);
    }

    // --- Term A^dag: result(t) = C (w0(t) + sum_m A_m R_m w0(t)), w0(t) = W^dag(t) xi ---
    // assemble the w0 block (scr.d_B); R_m w0(t) for all m via solve_multishift_block; reduce with the
    // residues A[m] (block_reduce_poles, handling the (col,pole) layout) -> Sw0(t); result = C Sw0(t).
    for(int t=0;t<NSTACK;t++){
      coo_WH[t](kop.d_tmp1, d_xi);                                   // w0(t) = W^dag(t) xi
      CUDA_CHECK(cudaDeviceSynchronize());
      CUDA_CHECK(cudaMemcpy(scr.d_B + (size_t)t*N, kop.d_tmp1, N*CD, D2D));
    }
    blk.solve_multishift_block(blk.d_Zblk, scr.d_B, Comp::TOL_INNER); // d_Zblk[(c,m)] = R_m w0(c)
    blk.load_residues();                                             // d_alm = A[1..size-1]
    block_reduce_poles<N,NSTACK><<<g,NThreadsPerBlock>>>(scr.d_V, scr.d_B, blk.d_alm, blk.d_Zblk, D.size-1); // Sw0 = w0 + sum A_m R_m w0
    CUDA_CHECK(cudaMemset(d_result_block, 0, Ns*CD));
    axpy_uniform<N><<<g,NThreadsPerBlock>>>(d_result_block, cplx(D.C), scr.d_V, d_result_block, NSTACK);     // result = C Sw0
    CUDA_CHECK(cudaDeviceSynchronize());

    // --- Term B^dag: result(t) -= C A_m R_m (W^dag(t) X - X^dag W(t)) Y_m ---
    for(int m=1;m<D.size;m++){
      // X Y_m = (1/lmax) M_DW Y_m  (t-independent)
      { MatPoly X; X.push_back(cplx(1.0/D.lambda_max), {&D.M_DW}); X.on_gpu<N>(scr.d_XZ, kop.d_Zs[m]); }
      // r_m(t) = W^dag(t) X Y_m - X^dag W(t) Y_m
      for(int t=0;t<NSTACK;t++){
        coo_WH[t](kop.d_tmp1, scr.d_XZ);                             // a = W^dag(t) X Y_m
        CUDA_CHECK(cudaDeviceSynchronize());
        coo_W[t](kop.d_tmp2, kop.d_Zs[m]);                           // W(t) Y_m
        CUDA_CHECK(cudaDeviceSynchronize());
        { MatPoly XH; XH.push_back(cplx(1.0/D.lambda_max), {&D.M_DWH}); XH.on_gpu<N>(kop.d_tmp3, kop.d_tmp2); } // b = X^dag W(t) Y_m
        Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(scr.d_B + (size_t)t*N, -1.0, kop.d_tmp3, kop.d_tmp1); // r_m(t) = a - b
      }
      CUDA_CHECK(cudaDeviceSynchronize());
      blk.solve_shift_block(scr.d_V, scr.d_B, sigma[m-1], Comp::TOL_INNER);                                  // s_m(t) = R_m r_m(t)
      axpy_uniform<N><<<g,NThreadsPerBlock>>>(d_result_block, cplx(-D.C*D.A[m]), scr.d_V, d_result_block, NSTACK); // -= C A_m s_m(t)
      CUDA_CHECK(cudaDeviceSynchronize());
    }
  }
};
