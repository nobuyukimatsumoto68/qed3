#pragma once

// overlap_force_Nfblock_claude.h -- BlockedForce<N,NSTACK,Op>: the Nf-flavor BLOCK form of the HMC
// fermion force. It reimplements the serial GRAD_L4 path (OverlapWMass::precalc_grad_deviceAsyncLaunch_ms
// :719 + grad_deviceAsyncLaunch_l4 :940, incl. the diagonal-mass Sec 5' fix) over an extra NSTACK
// (=Nf/2) axis, so the total force sum_f grad_f(link) is computed with ALL flavors packed. It does NOT
// call the serial routines -- those stay as the bit-for-tolerance reference and the non-block path.
//
// Reuse: the two multishift solves in precalc (Z_m = R_m eta, Y_m = R_m X^dag eta, seed
// A=(1/lmax^2) M_DWH M_DW, shifts -k^2/cp[m]) are EXACTLY BlockedMat::solve_multishift_block -> call it
// with the NSTACK eta-block; the [flavor*npole+pole] output layout is what the per-link grad kernels
// want (NO gather). Per-link work uses the runtime-ncol block kernels (mult_coo_block / block_dot /
// link_matvec_block) at ncol = NSTACK*npole, + block_reduce_poles for the m=0 combine. The action pool's
// multishift scratch is BORROWED (force runs after update_eta, never concurrent) so this owns only the
// force OUTPUT blocks.
//
// Reference (algorithm): B. Jegerlehner, hep-lat/9612014 (multishift). Include AFTER blocked_mat_claude.h.

#include <vector>

// scale a block IN PLACE by a uniform complex scalar: d_x_c *= a, over ncol columns of length N.
template<Idx N> __global__
void scale_uniform_blk( CuC* d_x, const CuC a, const int ncol ){
  const Idx gid = (Idx)blockIdx.x*blockDim.x + threadIdx.x;
  const Idx total = (Idx)N*ncol;
  if(gid < total) d_x[gid] = a*d_x[gid];
}

template<Idx N, int NSTACK, class Op>
struct BlockedForce {
  const Op& D;
  BlockedMat<N,NSTACK,Op>& blk;             // for solve_multishift_block (borrows its pool scratch)
  const int npole;                          // = D.size - 1

  // ---- force OUTPUT blocks (N*NSTACK*npole), layout [ (f*npole+m)*N + i ] ----
  CuC *d_Zg, *d_Yg, *d_XZg, *d_XYg, *d_CY, *d_CZ, *d_Wg;
  // ---- N*NSTACK work vectors ----
  CuC *d_Ys0, *d_eta_bra, *d_Mmeta, *d_comb, *d_combL, *d_tmpN;
  // ---- per-link single-link COO entries (raw; GRAD_L4 skips do_it) ----
  static constexpr int MAXENT = 256;
  Idx *d_ent_i, *d_ent_j;  CuC *d_ent_v;
  // ---- dots + residues ----
  CuC *d_dotA, *d_dotB, *d_innerF;          // NSTACK*npole, NSTACK*npole, NSTACK
  double *d_Ares;                           // npole  (Zolotarev residues A[1..size-1])

  static int ng(size_t tot){ return (int)(( tot + NThreadsPerBlock - 1 )/NThreadsPerBlock); }

  BlockedForce()=delete;

  explicit BlockedForce( const Op& D_, BlockedMat<N,NSTACK,Op>& blk_ )
    : D(D_), blk(blk_), npole(D_.size-1)
  {
    const size_t Ns  = (size_t)N*NSTACK;
    const size_t Nsp = (size_t)N*NSTACK*npole;
    CUDA_CHECK(cudaMalloc(&d_Zg,  Nsp*CD));  CUDA_CHECK(cudaMalloc(&d_Yg,  Nsp*CD));
    CUDA_CHECK(cudaMalloc(&d_XZg, Nsp*CD));  CUDA_CHECK(cudaMalloc(&d_XYg, Nsp*CD));
    CUDA_CHECK(cudaMalloc(&d_CY,  Nsp*CD));  CUDA_CHECK(cudaMalloc(&d_CZ,  Nsp*CD));
    CUDA_CHECK(cudaMalloc(&d_Wg,  Nsp*CD));
    CUDA_CHECK(cudaMalloc(&d_Ys0,     Ns*CD));  CUDA_CHECK(cudaMalloc(&d_eta_bra, Ns*CD));
    CUDA_CHECK(cudaMalloc(&d_Mmeta,   Ns*CD));  CUDA_CHECK(cudaMalloc(&d_comb,    Ns*CD));
    CUDA_CHECK(cudaMalloc(&d_combL,   Ns*CD));  CUDA_CHECK(cudaMalloc(&d_tmpN,    Ns*CD));
    CUDA_CHECK(cudaMalloc(&d_ent_i, MAXENT*sizeof(Idx)));
    CUDA_CHECK(cudaMalloc(&d_ent_j, MAXENT*sizeof(Idx)));
    CUDA_CHECK(cudaMalloc(&d_ent_v, MAXENT*CD));
    CUDA_CHECK(cudaMalloc(&d_dotA,   (size_t)NSTACK*npole*CD));
    CUDA_CHECK(cudaMalloc(&d_dotB,   (size_t)NSTACK*npole*CD));
    CUDA_CHECK(cudaMalloc(&d_innerF, (size_t)NSTACK*CD));
    CUDA_CHECK(cudaMalloc(&d_Ares,   (size_t)npole*DB));
    CUDA_CHECK(cudaDeviceSynchronize());
  }

  ~BlockedForce(){
    CUDA_CHECK(cudaFree(d_Zg));  CUDA_CHECK(cudaFree(d_Yg));  CUDA_CHECK(cudaFree(d_XZg));
    CUDA_CHECK(cudaFree(d_XYg)); CUDA_CHECK(cudaFree(d_CY));  CUDA_CHECK(cudaFree(d_CZ));
    CUDA_CHECK(cudaFree(d_Wg));
    CUDA_CHECK(cudaFree(d_Ys0)); CUDA_CHECK(cudaFree(d_eta_bra)); CUDA_CHECK(cudaFree(d_Mmeta));
    CUDA_CHECK(cudaFree(d_comb)); CUDA_CHECK(cudaFree(d_combL)); CUDA_CHECK(cudaFree(d_tmpN));
    CUDA_CHECK(cudaFree(d_ent_i)); CUDA_CHECK(cudaFree(d_ent_j)); CUDA_CHECK(cudaFree(d_ent_v));
    CUDA_CHECK(cudaFree(d_dotA)); CUDA_CHECK(cudaFree(d_dotB)); CUDA_CHECK(cudaFree(d_innerF));
    CUDA_CHECK(cudaFree(d_Ares));
  }

  BlockedForce(const BlockedForce&)=delete;
  BlockedForce& operator=(const BlockedForce&)=delete;

  size_t bytes() const {
    const size_t Nsp = (size_t)N*NSTACK*npole;
    const size_t Ns  = (size_t)N*NSTACK;
    return ( (size_t)7*Nsp + (size_t)6*Ns )*CD;
  }

  void load_residues() const {
    std::vector<double> Ah(npole);
    for(int m=1; m<D.size; m++) Ah[m-1]=D.A[m];
    CUDA_CHECK(cudaMemcpy(d_Ares, Ah.data(), (size_t)npole*DB, H2D));
  }

  // ---- precalc (link-independent): fill d_Zg, d_Yg(hat), d_XZg, d_XYg, d_eta_bra from d_eta_blk ----
  template<typename Gauge>
  void precalc( const Gauge& U, const CuC* d_eta_blk ) {
    const size_t Ns  = (size_t)N*NSTACK;
    const size_t Nsp = (size_t)N*NSTACK*npole;
    const int gN  = ng(Ns);
    const int gNp = ng(Nsp);
    const CuC inv_lmax = cplx(1.0/D.lambda_max);

    // Ys0 = (1/lmax) M_DWH eta
    mult_coo_block<CuC,N><<<gN,NThreadsPerBlock>>>(d_Ys0, d_eta_blk, D.M_DWH.val, D.M_DWH.cols, D.M_DWH.rows, NSTACK);
    scale_uniform_blk<N><<<gN,NThreadsPerBlock>>>(d_Ys0, inv_lmax, NSTACK);

    // Z_m = R_m eta ; Y_m = R_m Ys0   (block multishift over the eta-block -- borrows blk's pool)
    blk.solve_multishift_block(d_Zg, d_eta_blk, Comp::TOL_INNER);
    blk.solve_multishift_block(d_Yg, d_Ys0,     Comp::TOL_INNER);

    // massive (diagonal m_L): bra = (1+m_L)eta ; hatY_m = Y_m + mass_coeff * R_m X^dag(M_mass eta)
    if( D.mass != Complex(0.0,0.0) ){
      mult_coo_block<CuC,N><<<gN,NThreadsPerBlock>>>(d_Mmeta, d_eta_blk, D.M_mass.d_val, D.M_mass.d_cols, D.M_mass.d_rows, NSTACK); // M_mass eta
      axpy_uniform<N><<<gN,NThreadsPerBlock>>>(d_eta_bra, cplx(D.mass_coeff), d_Mmeta, d_eta_blk, NSTACK);                          // (1+m_L)eta
      mult_coo_block<CuC,N><<<gN,NThreadsPerBlock>>>(d_tmpN, d_Mmeta, D.M_DWH.val, D.M_DWH.cols, D.M_DWH.rows, NSTACK);            // M_DWH(M_mass eta)
      scale_uniform_blk<N><<<gN,NThreadsPerBlock>>>(d_tmpN, inv_lmax, NSTACK);                                                     // X^dag(M_mass eta)
      blk.solve_multishift_block(d_Wg, d_tmpN, Comp::TOL_INNER);                                                                   // W_m
      axpy_uniform<N><<<gNp,NThreadsPerBlock>>>(d_Yg, cplx(D.mass_coeff), d_Wg, d_Yg, NSTACK*npole);                              // hatY = Y + mass_coeff W
    }
    else {
      CUDA_CHECK(cudaMemcpy(d_eta_bra, d_eta_blk, Ns*CD, D2D));
    }

    // XZg = (1/lmax) M_DW Zg ; XYg = (1/lmax) M_DW hatYg   (block over NSTACK*npole columns)
    mult_coo_block<CuC,N><<<gNp,NThreadsPerBlock>>>(d_XZg, d_Zg, D.M_DW.val, D.M_DW.cols, D.M_DW.rows, NSTACK*npole);
    scale_uniform_blk<N><<<gNp,NThreadsPerBlock>>>(d_XZg, inv_lmax, NSTACK*npole);
    mult_coo_block<CuC,N><<<gNp,NThreadsPerBlock>>>(d_XYg, d_Yg, D.M_DW.val, D.M_DW.cols, D.M_DW.rows, NSTACK*npole);
    scale_uniform_blk<N><<<gNp,NThreadsPerBlock>>>(d_XYg, inv_lmax, NSTACK*npole);

    load_residues();
    CUDA_CHECK(cudaDeviceSynchronize());
  }

  // ---- per-link total force = sum_f grad_f(link)  (GRAD_L4 block form) ----
  template<typename Link, typename Gauge>
  double grad_block( const Link& link, const Gauge& U, const CuC* d_eta_blk ) const {
    const size_t Ns  = (size_t)N*NSTACK;
    const size_t Nsp = (size_t)N*NSTACK*npole;

    COO<N> coo;
    D.DW.d_coo_format(coo.en, U, link);          // host entries only; NO do_it
    const int nent = (int)coo.en.size();
    assert(nent <= MAXENT);
    std::vector<Idx> hi(nent), hj(nent);  std::vector<CuC> hv(nent);
    for(int k=0;k<nent;k++){ hi[k]=coo.en[k].i; hj[k]=coo.en[k].j; hv[k]=coo.en[k].v; }
    CUDA_CHECK(cudaMemcpy(d_ent_i, hi.data(), (size_t)nent*sizeof(Idx), H2D));
    CUDA_CHECK(cudaMemcpy(d_ent_j, hj.data(), (size_t)nent*sizeof(Idx), H2D));
    CUDA_CHECK(cudaMemcpy(d_ent_v, hv.data(), (size_t)nent*CD,          H2D));

    const int gL  = ng((size_t)nent*NSTACK*npole);
    const int gL0 = ng((size_t)nent*NSTACK);

    // CY = link . hatYg ; CZ = link . Zg   (over NSTACK*npole cols; out zeroed first)
    CUDA_CHECK(cudaMemset(d_CY, 0, Nsp*CD));
    CUDA_CHECK(cudaMemset(d_CZ, 0, Nsp*CD));
    link_matvec_block<N><<<gL,NThreadsPerBlock>>>(d_CY, d_Yg, d_ent_i, d_ent_j, d_ent_v, nent, NSTACK*npole);
    link_matvec_block<N><<<gL,NThreadsPerBlock>>>(d_CZ, d_Zg, d_ent_i, d_ent_j, d_ent_v, nent, NSTACK*npole);
    // a_{f,m} = <CY, XZg> ; b_{f,m} = <XYg, CZ>
    block_dot<N><<<NSTACK*npole,NThreadsPerBlock>>>(d_dotA, d_CY, d_XZg, NSTACK*npole);
    block_dot<N><<<NSTACK*npole,NThreadsPerBlock>>>(d_dotB, d_XYg, d_CZ, NSTACK*npole);

    // m=0 term: comb_f = eta_f + sum_m A[m] Z_{f,m} ; innerF_f = <eta_bra_f, link . comb_f>
    block_reduce_poles<N,NSTACK><<<ng(Ns),NThreadsPerBlock>>>(d_comb, d_eta_blk, d_Ares, d_Zg, npole);
    CUDA_CHECK(cudaMemset(d_combL, 0, Ns*CD));
    link_matvec_block<N><<<gL0,NThreadsPerBlock>>>(d_combL, d_comb, d_ent_i, d_ent_j, d_ent_v, nent, NSTACK);
    block_dot<N><<<NSTACK,NThreadsPerBlock>>>(d_innerF, d_eta_bra, d_combL, NSTACK);
    CUDA_CHECK(cudaDeviceSynchronize());

    std::vector<CuC> hA(NSTACK*npole), hB(NSTACK*npole), hI(NSTACK);
    CUDA_CHECK(cudaMemcpy(hA.data(), d_dotA,   (size_t)NSTACK*npole*CD, D2H));
    CUDA_CHECK(cudaMemcpy(hB.data(), d_dotB,   (size_t)NSTACK*npole*CD, D2H));
    CUDA_CHECK(cudaMemcpy(hI.data(), d_innerF, (size_t)NSTACK*CD,       D2H));

    double res = 0.0;
    for(int f=0; f<NSTACK; f++){
      for(int p=0; p<npole; p++){
        const int col = f*npole + p;
        res += D.A[p+1]*( -cuCreal(hA[col]) - cuCreal(hB[col]) );   // A[p+1] = residue of pole p
      }
      res += cuCreal(hI[f]);
    }
    res *= -2.0*D.C/D.lambda_max;
    return res;
  }

  // ---- total force over all links -> dSf (drives the same link loop as Force::compute) ----
  template<typename Force, typename Gauge>
  void compute( Force& dSf, const Gauge& U, const CuC* d_eta_blk ) {
    precalc( U, d_eta_blk );
    for(int s=0; s<Force::Nt; s++)
      for(Idx ell=0; ell<dSf.lattice.n_links; ell++)
        dSf.sp(s, ell) = grad_block( std::make_pair(s, dSf.lattice.links[ell]), U, d_eta_blk );
    for(int s=0; s<Force::Nt; s++)
      for(Idx ix=0; ix<dSf.lattice.n_sites; ix++)
        dSf.tp(s, ix) = grad_block( std::make_pair(s, ix), U, d_eta_blk );
  }

  // ========================= Term B (bilinear / external-bra) block force =========================
  // Hasenbusch ratio-frame cross term: force = -d[2 Re<phi | D_ov | eta>]/dU, K = dD_ov/dU (conserved
  // current), MASS-FREE. Block form of the serial OverlapWMass::precalc_grad_bilinear_deviceAsyncLaunch_ms
  // (overlap_wmass_claude.h:831) + grad_deviceAsyncLaunch_l4 (:1037), packed over the NSTACK (=Nf/2) axis.
  // Independent bra != ket: bra = phi-block, ket = eta-block, both N*NSTACK. Serial->block map:
  //   d_Ys[0]=(1/lmax)D_W^dag phi     -> d_Ys0 = (1/lmax) M_DWH bra
  //   Z_m = R_m (ket=eta)             -> d_Zg  = solve_multishift_block(ket)
  //   Y_m = R_m (X^dag bra)           -> d_Yg  = solve_multishift_block(Ys0)
  //   d_eta_bra = phi  (BARE, no (1+m_L) fold; Term B contracts the raw bra against K -- NO mass branch)
  //   X Z_m / X Y_m                   -> d_XZg / d_XYg
  // The resolvent R_m is the SAME mass-free seed A=(1/lmax^2) M_DWH M_DW, shifts -k^2/cp[m]
  // (B. Jegerlehner, hep-lat/9612014) that Term A uses -- so solve_multishift_block is reused verbatim.
  // Reuses the SAME owned output blocks (d_Zg/d_Yg/d_XZg/d_XYg/d_eta_bra) as Term A: precalc_bilinear
  // OVERWRITES them, so a bilinear precalc + grad must not interleave with a Term A precalc + grad on the
  // same BlockedForce (the caller does the two frames sequentially). NO new buffers are needed -- bra and
  // ket both land in the existing N*NSTACK / N*NSTACK*npole scratch (d_Ys0 holds X^dag bra; d_Zg the ket
  // resolvent; d_Yg the bra resolvent; d_eta_bra the bare bra).

  // ---- precalc_bilinear (link-independent): fill d_Zg(ket), d_Yg(bra), d_XZg, d_XYg, d_eta_bra(=bra) ----
  template<typename Gauge>
  void precalc_bilinear( const Gauge& U, const CuC* d_bra_blk, const CuC* d_ket_blk ) {
    const size_t Ns  = (size_t)N*NSTACK;
    const int gN  = ng(Ns);
    const int gNp = ng((size_t)N*NSTACK*npole);
    const CuC inv_lmax = cplx(1.0/D.lambda_max);

    // Ys0 = (1/lmax) M_DWH bra   (= X^dag bra, the Y-side is built from the BRA=phi)
    mult_coo_block<CuC,N><<<gN,NThreadsPerBlock>>>(d_Ys0, d_bra_blk, D.M_DWH.val, D.M_DWH.cols, D.M_DWH.rows, NSTACK);
    scale_uniform_blk<N><<<gN,NThreadsPerBlock>>>(d_Ys0, inv_lmax, NSTACK);

    // Z_m = R_m ket ; Y_m = R_m Ys0   (block multishift; borrows blk's pool -- Jegerlehner hep-lat/9612014)
    blk.solve_multishift_block(d_Zg, d_ket_blk, Comp::TOL_INNER);
    blk.solve_multishift_block(d_Yg, d_Ys0,     Comp::TOL_INNER);

    // BARE external bra: d_eta_bra = bra (Term B carries no (1+m_L) fold; mass-free, NO mass branch).
    CUDA_CHECK(cudaMemcpy(d_eta_bra, d_bra_blk, Ns*CD, D2D));

    // XZg = (1/lmax) M_DW Zg ; XYg = (1/lmax) M_DW Yg   (block over NSTACK*npole columns)
    mult_coo_block<CuC,N><<<gNp,NThreadsPerBlock>>>(d_XZg, d_Zg, D.M_DW.val, D.M_DW.cols, D.M_DW.rows, NSTACK*npole);
    scale_uniform_blk<N><<<gNp,NThreadsPerBlock>>>(d_XZg, inv_lmax, NSTACK*npole);
    mult_coo_block<CuC,N><<<gNp,NThreadsPerBlock>>>(d_XYg, d_Yg, D.M_DW.val, D.M_DW.cols, D.M_DW.rows, NSTACK*npole);
    scale_uniform_blk<N><<<gNp,NThreadsPerBlock>>>(d_XYg, inv_lmax, NSTACK*npole);

    load_residues();
    CUDA_CHECK(cudaDeviceSynchronize());
  }

  // ---- per-link block bilinear grad = sum_f -2(C/lmax) Re<bra_f | K | ket_f> ----
  // After precalc_bilinear the owned blocks encode (bra,ket): d_eta_bra = bra, d_Zg = R_m ket,
  // d_Yg = hat Y_m from bra, d_XZg/d_XYg the X-applies. grad_block's m=0 combine takes the KET (its
  // d_eta_blk argument) and its bra is d_eta_bra -- so grad_block(link, U, ket) yields EXACTLY the serial
  // bilinear grad_deviceAsyncLaunch_l4(link, U, ket). No separate per-link kernel path is needed.
  template<typename Link, typename Gauge>
  double grad_bilinear_block( const Link& link, const Gauge& U, const CuC* d_ket_blk ) const {
    return grad_block( link, U, d_ket_blk );
  }

  // ---- total bilinear force over all links -> dSf (mirror compute) ----
  // Fills dSf with the RAW bilinear force sum_f -2(C/lmax) Re<phi_f|K|eta_f>; the caller subtracts it
  // (Term B enters the ratio-frame assembly with a minus, serial: dSf += -1.0*ftmp).
  template<typename Force, typename Gauge>
  void compute_bilinear( Force& dSf, const Gauge& U, const CuC* d_bra_blk, const CuC* d_ket_blk ) {
    precalc_bilinear( U, d_bra_blk, d_ket_blk );
    for(int s=0; s<Force::Nt; s++)
      for(Idx ell=0; ell<dSf.lattice.n_links; ell++)
        dSf.sp(s, ell) = grad_bilinear_block( std::make_pair(s, dSf.lattice.links[ell]), U, d_ket_blk );
    for(int s=0; s<Force::Nt; s++)
      for(Idx ix=0; ix<dSf.lattice.n_sites; ix++)
        dSf.tp(s, ix) = grad_bilinear_block( std::make_pair(s, ix), U, d_ket_blk );
  }
};
