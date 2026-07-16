#pragma once

// pseudofermion_parallelized_claude.h -- BLOCK (mrhs) pseudo-fermion manager.
// Packs the Nf/2 = NSTACK flavors' phi/eta as ONE N*NSTACK device block and drives the ACTION
// inversions (heat-bath gen, update_eta, S) through the validated BlockedMat (blocked_mat_claude.h),
// replacing the NSTACK serial PseudoFermion solves of the integrator's per-pf loop. The FORCE stays
// SERIAL per flavor -- read column f via eta_col(f); a block force is a later step (eta_col is the seam).
//
// Operator equivalence with the serial PseudoFermion (pseudofermion_claude.h):
//   gen        : phi_c = D^dag_ov xi_c  (block adj == D.adj_deviceAsyncLaunch_ms per column)
//   update_eta : eta   = (D^dag D)^{-1} phi   (block solve_sq == Op_DHD.solve, TOL_OUTER)
//   S          : sum_c Re( phi_c^dag eta_c )  (MatPoly::dot is a plain cublasZdotc)
// Column fill order in gen matches the serial for(pf) pf->gen RNG draw order, so a given start
// config+RNG reproduces the serial phi (hence a comparable trajectory).
//
// Reference (algorithm): B. Jegerlehner, "Krylov space solvers for shifted linear systems",
// hep-lat/9612014 (inner multishift) + standard block CG (outer). Include AFTER blocked_mat_claude.h.

#include <cmath>

template<typename Fermion, int NSTACK>
struct PseudoFermionBlock {
  static constexpr Idx N = Comp::N;

  Fermion& D;
  BlockMemPool<N,NSTACK> pool;          // OWNS the block scratch (pole blocks needed for solve_sq/adj)
  BlockedMat<N,NSTACK,Fermion> blk;     // binds the pool; provides adj / solve_sq

  CuC* d_phi_blk;    // N*NSTACK : column c = flavor c heat-bath source
  CuC* d_eta_blk;    // N*NSTACK : column c = (D^dag D)^{-1} phi_c
  CuC* d_xi_blk;     // N*NSTACK : gen scratch (device)
  Complex* xi_host;  // N*NSTACK : gen scratch (host, pinned)

  PseudoFermionBlock()=delete;

  explicit PseudoFermionBlock( Fermion& D_ )
    : D(D_)
    , pool(D_.size-1, true)             // npole = D.size-1 ; with_pole_blocks = true
    , blk(D_, pool)
  {
    CUDA_CHECK(cudaMalloc(&d_phi_blk, (size_t)N*NSTACK*CD));
    CUDA_CHECK(cudaMalloc(&d_eta_blk, (size_t)N*NSTACK*CD));
    CUDA_CHECK(cudaMalloc(&d_xi_blk,  (size_t)N*NSTACK*CD));
    CUDA_CHECK(cudaMallocHost(&xi_host, (size_t)N*NSTACK*CD));
  }

  ~PseudoFermionBlock(){
    CUDA_CHECK(cudaFree(d_phi_blk));
    CUDA_CHECK(cudaFree(d_eta_blk));
    CUDA_CHECK(cudaFree(d_xi_blk));
    CUDA_CHECK(cudaFreeHost(xi_host));
  }

  PseudoFermionBlock(const PseudoFermionBlock&)=delete;
  PseudoFermionBlock& operator=(const PseudoFermionBlock&)=delete;

  // heat-bath: xi ~ Gaussian (column-by-column in flavor order), phi = D^dag_ov xi, then update_eta.
  template<class Rng>
  void gen( Rng& rng ) {
    for(int c=0; c<NSTACK; c++) rng.fill_gaussian( xi_host + (size_t)c*N );
    CUDA_CHECK(cudaMemcpy(d_xi_blk, reinterpret_cast<CuC*>(xi_host), (size_t)N*NSTACK*CD, H2D));
    blk.adj( d_phi_blk, d_xi_blk );                    // phi = D^dag_ov xi  (block multishift)
    update_eta();
  }

  // outer block CG over D^dag D : eta = (D^dag D)^{-1} phi for all NSTACK columns at once.
  inline void update_eta() { blk.solve_sq( d_eta_blk, d_phi_blk, Comp::TOL_OUTER ); }

  // S = sum_c Re( phi_c^dag eta_c )  (per-flavor pseudo-fermion action, summed over the block).
  double S() const {
    double res = 0.0;
    for(int c=0; c<NSTACK; c++){
      CuC d;
      CUBLAS_CHECK(cublasZdotc(blk.handle, N,
                               d_phi_blk + (size_t)c*N, 1,
                               d_eta_blk + (size_t)c*N, 1, &d));
      res += real(d);
    }
    return res;
  }

  // FORCE seam: pointer to column f of eta (the serial per-flavor force reads this).
  inline CuC*       eta_col(int f)       { return d_eta_blk + (size_t)f*N; }
  inline const CuC* eta_col(int f) const { return d_eta_blk + (size_t)f*N; }

  // Device-memory footprint of the block scratch (startup device-limit check). Counts the BlockMemPool
  // CuC arrays (13 of N*NSTACK + 3 of N*NSTACK*npole) plus this struct's own phi/eta/xi (3 of N*NSTACK);
  // the small NSTACK / NSTACK*npole scalar arrays are negligible.
  static size_t block_bytes(int npole) {
    const size_t Ns  = (size_t)N*NSTACK;
    const size_t Nsp = (size_t)N*NSTACK*npole;
    return ( (size_t)16*Ns + (size_t)3*Nsp ) * CD;
  }
};
