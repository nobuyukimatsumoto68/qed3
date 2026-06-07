#pragma once

// multishift_block_kernels_claude.h  (C6a -- mrhs block leaf kernels)
//
// Block (multi-RHS) leaf kernels for the multi-shift CG, widening the single-RHS kernels in
// multishift_kernels_claude.h / gpu_header.h over an extra "stack" axis of NSTACK independent
// right-hand sides. The occupancy lever: the seed matvec at N=3072 fills only ~13 of 80 SMs
// (Compute(SM) 0.8%, DRAM 13%, achieved occupancy 12% -- see the profile gate); batching NSTACK
// RHS into one wider launch raises occupancy toward saturation. Multi-shift (1 Krylov pass for all
// poles) and multi-RHS (this) are orthogonal and compose.
//
// NSTACK is a COMPILE-TIME template parameter (sourced from a Comp:: constexpr per .cu -- its value
// depends on L / hardware / jj-vs-hmc). npole (= size-1, the Zolotarev pole count) stays a RUNTIME
// arg where it appears, since it is set at operator-construction time.
//
// Reference (algorithm): B. Jegerlehner, "Krylov space solvers for shifted linear systems",
// hep-lat/9612014 -- the multi-shift CG these RHS are batched onto.
//
// LAYOUT (column-major; mirrors the single-RHS [m*N+i] of multishift_kernels_claude.h):
//   - seed/outer vectors (p0, r, q): length N*NSTACK, column c contiguous -> index [c*N + i],
//     c in [0,NSTACK), i in [0,N).
//   - per-(rhs,pole) blocks (X, p_m): length N*NSTACK*npole, "column" col = c*npole + m contiguous
//     -> index [(c*npole + m)*N + i] = [col*N + i].
//   - per-column scalars: seed-level a[c] (length NSTACK); per-(rhs,pole) a[col] (length
//     NSTACK*npole). Each thread reads ITS column's scalar by col = gid/N (N compile-time), and
//     recovers c = col/npole, m = col%npole when needed.
// NSTACK=1 reduces every kernel to its single-RHS sibling bit-for-bit (validation, C6b).
//
// Must be included AFTER gpu_header.h (CuC, Idx, the double*CuC / CuC+CuC operators, cplx).
// Grids are sized by the CALLER from the actual work size (N*NSTACK or N*NSTACK*npole); launch
// thread-block size is the global 256 for now (tuning deferred -- see the plan's launch knob).

// Block CSR matvec: res_c = A v_c for every RHS column c, sharing the ONE CSR structure
// (val/cols/rows) across columns. Thread gid -> (c = gid/N, i = gid%N); reads/writes column c.
// res[c*N + i] = sum_{k in row i} csr[k] * v[c*N + cols[k]].
template<typename T, Idx N, int NSTACK> __global__
void mult_block( T* res,
                 const T* v,
                 const T* v_csr,
                 const Idx* cols,
                 const Idx* rows ){
  const Idx gid = (Idx)blockIdx.x*blockDim.x + threadIdx.x;
  const Idx total = (Idx)N*NSTACK;
  if(gid < total){
    const Idx c    = gid / N;            // RHS column
    const Idx i     = gid - c*N;         // lattice index (gid % N)
    const Idx base  = c*N;               // start of column c
    res[gid] = cplx(0.0);
    const int row_start = rows[i];
    const int row_end   = rows[i+1];
    for(int jj=row_start; jj<row_end; jj++) res[gid] = res[gid] + v_csr[jj] * v[ base + cols[jj] ];
  }
}

// Block AXPY with a PER-COLUMN (real) scalar: res_c = a[c] * x_c + y_c, over ncol columns of
// length N (block length N*ncol). The scalar a is a length-ncol double array indexed by c = gid/N.
// Covers the seed-level updates q += sig0_c p0 (a=sig0), r -= al_c q (a=-al), p0 = r + bet_c p0
// (x=p0, y=r, a=bet). ncol is RUNTIME: ncol = NSTACK for seed-level, NSTACK*npole for per-(rhs,pole).
template<Idx N> __global__
void axpy_col( CuC* d_res, const double* d_a, const CuC* d_x, const CuC* d_y, const int ncol ){
  const Idx gid = (Idx)blockIdx.x*blockDim.x + threadIdx.x;
  const Idx total = (Idx)N*ncol;
  if(gid < total){
    const Idx col = gid / N;
    d_res[gid] = d_a[col]*d_x[gid] + d_y[gid];
  }
}

// Block AXPY with a PER-COLUMN COMPLEX scalar: res_c = a[c]*x_c + y_c, over ncol columns of length
// N. a is a length-ncol CuC array indexed by c = gid/N. Block analog of Taxpy<CuC,N>(res, CuC a,x,y)
// used by the outer block CG (C6d), whose al = mu/gam is complex (matches MatPoly::solve).
template<Idx N> __global__
void axpy_col_c( CuC* d_res, const CuC* d_a, const CuC* d_x, const CuC* d_y, const int ncol ){
  const Idx gid = (Idx)blockIdx.x*blockDim.x + threadIdx.x;
  const Idx total = (Idx)N*ncol;
  if(gid < total){
    const Idx col = gid / N;
    d_res[gid] = d_a[col]*d_x[gid] + d_y[gid];
  }
}

// Block AXPBY with TWO per-column (real) scalars: res_c = a[c]*x_c + b[c]*y_c, over ncol columns
// of length N. Used for the seed matvec combine q = coeff_seed * (M1 M0 p0) + sig0 * p0 (a=coeff,
// x=M1 M0 p0; b=sig0, y=p0). a, b are length-ncol double arrays indexed by c = gid/N (ncol RUNTIME).
template<Idx N> __global__
void axpby_col( CuC* d_res, const double* d_a, const CuC* d_x,
                const double* d_b, const CuC* d_y, const int ncol ){
  const Idx gid = (Idx)blockIdx.x*blockDim.x + threadIdx.x;
  const Idx total = (Idx)N*ncol;
  if(gid < total){
    const Idx col = gid / N;
    d_res[gid] = d_a[col]*d_x[gid] + d_b[col]*d_y[gid];
  }
}

// Block AXPY with a UNIFORM (complex) scalar: d_res = a * d_x + d_y, over ncol columns of length N
// (block length N*ncol; ncol RUNTIME). Block analog of Taxpy_gen<CuC,CuC,N>(res,a,x,y) and (with a =
// cplx(double)) of Taxpy_gen<CuC,double,N>. Used by the block overlap methods (C6c) for the C, mass,
// (1/lambda_max) folds -- scalar uniform across columns, passed by value.
template<Idx N> __global__
void axpy_uniform( CuC* d_res, const CuC a, const CuC* d_x, const CuC* d_y, const int ncol ){
  const Idx gid = (Idx)blockIdx.x*blockDim.x + threadIdx.x;
  const Idx total = (Idx)N*ncol;
  if(gid < total) d_res[gid] = a*d_x[gid] + d_y[gid];
}

// Per-pole residue fold: d_out_c = d_base_c + sum_{m=0}^{npole-1} A[m] * Z_{c,m}, over N*NSTACK
// (c=gid/N, i=gid%N). A is a length-npole device array (the Zolotarev residues A[1..size-1]); Z is
// the per-(rhs,pole) block [(c*npole+m)*N+i]. Block analog of the single-RHS loop
// d_Zs0 += A[m] Z_m (mult_ms/adj_ms). base MAY alias out (in-place; each thread touches its own gid).
// Sum order m=0..npole-1 matches the single-RHS accumulation -> bit-identical.
template<Idx N, int NSTACK> __global__
void block_reduce_poles( CuC* d_out, const CuC* d_base, const double* d_A,
                         const CuC* d_Z, const int npole ){
  const Idx gid = (Idx)blockIdx.x*blockDim.x + threadIdx.x;
  const Idx total = (Idx)N*NSTACK;
  if(gid < total){
    const Idx c = gid / N;
    const Idx i = gid - c*N;
    CuC acc = d_base[gid];
    for(int m=0;m<npole;m++) acc = acc + d_A[m]*d_Z[(size_t)(c*npole+m)*N + i];
    d_out[gid] = acc;
  }
}

// Block x-update over N*NSTACK*npole: X_{c,m} += alm_{c,m} * p_{c,m} for every (c,m), site i.
// col = gid/N = c*npole + m indexes both the block column and the per-(rhs,pole) scalar alm.
template<Idx N, int NSTACK> __global__
void multishift_x_update_block( CuC* d_x, const double* d_alm, const CuC* d_p, const int npole ){
  const Idx gid = (Idx)blockIdx.x*blockDim.x + threadIdx.x;
  const Idx total = (Idx)N*NSTACK*npole;
  if(gid < total){
    const Idx col = gid / N;             // = c*npole + m
    d_x[gid] = d_x[gid] + d_alm[col]*d_p[gid];
  }
}

// Block p-update over N*NSTACK*npole: p_{c,m} = zeta_{c,m} * r_c + betm_{c,m} * p_{c,m}.
// The per-RHS seed residual r (length N*NSTACK, layout [c*N+i]) is broadcast over the npole poles
// of the same RHS c. col = gid/N = c*npole+m -> c = col/npole; r index = c*N + (gid%N).
template<Idx N, int NSTACK> __global__
void multishift_p_update_block( CuC* d_p, const double* d_zeta, const CuC* d_r,
                                const double* d_betm, const int npole ){
  const Idx gid = (Idx)blockIdx.x*blockDim.x + threadIdx.x;
  const Idx total = (Idx)N*NSTACK*npole;
  if(gid < total){
    const Idx col = gid / N;             // = c*npole + m
    const Idx i   = gid - col*N;         // lattice index (gid % N)
    const Idx c   = col / npole;         // RHS column
    d_p[gid] = d_zeta[col]*d_r[c*N + i] + d_betm[col]*d_p[gid];
  }
}

// ===== L2 (HMC force opt): block-COO matvec + block-dot with RUNTIME ncol =====
// mult_coo_block: apply a SINGLE CSR (the per-link grad COO) to ncol RHS at once (block length N*ncol);
// runtime ncol = npole = size-1 (the Zolotarev poles). Same body as mult_block but ncol is a kernel arg.
template<typename T, Idx N> __global__
void mult_coo_block( T* res, const T* v, const T* v_csr, const Idx* cols, const Idx* rows, const int ncol ){
  const Idx gid = (Idx)blockIdx.x*blockDim.x + threadIdx.x;
  const Idx total = (Idx)N*ncol;
  if(gid < total){
    const Idx c    = gid / N;            // RHS column
    const Idx i    = gid - c*N;          // lattice index
    const Idx base = c*N;
    res[gid] = cplx(0.0);
    const int row_start = rows[i];
    const int row_end   = rows[i+1];
    for(int jj=row_start; jj<row_end; jj++) res[gid] = res[gid] + v_csr[jj] * v[ base + cols[jj] ];
  }
}

// block_dot: per-column conjugate dot  d_out[c] = sum_i conj(a[c*N+i]) b[c*N+i]  (= cublasZdotc(a_c,b_c)),
// for c = 0..ncol-1. ONE thread block per column reduces N elements (shared-mem reduction). Launch with
// ncol blocks x NThreadsPerBlock threads (NThreadsPerBlock <= 256). Replaces the per-pole streamed
// cublasZdotc + host sync in grad with ONE kernel + ONE memcpy (the L2 host-sync killer).
template<Idx N> __global__
void block_dot( CuC* d_out, const CuC* d_a, const CuC* d_b, const int ncol ){
  const int c = blockIdx.x;
  if(c >= ncol) return;
  __shared__ double sh_re[256];
  __shared__ double sh_im[256];
  const Idx base = (Idx)c*N;
  double re=0.0, im=0.0;
  for(Idx i=threadIdx.x; i<N; i+=blockDim.x){
    const CuC a = d_a[base+i], b = d_b[base+i];
    re += cuCreal(a)*cuCreal(b) + cuCimag(a)*cuCimag(b);   // Re conj(a)b
    im += cuCreal(a)*cuCimag(b) - cuCimag(a)*cuCreal(b);   // Im conj(a)b
  }
  sh_re[threadIdx.x]=re; sh_im[threadIdx.x]=im;
  __syncthreads();
  for(int s=blockDim.x/2; s>0; s>>=1){
    if(threadIdx.x<s){ sh_re[threadIdx.x]+=sh_re[threadIdx.x+s]; sh_im[threadIdx.x]+=sh_im[threadIdx.x+s]; }
    __syncthreads();
  }
  if(threadIdx.x==0) d_out[c] = make_cuDoubleComplex(sh_re[0], sh_im[0]);
}

// L4 (HMC force): single-link matvec from the RAW COO entries -- SKIPS coo.do_it() (its O(N) CSR-rows
// loop + 3 cudaMalloc, ~11% of grad_l2). out[c*N+ei[k]] += ev[k] * in[c*N+ej[k]] over nent entries x ncol
// RHS; out MUST be zeroed first (memset). atomicAdd: a single-link COO's entries are on a few distinct
// rows so contention is ~nil (=> deterministic); atomics just keep it correct if any row repeats.
__device__ inline void atomicAddCuC( CuC* a, const CuC b ){
  atomicAdd( &(reinterpret_cast<double*>(a)[0]), cuCreal(b) );
  atomicAdd( &(reinterpret_cast<double*>(a)[1]), cuCimag(b) );
}
template<Idx N> __global__
void link_matvec_block( CuC* out, const CuC* in, const Idx* ei, const Idx* ej, const CuC* ev,
                        const int nent, const int ncol ){
  const int gid = blockIdx.x*blockDim.x + threadIdx.x;
  if(gid < nent*ncol){
    const int c = gid / nent;
    const int k = gid - c*nent;
    atomicAddCuC( &out[(size_t)c*N + ei[k]], ev[k] * in[(size_t)c*N + ej[k]] );
  }
}
