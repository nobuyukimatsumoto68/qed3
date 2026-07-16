#pragma once

// sparse_matrix_mixed_claude.h -- fp32 ("mixed inner precision") companions to the fp64
// sparse CSR matvec in gpu_header.h / sparse_matrix_claude.h. Chunk 0 of the mixed-precision
// inner-solver plan (mixedprec_impl_plan_claude.md): supplies the low-precision matvec + the
// fp64<->fp32 cast kernels so the reliable-update CG (Chunk 1+) can run its Krylov iteration in
// single precision while the true residual is corrected in double.
//
// Precision policy: the INNER Krylov type is cuFloatComplex (CuCf); the OUTER accumulate/residual
// stays cuDoubleComplex (CuC). Kept as a separate header (not an edit of the fp64 versions) per NM's
// preserve-original convention. Algorithm ref: M.A. Clark et al., arXiv:0911.3191.

#include <cuComplex.h>

using CuCf = cuFloatComplex;
#define CDf (sizeof(cuFloatComplex))

// ---- minimal cuFloatComplex arithmetic (mirror of gpu_header.h's CuC operators) --------------
__device__ __host__ inline CuCf cplxf(const float c) { return make_cuFloatComplex(c, 0.0f); }
__device__ __host__ inline CuCf operator+(CuCf a, CuCf b) { return cuCaddf(a,b); }
__device__ __host__ inline CuCf operator-(CuCf a, CuCf b) { return cuCsubf(a,b); }
__device__ __host__ inline CuCf operator*(CuCf a, CuCf b) { return cuCmulf(a,b); }
__device__ __host__ inline CuCf operator*(CuCf a, float b) { return cuCmulf(a, cplxf(b)); }
__device__ __host__ inline CuCf operator*(float a, CuCf b) { return cuCmulf(cplxf(a), b); }

// ---- fp64 <-> fp32 element-wise casts (n = element count) -------------------------------------
__global__ void cast_z2c( CuCf* out, const CuC* in, const Idx n ){
  Idx i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<n) out[i] = make_cuFloatComplex( (float)cuCreal(in[i]), (float)cuCimag(in[i]) );
}

__global__ void cast_c2z( CuC* out, const CuCf* in, const Idx n ){
  Idx i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<n) out[i] = make_cuDoubleComplex( (double)cuCrealf(in[i]), (double)cuCimagf(in[i]) );
}

// host helpers: launch a cast over n elements (generic block count, not the N-sized NBlocks macro).
inline void cast_z2c_launch( CuCf* d_out, const CuC* d_in, const Idx n ){
  const int nb = (n + NThreadsPerBlock)/NThreadsPerBlock;
  cast_z2c<<<nb, NThreadsPerBlock>>>(d_out, d_in, n);
}
inline void cast_c2z_launch( CuC* d_out, const CuCf* d_in, const Idx n ){
  const int nb = (n + NThreadsPerBlock)/NThreadsPerBlock;
  cast_c2z<<<nb, NThreadsPerBlock>>>(d_out, d_in, n);
}

// ---- fp32 CSR matvec (single-precision copy of gpu_header.h `mult<T,N>`) -----------------------
template<Idx N> __global__
void mult_f( CuCf* res,
             const CuCf* v,
             const CuCf* v_csr,
             const Idx* cols,
             const Idx* rows ){
  Idx i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<N) {
    CuCf s = cplxf(0.0f);
    const int row_start = rows[i];
    const int row_end   = rows[i+1];
    for(int jj=row_start; jj<row_end; jj++) s = s + v_csr[jj] * v[ cols[jj] ];
    res[i] = s;
  }
}

// ---- fp32 vector primitives for the inner Krylov (mirror gpu_header.h Taxpy) -----------------
// d_res = a * d_x + d_y   (fp32 complex axpy; in-place safe when d_res aliases d_x or d_y).
template<Idx N> __global__
void Taxpy_f(CuCf* d_res, const CuCf a, const CuCf* d_x, const CuCf* d_y){
  Idx i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<N) d_res[i] = a * d_x[i] + d_y[i];
}

// out = sc * Av + sh * v   (fused scale+shift; used to add the Zolotarev pole shift sh and the
// (1/lambda_max^2) scale sc to the raw D_W^dag D_W apply Av). fp32 and fp64 variants.
template<Idx N> __global__
void axpby_shift_f(CuCf* out, const CuCf* Av, const CuCf* v, const float sc, const float sh){
  Idx i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<N) out[i] = sc * Av[i] + sh * v[i];
}
template<Idx N> __global__
void axpby_shift_d(CuC* out, const CuC* Av, const CuC* v, const double sc, const double sh){
  Idx i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<N) out[i] = sc * Av[i] + sh * v[i];
}

// fp32 CSR operator mirroring CSR<N> (sparse_matrix_claude.h). Does NOT inherit LinOp (that
// interface is fp64-typed). Holds a device fp32 val array cast from a fp64 CSR; shares cols/rows
// pointers with the source fp64 CSR (same sparsity pattern -- integer indices are precision-free).
//   build(d_val_fp64, d_cols, d_rows, nnz): (re)cast the fp32 val from the fp64 source.
// Call build() after every gauge update (in the FORCE, every MD step) -- a single cheap pass.
template<Idx N>
struct CSRf {
  CuCf* val = nullptr;         // owned fp32 values (nnz)
  const Idx* cols = nullptr;   // borrowed (shared with the fp64 CSR)
  const Idx* rows = nullptr;   // borrowed
  Idx nnz = 0;
  bool is_set = false;

  ~CSRf(){ if(is_set) CUDA_CHECK(cudaFree(val)); }

  // associate the (fixed) sparsity pattern and allocate the fp32 val buffer.
  void associate( const Idx* d_cols, const Idx* d_rows, const Idx nnz_ ){
    cols = d_cols;
    rows = d_rows;
    nnz  = nnz_;
    CUDA_CHECK(cudaMalloc(&val, (size_t)nnz*CDf));
    is_set = true;
  }

  // (re)cast the fp32 values from a fp64 source val array (same nnz ordering).
  void cast_from( const CuC* d_val_fp64 ) const {
    cast_z2c_launch( val, d_val_fp64, nnz );
  }

  // res = A v  (fp32), single CSR apply.
  void apply( CuCf* d_res, const CuCf* d_v ) const {
    mult_f<N><<<NBlocks, NThreadsPerBlock>>>(d_res, d_v, val, cols, rows);
  }
};
