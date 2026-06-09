// test_multishift_block_kernels_claude.cu
// C6a (device): validate the mrhs block leaf kernels (multishift_block_kernels_claude.h) against
// their single-RHS siblings, column by column. Self-contained: a random CSR matrix and random
// vectors -- no lattice setup -- since C6a only exercises the block index arithmetic.
//
// Block layouts (column-major): seed/outer [c*N+i] (c in [0,nstack)); per-(rhs,pole) [(c*npole+m)*N+i]
// (col = c*npole+m). Each block kernel reduces, on the contiguous column range of RHS c, to its
// single-RHS sibling run on the sub-block -> the results must be BIT-IDENTICAL (same FP ops).
// PASS = max|block - per-column-single| < 1e-12 for all four kernels.
//
//   mult_block               vs  mult                 (CSR matvec)
//   axpy_col                 vs  Taxpy(cplx(a[c]),..)  (per-column AXPY)
//   multishift_x_update_block vs multishift_x_update   (x_m += alm p_m)
//   multishift_p_update_block vs multishift_p_update   (p_m = zeta r + betm p_m)
//
// Compile: handled by the both_3d Makefile (auto-picks every *.cu).
// Run: ./test_multishift_block_kernels_claude.o

#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <vector>
#include <random>
#include <complex>
#include <algorithm>

#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cublas_v2.h>

using Idx = std::int32_t;
using CuC = cuDoubleComplex;

#include "gpu_header.h"
#include "multishift_kernels_claude.h"
#include "multishift_block_kernels_claude.h"

// test sizes (small; only index math matters)
static constexpr Idx N      = 64;
static constexpr int NSTACK = 3;
static constexpr int NPOLE  = 5;
static constexpr int NNZ_PER_ROW = 5;

static inline int ngrid(Idx total){ return (int)((total + NThreadsPerBlock - 1)/NThreadsPerBlock); }

// max | a - b | over n complex entries (host arrays)
static double maxdiff(const std::vector<CuC>& a, const std::vector<CuC>& b){
  double d = 0.0;
  for(size_t i=0;i<a.size();++i){
    const double dr = a[i].x - b[i].x, di = a[i].y - b[i].y;
    d = std::max(d, std::sqrt(dr*dr + di*di));
  }
  return d;
}

int main(){
  std::mt19937 rng(20260606);
  std::uniform_real_distribution<double> U(-1.0, 1.0);
  std::uniform_int_distribution<Idx> Ucol(0, N-1);

  // ---- random CSR (N x N, ~NNZ_PER_ROW distinct cols/row), complex values ----
  std::vector<Idx> rows(N+1, 0), cols;
  std::vector<CuC> val;
  for(Idx i=0;i<N;++i){
    std::vector<Idx> rc;
    for(int k=0;k<NNZ_PER_ROW;++k){ Idx c=Ucol(rng); if(std::find(rc.begin(),rc.end(),c)==rc.end()) rc.push_back(c); }
    std::sort(rc.begin(), rc.end());
    for(Idx c: rc){ cols.push_back(c); val.push_back( make_cuDoubleComplex(U(rng),U(rng)) ); }
    rows[i+1] = (Idx)cols.size();
  }
  const Idx nnz = (Idx)cols.size();

  Idx *d_rows,*d_cols; CuC *d_val;
  CUDA_CHECK(cudaMalloc(&d_rows,(N+1)*sizeof(Idx)));
  CUDA_CHECK(cudaMalloc(&d_cols, nnz*sizeof(Idx)));
  CUDA_CHECK(cudaMalloc(&d_val,  nnz*CD));
  CUDA_CHECK(cudaMemcpy(d_rows, rows.data(), (N+1)*sizeof(Idx), H2D));
  CUDA_CHECK(cudaMemcpy(d_cols, cols.data(), nnz*sizeof(Idx), H2D));
  CUDA_CHECK(cudaMemcpy(d_val,  val.data(),  nnz*CD, H2D));

  auto rand_block = [&](Idx n){ std::vector<CuC> h(n); for(Idx i=0;i<n;++i) h[i]=make_cuDoubleComplex(U(rng),U(rng)); return h; };
  bool all_pass = true;
  std::cout << std::scientific << std::setprecision(3);
  std::cout << "# C6a block leaf kernels   N="<<N<<" nstack="<<NSTACK<<" npole="<<NPOLE<<"\n";

  // ======================= Test 1: mult_block vs mult =======================
  {
    const Idx ncol = NSTACK, len = N*ncol;
    auto hv = rand_block(len);
    CuC *d_v,*d_blk,*d_ref;
    CUDA_CHECK(cudaMalloc(&d_v,len*CD)); CUDA_CHECK(cudaMalloc(&d_blk,len*CD)); CUDA_CHECK(cudaMalloc(&d_ref,len*CD));
    CUDA_CHECK(cudaMemcpy(d_v,hv.data(),len*CD,H2D));

    mult_block<CuC,N,NSTACK><<<ngrid(len),NThreadsPerBlock>>>(d_blk, d_v, d_val, d_cols, d_rows);
    for(Idx c=0;c<ncol;++c) mult<CuC,N><<<ngrid(N),NThreadsPerBlock>>>(d_ref + c*N, d_v + c*N, d_val, d_cols, d_rows);
    CUDA_CHECK(cudaDeviceSynchronize());

    std::vector<CuC> hblk(len),href(len);
    CUDA_CHECK(cudaMemcpy(hblk.data(),d_blk,len*CD,D2H)); CUDA_CHECK(cudaMemcpy(href.data(),d_ref,len*CD,D2H));
    const double d = maxdiff(hblk,href); const bool ok = d<1e-12; all_pass &= ok;
    std::cout << "# mult_block               : max|diff|="<<d<<"  ("<<(ok?"PASS":"FAIL")<<")\n";
    cudaFree(d_v); cudaFree(d_blk); cudaFree(d_ref);
  }

  // ======================= Test 2: axpy_col vs Taxpy =======================
  {
    const Idx ncol = NSTACK, len = N*ncol;
    auto hx = rand_block(len), hy = rand_block(len);
    std::vector<double> ha(ncol); for(Idx c=0;c<ncol;++c) ha[c]=U(rng);
    CuC *d_x,*d_y,*d_blk,*d_ref; double* d_a;
    CUDA_CHECK(cudaMalloc(&d_x,len*CD)); CUDA_CHECK(cudaMalloc(&d_y,len*CD));
    CUDA_CHECK(cudaMalloc(&d_blk,len*CD)); CUDA_CHECK(cudaMalloc(&d_ref,len*CD));
    CUDA_CHECK(cudaMalloc(&d_a,ncol*DB));
    CUDA_CHECK(cudaMemcpy(d_x,hx.data(),len*CD,H2D)); CUDA_CHECK(cudaMemcpy(d_y,hy.data(),len*CD,H2D));
    CUDA_CHECK(cudaMemcpy(d_a,ha.data(),ncol*DB,H2D));

    axpy_col<N><<<ngrid(len),NThreadsPerBlock>>>(d_blk, d_a, d_x, d_y, ncol);
    for(Idx c=0;c<ncol;++c) Taxpy<CuC,N><<<ngrid(N),NThreadsPerBlock>>>(d_ref + c*N, cplx(ha[c]), d_x + c*N, d_y + c*N);
    CUDA_CHECK(cudaDeviceSynchronize());

    std::vector<CuC> hblk(len),href(len);
    CUDA_CHECK(cudaMemcpy(hblk.data(),d_blk,len*CD,D2H)); CUDA_CHECK(cudaMemcpy(href.data(),d_ref,len*CD,D2H));
    const double d = maxdiff(hblk,href); const bool ok = d<1e-12; all_pass &= ok;
    std::cout << "# axpy_col                 : max|diff|="<<d<<"  ("<<(ok?"PASS":"FAIL")<<")\n";
    cudaFree(d_x);cudaFree(d_y);cudaFree(d_blk);cudaFree(d_ref);cudaFree(d_a);
  }

  // ============== Test 3: multishift_x_update_block vs multishift_x_update ==============
  {
    const Idx ncol = (Idx)NSTACK*NPOLE, len = N*ncol;
    auto hx = rand_block(len), hp = rand_block(len);
    std::vector<double> halm(ncol); for(Idx col=0;col<ncol;++col) halm[col]=U(rng);
    CuC *d_xb,*d_xr,*d_p; double* d_alm;
    CUDA_CHECK(cudaMalloc(&d_xb,len*CD)); CUDA_CHECK(cudaMalloc(&d_xr,len*CD)); CUDA_CHECK(cudaMalloc(&d_p,len*CD));
    CUDA_CHECK(cudaMalloc(&d_alm,ncol*DB));
    CUDA_CHECK(cudaMemcpy(d_xb,hx.data(),len*CD,H2D)); CUDA_CHECK(cudaMemcpy(d_xr,hx.data(),len*CD,H2D));
    CUDA_CHECK(cudaMemcpy(d_p,hp.data(),len*CD,H2D)); CUDA_CHECK(cudaMemcpy(d_alm,halm.data(),ncol*DB,H2D));

    multishift_x_update_block<N,NSTACK><<<ngrid(len),NThreadsPerBlock>>>(d_xb, d_alm, d_p, NPOLE);
    for(int c=0;c<NSTACK;++c)
      multishift_x_update<N><<<ngrid((Idx)N*NPOLE),NThreadsPerBlock>>>(d_xr + (Idx)c*NPOLE*N, d_alm + c*NPOLE, d_p + (Idx)c*NPOLE*N, NPOLE);
    CUDA_CHECK(cudaDeviceSynchronize());

    std::vector<CuC> hblk(len),href(len);
    CUDA_CHECK(cudaMemcpy(hblk.data(),d_xb,len*CD,D2H)); CUDA_CHECK(cudaMemcpy(href.data(),d_xr,len*CD,D2H));
    const double d = maxdiff(hblk,href); const bool ok = d<1e-12; all_pass &= ok;
    std::cout << "# multishift_x_update_block: max|diff|="<<d<<"  ("<<(ok?"PASS":"FAIL")<<")\n";
    cudaFree(d_xb);cudaFree(d_xr);cudaFree(d_p);cudaFree(d_alm);
  }

  // ============== Test 4: multishift_p_update_block vs multishift_p_update ==============
  {
    const Idx ncol = (Idx)NSTACK*NPOLE, len = N*ncol, rlen = (Idx)N*NSTACK;
    auto hp = rand_block(len), hr = rand_block(rlen);
    std::vector<double> hzeta(ncol), hbetm(ncol);
    for(Idx col=0;col<ncol;++col){ hzeta[col]=U(rng); hbetm[col]=U(rng); }
    CuC *d_pb,*d_pr,*d_r; double *d_zeta,*d_betm;
    CUDA_CHECK(cudaMalloc(&d_pb,len*CD)); CUDA_CHECK(cudaMalloc(&d_pr,len*CD)); CUDA_CHECK(cudaMalloc(&d_r,rlen*CD));
    CUDA_CHECK(cudaMalloc(&d_zeta,ncol*DB)); CUDA_CHECK(cudaMalloc(&d_betm,ncol*DB));
    CUDA_CHECK(cudaMemcpy(d_pb,hp.data(),len*CD,H2D)); CUDA_CHECK(cudaMemcpy(d_pr,hp.data(),len*CD,H2D));
    CUDA_CHECK(cudaMemcpy(d_r,hr.data(),rlen*CD,H2D));
    CUDA_CHECK(cudaMemcpy(d_zeta,hzeta.data(),ncol*DB,H2D)); CUDA_CHECK(cudaMemcpy(d_betm,hbetm.data(),ncol*DB,H2D));

    multishift_p_update_block<N,NSTACK><<<ngrid(len),NThreadsPerBlock>>>(d_pb, d_zeta, d_r, d_betm, NPOLE);
    for(int c=0;c<NSTACK;++c)
      multishift_p_update<N><<<ngrid((Idx)N*NPOLE),NThreadsPerBlock>>>(d_pr + (Idx)c*NPOLE*N, d_zeta + c*NPOLE, d_r + (Idx)c*N, d_betm + c*NPOLE, NPOLE);
    CUDA_CHECK(cudaDeviceSynchronize());

    std::vector<CuC> hblk(len),href(len);
    CUDA_CHECK(cudaMemcpy(hblk.data(),d_pb,len*CD,D2H)); CUDA_CHECK(cudaMemcpy(href.data(),d_pr,len*CD,D2H));
    const double d = maxdiff(hblk,href); const bool ok = d<1e-12; all_pass &= ok;
    std::cout << "# multishift_p_update_block: max|diff|="<<d<<"  ("<<(ok?"PASS":"FAIL")<<")\n";
    cudaFree(d_pb);cudaFree(d_pr);cudaFree(d_r);cudaFree(d_zeta);cudaFree(d_betm);
  }

  cudaFree(d_rows);cudaFree(d_cols);cudaFree(d_val);
  std::cout << "# C6a RESULT: " << (all_pass ? "ALL PASS" : "FAIL") << "\n";
  return all_pass ? 0 : 1;
}
