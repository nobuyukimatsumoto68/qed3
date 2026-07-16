// test_matvec_mixed_claude.cu
// _claude: Chunk-0 GATE of the mixed-precision inner-solver plan (mixedprec_impl_plan_claude.md).
// Times the raw D_W^dag D_W matvec (the hot op of the overlap inner CG) in fp64 vs fp32 at a given
// L, and reports the fp32-vs-fp64 rel error + speedup. Purpose: confirm the fp32 bandwidth win is
// worthwhile before building the reliable-update CG machinery (small-N L1 may be latency-bound).
//
// D_W^dag D_W v = M_DWH ( M_DW v ) -- the seed matrix of the Zolotarev multishift, minus the
// (1/lambda_max^2) scale (irrelevant to timing). fp32 path casts the CSR val once, then applies the
// two fp32 CSR matvecs on fp32 vectors (exactly what an fp32 Krylov iteration does).
//
// Build: mirror a hmc_* target; run per L via -DLREF={1,2,4} on a GPU. Handoff:
// tmp_matvec_mixed_build_claude.sh.

#include <typeinfo>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <filesystem>
#include <cstdint>
#include <complex>
#include <array>
#include <vector>
#include <map>
#include <random>
#include <Eigen/Dense>

using Double = double;
using Idx = std::int32_t;
using Complex = std::complex<double>;
using Face = std::vector<Idx>;
using MS=Eigen::Matrix2cd;
using VD=Eigen::Vector2d;
using VE=Eigen::Vector3d;
using VC=Eigen::VectorXcd;

static constexpr int NS = 2;
static constexpr int DIM = 2;
static constexpr Complex I = Complex(0.0, 1.0);

namespace Comp{
  constexpr bool is_compact=false;
  constexpr int NPARALLEL_DUPDATE=1;
  constexpr int NPARALLEL=NPARALLEL_DUPDATE;
  constexpr int NSTREAMS=NPARALLEL_DUPDATE;
  constexpr int NPARALLEL_GAUGE=16;
  constexpr int NPARALLEL_SORT=16;

#ifndef LREF
#define LREF 1
#endif
  constexpr int N_REFINE=LREF;
  constexpr int NS=2;
#ifndef NTIME
#define NTIME 16
#endif
  constexpr int Nt=NTIME;      // matvec cost scales with N = 2*N_SITES*Nt; use a production-ish Nt

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;
  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-8;
}

const std::string dir = "../../geometry/data/";
#include "../../geometry/geodesic.h"

#include "timer.h"
#include "s2n_simp.h"
#include "rng.h"

#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cublas_api.h>
#include <cusolverDn.h>
using CuC = cuDoubleComplex;
#include "gpu_header.h"

#include "valence.h"
#include "gauge_ext.h"
#include "action_ext.h"

#include "sparse_matrix_claude.h"
#include "dirac_simp.h"
#include "dirac_ext.h"
#include "sparse_dirac_claude.h"
#include "matpoly_claude.h"
#include "includes/overlap_wmass_claude.h"
#include "includes/sparse_matrix_mixed_claude.h"   // fp32 CSR matvec + casts (Chunk 0)

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(6);

  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();

  constexpr Idx N = Comp::N;
  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;

  Base base(Comp::N_REFINE);
  const double M5 = -1.0, at = 0.2, nu0 = 1.0;
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);

  using Rng = ParallelRngExt<Base, Comp::Nt>;
  Rng rng(base);
  Gauge U(base);
  U.gaussian(rng, 0.3);

  const int npole = 21;
  OverlapWMass<WilsonDirac> Dov(DW, Complex(0.0,0.0), npole);
  Dov.update(U);   // fills M_DW / M_DWH Wilson CSR values (and lambda_max)

  const Idx nnz = Dov.d_DW.len;
  std::cout << "# L=" << Comp::N_REFINE << "  N=" << N << "  Nt=" << Comp::Nt
            << "  nnz(D_W)=" << nnz << std::endl;

  // ----- fp32 mirrors of M_DW (val=d_val, cols=d_cols, rows=d_rows) and
  //       M_DWH (val=d_valH, cols=d_colsT, rows=d_rowsT) -----
  CSRf<N> dw_f, dwh_f;
  dw_f .associate( Dov.d_DW.d_cols,  Dov.d_DW.d_rows,  nnz );
  dwh_f.associate( Dov.d_DW.d_colsT, Dov.d_DW.d_rowsT, nnz );
  dw_f .cast_from( Dov.d_DW.d_val  );
  dwh_f.cast_from( Dov.d_DW.d_valH );

  // random source
  std::vector<Complex> xi(N);
  std::mt19937_64 gen(12345ull);
  std::normal_distribution<double> nd(0.0,1.0);
  for(Idx i=0;i<N;i++) xi[i] = Complex(nd(gen), nd(gen));

  CuC  *d_in=nullptr, *d_t64=nullptr, *d_out64=nullptr;
  CuCf *d_inf=nullptr, *d_tf=nullptr, *d_outf=nullptr;
  CUDA_CHECK(cudaMalloc(&d_in,    N*CD));
  CUDA_CHECK(cudaMalloc(&d_t64,   N*CD));
  CUDA_CHECK(cudaMalloc(&d_out64, N*CD));
  CUDA_CHECK(cudaMalloc(&d_inf,   N*CDf));
  CUDA_CHECK(cudaMalloc(&d_tf,    N*CDf));
  CUDA_CHECK(cudaMalloc(&d_outf,  N*CDf));
  CUDA_CHECK(cudaMemcpy(d_in, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));
  cast_z2c_launch(d_inf, d_in, N);   // fp32 input (kept in fp32 across the loop, as in an fp32 CG)

  // ===== fp64 D_W^dag D_W apply: out = M_DWH ( M_DW in ) =====
  auto dhdh_fp64 = [&](){
    Dov.M_DW ( d_t64,   d_in  );
    Dov.M_DWH( d_out64, d_t64 );
  };
  // ===== fp32 D_W^dag D_W apply: outf = M_DWH_f ( M_DW_f inf ) =====
  auto dhdh_fp32 = [&](){
    dw_f .apply( d_tf,   d_inf );
    dwh_f.apply( d_outf, d_tf  );
  };

  // ----- accuracy: fp32 result vs fp64 result (relative L2) -----
  dhdh_fp64();
  dhdh_fp32();
  CUDA_CHECK(cudaDeviceSynchronize());
  std::vector<Complex> o64(N), o32(N);
  CuC* d_o32_as64=nullptr;  CUDA_CHECK(cudaMalloc(&d_o32_as64, N*CD));
  cast_c2z_launch(d_o32_as64, d_outf, N);
  CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(o64.data()), d_out64,    N*CD, D2H));
  CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(o32.data()), d_o32_as64, N*CD, D2H));
  double nd2=0.0, nb2=0.0;
  for(Idx i=0;i<N;i++){ nd2 += std::norm(o32[i]-o64[i]); nb2 += std::norm(o64[i]); }
  const double relerr = std::sqrt(nd2)/std::max(std::sqrt(nb2),1.0e-300);
  std::cout << "# fp32-vs-fp64 matvec rel err = " << relerr << std::endl;

  // ----- timing (cudaEvent), NREP applies after a warmup -----
  const int NWARM = 20;
  const int NREP  = 400;
  cudaEvent_t e0, e1;
  CUDA_CHECK(cudaEventCreate(&e0));
  CUDA_CHECK(cudaEventCreate(&e1));

  for(int r=0;r<NWARM;r++) dhdh_fp64();
  CUDA_CHECK(cudaDeviceSynchronize());
  CUDA_CHECK(cudaEventRecord(e0));
  for(int r=0;r<NREP;r++) dhdh_fp64();
  CUDA_CHECK(cudaEventRecord(e1));
  CUDA_CHECK(cudaEventSynchronize(e1));
  float ms64=0.0f;  CUDA_CHECK(cudaEventElapsedTime(&ms64, e0, e1));

  for(int r=0;r<NWARM;r++) dhdh_fp32();
  CUDA_CHECK(cudaDeviceSynchronize());
  CUDA_CHECK(cudaEventRecord(e0));
  for(int r=0;r<NREP;r++) dhdh_fp32();
  CUDA_CHECK(cudaEventRecord(e1));
  CUDA_CHECK(cudaEventSynchronize(e1));
  float ms32=0.0f;  CUDA_CHECK(cudaEventElapsedTime(&ms32, e0, e1));

  const double us64 = 1.0e3*ms64/NREP;   // microseconds per D_W^dag D_W apply
  const double us32 = 1.0e3*ms32/NREP;
  std::cout << "# fp64 matvec: " << us64 << " us/apply" << std::endl;
  std::cout << "# fp32 matvec: " << us32 << " us/apply" << std::endl;
  std::cout << "# SPEEDUP (fp64/fp32) = " << (us64/us32) << "x" << std::endl;

  CUDA_CHECK(cudaFree(d_in));   CUDA_CHECK(cudaFree(d_t64));  CUDA_CHECK(cudaFree(d_out64));
  CUDA_CHECK(cudaFree(d_inf));  CUDA_CHECK(cudaFree(d_tf));   CUDA_CHECK(cudaFree(d_outf));
  CUDA_CHECK(cudaFree(d_o32_as64));
  return 0;
}
