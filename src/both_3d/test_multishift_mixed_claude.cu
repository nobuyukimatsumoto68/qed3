// test_multishift_mixed_claude.cu
// _claude: Chunk-2 validation of the mixed-precision MULTISHIFT (matpoly_mixed_claude.h,
// mixedprec_impl_plan_claude.md). Solves (scale D_W^dag D_W + sigma_m) X_m = b for ALL Zolotarev
// poles m in ONE Krylov pass, with the fp64 reference MatPoly::solve_multishift and the mixed
// MixedSolver::solve_multishift (fp32 matvec -> crossover at tol_f -> fp64 tail to tol_d). Checks
// per-pole solution agreement + TRUE fp64 residual, then a tol_f scan (whole-solve iters + speedup).
//
// Build per L via -DLREF={1,2,4} on a GPU. Handoff: tmp_multishift_mixed_build_claude.sh.

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
  constexpr int Nt=NTIME;

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
#include "includes/sparse_matrix_mixed_claude.h"
#include "includes/matpoly_mixed_claude.h"

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

  Gauge U(base);
  if(argc>1){
    U.read(argv[1]);                                    // real thermalized config (GaugeExt raw-double dump)
    std::cout << "# config = " << argv[1] << std::endl;
  }
  else{
    using Rng = ParallelRngExt<Base, Comp::Nt>;
    Rng rng(base);
    U.gaussian(rng, 0.3);
    std::cout << "# config = gaussian(w=0.3)" << std::endl;
  }

  const int npole = 21;
  OverlapWMass<WilsonDirac> Dov(DW, Complex(0.0,0.0), npole);
  // frozen PRODUCTION window per L (frozen_window_claude.md): L1 (0.1,13), L2 (0.06,8), L4 (0.008,5)
  // -> seed conditioning matches real running (esp. L4's near-zero tail), not the default k=0.01.
  double lmin_fr=0.1, lmax_fr=13.0;
  if(Comp::N_REFINE==2){ lmin_fr=0.06; lmax_fr=8.0; }
  else if(Comp::N_REFINE==4){ lmin_fr=0.008; lmax_fr=5.0; }
  Dov.set_lambda(lmin_fr, lmax_fr);
  Dov.update(U);

  const double scale = 1.0/(Dov.lambda_max*Dov.lambda_max);
  const Idx nnz = Dov.d_DW.len;
  const int npl = Dov.size-1;                        // number of Zolotarev poles
  std::vector<double> sigma(npl);
  for(int m=1; m<Dov.size; m++) sigma[m-1] = -Dov.k*Dov.k/Dov.cp[m];
  std::cout << "# L=" << Comp::N_REFINE << "  N=" << N << "  nnz=" << nnz << "  npole=" << npl
            << "  lambda_max=" << Dov.lambda_max << std::endl;

  // random RHS (single, shared by all poles)
  std::vector<Complex> bvec(N);
  std::mt19937_64 gen(2024ull);
  std::normal_distribution<double> nd(0.0,1.0);
  for(Idx i=0;i<N;i++) bvec[i] = Complex(nd(gen), nd(gen));

  CuC *d_b=nullptr, *d_Xref=nullptr, *d_Xmix=nullptr, *d_chk=nullptr;
  CUDA_CHECK(cudaMalloc(&d_b,    N*CD));
  CUDA_CHECK(cudaMalloc(&d_Xref, (size_t)N*npl*CD));
  CUDA_CHECK(cudaMalloc(&d_Xmix, (size_t)N*npl*CD));
  CUDA_CHECK(cudaMalloc(&d_chk,  N*CD));
  CUDA_CHECK(cudaMemcpy(d_b, reinterpret_cast<const CuC*>(bvec.data()), N*CD, H2D));

  double bnorm2_host=0.0;
  for(const Complex& z : bvec) bnorm2_host += std::norm(z);
  const double bnorm = std::sqrt(bnorm2_host);

  MixedSolver<N> ms(Dov.M_DW, Dov.M_DWH, scale, nnz);
  cublasHandle_t hb;  CUBLAS_CHECK(cublasCreate(&hb));

  const double tol_d  = Comp::TOL_INNER;   // 1e-9
  const double tol_f0 = 1.0e-5;            // default fp32->fp64 crossover
  const int    freq0  = 10;                // default reliable-update cadence

  // ---- reference fp64 multishift ----
  MatPoly Aseed;
  Aseed.push_back( cplx(scale), {&Dov.M_DW, &Dov.M_DWH} );
  reset_cg_iters();
  Aseed.solve_multishift<N>(d_Xref, d_b, sigma.data(), npl, tol_d);
  const unsigned long long iters_ref = get_cg_iters();

  // ---- mixed multishift (default tol_f/freq) ----
  ms.solve_multishift(d_Xmix, d_b, sigma.data(), npl, tol_d, tol_f0, freq0);

  std::cout << "# fp64 ref multishift : Krylov iters = " << iters_ref << std::endl;
  std::cout << "# mixed multishift (tol_f=" << tol_f0 << ", freq=" << freq0 << ", n_nested=1): fp32it=" << ms.last_f
            << " fp64it=" << ms.last_d << " RUs=" << ms.last_stage
            << " clean_lo=" << ms.last_clean_lo << " clean_hi=" << ms.last_clean_hi << std::endl;

  // ---- per-pole agreement + TRUE fp64 residual ----
  std::vector<Complex> xr(N), xm(N);
  double max_reldiff=0.0, max_res=0.0;
  for(int m=0; m<npl; m++){
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(xr.data()), d_Xref+(size_t)m*N, N*CD, D2H));
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(xm.data()), d_Xmix+(size_t)m*N, N*CD, D2H));
    double nd2=0.0, nb2=0.0;
    for(Idx i=0;i<N;i++){ nd2+=std::norm(xm[i]-xr[i]); nb2+=std::norm(xr[i]); }
    const double reld = std::sqrt(nd2)/std::max(std::sqrt(nb2),1e-300);

    ms.applyA_fp64(d_chk, d_Xmix+(size_t)m*N, sigma[m]);
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_chk, -1.0, d_b, d_chk);
    double rn=0.0;  CUBLAS_CHECK(cublasDznrm2(hb, N, d_chk, 1, &rn));
    const double tres = rn/std::max(bnorm,1e-300);

    if(reld>max_reldiff) max_reldiff=reld;
    if(tres>max_res)     max_res=tres;
  }
  std::cout << "# max per-pole reldiff (mix vs ref) = " << max_reldiff << std::endl;
  std::cout << "# max per-pole TRUE fp64 residual   = " << max_res << "   (tol_d=" << tol_d << ")" << std::endl;
  const bool ok = (max_reldiff < 1e-6) && (max_res < 5.0*tol_d);
  std::cout << "# -> " << (ok ? "PASS" : "FAIL") << std::endl;

  // ---- tol_f scan (whole multishift): iters + max true residual + speedup ----
  cudaEvent_t e0,e1;  CUDA_CHECK(cudaEventCreate(&e0));  CUDA_CHECK(cudaEventCreate(&e1));
  const int NREP=5;   // nt128 real configs are ~8x bigger than the nt16 synthetic test
  Aseed.solve_multishift<N>(d_Xref, d_b, sigma.data(), npl, tol_d);  // warm
  CUDA_CHECK(cudaDeviceSynchronize());
  CUDA_CHECK(cudaEventRecord(e0));
  for(int r=0;r<NREP;r++) Aseed.solve_multishift<N>(d_Xref, d_b, sigma.data(), npl, tol_d);
  CUDA_CHECK(cudaEventRecord(e1));  CUDA_CHECK(cudaEventSynchronize(e1));
  float ms_ref=0.0f;  CUDA_CHECK(cudaEventElapsedTime(&ms_ref, e0, e1));

  std::cout << "\n# ===== (tol_f x n_nested) scan (mixed multishift, tol_d=" << tol_d << ", freq=" << freq0 << ") =====" << std::endl;
  std::cout << "#   fp64 ref multishift: " << (ms_ref/NREP) << " ms/solve   (nn=0 = pure fp64 cleanup)" << std::endl;
  std::cout << "#   tol_f    nn  fp32it fp64it RUs cln_lo cln_hi   max_true_res    ms/solve   speedup" << std::endl;
  const double tol_f_vals[] = { 1.0e-4, 1.0e-5, 1.0e-6 };
  const int    nn_vals[]    = { 0, 1, 2 };
  for(double tf : tol_f_vals){
    for(int nn : nn_vals){
      ms.solve_multishift(d_Xmix, d_b, sigma.data(), npl, tol_d, tf, freq0, true, nn);   // warm
      CUDA_CHECK(cudaDeviceSynchronize());
      CUDA_CHECK(cudaEventRecord(e0));
      for(int r=0;r<NREP;r++) ms.solve_multishift(d_Xmix, d_b, sigma.data(), npl, tol_d, tf, freq0, true, nn);
      CUDA_CHECK(cudaEventRecord(e1));  CUDA_CHECK(cudaEventSynchronize(e1));
      float ms_mix=0.0f;  CUDA_CHECK(cudaEventElapsedTime(&ms_mix, e0, e1));

      double mres=0.0;
      for(int m=0; m<npl; m++){
        ms.applyA_fp64(d_chk, d_Xmix+(size_t)m*N, sigma[m]);
        Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_chk, -1.0, d_b, d_chk);
        double rn=0.0;  CUBLAS_CHECK(cublasDznrm2(hb, N, d_chk, 1, &rn));
        const double tres = rn/std::max(bnorm,1e-300);
        if(tres>mres) mres=tres;
      }
      std::cout << "#   " << std::setw(8) << tf
                << std::setw(4) << nn
                << std::setw(8) << ms.last_f
                << std::setw(7) << ms.last_d
                << std::setw(4) << ms.last_stage
                << std::setw(7) << ms.last_clean_lo
                << std::setw(7) << ms.last_clean_hi
                << "   " << std::setw(10) << mres
                << "  " << std::setw(10) << (ms_mix/NREP)
                << "  " << std::setw(8) << (ms_ref/ms_mix) << "x" << std::endl;
    }
  }

  std::cout << "\n# ============ " << (ok ? "PASS" : "FAIL") << " ============" << std::endl;

  CUBLAS_CHECK(cublasDestroy(hb));
  CUDA_CHECK(cudaFree(d_b)); CUDA_CHECK(cudaFree(d_Xref)); CUDA_CHECK(cudaFree(d_Xmix)); CUDA_CHECK(cudaFree(d_chk));
  return ok ? 0 : 1;
}
