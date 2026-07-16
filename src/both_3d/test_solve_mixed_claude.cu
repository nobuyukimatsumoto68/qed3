// test_solve_mixed_claude.cu
// _claude: Chunk-1 validation of the single-shift mixed-precision reliable-update CG
// (matpoly_mixed_claude.h, mixedprec_impl_plan_claude.md). Solves the shifted seed system
//   (scale * D_W^dag D_W + sigma) x = b ,  scale = 1/lambda_max^2 ,  sigma = -k^2 / cp[m]
// (an actual Zolotarev pole of the overlap) with BOTH the reference fp64 MatPoly::solve and the
// mixed MixedSolver, and checks: (a) the two solutions agree to ~tol, (b) the mixed solution's TRUE
// fp64 residual meets tol (reliable update reaches full precision), (c) reports iteration counts +
// wall-time speedup. Tests the seed (smallest, hardest shift) + a mid + the largest shift.
//
// Build per L via -DLREF={1,2,4} on a GPU. Handoff: tmp_solve_mixed_build_claude.sh.

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

  using Rng = ParallelRngExt<Base, Comp::Nt>;
  Rng rng(base);
  Gauge U(base);
  U.gaussian(rng, 0.3);

  const int npole = 21;
  OverlapWMass<WilsonDirac> Dov(DW, Complex(0.0,0.0), npole);
  Dov.update(U);

  const double scale = 1.0/(Dov.lambda_max*Dov.lambda_max);
  const Idx nnz = Dov.d_DW.len;
  std::cout << "# L=" << Comp::N_REFINE << "  N=" << N << "  nnz=" << nnz
            << "  lambda_max=" << Dov.lambda_max << "  k=" << Dov.k << std::endl;

  // Zolotarev pole shifts sigma_m = -k^2 / cp[m], m = 1..size-1 (as in mult_..._ms).
  std::vector<double> sigma(Dov.size-1);
  for(int m=1; m<Dov.size; m++) sigma[m-1] = -Dov.k*Dov.k/Dov.cp[m];
  std::sort(sigma.begin(), sigma.end());
  const int npl = sigma.size();
  const std::array<int,3> pick = { 0, npl/2, npl-1 };   // seed (smallest) / mid / largest

  // random RHS
  std::vector<Complex> bvec(N);
  std::mt19937_64 gen(2024ull);
  std::normal_distribution<double> nd(0.0,1.0);
  for(Idx i=0;i<N;i++) bvec[i] = Complex(nd(gen), nd(gen));

  CuC *d_b=nullptr, *d_xref=nullptr, *d_xmix=nullptr, *d_chk=nullptr;
  CUDA_CHECK(cudaMalloc(&d_b,    N*CD));
  CUDA_CHECK(cudaMalloc(&d_xref, N*CD));
  CUDA_CHECK(cudaMalloc(&d_xmix, N*CD));
  CUDA_CHECK(cudaMalloc(&d_chk,  N*CD));
  CUDA_CHECK(cudaMemcpy(d_b, reinterpret_cast<const CuC*>(bvec.data()), N*CD, H2D));

  double bnorm2_host=0.0;
  for(const Complex& z : bvec) bnorm2_host += std::norm(z);
  const double bnorm = std::sqrt(bnorm2_host);

  MixedSolver<N> ms(Dov.M_DW, Dov.M_DWH, scale, nnz);

  cublasHandle_t hb;  CUBLAS_CHECK(cublasCreate(&hb));

  const double tol_d  = Comp::TOL_INNER;   // 1e-9 physics target (fp64 finish)
  const double tol_f0 = 1.0e-5;            // default fp32 stage tol for the correctness pass
  cudaEvent_t e0,e1;  CUDA_CHECK(cudaEventCreate(&e0));  CUDA_CHECK(cudaEventCreate(&e1));

  bool all_ok = true;
  for(int idx : pick){
    const double sig = sigma[idx];
    std::cout << "\n# ===== pole " << idx << "/" << (npl-1) << "  sigma = " << sig << " =====" << std::endl;

    // ---- reference fp64 solve: (scale D_W^dag D_W + sigma) x = b ----
    MatPoly Aref;
    Aref.push_back( cplx(scale), {&Dov.M_DW, &Dov.M_DWH} );
    Aref.push_back( cplx(sig),   {} );
    reset_cg_iters();
    Aref.solve<N>(d_xref, d_b, tol_d);
    const unsigned long long iters_ref = get_cg_iters();

    // ---- mixed solve (default tol_f) ----
    ms.solve(d_xmix, d_b, sig, tol_f0, tol_d);

    // ---- agreement of the two solutions ----
    std::vector<Complex> xr(N), xm(N);
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(xr.data()), d_xref, N*CD, D2H));
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(xm.data()), d_xmix, N*CD, D2H));
    double nd2=0.0, nb2=0.0;
    for(Idx i=0;i<N;i++){ nd2+=std::norm(xm[i]-xr[i]); nb2+=std::norm(xr[i]); }
    const double sol_reldiff = std::sqrt(nd2)/std::max(std::sqrt(nb2),1e-300);

    // ---- TRUE fp64 residual of the mixed solution: ||(A+sigma)x - b|| / ||b|| ----
    ms.applyA_fp64(d_chk, d_xmix, sig);
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_chk, -1.0, d_b, d_chk); // (A+sig)x - b
    double rn=0.0;  CUBLAS_CHECK(cublasDznrm2(hb, N, d_chk, 1, &rn));
    const double true_res = rn/std::max(bnorm,1e-300);

    std::cout << "#   fp64 ref : iters=" << iters_ref << std::endl;
    std::cout << "#   mixed (tol_f=" << tol_f0 << "): stages=" << ms.last_stage << " fp32it=" << ms.last_f
              << " fp64it=" << ms.last_d << std::endl;
    std::cout << "#   sol reldiff (mix vs ref) = " << sol_reldiff << std::endl;
    std::cout << "#   mixed TRUE fp64 residual = " << true_res << "   (tol_d=" << tol_d << ")" << std::endl;

    const bool ok = (sol_reldiff < 1e-6) && (true_res < 5.0*tol_d);
    std::cout << "#   -> " << (ok ? "PASS" : "FAIL") << std::endl;
    all_ok = all_ok && ok;
  }

  // ================= tol_f SCAN on the seed pole (the count-driving, hardest shift) ==============
  {
    const double sig = sigma[pick[0]];
    std::cout << "\n# ===== tol_f scan  (seed pole sigma=" << sig << ",  tol_d=" << tol_d << ") =====" << std::endl;

    MatPoly Aref;
    Aref.push_back( cplx(scale), {&Dov.M_DW, &Dov.M_DWH} );
    Aref.push_back( cplx(sig),   {} );
    const int NREP=20;
    Aref.solve<N>(d_xref, d_b, tol_d);   // warm
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaEventRecord(e0));
    for(int r=0;r<NREP;r++) Aref.solve<N>(d_xref, d_b, tol_d);
    CUDA_CHECK(cudaEventRecord(e1));  CUDA_CHECK(cudaEventSynchronize(e1));
    float ms_ref=0.0f;  CUDA_CHECK(cudaEventElapsedTime(&ms_ref, e0, e1));
    std::cout << "#   fp64 ref: " << (ms_ref/NREP) << " ms/solve" << std::endl;
    std::cout << "#   tol_f     stages  fp32it  fp64it   true_res     ms/solve    speedup" << std::endl;

    const double tol_f_vals[] = { 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7 };
    for(double tf : tol_f_vals){
      ms.solve(d_xmix, d_b, sig, tf, tol_d);   // warm
      CUDA_CHECK(cudaDeviceSynchronize());
      CUDA_CHECK(cudaEventRecord(e0));
      for(int r=0;r<NREP;r++) ms.solve(d_xmix, d_b, sig, tf, tol_d);
      CUDA_CHECK(cudaEventRecord(e1));  CUDA_CHECK(cudaEventSynchronize(e1));
      float ms_mix=0.0f;  CUDA_CHECK(cudaEventElapsedTime(&ms_mix, e0, e1));

      ms.applyA_fp64(d_chk, d_xmix, sig);
      Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_chk, -1.0, d_b, d_chk);
      double rns=0.0;  CUBLAS_CHECK(cublasDznrm2(hb, N, d_chk, 1, &rns));
      const double tres = rns/std::max(bnorm,1e-300);

      std::cout << "#   " << std::setw(8) << tf
                << std::setw(8) << ms.last_stage
                << std::setw(8) << ms.last_f
                << std::setw(8) << ms.last_d
                << "   " << std::setw(10) << tres
                << "  " << std::setw(10) << (ms_mix/NREP)
                << "  " << std::setw(8) << (ms_ref/ms_mix) << "x" << std::endl;
    }
  }

  std::cout << "\n# ============ " << (all_ok ? "ALL PASS" : "SOME FAIL") << " ============" << std::endl;

  CUBLAS_CHECK(cublasDestroy(hb));
  CUDA_CHECK(cudaFree(d_b)); CUDA_CHECK(cudaFree(d_xref)); CUDA_CHECK(cudaFree(d_xmix));
  CUDA_CHECK(cudaFree(d_chk));
  return all_ok ? 0 : 1;
}
