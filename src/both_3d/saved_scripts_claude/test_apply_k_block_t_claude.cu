// test_apply_k_block_t_claude.cu
// C6f-b1 (device): validate + BENCHMARK the t-blocked sink K apply
// (ConservedCurrentBlockT::apply_k_block_t, which computes K(t,fixed) xi for all t in one pass) against
// the per-t reference ConservedCurrent::apply_k_ms. By construction the only non-bit-identical step is
// the Term B inner solve (solve_shift_block vs MatPoly::solve, ~1e-11 per column, C6f-b0); Step 1
// (shared multishift) + Term A + the COO/Wilson matvecs are identical arithmetic. So block column t
// equals apply_k_ms(t,fixed) to ~1e-11.
//   block: wrap.apply_k_block_t<NSTACK=Nt>(d_block, d_xi, U, fixed)   (fixed = site n / link lk)
//   ref:   for each t, kop.apply_k_ms(d_ref, d_xi, U, {t, fixed})
// Tested TEMPORAL (site n=0) and SPATIAL (link base.links[0]). PASS = max|block-perT| < 1e-8.
// Free field (U=1). Reports wall-time block vs Nt-serial apply_k_ms + speedup (the C6f lever).
//
// Reference (algorithm): conserved current K eq.(3.34); inner shifted solves = B. Jegerlehner,
// hep-lat/9612014. Compile: both_3d Makefile. Run: ./test_apply_k_block_t_claude.o [--mass-re x --mass-im y]

#include <typeinfo>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <complex>
#include <array>
#include <vector>
#include <map>
#include <string>
#include <random>
#include <Eigen/Dense>

using Double = double;
using Idx = std::int32_t;
using Complex = std::complex<double>;

using Link = std::array<Idx,2>;
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
  constexpr int NSTREAMS=4;
  constexpr int NPARALLEL_GAUGE=NPARALLEL_DUPDATE;
  constexpr int NPARALLEL_SORT=NPARALLEL_DUPDATE;

  constexpr int N_REFINE=1;
  constexpr int NS=2;
  constexpr int Nt=128;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;

  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-5;
}

const std::string dir = "../../geometry/data/";

#include "timer.h"
#include "s2n_simp.h"
#include "s2n_dual.h"
#include "rng.h"

#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cublas_api.h>
#include <cusolverDn.h>
using CuC = cuDoubleComplex;
#include "gpu_header.h"

#include "valence_claude.h"
#include "gauge_ext.h"
#include "action_ext.h"

#include "sparse_matrix.h"
#include "dirac_base.h"
#include "dirac_simp.h"
#include "dirac_dual.h"
#include "dirac_ext.h"

#include "sparse_dirac.h"
#include "matpoly_claude.h"
#include "overlap_wmass_claude.h"
#include "blocked_mat_claude.h"             // BlockedMat + BlockMemPool + solve_shift_block
#include "conserved_current_claude.h"       // ConservedCurrent (apply_k_ms reference)
#include "conserved_current_block_claude.h" // ConservedCurrentBlockT::apply_k_block_t (C6f-b1)

#include <getopt.h>

static double cmp_col(const CuC* d_blk_col, const CuC* d_ref, Idx n){
  std::vector<Complex> a(n), b(n);
  CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(a.data()), d_blk_col, n*CD, D2H));
  CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(b.data()), d_ref,     n*CD, D2H));
  double d=0.0; for(Idx i=0;i<n;++i) d=std::max(d, std::abs(a[i]-b[i]));
  return d;
}

int main(int argc, char* argv[]){
  double mass_re=0.0, mass_im=0.0;
  static struct option long_opts[] = {
    {"mass-re", required_argument, nullptr, 'r'},
    {"mass-im", required_argument, nullptr, 'i'},
    {nullptr,0,nullptr,0}
  };
  int opt, idx;
  while((opt=getopt_long(argc,argv,"r:i:",long_opts,&idx))!=-1){
    if(opt=='r') mass_re=std::stod(optarg);
    else if(opt=='i') mass_im=std::stod(optarg);
  }
  const Complex valence_mass(mass_re, mass_im);

  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();

  cudaDeviceProp dp; cudaGetDeviceProperties(&dp,0);
  std::cout << "# dev = " << dp.name << std::endl;
  CUDA_CHECK(cudaSetDevice(0));

  constexpr Idx N = Comp::N;

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;
  using Rng=ParallelRngExt<Base,Comp::Nt>;
  using Fermion=OverlapWMass<WilsonDirac>;

  Base base(Comp::N_REFINE);
  const double M5=-1.0, at=0.2, nu1=1.0;
  if(Comp::Nt!=1) assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu1);

  Gauge U(base);            // free field U=1
  Rng rng(base, 1234);

  Fermion Dm(DW, valence_mass, 21);
  Dm.update(U);
  std::cout << "# D_m updated: lambda_min/max="<<Dm.lambda_min<<"/"<<Dm.lambda_max
            << "  size="<<Dm.size << std::endl;

  ConservedCurrent<Fermion,Gauge> kop(Dm);
  ConservedCurrentBlockT<Fermion,Gauge,Comp::Nt> wrap(kop);

  // random source xi (any xi works -- block vs per-t use the SAME xi)
  CuC *d_xi, *d_block, *d_ref;
  CUDA_CHECK(cudaMalloc(&d_xi,    (size_t)N*CD));
  CUDA_CHECK(cudaMalloc(&d_block, (size_t)N*Comp::Nt*CD));
  CUDA_CHECK(cudaMalloc(&d_ref,   (size_t)N*CD));
  { FermionVector xi; xi.fill_z2_source(rng);
    CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<CuC*>(xi.field), N*CD, H2D)); }

  std::cout << std::scientific << std::setprecision(3);
  bool all_pass = true;

  // ---- TEMPORAL: fixed dual site n=0, el(t) = {t, 0} ----
  {
    const Idx n = 0;
    CUDA_CHECK(cudaDeviceSynchronize());
    auto t0=std::chrono::steady_clock::now();
    wrap.apply_k_block_t(d_block, d_xi, U, n);
    CUDA_CHECK(cudaDeviceSynchronize());
    const double t_blk=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();

    double worst=0.0, t_ref=0.0;
    for(int t=0;t<Comp::Nt;t++){
      auto tr=std::chrono::steady_clock::now();
      kop.apply_k_ms(d_ref, d_xi, U, std::pair<int,Idx>{t, n});
      CUDA_CHECK(cudaDeviceSynchronize());
      t_ref += std::chrono::duration<double>(std::chrono::steady_clock::now()-tr).count();
      worst = std::max(worst, cmp_col(d_block + (size_t)t*N, d_ref, N));
    }
    const bool ok = worst<1e-8; all_pass &= ok;
    std::cout << "# TEMPORAL n=0   max|block-perT|="<<worst<<"  ("<<(ok?"PASS":"FAIL")<<")"
              << "   block="<<t_blk<<"s  serial(Nt="<<Comp::Nt<<")="<<t_ref<<"s  speedup="
              << (t_blk>0? t_ref/t_blk : 0.0) <<"x\n";
  }

  // ---- SPATIAL: fixed link base.links[0], el(t) = {t, lk} ----
  {
    const BaseLink lk = base.links[0];
    CUDA_CHECK(cudaDeviceSynchronize());
    auto t0=std::chrono::steady_clock::now();
    wrap.apply_k_block_t(d_block, d_xi, U, lk);
    CUDA_CHECK(cudaDeviceSynchronize());
    const double t_blk=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();

    double worst=0.0, t_ref=0.0;
    for(int t=0;t<Comp::Nt;t++){
      auto tr=std::chrono::steady_clock::now();
      kop.apply_k_ms(d_ref, d_xi, U, std::pair<int,BaseLink>{t, lk});
      CUDA_CHECK(cudaDeviceSynchronize());
      t_ref += std::chrono::duration<double>(std::chrono::steady_clock::now()-tr).count();
      worst = std::max(worst, cmp_col(d_block + (size_t)t*N, d_ref, N));
    }
    const bool ok = worst<1e-8; all_pass &= ok;
    std::cout << "# SPATIAL lk={"<<lk[0]<<","<<lk[1]<<"}  max|block-perT|="<<worst<<"  ("<<(ok?"PASS":"FAIL")<<")"
              << "   block="<<t_blk<<"s  serial(Nt="<<Comp::Nt<<")="<<t_ref<<"s  speedup="
              << (t_blk>0? t_ref/t_blk : 0.0) <<"x\n";
  }

  // ---- K^dag TEMPORAL: fixed dual site n=0 ----
  {
    const Idx n = 0;
    CUDA_CHECK(cudaDeviceSynchronize());
    auto t0=std::chrono::steady_clock::now();
    wrap.apply_k_dag_block_t(d_block, d_xi, U, n);
    CUDA_CHECK(cudaDeviceSynchronize());
    const double t_blk=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();

    double worst=0.0, t_ref=0.0;
    for(int t=0;t<Comp::Nt;t++){
      auto tr=std::chrono::steady_clock::now();
      kop.apply_k_dag_ms(d_ref, d_xi, U, std::pair<int,Idx>{t, n});
      CUDA_CHECK(cudaDeviceSynchronize());
      t_ref += std::chrono::duration<double>(std::chrono::steady_clock::now()-tr).count();
      worst = std::max(worst, cmp_col(d_block + (size_t)t*N, d_ref, N));
    }
    const bool ok = worst<1e-8; all_pass &= ok;
    std::cout << "# K^dag TEMPORAL n=0   max|block-perT|="<<worst<<"  ("<<(ok?"PASS":"FAIL")<<")"
              << "   block="<<t_blk<<"s  serial(Nt="<<Comp::Nt<<")="<<t_ref<<"s  speedup="
              << (t_blk>0? t_ref/t_blk : 0.0) <<"x\n";
  }

  // ---- K^dag SPATIAL: fixed link base.links[0] ----
  {
    const BaseLink lk = base.links[0];
    CUDA_CHECK(cudaDeviceSynchronize());
    auto t0=std::chrono::steady_clock::now();
    wrap.apply_k_dag_block_t(d_block, d_xi, U, lk);
    CUDA_CHECK(cudaDeviceSynchronize());
    const double t_blk=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();

    double worst=0.0, t_ref=0.0;
    for(int t=0;t<Comp::Nt;t++){
      auto tr=std::chrono::steady_clock::now();
      kop.apply_k_dag_ms(d_ref, d_xi, U, std::pair<int,BaseLink>{t, lk});
      CUDA_CHECK(cudaDeviceSynchronize());
      t_ref += std::chrono::duration<double>(std::chrono::steady_clock::now()-tr).count();
      worst = std::max(worst, cmp_col(d_block + (size_t)t*N, d_ref, N));
    }
    const bool ok = worst<1e-8; all_pass &= ok;
    std::cout << "# K^dag SPATIAL lk={"<<lk[0]<<","<<lk[1]<<"}  max|block-perT|="<<worst<<"  ("<<(ok?"PASS":"FAIL")<<")"
              << "   block="<<t_blk<<"s  serial(Nt="<<Comp::Nt<<")="<<t_ref<<"s  speedup="
              << (t_blk>0? t_ref/t_blk : 0.0) <<"x\n";
  }

  CUDA_CHECK(cudaFree(d_xi)); CUDA_CHECK(cudaFree(d_block)); CUDA_CHECK(cudaFree(d_ref));
  std::cout << "# C6f-b1 RESULT: " << (all_pass ? "ALL PASS" : "FAIL") << std::endl;
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
  return all_pass ? 0 : 1;
}
