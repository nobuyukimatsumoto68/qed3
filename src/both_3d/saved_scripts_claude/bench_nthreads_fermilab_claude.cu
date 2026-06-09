// bench_nthreads_claude.cu
// Benchmark NThreadsPerBlock for the two hot code paths:
//   Part 1 (HMC): precalc_grad + full N_LINKS * grad_l4 sweep (GRAD_L4 path, free field U=1)
//   Part 2 (jj):  BlockedMat::solve_sq with NSTACK=N_SITES (tp, the smaller batch used in jj source loop)
//
// Compiled three times with -DNThreadsPerBlock=256/512/1024 (gpu_header.h now #ifndef-guarded).
// Each binary prints its NThreadsPerBlock and median wall time for each part.
// Run from benchmark/ (outputs go to benchmark/<jobname>_<jobid>.out via SLURM).

#define GRAD_L4   // same flag as production hmc_fermilab_wmass_claude.cu

#include <typeinfo>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <filesystem>
#include <thread>
#include <chrono>
#include <cstdint>
#include <complex>
#include <array>
#include <vector>
#include <map>
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

#define InfoDelta

namespace Comp{
  constexpr bool is_compact=false;
  constexpr int NPARALLEL_DUPDATE=4;
  constexpr int NPARALLEL=NPARALLEL_DUPDATE;
  constexpr int NSTREAMS=NPARALLEL_DUPDATE;
  constexpr int NPARALLEL_GAUGE=16;
  constexpr int NPARALLEL_SORT=16;
  constexpr int N_REFINE=1;
  constexpr int NS=2;
  constexpr int Nt=128;
  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;
  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;
  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-8;
}

const std::string dir = "/project/affine/nmatsum/qed3/geometry/data/";
#include "/project/affine/nmatsum/qed3/geometry/geodesic.h"

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
#include "sparse_dirac.h"
#include "matpoly_claude.h"
#include "includes/overlap_wmass_claude.h"
#include "blocked_mat_claude.h"

static double wall_sec(){
  struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t);
  return t.tv_sec + 1e-9*t.tv_nsec;
}

int main(){
  std::cout << std::scientific << std::setprecision(6);
  std::cout << "# NThreadsPerBlock = " << NThreadsPerBlock << std::endl;
  std::cout << "# N = " << Comp::N << "  N_LINKS = " << Comp::N_LINKS << std::endl;

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();

  constexpr Idx N = Comp::N;
  constexpr int Nt = Comp::Nt;
  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base,DiracS2Simp>;
  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  using Rng=ParallelRngExt<Base,Nt>;
  using Fermion=OverlapWMass<WilsonDirac>;

  Base base(Comp::N_REFINE);
  const double M5=-1.0, at=0.2, nu0=1.0;
  WilsonDirac DW(base,0.0,1.0,M5,at,nu0);

  Gauge U(base);   // default = free field (U=1)
  Rng rng(base);

  const Complex mass(0.1, 0.0);
  Fermion D(DW, mass, 21);
  D.update(U);
  std::cout << "# lambda_min/max = " << D.lambda_min << " / " << D.lambda_max
            << "  delta = " << D.Delta() << std::endl;

  // pseudofermion: gaussian random eta on host, then push to device
  CuC* d_eta;
  CUDA_CHECK(cudaMalloc(&d_eta, (size_t)N*sizeof(CuC)));
  {
    std::vector<CuC> h(N);
    std::mt19937_64 eng(42);
    std::normal_distribution<double> nd(0.0,1.0);
    for(auto& c : h) c=make_cuDoubleComplex(nd(eng),nd(eng));
    CUDA_CHECK(cudaMemcpy(d_eta,h.data(),(size_t)N*sizeof(CuC),cudaMemcpyHostToDevice));
  }

  // -------------------------------------------------------
  // Part 1: HMC force  (precalc_grad + full GaugeExt::compute, GRAD_L4 path)
  // GaugeExt::compute sweeps all spatial links (pair<int,BaseLink>) and temporal
  // sites (pair<int,Idx>) -- the same per-link grad the production force does.
  // -------------------------------------------------------
  std::cout << "\n# --- Part 1: HMC force (GRAD_L4) ---" << std::endl;
  using Force=GaugeExt<Base,Nt,Comp::is_compact>;
  Force pi(base);
  const int N_WARM_HMC = 3;
  const int N_BENCH_HMC = 10;
  std::vector<double> hmc_times;

  for(int iter=0; iter < N_WARM_HMC + N_BENCH_HMC; iter++){
    CUDA_CHECK(cudaDeviceSynchronize());
    double t0 = wall_sec();
    D.precalc_grad_deviceAsyncLaunch(U, d_eta);
    pi.compute(U, d_eta, D);
    CUDA_CHECK(cudaDeviceSynchronize());
    double t1 = wall_sec();
    if(iter >= N_WARM_HMC) hmc_times.push_back(t1-t0);
  }
  std::sort(hmc_times.begin(), hmc_times.end());
  double hmc_med = hmc_times[N_BENCH_HMC/2];
  std::cout << "# HMC force: min=" << hmc_times.front() << "s  median=" << hmc_med
            << "s  max=" << hmc_times.back() << "s" << std::endl;

  // -------------------------------------------------------
  // Part 2: jj source solve  (BlockedMat::solve_sq, NSTACK=N_SITES=12)
  // -------------------------------------------------------
  std::cout << "\n# --- Part 2: jj solve_sq (NSTACK=N_SITES=" << Comp::N_SITES << ") ---" << std::endl;
  constexpr int NSTACK_TP = Comp::N_SITES;
  BlockedMat<N,NSTACK_TP,Fermion> blk_tp(D);

  // random RHS block on device
  CuC *d_B_tp, *d_X_tp;
  const size_t Ns_tp = (size_t)N*NSTACK_TP;
  CUDA_CHECK(cudaMalloc(&d_B_tp, Ns_tp*sizeof(CuC)));
  CUDA_CHECK(cudaMalloc(&d_X_tp, Ns_tp*sizeof(CuC)));
  {
    std::vector<CuC> h(Ns_tp);
    std::mt19937_64 eng(99);
    std::normal_distribution<double> nd(0.0,1.0);
    for(auto& c : h) c=make_cuDoubleComplex(nd(eng),nd(eng));
    CUDA_CHECK(cudaMemcpy(d_B_tp,h.data(),Ns_tp*sizeof(CuC),cudaMemcpyHostToDevice));
  }

  const int N_WARM_JJ = 1;
  const int N_BENCH_JJ = 3;
  std::vector<double> jj_times;

  for(int iter=0; iter < N_WARM_JJ + N_BENCH_JJ; iter++){
    CUDA_CHECK(cudaDeviceSynchronize());
    double t0 = wall_sec();
    blk_tp.solve_sq(d_X_tp, d_B_tp, Comp::TOL_OUTER);
    CUDA_CHECK(cudaDeviceSynchronize());
    double t1 = wall_sec();
    if(iter >= N_WARM_JJ) jj_times.push_back(t1-t0);
    std::cout << "# jj solve_sq iter " << iter << ": " << (t1-t0) << "s" << std::endl;
  }
  std::sort(jj_times.begin(), jj_times.end());
  double jj_med = jj_times[N_BENCH_JJ/2];
  std::cout << "# jj solve_sq: min=" << jj_times.front() << "s  median=" << jj_med
            << "s  max=" << jj_times.back() << std::endl;

  // -------------------------------------------------------
  // Summary
  // -------------------------------------------------------
  std::cout << "\n# === SUMMARY ===" << std::endl;
  std::cout << "# NThreadsPerBlock=" << NThreadsPerBlock
            << "  HMC_force_median=" << hmc_med << "s"
            << "  jj_solve_sq_median=" << jj_med << "s" << std::endl;

  CUDA_CHECK(cudaFree(d_eta));
  CUDA_CHECK(cudaFree(d_B_tp));
  CUDA_CHECK(cudaFree(d_X_tp));
  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();
  return 0;
}
