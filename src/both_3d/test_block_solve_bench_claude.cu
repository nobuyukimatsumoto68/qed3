// test_block_solve_bench_claude.cu -- microbenchmark: BLOCK solve (BlockedMat mrhs) vs SERIAL single solves,
// for the outer op_Dmsq = D_m^dag D_m (the dominant cost of the conn/disc drivers).  Sweeps NBLOCK in
// {1,8,16,32,64}, on a REAL massless config, for a given L (-DN_REFINE_CLI).  Reports per-RHS time, speedup,
// pool MB, and a correctness check (max |X_block - X_serial|, should be ~tol -- bit-identical per column).
// Decides: best NBLOCK, block-vs-MPS.  Handoff: run per L via tmp; read the log.  nhit-blocking deferred
// (same machinery, more tiles).  Plan: blocking_impl_plan_claude.md.  Op setup MATCHES the conn driver
// (Fermion(DW, mass=0, 11) + update(U); npole = size-1).
#include <typeinfo>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <filesystem>
#include <memory>
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

#ifndef N_REFINE_CLI
#define N_REFINE_CLI 1
#endif
  constexpr int N_REFINE=N_REFINE_CLI;
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

#include "sparse_dirac_claude.h"
#include "matpoly_claude.h"

#include "overlap_wmass_claude.h"
#include "blocked_mat_claude.h"          // BlockedMat<N,NSTACK,Op> + BlockMemPool (mrhs block solves)

#include <getopt.h>

using Base=S2Simp;
using WilsonDirac=DiracExt<Base, DiracS2Simp>;
using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;
using Rng=ParallelRngExt<Base,Comp::Nt>;
using Fermion=OverlapWMass<WilsonDirac>;

// fill a column-blocked host buffer B[c*N + i] with per-element Z2 (independent per column)
static void fill_block_z2(std::vector<Complex>& B, std::mt19937& gen){
  std::uniform_int_distribution<int> bit(0, 1);
  const double s = 1.0/std::sqrt(2.0);
  for(size_t i=0; i<B.size(); i++){
    const double re = bit(gen) ? s : -s;
    const double im = bit(gen) ? s : -s;
    B[i] = Complex(re, im);
  }
}

// one NBLOCK point: time SERIAL (NB single op_Dmsq.solve) vs BLOCK (one solve_sq_from_cpu of width NB).
template<int NB>
static void run_bench(Fermion& Dm, MatPoly& op_Dmsq, double tol, std::mt19937& gen){
  constexpr Idx N = Comp::N;
  const int npole = Dm.size - 1;
  const double pool_MB = (3.0*npole + 15.0) * (double)N * NB * 16.0 / 1.0e6;   // ~ 3 pole arrays + ~13 work + slack

  std::vector<Complex> B((size_t)N*NB);
  std::vector<Complex> Xser((size_t)N*NB);
  std::vector<Complex> Xblk((size_t)N*NB);
  fill_block_z2(B, gen);

  BlockedMat<N,NB,Fermion> blk(Dm);   // owns its pool (RAII); freed at scope exit

  // ---- warm-up block solve (allocations / caches) ----
  Xblk = B;
  blk.solve_sq_from_cpu(Xblk.data(), Xblk.data(), tol);

  // ---- SERIAL: NB single-vector op_Dmsq solves ----
  CUDA_CHECK(cudaDeviceSynchronize());
  const auto s0 = std::chrono::steady_clock::now();
  for(int c=0; c<NB; c++){
    op_Dmsq.solve<N>(Xser.data() + (size_t)c*N, B.data() + (size_t)c*N, tol);
  }
  CUDA_CHECK(cudaDeviceSynchronize());
  const double t_ser = std::chrono::duration<double>(std::chrono::steady_clock::now()-s0).count();

  // ---- BLOCK: one solve_sq of width NB (in-place; re-seed from B) ----
  Xblk = B;
  CUDA_CHECK(cudaDeviceSynchronize());
  const auto b0 = std::chrono::steady_clock::now();
  blk.solve_sq_from_cpu(Xblk.data(), Xblk.data(), tol);
  CUDA_CHECK(cudaDeviceSynchronize());
  const double t_blk = std::chrono::duration<double>(std::chrono::steady_clock::now()-b0).count();

  // ---- correctness: max column-wise |X_block - X_serial| (should be ~tol, bit-identical per column) ----
  double maxdiff = 0.0;
  for(size_t i=0; i<(size_t)N*NB; i++){
    maxdiff = std::max(maxdiff, std::abs(Xblk[i] - Xser[i]));
  }

  const double per_ser = t_ser / NB;
  const double per_blk = t_blk / NB;
  const double speedup = t_ser / t_blk;
  std::cout << "NB="<<std::setw(3)<<NB
            << "  serial="<<std::setw(9)<<t_ser<<" s"
            << "  block="<<std::setw(9)<<t_blk<<" s"
            << "  per-RHS: ser="<<per_ser<<" blk="<<per_blk<<" s"
            << "  SPEEDUP="<<std::fixed<<std::setprecision(2)<<speedup<<std::scientific
            << "  pool~"<<std::setprecision(0)<<std::fixed<<pool_MB<<" MB"<<std::scientific
            << "  maxdiff="<<std::setprecision(3)<<maxdiff
            << std::setprecision(6) << std::endl;
}

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(6);

  std::string ens_dir = "";
  int k = 0;
  double tol = Comp::TOL_OUTER;
  int maxnb = 1000;                 // only run NB <= maxnb (skip the expensive tail; default = all)
  static struct option long_opts[] = {
    {"ens-dir", required_argument, nullptr, 'e'},
    {"k",       required_argument, nullptr, 'k'},
    {"tol",     required_argument, nullptr, 't'},
    {"maxnb",   required_argument, nullptr, 'm'},
    {nullptr, 0, nullptr, 0}
  };
  int opt, idx;
  while((opt = getopt_long(argc, argv, "e:k:t:m:", long_opts, &idx)) != -1){
    switch(opt){
    case 'e': ens_dir = optarg; break;
    case 'k': k = std::stoi(optarg); break;
    case 't': tol = std::stod(optarg); break;
    case 'm': maxnb = std::stoi(optarg); break;
    default: break;
    }
  }
  const bool free_field = ens_dir.empty();
  std::cout << "# L(N_REFINE)="<<Comp::N_REFINE<<" N="<<Comp::N<<" Nt="<<Comp::Nt
            << " ens_dir="<<(free_field?std::string("<free>"):ens_dir)<<" k="<<k<<" tol="<<tol<<std::endl;

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();
  int device;
  CUDA_CHECK(cudaGetDeviceCount(&device));
  cudaDeviceProp dp;
  cudaGetDeviceProperties(&dp, 0);
  std::cout << "# dev = " << dp.name << "  ("
            << (dp.totalGlobalMem/1048576) << " MB, memBW "
            << (2.0*dp.memoryClockRate*(dp.memoryBusWidth/8)/1.0e6) << " GB/s)" << std::endl;
  CUDA_CHECK(cudaSetDevice(0));

  Base base(Comp::N_REFINE);
  const double M5 = -1.0;
  const double at = 0.2;
  WilsonDirac DW(base, 0.0, 1.0, M5, at, 1.0);
  Gauge U(base);
  Rng rng(base, 1234);

  Fermion Dm(DW, Complex(0.0, 0.0), 11);   // MASSLESS, pole count 11 (matches conn driver)
  if(!free_field){
    const std::string str_lat = ens_dir + "ckpoint_lat." + std::to_string(k);
    assert(std::filesystem::exists(str_lat) && "ckpoint_lat.k not found");
    U.read(str_lat);
  }
  Dm.update(U);
  std::cout << "# Dm.update: lambda_min/max="<<Dm.lambda_min<<"/"<<Dm.lambda_max
            << "  size="<<Dm.size<<" (npole="<<Dm.size-1<<")"<<std::endl;

  auto f_Dmsq = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms, &Dm, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_Dmsq(f_Dmsq);
  MatPoly op_Dmsq;
  op_Dmsq.push_back(cplx(1.0), {&M_Dmsq});

  std::mt19937 gen(20260714u);
  std::cout << "# --- serial vs block solve_sq (tol="<<tol<<") ---" << std::endl;
  if(maxnb >= 1)  run_bench<1> (Dm, op_Dmsq, tol, gen);
  if(maxnb >= 2)  run_bench<2> (Dm, op_Dmsq, tol, gen);
  if(maxnb >= 4)  run_bench<4> (Dm, op_Dmsq, tol, gen);
  if(maxnb >= 8)  run_bench<8> (Dm, op_Dmsq, tol, gen);
  if(maxnb >= 16) run_bench<16>(Dm, op_Dmsq, tol, gen);
  if(maxnb >= 32) run_bench<32>(Dm, op_Dmsq, tol, gen);
  if(maxnb >= 64) run_bench<64>(Dm, op_Dmsq, tol, gen);

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();
  return 0;
}
