// test_block_outer_solve_claude.cu
// C6d (device): validate + BENCHMARK the block outer CG (MatPoly::solve_block_cg over the block
// operator DDH_..._ms_block = op_Dmsq) against nstack single op_Dmsq solves (MatPoly::solve over a
// LinOpWrapper of DDH_..._ms). By construction (block apply == per-column single apply bit-for-bit,
// C6c; per-column freeze) the block solution column c equals the single solve on b_c BIT-FOR-BIT.
//   block:  op.solve_block_cg<N,NSTACK>(apply_block, d_Xblk, d_B, TOL_OUTER)
//   ref:    for each c, op.solve<N>(d_ref, d_B + c*N, TOL_OUTER)   (op = MatPoly{DDH_ms wrapper})
// NSTACK = 1 (reproduces solve) / 4 / N_SITES. PASS = max|block-single| < 1e-10. Free field (U=1).
// Also reports wall-time block vs nstack-serial + speedup (the mrhs win).
//
// Compile: handled by the both_3d Makefile. Run: ./test_block_outer_solve_claude.o [--mass-re x --mass-im y]

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
#include "blocked_mat_claude.h"   // C7: BlockedMat<N,NSTACK,Op>

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
  const int npole = Dm.size - 1;

  // op_Dmsq reference: MatPoly with one LinOp = (DDH = (D+m)^dag(D+m)) via the multishift _ms apply
  // variadic tail absorbs LinOpWrapper::Async's extra stream arg (only operator() is used here)
  auto ddh_ms = [&](CuC* res, const CuC* v, auto&&...){ Dm.DDH_deviceAsyncLaunch_ms(res, v); };
  LinOpWrapper<decltype(ddh_ms)> wrap(ddh_ms);
  MatPoly op; op.push_back(cplx(1.0), {&wrap});

  std::cout << std::scientific << std::setprecision(3);
  bool all_pass = true;

  auto run = [&]<int NSTACK>(){
    const size_t Ns = (size_t)N*NSTACK;
    CuC *d_B,*d_Xblk,*d_ref;
    CUDA_CHECK(cudaMalloc(&d_B,   Ns*CD));
    CUDA_CHECK(cudaMalloc(&d_Xblk,Ns*CD));
    CUDA_CHECK(cudaMalloc(&d_ref, (size_t)N*CD));
    for(int c=0;c<NSTACK;c++){ FermionVector eta; eta.fill_z2_source(rng);
      CUDA_CHECK(cudaMemcpy(d_B + (size_t)c*N, reinterpret_cast<CuC*>(eta.field), N*CD, H2D)); }

    BlockedMat<N,NSTACK,Fermion> blk(Dm);   // owns block scratch (RAII), wraps Dm

    // ---- block solve (one call, NSTACK RHS) ----
    CUDA_CHECK(cudaDeviceSynchronize());
    auto t0=std::chrono::steady_clock::now();
    blk.solve_sq(d_Xblk, d_B, Comp::TOL_OUTER);
    CUDA_CHECK(cudaDeviceSynchronize());
    const double t_blk=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();

    // ---- NSTACK single op_Dmsq solves (reference + serial timing) ----
    double worst=0.0, t_ref=0.0;
    for(int c=0;c<NSTACK;c++){
      auto tr=std::chrono::steady_clock::now();
      op.solve<N>(d_ref, d_B + (size_t)c*N, Comp::TOL_OUTER);
      CUDA_CHECK(cudaDeviceSynchronize());
      t_ref += std::chrono::duration<double>(std::chrono::steady_clock::now()-tr).count();
      worst = std::max(worst, cmp_col(d_Xblk + (size_t)c*N, d_ref, N));
    }
    const bool ok = worst<1e-10; all_pass &= ok;
    std::cout << "# nstack="<<std::setw(2)<<NSTACK<<"  max|block-single|="<<worst<<"  ("<<(ok?"PASS":"FAIL")<<")"
              << "   block="<<t_blk<<"s  serial("<<NSTACK<<")="<<t_ref<<"s  speedup="
              << (t_blk>0? t_ref/t_blk : 0.0) <<"x\n";
    CUDA_CHECK(cudaFree(d_B)); CUDA_CHECK(cudaFree(d_Xblk)); CUDA_CHECK(cudaFree(d_ref));
  };

  run.operator()<1>();
  run.operator()<4>();
  run.operator()<Comp::N_SITES>();

  std::cout << "# C6d RESULT: " << (all_pass ? "ALL PASS" : "FAIL") << std::endl;
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
  return all_pass ? 0 : 1;
}
