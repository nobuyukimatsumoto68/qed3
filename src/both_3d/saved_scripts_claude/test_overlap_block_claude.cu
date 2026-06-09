// test_overlap_block_claude.cu
// C6c (device): validate the mrhs BLOCK overlap ops (mult/adj/DDH _deviceAsyncLaunch_ms_block)
// against the single-RHS _ms siblings, column by column, on the REAL operator. By construction
// (per-column independence + the block solve == per-column single solve bit-for-bit, C6b) the block
// column c must equal the single _ms applied to b_c BIT-FOR-BIT.
//
//   block:  Dm.{mult,adj,DDH}_deviceAsyncLaunch_ms_block<NSTACK>(d_blk, d_B)
//   ref:    for each c, Dm.{mult,adj,DDH}_deviceAsyncLaunch_ms(d_ref, d_B + c*N)
// NSTACK = 1 (reproduces _ms) / 4 / N_SITES. PASS = max|block-single| < 1e-10 for all three ops.
// Free field (U=1): the operator math is gauge-independent.
//
// Compile: handled by the both_3d Makefile. Run: ./test_overlap_block_claude.o [--mass-re x --mass-im y]

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

// max | block_column_c - single | over N complex entries
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

  // preallocate the thread-owned block scratch once, sized for the largest batch (N_SITES)

  std::cout << std::scientific << std::setprecision(3);
  bool all_pass = true;

  auto run = [&]<int NSTACK>(){
    const size_t Ns = (size_t)N*NSTACK;
    CuC *d_B,*d_blk,*d_ref;
    CUDA_CHECK(cudaMalloc(&d_B,   Ns*CD));
    CUDA_CHECK(cudaMalloc(&d_blk, Ns*CD));
    CUDA_CHECK(cudaMalloc(&d_ref, (size_t)N*CD));
    for(int c=0;c<NSTACK;c++){ FermionVector eta; eta.fill_z2_source(rng);
      CUDA_CHECK(cudaMemcpy(d_B + (size_t)c*N, reinterpret_cast<CuC*>(eta.field), N*CD, H2D)); }

    BlockedMat<N,NSTACK,Fermion> blk(Dm);   // owns block scratch (RAII), wraps Dm

    // ---- mult ----
    blk.mult(d_blk, d_B); CUDA_CHECK(cudaDeviceSynchronize());
    double w_mult=0.0;
    for(int c=0;c<NSTACK;c++){ Dm.mult_deviceAsyncLaunch_ms(d_ref, d_B + (size_t)c*N); CUDA_CHECK(cudaDeviceSynchronize());
      w_mult=std::max(w_mult, cmp_col(d_blk + (size_t)c*N, d_ref, N)); }
    // ---- adj ----
    blk.adj(d_blk, d_B); CUDA_CHECK(cudaDeviceSynchronize());
    double w_adj=0.0;
    for(int c=0;c<NSTACK;c++){ Dm.adj_deviceAsyncLaunch_ms(d_ref, d_B + (size_t)c*N); CUDA_CHECK(cudaDeviceSynchronize());
      w_adj=std::max(w_adj, cmp_col(d_blk + (size_t)c*N, d_ref, N)); }
    // ---- DDH ----
    blk.DDH(d_blk, d_B); CUDA_CHECK(cudaDeviceSynchronize());
    double w_ddh=0.0;
    for(int c=0;c<NSTACK;c++){ Dm.DDH_deviceAsyncLaunch_ms(d_ref, d_B + (size_t)c*N); CUDA_CHECK(cudaDeviceSynchronize());
      w_ddh=std::max(w_ddh, cmp_col(d_blk + (size_t)c*N, d_ref, N)); }

    const bool ok = (w_mult<1e-10)&&(w_adj<1e-10)&&(w_ddh<1e-10); all_pass &= ok;
    std::cout << "# nstack="<<std::setw(2)<<NSTACK
              << "  mult="<<w_mult<<" adj="<<w_adj<<" DDH="<<w_ddh<<"  ("<<(ok?"PASS":"FAIL")<<")\n";
    CUDA_CHECK(cudaFree(d_B)); CUDA_CHECK(cudaFree(d_blk)); CUDA_CHECK(cudaFree(d_ref));
  };

  run.operator()<1>();
  run.operator()<4>();
  run.operator()<Comp::N_SITES>();

  std::cout << "# C6c RESULT: " << (all_pass ? "ALL PASS" : "FAIL") << std::endl;
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
  return all_pass ? 0 : 1;
}
