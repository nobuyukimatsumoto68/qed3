// test_overlap_multishift_claude.cu
// C3 (device): validate the multi-shift overlap variants against the originals.
//   D_ov v        : mult_deviceAsyncLaunch        vs mult_deviceAsyncLaunch_ms
//   D_ov^dag v     : adj_deviceAsyncLaunch         vs adj_deviceAsyncLaunch_ms
// Both inner paths solve to TOL_INNER, so they agree to ~TOL_INNER (not machine eps).
// PASS = relative difference below 1e-6 for both. Free field (U=1).
//
// Compile: handled by the both_3d Makefile. Run: ./test_overlap_multishift_claude.o

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

#include <getopt.h>

// max abs diff and relative diff (vs ref) between two device vectors of length N
static void compare(const char* tag, const CuC* d_a, const CuC* d_ref){
  constexpr Idx N = Comp::N;
  std::vector<Complex> a(N), r(N);
  CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(a.data()), d_a,   N*CD, D2H));
  CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(r.data()), d_ref, N*CD, D2H));
  double diff=0.0, rn=0.0;
  for(Idx i=0;i<N;i++){ diff=std::max(diff, std::abs(a[i]-r[i])); rn+=std::norm(r[i]); }
  rn=std::sqrt(rn);
  const double rel = (rn>0?diff/rn:diff);
  std::cout << "#  "<<tag<<": max|diff|="<<diff<<"  rel="<<rel
            << "  ("<<(rel<1.0e-6?"PASS":"FAIL")<<")" << std::endl;
}

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
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
  int device; CUDA_CHECK(cudaGetDeviceCount(&device));
  cudaDeviceProp dp; cudaGetDeviceProperties(&dp,0);
  std::cout << "# dev = " << dp.name << "  mass=("<<mass_re<<","<<mass_im<<")"<<std::endl;
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

  FermionVector eta; eta.fill_z2_source(rng);
  CuC *d_xi, *d_ref, *d_ms;
  CUDA_CHECK(cudaMalloc(&d_xi,  N*CD));
  CUDA_CHECK(cudaMalloc(&d_ref, N*CD));
  CUDA_CHECK(cudaMalloc(&d_ms,  N*CD));
  CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<CuC*>(eta.field), N*CD, H2D));

  auto timeit = [&](const char* tag, auto&& fn)->double{
    CUDA_CHECK(cudaDeviceSynchronize());
    auto t0=std::chrono::steady_clock::now();
    fn();
    CUDA_CHECK(cudaDeviceSynchronize());
    const double s=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    std::cout << "#  "<<tag<<" : "<<s<<" s"<<std::endl;
    return s;
  };

  std::cout << "\n# === D_ov  (mult vs mult_ms) ===" << std::endl;
  const double tm0 = timeit("mult     (pole loop)", [&]{ Dm.mult_deviceAsyncLaunch   (d_ref, d_xi); });
  const double tm1 = timeit("mult_ms  (multishift)", [&]{ Dm.mult_deviceAsyncLaunch_ms(d_ms,  d_xi); });
  compare("D_ov", d_ms, d_ref);
  std::cout << "#  speedup (loop/ms) = "<<(tm1>0?tm0/tm1:0.0)<<"x"<<std::endl;

  std::cout << "\n# === D_ov^dag  (adj vs adj_ms) ===" << std::endl;
  const double ta0 = timeit("adj      (pole loop)", [&]{ Dm.adj_deviceAsyncLaunch   (d_ref, d_xi); });
  const double ta1 = timeit("adj_ms   (multishift)", [&]{ Dm.adj_deviceAsyncLaunch_ms(d_ms,  d_xi); });
  compare("D_ov^dag", d_ms, d_ref);
  std::cout << "#  speedup (loop/ms) = "<<(ta1>0?ta0/ta1:0.0)<<"x"<<std::endl;

  CUDA_CHECK(cudaFree(d_xi)); CUDA_CHECK(cudaFree(d_ref)); CUDA_CHECK(cudaFree(d_ms));
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
  return 0;
}
