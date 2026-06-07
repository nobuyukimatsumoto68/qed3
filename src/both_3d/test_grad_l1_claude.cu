// test_grad_l1_claude.cu
// L1 (HMC force opt, hmc_force_opt_impl_plan_claude.md): validate + BENCHMARK
// OverlapWMass::grad_deviceAsyncLaunch_l1 (X Z_m / X Y_m hoisted into precalc_grad) against the original
// grad_deviceAsyncLaunch, per gauge link. The L1 variant is the SAME math with the link-independent X
// applies precomputed once, so the per-link force value must match to ~machine.
//   ref: D.grad_deviceAsyncLaunch   (link, U, d_eta)
//   L1 : D.grad_deviceAsyncLaunch_l1(link, U, d_eta)   (reads d_XZpre/d_XYpre from precalc_grad)
// Sweeps NS_TEST timeslices x (spatial links + temporal sites). PASS = max|ref-L1| < 1e-10. Reports
// wall-time ref vs L1 over the sweep + speedup (the force gain). Free field (U=1).
// NB: compile WITHOUT -DGRAD_L1 (this test calls BOTH methods explicitly; -DGRAD_L1 would alias the ref).
//
// Reference (algorithm): overlap force = derivative of the Zolotarev sign approximation; inner shifted
// solves = B. Jegerlehner hep-lat/9612014. Compile: both_3d Makefile. Run: ./test_grad_l1_claude.o [--mass-re x --mass-im y]

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
using BaseLink = std::array<int,2>;   // matches GaugeExt::BaseLink (gauge_ext.h:5)
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
#include "overlap_wmass_claude.h"   // grad_deviceAsyncLaunch + _l1 (L1)

#include <getopt.h>

int main(int argc, char* argv[]){
  double mass_re=0.1, mass_im=0.0;   // non-zero default to exercise the (1+conj(mass)) force factor
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
            << "  size="<<Dm.size << "  mass=("<<mass_re<<","<<mass_im<<")" << std::endl;

  // pseudofermion source eta (any eta works -- ref and L1 use the same one)
  CuC* d_eta; CUDA_CHECK(cudaMalloc(&d_eta, (size_t)N*CD));
  { FermionVector eta; eta.fill_z2_source(rng);
    CUDA_CHECK(cudaMemcpy(d_eta, reinterpret_cast<CuC*>(eta.field), N*CD, H2D)); }

  // precalc once: fills d_Zs/d_Ys AND (L1) d_XZpre/d_XYpre. grad / grad_l1 both consume these.
  Dm.precalc_grad_deviceAsyncLaunch(U, d_eta);
  CUDA_CHECK(cudaDeviceSynchronize());

  const int n_sites = static_cast<int>(base.n_sites);
  const int n_links = static_cast<int>(base.links.size());
  const int NS_TEST = 8;   // timeslices to sweep (per-link cost is s-independent; 8 x (links+sites) is plenty)

  // ---- reference sweep (original grad) ----
  std::vector<double> ref;
  CUDA_CHECK(cudaDeviceSynchronize());
  auto t0=std::chrono::steady_clock::now();
  for(int s=0;s<NS_TEST;s++){
    for(int ell=0;ell<n_links;ell++){ const BaseLink lk{base.links[ell][0], base.links[ell][1]};
      ref.push_back( Dm.grad_deviceAsyncLaunch( std::pair<int,BaseLink>(s, lk), U, d_eta ) ); }
    for(int ix=0;ix<n_sites;ix++)    ref.push_back( Dm.grad_deviceAsyncLaunch( std::pair<int,Idx>(s, ix),            U, d_eta ) );
  }
  CUDA_CHECK(cudaDeviceSynchronize());
  const double t_ref=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();

  // ---- L1 sweep ----
  std::vector<double> l1;
  CUDA_CHECK(cudaDeviceSynchronize());
  t0=std::chrono::steady_clock::now();
  for(int s=0;s<NS_TEST;s++){
    for(int ell=0;ell<n_links;ell++){ const BaseLink lk{base.links[ell][0], base.links[ell][1]};
      l1.push_back( Dm.grad_deviceAsyncLaunch_l1( std::pair<int,BaseLink>(s, lk), U, d_eta ) ); }
    for(int ix=0;ix<n_sites;ix++)    l1.push_back( Dm.grad_deviceAsyncLaunch_l1( std::pair<int,Idx>(s, ix),            U, d_eta ) );
  }
  CUDA_CHECK(cudaDeviceSynchronize());
  const double t_l1=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();

  // ---- COO-build-only timing: isolate the per-link do_it cost (O(N) rows loop + 3 cudaMalloc).
  //      L2 (block poles) does NOT touch this; L4 (single-link matvec / cache) does. Tells us the ceiling.
  CUDA_CHECK(cudaDeviceSynchronize());
  t0=std::chrono::steady_clock::now();
  for(int s=0;s<NS_TEST;s++){
    for(int ell=0;ell<n_links;ell++){ const BaseLink lk{base.links[ell][0], base.links[ell][1]};
      COO<N> coo; DW.d_coo_format(coo.en, U, std::pair<int,BaseLink>(s, lk)); coo.do_it(); }
    for(int ix=0;ix<n_sites;ix++){ COO<N> coo; DW.d_coo_format(coo.en, U, std::pair<int,Idx>(s, ix)); coo.do_it(); }
  }
  CUDA_CHECK(cudaDeviceSynchronize());
  const double t_coo=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();

  // ---- compare ----
  double maxabs=0.0, maxrel=0.0;
  for(size_t i=0;i<ref.size();i++){
    const double d=std::abs(ref[i]-l1[i]);
    maxabs=std::max(maxabs,d);
    const double den=std::abs(ref[i]); if(den>1e-300) maxrel=std::max(maxrel, d/den);
  }
  const int nlink = (int)ref.size();
  const bool ok = maxabs<1e-10;
  std::cout << std::scientific << std::setprecision(3);
  std::cout << "# nlinks="<<nlink<<" (NS_TEST="<<NS_TEST<<" x ("<<n_links<<" sp + "<<n_sites<<" tp))\n";
  std::cout << "# max|ref-L1|="<<maxabs<<"  max rel="<<maxrel<<"  ("<<(ok?"PASS":"FAIL")<<")\n";
  std::cout << "# force sweep:  grad(ref)="<<t_ref<<"s  grad_l1="<<t_l1<<"s  speedup="<<(t_l1>0?t_ref/t_l1:0.0)<<"x\n";
  std::cout << "# COO-build only="<<t_coo<<"s  ("<<std::fixed<<std::setprecision(1)<<(t_l1>0?100.0*t_coo/t_l1:0.0)
            <<"% of grad_l1; remainder = pole-loop + post)"<<std::scientific<<"\n";
  std::cout << "#   => " << (t_coo>0.5*t_l1
              ? "COO BUILD DOMINATES -> L4 (single-link matvec / cache do_it) is the lever"
              : "pole loop dominates -> L2 (block poles) is the lever") << "\n";
  std::cout << "# L1 RESULT: " << (ok ? "PASS" : "FAIL") << std::endl;

  CUDA_CHECK(cudaFree(d_eta));
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
  return ok ? 0 : 1;
}
