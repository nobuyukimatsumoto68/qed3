// test_multishift_gpu_claude.cu
// C2 (device): validate MatPoly::solve_multishift against independent per-pole CG
// solves on the REAL Zolotarev inner-pole system. The seed matrix is the exact
// production operator  A = (1/lambda_max^2) D_W^dag D_W  and the shifts are the
// actual poles  sigma_m = -k^2/c'_m > 0  (m = 1..size-1), as used inside
// OverlapWMass::mult_deviceAsyncLaunch (overlap_wmass_claude.h:315-322).
//
//   block:  A.solve_multishift(d_Xblock, d_xi, sigma, npole, TOL_INNER)   (one pass)
//   ref:    for each m, (A + sigma_m) X_m = d_xi via foreground solve()    (independent)
// PASS = block column m matches the reference X_m to ~TOL_INNER for every pole.
// Free field (U=1): the solver math is gauge-independent.
//
// Compile: handled by the both_3d Makefile (auto-picks every *.cu).
// Run: ./test_multishift_gpu_claude.o [--mass-re x] [--mass-im y]

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
  Dm.update(U);             // sets lambda_min/max (and, via mult path, k/cp)
  std::cout << "# D_m updated: lambda_min/max="<<Dm.lambda_min<<"/"<<Dm.lambda_max
            << "  size="<<Dm.size<<"  k="<<Dm.k << std::endl;

  // poles exactly as in mult_deviceAsyncLaunch: shift = -k^2/cp[m], m=1..size-1
  const int npole = Dm.size - 1;
  std::vector<double> sigma(npole);
  for(int m=1;m<Dm.size;m++) sigma[m-1] = -Dm.k*Dm.k/Dm.cp[m];
  { double smin=sigma[0], smax=sigma[0];
    for(int j=0;j<npole;j++){ smin=std::min(smin,sigma[j]); smax=std::max(smax,sigma[j]); }
    std::cout << "# npole="<<npole<<"  sigma in ["<<smin<<", "<<smax<<"]"<<std::endl; }

  // seed operator A = (1/lambda_max^2) D_W^dag D_W (no shift)
  const double inv_lmax2 = 1.0/(Dm.lambda_max*Dm.lambda_max);
  MatPoly Aseed; Aseed.push_back(cplx(inv_lmax2), {&Dm.M_DW, &Dm.M_DWH});

  // random source xi on device
  FermionVector eta; eta.fill_z2_source(rng);
  CuC *d_xi, *d_Xblock, *d_Xref;
  CUDA_CHECK(cudaMalloc(&d_xi,     N*CD));
  CUDA_CHECK(cudaMalloc(&d_Xblock, (size_t)N*npole*CD));
  CUDA_CHECK(cudaMalloc(&d_Xref,   N*CD));
  CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<CuC*>(eta.field), N*CD, H2D));

  // ---- block multi-shift solve (one pass) ----
  CUDA_CHECK(cudaDeviceSynchronize());
  auto tb0=std::chrono::steady_clock::now();
  Aseed.solve_multishift<N>(d_Xblock, d_xi, sigma.data(), npole, Comp::TOL_INNER);
  CUDA_CHECK(cudaDeviceSynchronize());
  const double t_block=std::chrono::duration<double>(std::chrono::steady_clock::now()-tb0).count();

  // ---- per-pole independent reference solves ----
  std::vector<Complex> Xb(N), Xr(N), Axi(N), xi_h(N);
  CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(xi_h.data()), d_xi, N*CD, D2H));
  const double xin = [&]{ double s=0; for(Idx i=0;i<N;i++) s+=std::norm(xi_h[i]); return std::sqrt(s); }();

  double worst_diff=0.0, worst_rel=0.0, worst_res=0.0;
  auto tr0=std::chrono::steady_clock::now();
  for(int j=0;j<npole;j++){
    MatPoly Am;
    Am.push_back(cplx(inv_lmax2), {&Dm.M_DW, &Dm.M_DWH});
    Am.push_back(cplx(sigma[j]), {});               // + sigma_j I
    CUDA_CHECK(cudaMemset(d_Xref, 0, N*CD));
    Am.solve<N>(d_Xref, d_xi, Comp::TOL_INNER);

    // pull block column j and reference, compare
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(Xb.data()), d_Xblock + (size_t)j*N, N*CD, D2H));
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(Xr.data()), d_Xref, N*CD, D2H));
    double diff=0.0, refn=0.0;
    for(Idx i=0;i<N;i++){ diff=std::max(diff, std::abs(Xb[i]-Xr[i])); refn+=std::norm(Xr[i]); }
    refn=std::sqrt(refn);
    const double rel = (refn>0?diff/refn:diff);

    // residual ||(A+sigma_j) Xblock_j - xi|| / ||xi||  (block solution, independent check)
    CUDA_CHECK(cudaMemcpy(d_Xref, d_Xblock + (size_t)j*N, N*CD, D2D));   // reuse d_Xref as scratch
    CuC* d_AX; CUDA_CHECK(cudaMalloc(&d_AX, N*CD));
    Am.on_gpu<N>(d_AX, d_Xref);
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(Axi.data()), d_AX, N*CD, D2H));
    CUDA_CHECK(cudaFree(d_AX));
    double rn=0.0; for(Idx i=0;i<N;i++) rn+=std::norm(Axi[i]-xi_h[i]);
    const double res=std::sqrt(rn)/(xin>0?xin:1.0);

    worst_diff=std::max(worst_diff,diff); worst_rel=std::max(worst_rel,rel); worst_res=std::max(worst_res,res);
    std::cout << "#  pole "<<std::setw(2)<<j<<"  sigma="<<sigma[j]
              << "  max|dX|="<<diff<<"  rel="<<rel<<"  resid(block)="<<res << std::endl;
  }
  const double t_ref=std::chrono::duration<double>(std::chrono::steady_clock::now()-tr0).count();

  std::cout << "# WORST  diff="<<worst_diff<<"  rel="<<worst_rel<<"  resid(block)="<<worst_res<<std::endl;
  const bool pass = (worst_rel<1.0e-6) && (worst_res<1.0e-5);
  std::cout << "# C2 solve_multishift: " << (pass?"PASS":"FAIL") << std::endl;
  std::cout << "# timing: block(one pass)="<<t_block<<" s   per-pole-loop(serial,"<<npole
            <<" solves)="<<t_ref<<" s   (note: production overlaps 4 streams)"<<std::endl;

  CUDA_CHECK(cudaFree(d_xi)); CUDA_CHECK(cudaFree(d_Xblock)); CUDA_CHECK(cudaFree(d_Xref));
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
  return pass?0:1;
}
