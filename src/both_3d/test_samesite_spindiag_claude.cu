// test_samesite_spindiag_claude.cu
// -----------------------------------------------------------------------------
// Probe: is the SAME-spatial-site overlap propagator block  B(t) = [D_m^{-1}]_{(n0,t),(n0,t0)}
// SPIN-DIAGONAL?  If yes (off-diagonal ~ 0), then sigma_1 B sigma_1 = sigma_2 B sigma_2 and the LOCAL
// current channels s1 == s2 exactly (the property seen in free field).  We measured s1==s2 to 1e-16 at
// U=1; this checks whether it SURVIVES a non-trivial gauge field (a Gaussian U).
//
// Method (cheap, 2 solves per gauge case): point sources e at (n0,t0,spin 0) and (n0,t0,spin 1), solve
//   col_j = D_m^{-1} e_j;  then B(t)_{i,j} = col_j(t,n0,i) = [D_m^{-1}]_{(n0,t,i),(n0,t0,j)}.
//   Report max_t |off-diag(B)| / max_t |diag(B)| for  U=1  vs  Gaussian U (a few widths).
// Expectation (see jj_corr_dilute discussion): off-diag ~ 1e-16 at U=1 (per-site isotropy), but O(signal)
//   for Gaussian U (gauge links break the isotropy) -> s1 != s2 per config for interacting.
// Reuses the jj_corr_dilute operator setup (OverlapWMass + multishift _ms + MatPoly op_DmH/op_Dmsq).
// -----------------------------------------------------------------------------

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
using MS=Eigen::Matrix2cd; using VD=Eigen::Vector2d; using VE=Eigen::Vector3d; using VC=Eigen::VectorXcd;
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
  const double TOL_OUTER=1.0e-8;
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

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(6);
  int device; CUDA_CHECK(cudaGetDeviceCount(&device)); CUDA_CHECK(cudaSetDevice(0));
  cudaDeviceProp dp; cudaGetDeviceProperties(&dp,0); std::cout<<"# dev = "<<dp.name<<std::endl;

  constexpr Idx N = Comp::N;
  constexpr int Nt = Comp::Nt;
  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  using Rng=ParallelRngExt<Base,Nt>;
  using Fermion=OverlapWMass<WilsonDirac>;

  Base base(Comp::N_REFINE);
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();
  const double r=1.0, M5=-1.0, at=0.2, nu1=1.0;
  WilsonDirac DW(base, 0.0, r, M5, at, nu1);
  Gauge U(base);                                  // default = U=1 (zero gauge phases)
  Rng rng(base, 1234);
  Fermion Dm(DW, Complex(0.0), 11);               // massless overlap D_ov (npole=11)

  // X^{-1} b = op_DmH (RHS X^dag b) + op_Dmsq (CG) -- the multishift _ms entry points.
  auto f_Dm   = std::bind(&Fermion::mult_deviceAsyncLaunch_ms, &Dm, std::placeholders::_1, std::placeholders::_2);
  auto f_DmH  = std::bind(&Fermion::adj_deviceAsyncLaunch_ms,  &Dm, std::placeholders::_1, std::placeholders::_2);
  auto f_Dmsq = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms,  &Dm, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_Dm(f_Dm), M_DmH(f_DmH), M_Dmsq(f_Dmsq);
  MatPoly op_DmH;  op_DmH.push_back(cplx(1.0), {&M_DmH});
  MatPoly op_Dmsq; op_Dmsq.push_back(cplx(1.0), {&M_Dmsq});

  const Idx  n0 = 0;     // test spatial site
  const int  t0 = 0;     // source timeslice
  FermionVector e, tmp, col0, col1;

  // gauge cases: 0 = free (U=1), then a few Gaussian-U widths.
  const std::vector<double> widths = {0.0, 0.1, 0.3, 1.0};
  for(double w : widths){
    if(w==0.0){ /* keep U=1 */ } else { U.gaussian(rng, w); }   // Gaussian gauge field (gauge_ext.h:253)
    Dm.update(U);

    e.set_pt_source(t0, n0, 0); op_DmH.from_cpu<N>(tmp.field, e.field); op_Dmsq.solve<N>(col0.field, tmp.field, Comp::TOL_OUTER);
    e.set_pt_source(t0, n0, 1); op_DmH.from_cpu<N>(tmp.field, e.field); op_Dmsq.solve<N>(col1.field, tmp.field, Comp::TOL_OUTER);

    // B(t)_{i,j} = col_j(t,n0,i) = [D_m^{-1}]_{(n0,t,i),(n0,t0,j)}
    double maxdiag=0.0, maxoff=0.0;
    for(int t=0;t<Nt;t++){
      const Complex B00=col0(t,n0,0), B11=col1(t,n0,1);   // diagonal (spin-preserving)
      const Complex B01=col1(t,n0,0), B10=col0(t,n0,1);   // off-diagonal (spin-flip)
      maxdiag=std::max(maxdiag, std::max(std::abs(B00),std::abs(B11)));
      maxoff =std::max(maxoff,  std::max(std::abs(B01),std::abs(B10)));
    }
    std::cout << "# gauge "<<(w==0.0?std::string("U=1 (free)"):("Gaussian w="+std::to_string(w)))
              << " : max|offdiag|="<<maxoff<<"  max|diag|="<<maxdiag
              << "  off/diag="<<(maxoff/maxdiag)<<std::endl;
    for(int t : {1,2,4}){
      std::cout << "#   t="<<t<<"  diag(B00="<<std::abs(col0(t,n0,0))<<",B11="<<std::abs(col1(t,n0,1))
                << ")  offdiag(B01="<<std::abs(col1(t,n0,0))<<",B10="<<std::abs(col0(t,n0,1))<<")"<<std::endl;
    }
  }

  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
  return 0;
}
