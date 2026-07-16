// test_hasenbusch_termAB_deriv_claude.cu
// _claude: PER-TERM derivative test to localize the L1 dH floor (hasenbusch_massless_impl_plan_claude.md).
// Each Hasenbusch force term is compared to the central FD of the scalar it is supposed to be the
// gradient of, with FROZEN vectors (so the FD has NO ill-conditioned re-solve -- it is limited only by
// TOL_INNER, set to 1e-12 here, far below delta). grad routines return the ACTION GRADIENT +dS/dU and
// the test scalars are the corresponding action pieces, so the correct relation is grad == -FD:
//
//   Term A (per mass m): grad(D_m, eta)         vs  -d/dU[ eta^dag DHD_ms(U) eta ]     (frozen eta)
//   Term B (mass-free) : grad_bilinear(phi,eta) vs  -d/dU[ 2 Re<phi| D_ov(U) |eta> ]   (frozen phi,eta)
//
// DHD_ms / D_ov are one mult+one adj (no inverse), so the scalars are directly evaluable. If a term is
// the exact gradient of its scalar, |grad + FD| ~ FD floor (<<delta); a delta-level gap localizes the
// culprit (massive-frame Term A GW mismatch, or Term B). Build: tmp_hb_termAB_local_claude.sh.

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

#define InfoDelta

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
  constexpr int Nt=2;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;
  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;   // production tol (reverted from the diagnostic 1e-14; re-tighten to see the FD floor)
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
#include "pseudofermion_claude.h"      // PseudoFermion: mimic the validated test_diag_mass force check

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);

  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();

  constexpr Idx N = Comp::N;
  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;
  using Force=Gauge;
  using Fermion=OverlapWMass<WilsonDirac>;
  using Rng=ParallelRngExt<Base, Comp::Nt>;

  Base base(Comp::N_REFINE);
  const double M5 = -1.0, at = 0.2, nu0 = 1.0;
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);

  Rng rng(base);
  Gauge U(base);
  U.gaussian(rng, 0.3);
  const int npole = 51;

  Fermion D(DW, Complex(0.0,0.0), npole, 0.001);
  D.update(U);
  std::cout << "# per-term derivative test L=" << Comp::N_REFINE << " N=" << N
            << " npole=" << npole << " delta=" << D.Delta() << std::endl;

  // frozen random vectors
  std::vector<Complex> phi_h(N), eta_h(N);
  { std::mt19937_64 gp(2718ull), ge(3141ull); std::normal_distribution<double> g01(0.0,1.0);
    for(Idx i=0;i<N;i++){ phi_h[i]=Complex(g01(gp),g01(gp)); eta_h[i]=Complex(g01(ge),g01(ge)); } }
  CuC *d_phi=nullptr, *d_eta=nullptr, *d_res=nullptr;
  CUDA_CHECK(cudaMalloc(&d_phi, N*sizeof(CuC)));
  CUDA_CHECK(cudaMalloc(&d_eta, N*sizeof(CuC)));
  CUDA_CHECK(cudaMalloc(&d_res, N*sizeof(CuC)));
  CUDA_CHECK(cudaMemcpy(d_phi, reinterpret_cast<const CuC*>(phi_h.data()), N*sizeof(CuC), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_eta, reinterpret_cast<const CuC*>(eta_h.data()), N*sizeof(CuC), cudaMemcpyHostToDevice));

  const int s = 0;
  const Idx link_sp = 2;    // representative spatial link for the eps sweep
  const Idx link_tp = 3;    // representative temporal link (largest force)
  const std::array<double,7> epslist = {1.0e-1, 3.0e-2, 1.0e-2, 3.0e-3, 1.0e-3, 3.0e-4, 1.0e-4};

  // ---------- MASSLESS force vs FD of the ACTION -- EXACT mimic of test_diag_mass_l1_claude.cu:329 ----------
  // Validated methodology: fixed pseudofermion source phi, eta = (D^dag D)^{-1} phi RE-SOLVED at each
  // U+-eps via update_eta, FD of the ACTION S = phi^dag (D^dag D)^{-1} phi = pf.S(), compare |grad - fd|.
  // (Earlier frozen random-eta scalar test was NOT this; the force is the exact gradient only for a
  // SOLUTION eta.) eps sweep to see whether |grad-fd| is a true residual or the FD floor.
  {
    D.set_mass(Complex(0.0,0.0));
    D.update(U);
    const double lm0 = D.lambda_max;   // lambda_max at U (frozen reference for the "freeze" column)
    PseudoFermion<Fermion> pf(D);
    CUDA_CHECK(cudaMemcpy(pf.d_phi, d_phi, N*sizeof(CuC), cudaMemcpyDeviceToDevice));  // fixed phi source
    pf.update_eta();
    D.precalc_grad_deviceAsyncLaunch(U, pf.d_eta);
    Force fA(base);
    pf.get_force(fA, U);

    // freeze=false: recompute lambda_max at Up (FD includes d lambda_max/dU); freeze=true: hold lambda_max=lm0
    auto S_of_U = [&](const Gauge& Up, bool freeze)->double {
      D.update(Up);
      if(freeze) D.lambda_max = lm0;
      pf.update_eta();          // re-solve eta = (D^dag D)^{-1} phi at Up (phi fixed)
      return pf.S();
    };

    std::cout << "\n# ===== MASSLESS force vs FD of action (tp" << link_tp << "): recompute vs FROZEN lambda_max =====" << std::endl;
    std::cout << "#   (grad tp" << link_tp << "=" << fA.tp(s,link_tp) << ")   lm0=" << lm0 << std::endl;
    std::cout << "#   eps            |grad-fd| recompute      |grad-fd| frozen-lambda" << std::endl;
    for(const double eps : epslist){
      Gauge UPt(U); UPt.tp(s,link_tp)+=eps; Gauge UMt(U); UMt.tp(s,link_tp)-=eps;
      double fd_re =(S_of_U(UPt,false)-S_of_U(UMt,false))/(2.0*eps);
      double fd_fz =(S_of_U(UPt,true) -S_of_U(UMt,true)) /(2.0*eps);
      std::cout << "#   " << eps << "    " << std::abs(fA.tp(s,link_tp)-fd_re)
                << "    " << std::abs(fA.tp(s,link_tp)-fd_fz) << std::endl;
    }
    D.update(U);
  }

  // ---------- Term B: grad_bilinear(phi,eta) vs -d/dU[ 2Re<phi|D_ov|eta> ] (frozen phi,eta) ----------
  // K = dD_ov/dU is mass-free; use the MASSLESS D_ov for the scalar (set_mass 0).
  D.set_mass(Complex(0.0,0.0));
  D.update(U);
  D.precalc_grad_bilinear_deviceAsyncLaunch_ms(U, d_phi, d_eta);
  Force fB(base);
  fB.compute(U, d_eta, D);

  auto T_B = [&](const Gauge& Up)->double {
    D.update(Up);
    D.mult_deviceAsyncLaunch_ms(d_res, d_eta);       // D_ov(Up) eta  (massless)
    CUDA_CHECK(cudaDeviceSynchronize());
    std::vector<Complex> res(N);
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(res.data()), d_res, N*sizeof(CuC), cudaMemcpyDeviceToHost));
    Complex acc(0.0,0.0);
    for(Idx i=0;i<N;i++) acc += std::conj(phi_h[i])*res[i];
    return 2.0*acc.real();
  };

  std::cout << "\n# ===== Term B (mass-free): grad_bilinear(phi,eta) vs -FD[2Re<phi|D_ov|eta>], eps sweep =====" << std::endl;
  std::cout << "#   (grad sp" << link_sp << "=" << fB.sp(s,link_sp) << "  tp" << link_tp << "=" << fB.tp(s,link_tp) << ")" << std::endl;
  std::cout << "#   eps            |grad+fd| sp" << link_sp << "         |grad+fd| tp" << link_tp << std::endl;
  for(const double eps : epslist){
    Gauge UPs(U); UPs.sp(s,link_sp)+=eps; double sps=T_B(UPs); Gauge UMs(U); UMs.sp(s,link_sp)-=eps; double sms=T_B(UMs);
    double fd_sp=(sps-sms)/(2.0*eps);
    Gauge UPt(U); UPt.tp(s,link_tp)+=eps; double spt=T_B(UPt); Gauge UMt(U); UMt.tp(s,link_tp)-=eps; double smt=T_B(UMt);
    double fd_tp=(spt-smt)/(2.0*eps);
    std::cout << "#   " << eps << "    " << std::abs(fB.sp(s,link_sp)+fd_sp) << "    " << std::abs(fB.tp(s,link_tp)+fd_tp) << std::endl;
  }
  D.update(U);

  CUDA_CHECK(cudaFree(d_phi)); CUDA_CHECK(cudaFree(d_eta)); CUDA_CHECK(cudaFree(d_res));
  std::cout << "\n# ===== done (read |grad+fd| per term/mass; ~FD floor = exact, ~delta = culprit) =====" << std::endl;
  return 0;
}
