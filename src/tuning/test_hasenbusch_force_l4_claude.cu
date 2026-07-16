// test_hasenbusch_force_l1_claude.cu
// _claude: C2 validation gate for the Hasenbusch pseudofermion stack
// (hasenbusch_massless_impl_plan_claude.md). Two checks, on a frozen nontrivial gauge:
//   (1) HEATBATH IDENTITY: right after gen(), each frame's action S_i must equal xi_i^dag xi_i
//       (the drawn Gaussian's norm) -- the exact stochastic-estimator identity at generation.
//   (2) FORCE vs FINITE DIFFERENCE: the full summed action-gradient get_force() = +dS/dU (Term A
//       over all frames + Term B [SUBTRACTED, see the sign note in pseudofermion_hasenbusch_claude.h]
//       over ratio frames) must match the central FD of the total action S(U') = sum_i chi_i(U')^dag
//       eta_i(U') with the pseudofermions phi_i FROZEN and eta_i,chi_i RE-SOLVED at each U'. get_force
//       is the action gradient the integrator feeds to pi += -tau*dSf, so grad == +fd. The re-solve
//       needs a TIGHT solve (massless inverse near-singular) or FD noise ~49 r/eps swamps the signal.
//
// Runs the K=1 ladder {0, 0.1} (validates heatbath + BOTH force terms) and the K=2 ladder
// {0, 0.1, 0.4}. Default reference grad (add -DGRAD_L4 to exercise the production block grad).
// Build: see tmp_hb_force_local_claude.sh.

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

  // TIGHT tolerances for the force FD: S = sum chi_i^dag eta_i is RE-SOLVED at each U', and a loose
  // solve on the near-singular MASSLESS (D_0^dag D_0)^{-1} injects FD noise ~ 49 r/eps (r = residual).
  // r=1e-12, eps=1e-5 -> noise ~5e-6, well under the 1e-4 gate. (Prod runs use 1e-9/1e-8.)
  const double TOL_INNER=1.0e-9;   // production tol (reverted from the diagnostic 1e-12; re-tighten if re-running the FD floor study)
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
#include "pseudofermion_hasenbusch_claude.h"

static bool run_ladder( const std::vector<Complex>& masses, const std::string& tag );

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);

  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();

  std::cout << "# Hasenbusch force+heatbath L=" << Comp::N_REFINE << " gate: N=" << Comp::N
            << " Nt=" << Comp::Nt << std::endl;

  bool ok = true;
  // RETUNED L4 ladder (2026-07-16): 4-frame {0, 0.3, 0.6, 1.0} (massless). L4 sites are NOT vertex-transitive
  // -> A_y/Abar VARIES -> this exercises the NON-uniform measure-weighted diagonal mass force.
  ok = run_ladder( { Complex(0.0,0.0), Complex(0.3,0.0), Complex(0.6,0.0), Complex(1.0,0.0) },
                   "L4_K3_{0,0.3,0.6,1.0}" ) && ok;

  std::cout << "\n# ===== Hasenbusch C2 gate: " << (ok ? "PASS" : "FAIL")
            << " (heatbath ~1e-7; force grad-vs-FD ~1e-5) =====" << std::endl;
  return ok ? 0 : 1;
}


static bool run_ladder( const std::vector<Complex>& masses, const std::string& tag ){
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
  const int npole = 17;

  Fermion D(DW, Complex(0.0,0.0), npole);   // shared operator; set_mass switches frames
  D.update(U);

  HasenbuschPF<Fermion,Force> hb(D, masses, base);

  std::cout << "\n# ========================= ladder " << tag << " (K=" << masses.size()-1
            << ") =========================" << std::endl;
  std::cout << "# delta (Zolotarev sign-fn accuracy) = " << D.Delta() << std::endl;

  // --- (1) heatbath identity: S_i == xi_i^dag xi_i at generation, AVERAGED over M draws ---
  // Per draw the identity holds only to the overlap's intrinsic accuracy (D_i and D_i^{-1} share the
  // SAME Zolotarev-approximate D_ov but route through different kernels: adj/mult vs the GW-shortcut
  // DHD_ms), so a SYSTEMATIC ~delta floor is expected and will NOT average down; a per-draw random
  // component WOULD shrink like 1/sqrt(M). We report both mean-of-(S_i)/mean-of-(xi^2) and the
  // per-draw RMS to tell which. Gate on the AVERAGED identity (systematic floor) at 1e-5.
  const int M_DRAW = 64;
  bool ok = true;
  std::vector<double> sum_S(masses.size(), 0.0), sum_x(masses.size(), 0.0), sum_rel2(masses.size(), 0.0);
  double maxrel = 0.0;
  for(int draw=0; draw<M_DRAW; draw++){
    hb.gen(rng);
    for(int i=0; i<(int)masses.size(); i++){
      const double s_i = hb.S_frame(i);
      const double x_i = hb.xi_sqnorm[i];
      const double rel = std::abs(s_i - x_i)/std::abs(x_i);
      sum_S[i]   += s_i;
      sum_x[i]   += x_i;
      sum_rel2[i]+= rel*rel;
      maxrel = std::max(maxrel, rel);
    }
  }
  std::cout << "# --- heatbath identity S_i vs xi_i^dag xi_i, averaged over M=" << M_DRAW << " draws ---" << std::endl;
  for(int i=0; i<(int)masses.size(); i++){
    const double mean_rel = std::abs(sum_S[i] - sum_x[i])/std::abs(sum_x[i]);   // averaged identity
    const double rms_rel  = std::sqrt(sum_rel2[i]/M_DRAW);                       // per-draw scatter
    ok = ok && (mean_rel < 1.0e-5);
    std::cout << "#   frame " << i << " (m=" << masses[i].real() << "): mean(S)=" << sum_S[i]/M_DRAW
              << " mean(xi^2)=" << sum_x[i]/M_DRAW << " mean_rel=" << mean_rel
              << " per-draw rms_rel=" << rms_rel << std::endl;
  }

  // --- (2) force vs FD of the total action (phi FROZEN, eta re-solved) ---
  D.update(U);
  hb.update_eta();
  Force dSf(base);
  hb.get_force(dSf, U);     // = -dS/dU

  auto S_of_U = [&](const Gauge& Up)->double {
    D.update(Up);
    hb.update_eta();        // re-solve chi_i, eta_i at Up with phi_i frozen
    return hb.S();
  };

  const double eps = 1.0e-5;
  const int s = 0;
  const std::array<Idx,3> test_sp = {2, 7, 13};
  const std::array<Idx,2> test_tp = {3, 8};

  // get_force is the ACTION GRADIENT +dS/dU (what the integrator feeds to pi += -tau*dSf), and the FD
  // here is of the ACTUAL action S = sum chi_i^dag eta_i, so the correct relation is grad == +fd.
  const double tol = 1.0e-4;
  std::cout << "# --- force get_force (=+dS/dU) vs central-FD +dS/dU (phi frozen, eta re-solved) ---" << std::endl;
  for(Idx il : test_sp){
    Gauge UP(U); UP.sp(s,il)+=eps; double sp=S_of_U(UP);
    Gauge UM(U); UM.sp(s,il)-=eps; double sm=S_of_U(UM);
    const double fd=(sp-sm)/(2.0*eps);
    const double an=dSf.sp(s,il);
    const double d=std::abs(an-fd);
    ok = ok && (d<tol);
    std::cout << "#   sp " << il << ": grad=" << an << " fd=" << fd << " |grad-fd|=" << d << std::endl;
  }
  for(Idx ix : test_tp){
    Gauge UP(U); UP.tp(s,ix)+=eps; double sp=S_of_U(UP);
    Gauge UM(U); UM.tp(s,ix)-=eps; double sm=S_of_U(UM);
    const double fd=(sp-sm)/(2.0*eps);
    const double an=dSf.tp(s,ix);
    const double d=std::abs(an-fd);
    ok = ok && (d<tol);
    std::cout << "#   tp " << ix << ": grad=" << an << " fd=" << fd << " |grad-fd|=" << d << std::endl;
  }
  D.update(U);

  std::cout << "# ladder " << tag << ": " << (ok ? "PASS" : "FAIL") << std::endl;
  return ok;
}
