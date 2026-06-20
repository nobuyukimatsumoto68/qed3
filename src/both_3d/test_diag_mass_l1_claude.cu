// test_diag_mass_l1_claude.cu
// _claude: L=1 equivalence check for the measure-weighted diagonal mass.
// At L=1 the icosahedral symmetry makes m_L = diag(m A_y/abar_s) UNIFORM, so the
// production OverlapWMass (diagonal m_L, physical m) must reduce to the frozen scalar
// obsolete::OverlapWMass with scalar mass c = m * mean_dual_area / mean_ell, applied to
// the same random vector. Tests mult, adj, DHD for m = 0, 0.1 (real), 0.1 i (complex).
// See mass_diag_l1_task_claude.md / mass_measure_factor_impl_plan_claude.md.
//
// Build: add a Makefile.fnal target (mirror the hmc_* targets); run on a GPU.
// Pass: all relative differences ~1e-13..1e-14.

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
  constexpr int N_REFINE=LREF;   // L; build -DLREF=2 / -DLREF=4 to run the FD check at L=2,4
  constexpr int NS=2;
  constexpr int Nt=2;         // small; the mass check is t-independent

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;
  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-8;
}

// const std::string dir = "/project/affine/nmatsum/qed3/geometry/data/";   // remote (fnal) path
// #include "/project/affine/nmatsum/qed3/geometry/geodesic.h"
const std::string dir = "../../geometry/data/";                              // local path
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
// (no GRAD_L1/2/4 here -> default reference grad, which IS converted. Add -DGRAD_L4 on the
//  command line to exercise the block grad_l4 once it too is folded.)
#include "includes/overlap_wmass_obsolete_claude.h"   // obsolete::OverlapWMass (scalar mass)
#include "includes/overlap_wmass_claude.h"            // OverlapWMass (diagonal m_L)
#include "includes/blocked_mat_claude.h"              // BlockedMat: diagonal m_L block-apply consistency regression
#include "pseudofermion_claude.h"                     // PseudoFermion: gen/update_eta/S/get_force (force-vs-FD check)

// relative L2 difference ||a-b|| / max(||b||, tiny)
static double reldiff(const std::vector<Complex>& a, const std::vector<Complex>& b){
  double nd=0.0, nb=0.0;
  for(size_t i=0;i<a.size();i++){ nd += std::norm(a[i]-b[i]); nb += std::norm(b[i]); }
  return std::sqrt(nd) / std::max(std::sqrt(nb), 1.0e-300);
}

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);

  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();

  constexpr Idx N = Comp::N;
  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;

  Base base(Comp::N_REFINE);
  const double M5 = -1.0, at = 0.2, nu0 = 1.0;
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);

  // _claude: NONTRIVIAL gauge via the gaussian routine (gauge_ext.h:253). The L=1 operator
  // equivalence is gauge-independent (still passes), but the force-vs-FD and obsolete-vs-production
  // FORCE checks are only meaningful on a non-cold config (a cold gauge can hide force bugs by
  // special cancellation). Fixed seed -> reproducible. Width kept modest so the overlap stays
  // well-conditioned (lambda_min/lambda_max not tiny).
  using Rng = ParallelRngExt<Base, Comp::Nt>;
  Rng rng(base);
  Gauge U(base);
  U.gaussian(rng, 0.3);
  const int npole = 17;

  std::cout << "# L=" << Comp::N_REFINE << " check: N_SITES=" << Comp::N_SITES << " Nt=" << Comp::Nt << " N=" << N
            << " mean_dual_area=" << base.mean_dual_area << " mean_ell=" << base.mean_ell << std::endl;

  // fixed random source
  std::vector<Complex> xi(N);
  std::mt19937_64 gen(12345ull);
  std::normal_distribution<double> nd(0.0,1.0);
  for(Idx i=0;i<N;i++) xi[i] = Complex(nd(gen), nd(gen));

  // apply a device operator (void(CuC*,const CuC*)) to xi, return host result
  auto applyDev = [&](auto applyFn){
    CuC *d_in=nullptr, *d_out=nullptr;
    CUDA_CHECK(cudaMalloc(&d_in,  N*sizeof(CuC)));
    CUDA_CHECK(cudaMalloc(&d_out, N*sizeof(CuC)));
    CUDA_CHECK(cudaMemcpy(d_in, reinterpret_cast<const CuC*>(xi.data()), N*sizeof(CuC), cudaMemcpyHostToDevice));
    applyFn(d_out, d_in);
    CUDA_CHECK(cudaDeviceSynchronize());
    std::vector<Complex> out(N);
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(out.data()), d_out, N*sizeof(CuC), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaFree(d_in)); CUDA_CHECK(cudaFree(d_out));
    return out;
  };

  const std::array<Complex,3> masses = { Complex(0.0,0.0), Complex(0.1,0.0), Complex(0.0,0.1) };

  bool all_ok = true;

  // ===== operator equivalence: obsolete(scalar c) vs production(diagonal m_L) -- L=1 ONLY =====
  // (at L>1 the measure varies, so there is no scalar reference; force-vs-FD below is the L>1 validator)
  if constexpr (Comp::N_REFINE == 1)
  for(const Complex m : masses){
    const Complex c = m * (base.mean_dual_area / base.mean_ell);   // uniform L=1 m_L value
    std::cout << "\n# === m = " << m << "  (scalar c = m*Abar/abar_s = " << c << ") ===" << std::endl;

    obsolete::OverlapWMass<WilsonDirac> Dold(DW, c, npole);   Dold.update(U);   // scalar reference
    OverlapWMass<WilsonDirac>           Dnew(DW, m, npole);   Dnew.update(U);   // diagonal production

    // --- non-ms path ---
    auto o_mult = applyDev([&](CuC* o,const CuC* i){ Dold.mult_deviceAsyncLaunch(o,i); });
    auto n_mult = applyDev([&](CuC* o,const CuC* i){ Dnew.mult_deviceAsyncLaunch(o,i); });
    auto o_adj  = applyDev([&](CuC* o,const CuC* i){ Dold.adj_deviceAsyncLaunch (o,i); });
    auto n_adj  = applyDev([&](CuC* o,const CuC* i){ Dnew.adj_deviceAsyncLaunch (o,i); });
    auto o_dhd  = applyDev([&](CuC* o,const CuC* i){ Dold.DHD_deviceAsyncLaunch (o,i); });
    auto n_dhd  = applyDev([&](CuC* o,const CuC* i){ Dnew.DHD_deviceAsyncLaunch (o,i); });
    // --- multi-shift path ---
    auto o_multms = applyDev([&](CuC* o,const CuC* i){ Dold.mult_deviceAsyncLaunch_ms(o,i); });
    auto n_multms = applyDev([&](CuC* o,const CuC* i){ Dnew.mult_deviceAsyncLaunch_ms(o,i); });
    auto o_adjms  = applyDev([&](CuC* o,const CuC* i){ Dold.adj_deviceAsyncLaunch_ms (o,i); });
    auto n_adjms  = applyDev([&](CuC* o,const CuC* i){ Dnew.adj_deviceAsyncLaunch_ms (o,i); });
    auto o_dhdms  = applyDev([&](CuC* o,const CuC* i){ Dold.DHD_deviceAsyncLaunch_ms (o,i); });
    auto n_dhdms  = applyDev([&](CuC* o,const CuC* i){ Dnew.DHD_deviceAsyncLaunch_ms (o,i); });

    const double r_mult   = reldiff(n_mult,   o_mult);
    const double r_adj    = reldiff(n_adj,    o_adj);
    const double r_dhd    = reldiff(n_dhd,    o_dhd);
    const double r_multms = reldiff(n_multms, o_multms);
    const double r_adjms  = reldiff(n_adjms,  o_adjms);
    const double r_dhdms  = reldiff(n_dhdms,  o_dhdms);
    std::cout << "#   non-ms:  mult = " << r_mult   << "   adj = " << r_adj   << "   DHD = " << r_dhd   << std::endl;
    std::cout << "#   _ms   :  mult = " << r_multms << "   adj = " << r_adjms << "   DHD = " << r_dhdms << std::endl;
    const double tol = 1.0e-10;
    const bool ok = (r_mult<tol && r_adj<tol && r_dhd<tol &&
                     r_multms<tol && r_adjms<tol && r_dhdms<tol);
    std::cout << "#   -> " << (ok ? "PASS" : "FAIL") << std::endl;
    all_ok = all_ok && ok;
  }
  else std::cout << "\n# L>1: operator obsolete-vs-production skipped (no scalar reference; see force-vs-FD)." << std::endl;

  // ===== BlockedMat (NSTACK=1) vs MatPoly _ms : diagonal m_L consistency (ALL L) =====
  // Regression for the blocked_mat_claude.h scalar-mass bug (mass_measure_audit_handoff_claude.md):
  // at NSTACK=1 BlockedMat reproduces the single MatPoly op bit-for-bit, so once BlockedMat also
  // applies m_L, blk.{mult,adj,DDH} must match Dnew.{mult,adj,DHD}_ms to ~1e-13. Pre-fix this FAILS
  // (scalar mass vs diagonal m_L), even at L=1 (the ~0.4% Abar/abar_s factor).
  {
    std::cout << "\n# ===== BlockedMat(NSTACK=1) vs MatPoly _ms (diagonal m_L block consistency) =====" << std::endl;
    for(const Complex m : masses){
      OverlapWMass<WilsonDirac> Dnew(DW, m, npole);  Dnew.update(U);
      BlockedMat<N, 1, OverlapWMass<WilsonDirac>> blk(Dnew);
      auto b_mult = applyDev([&](CuC* o,const CuC* i){ blk.mult(o,i); });
      auto m_mult = applyDev([&](CuC* o,const CuC* i){ Dnew.mult_deviceAsyncLaunch_ms(o,i); });
      auto b_adj  = applyDev([&](CuC* o,const CuC* i){ blk.adj (o,i); });
      auto m_adj  = applyDev([&](CuC* o,const CuC* i){ Dnew.adj_deviceAsyncLaunch_ms (o,i); });
      auto b_ddh  = applyDev([&](CuC* o,const CuC* i){ blk.DDH (o,i); });
      auto m_ddh  = applyDev([&](CuC* o,const CuC* i){ Dnew.DHD_deviceAsyncLaunch_ms (o,i); });
      const double r_mult = reldiff(b_mult, m_mult);
      const double r_adj  = reldiff(b_adj,  m_adj);
      const double r_ddh  = reldiff(b_ddh,  m_ddh);
      const double tol = 1.0e-10;
      const bool ok = (r_mult<tol && r_adj<tol && r_ddh<tol);
      std::cout << "#   m = " << m << " : mult = " << r_mult << "  adj = " << r_adj << "  DDH = " << r_ddh
                << "  -> " << (ok ? "PASS" : "FAIL") << std::endl;
      all_ok = all_ok && ok;
    }
  }

  // ===== adjointness: <u|D_m v> == <D_m^dag u|v>  (audit check 3; ALL L) =====
  // D_m^dag must be the exact adjoint of D_m; a wrong/missing conj on the per-site complex mass_coeff
  // (m_L^dag = conj(mass_coeff) M_mass) breaks it. Most sensitive at imaginary m (m=0.1i): a wrong conj
  // there gives an O(|mass_coeff|) defect, vs the ~solver-tol (~1e-9) defect when correct. Two independent
  // random vectors u,v; both non-ms and _ms paths. Site-varying m_L at L>1 is the real target.
  {
    std::cout << "\n# ===== adjointness <u|D_m v> vs <D_m^dag u|v> (m_L self-adjoint consistency) =====" << std::endl;
    std::vector<Complex> uu(N), vv(N);
    std::mt19937_64 ga(2024ull);
    std::normal_distribution<double> g01(0.0,1.0);
    for(Idx i=0;i<N;i++){ uu[i]=Complex(g01(ga),g01(ga)); vv[i]=Complex(g01(ga),g01(ga)); }

    // apply a device op (void(CuC*,const CuC*)) to an arbitrary host vector -> host result
    auto applyTo = [&](const std::vector<Complex>& in, auto applyFn){
      CuC *d_in=nullptr, *d_out=nullptr;
      CUDA_CHECK(cudaMalloc(&d_in,  N*sizeof(CuC)));
      CUDA_CHECK(cudaMalloc(&d_out, N*sizeof(CuC)));
      CUDA_CHECK(cudaMemcpy(d_in, reinterpret_cast<const CuC*>(in.data()), N*sizeof(CuC), cudaMemcpyHostToDevice));
      applyFn(d_out, d_in);
      CUDA_CHECK(cudaDeviceSynchronize());
      std::vector<Complex> out(N);
      CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(out.data()), d_out, N*sizeof(CuC), cudaMemcpyDeviceToHost));
      CUDA_CHECK(cudaFree(d_in)); CUDA_CHECK(cudaFree(d_out));
      return out;
    };
    auto braket = [&](const std::vector<Complex>& a, const std::vector<Complex>& b){   // <a|b> = sum conj(a) b
      Complex s(0.0,0.0);
      for(Idx i=0;i<N;i++) s += std::conj(a[i])*b[i];
      return s;
    };

    for(const Complex m : masses){
      OverlapWMass<WilsonDirac> Dnew(DW, m, npole);  Dnew.update(U);
      // non-ms
      auto Dv = applyTo(vv, [&](CuC* o,const CuC* i){ Dnew.mult_deviceAsyncLaunch(o,i); });
      auto Hu = applyTo(uu, [&](CuC* o,const CuC* i){ Dnew.adj_deviceAsyncLaunch (o,i); });
      const Complex lhs  = braket(uu, Dv);     // <u|D_m v>
      const Complex rhs  = braket(Hu, vv);     // <D_m^dag u|v>
      const double rel   = std::abs(lhs-rhs)/std::max(std::abs(lhs),1.0e-300);
      // _ms
      auto Dvm = applyTo(vv, [&](CuC* o,const CuC* i){ Dnew.mult_deviceAsyncLaunch_ms(o,i); });
      auto Hum = applyTo(uu, [&](CuC* o,const CuC* i){ Dnew.adj_deviceAsyncLaunch_ms (o,i); });
      const Complex lhsm = braket(uu, Dvm);
      const Complex rhsm = braket(Hum, vv);
      const double relm  = std::abs(lhsm-rhsm)/std::max(std::abs(lhsm),1.0e-300);
      const double tol = 1.0e-7;   // solver-tol limited (TOL_INNER=1e-9); a wrong conj gives O(0.1)
      const bool ok = (rel<tol && relm<tol);
      std::cout << "#   m = " << m << " : |<u|Dv>-<Hu|v>|/|.| = " << rel << " (non-ms)  " << relm << " (_ms)"
                << "  -> " << (ok ? "PASS" : "FAIL") << std::endl;
      all_ok = all_ok && ok;
    }
  }

  // ===== HMC force checks (mimics hmc_fermilab_claude.cu:357-487) =====
  // (1) PRIMARY, machine precision: at L=1 m_L is uniform = c, so the diagonal-mass grad must
  //     equal the obsolete SCALAR grad on the same fixed phi -> compare force over ALL links.
  // (2) production force vs central finite-difference of S = phi^dag (D_m^dag D_m)^-1 phi (a few
  //     links). Build WITHOUT GRAD_L1/2/4 to exercise the default reference grad.
  {
    using Force = Gauge;
    const double eps = 1.0e-5;
    const int s = 0;
    const std::array<Idx,3> test_sp = {2, 7, 13};
    const std::array<Idx,2> test_tp = {3, 8};

    // one fixed Gaussian pseudofermion source, shared by obsolete + production
    std::vector<Complex> phi_h(N);
    { std::mt19937_64 g(777ull); std::normal_distribution<double> g01(0.0,1.0);
      for(Idx i=0;i<N;i++) phi_h[i]=Complex(g01(g),g01(g)); }

    std::cout << "\n# ===== HMC force: vs finite-diff (all L)"
              << (Comp::N_REFINE==1 ? " + obsolete-vs-production (L=1)" : "") << " =====" << std::endl;
    for(const Complex m : masses){
      OverlapWMass<WilsonDirac> Dnew(DW, m, npole);  Dnew.update(U);
      PseudoFermion<OverlapWMass<WilsonDirac>> pfn(Dnew);
      CUDA_CHECK(cudaMemcpy(pfn.d_phi, reinterpret_cast<const CuC*>(phi_h.data()), N*sizeof(CuC), cudaMemcpyHostToDevice));
      pfn.update_eta();
      Force fn(base);
      Dnew.precalc_grad_deviceAsyncLaunch(U, pfn.d_eta);  pfn.get_force(fn, U);

      // (1) L=1 only: obsolete(scalar c) vs production, over ALL links (machine precision)
      if constexpr (Comp::N_REFINE == 1) {
        const Complex c = m * (base.mean_dual_area / base.mean_ell);
        obsolete::OverlapWMass<WilsonDirac> Dold(DW, c, npole);  Dold.update(U);
        PseudoFermion<obsolete::OverlapWMass<WilsonDirac>> pfo(Dold);
        CUDA_CHECK(cudaMemcpy(pfo.d_phi, reinterpret_cast<const CuC*>(phi_h.data()), N*sizeof(CuC), cudaMemcpyHostToDevice));
        pfo.update_eta();
        Force fo(base);
        Dold.precalc_grad_deviceAsyncLaunch(U, pfo.d_eta);  pfo.get_force(fo, U);
        double maxd = 0.0;
        for(int t=0;t<Comp::Nt;t++){
          for(Idx il=0; il<base.n_links; il++) maxd = std::max(maxd, std::abs(fo.sp(t,il)-fn.sp(t,il)));
          for(Idx ix=0; ix<base.n_sites; ix++) maxd = std::max(maxd, std::abs(fo.tp(t,ix)-fn.tp(t,ix)));
        }
        const bool okF = (maxd < 1.0e-9);
        std::cout << "#  m = " << m << "  max|force_obsolete - force_production| = " << maxd
                  << "  -> " << (okF ? "PASS" : "FAIL") << std::endl;
        all_ok = all_ok && okF;
      } else {
        std::cout << "#  m = " << m << std::endl;
      }

      // (2) production force vs finite-difference of the action (a few links) -- ALL L
      auto fd_sp = [&](Idx il){
        Gauge UP(U); UP.sp(s,il)+=eps; Dnew.update(UP); pfn.update_eta(); double sp=pfn.S();
        Gauge UM(U); UM.sp(s,il)-=eps; Dnew.update(UM); pfn.update_eta(); double sm=pfn.S();
        return (sp-sm)/(2.0*eps);
      };
      auto fd_tp = [&](Idx ix){
        Gauge UP(U); UP.tp(s,ix)+=eps; Dnew.update(UP); pfn.update_eta(); double sp=pfn.S();
        Gauge UM(U); UM.tp(s,ix)-=eps; Dnew.update(UM); pfn.update_eta(); double sm=pfn.S();
        return (sp-sm)/(2.0*eps);
      };
      for(Idx il : test_sp){ double an=fn.sp(s,il), fd=fd_sp(il);
        std::cout << "#    sp " << il << ": grad=" << an << " fd=" << fd << " |d|=" << std::abs(an-fd) << std::endl; }
      for(Idx ix : test_tp){ double an=fn.tp(s,ix), fd=fd_tp(ix);
        std::cout << "#    tp " << ix << ": grad=" << an << " fd=" << fd << " |d|=" << std::abs(an-fd) << std::endl; }
      Dnew.update(U);   // restore
    }
  }

  std::cout << "\n# ===== operator + L=1 force check: " << (all_ok ? "ALL PASS" : "SOME FAIL")
            << "  (FD spot-check: grad vs fd above, tol ~1e-5) =====" << std::endl;
  return all_ok ? 0 : 1;
}
// ---------------------------------------------------------------------------
// Both gradient routines: the force checks call get_force -> grad_deviceAsyncLaunch,
// which dispatches by #ifdef. All grad variants (default, l1, l2, l4) are folded to the
// diagonal measure-weighted mass m_L -- each carries (1+m_L^*) into the force dot product
// via M_mass and std::conj(mass_coeff). Build TWICE to exercise both dispatch paths; both
// builds are expected to PASS:
//   nvcc ... test_diag_mass_l1_claude.cu                 # default reference grad (folded)
//   nvcc ... -DGRAD_L4 test_diag_mass_l1_claude.cu       # block grad_l4 (folded)
// ---------------------------------------------------------------------------
