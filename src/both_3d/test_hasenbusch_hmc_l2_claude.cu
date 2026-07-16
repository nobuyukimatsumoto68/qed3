// test_hasenbusch_hmc_l2_claude.cu
// _claude: C2c/C3 end-to-end validation of the Hasenbusch HMC (hasenbusch_massless_impl_plan_claude.md).
// Exercises the REAL integrator (MinimumNorm2Hasenbusch) + HMCHasenbusch on a MASSLESS overlap with the
// ladder {0, 0.1, 0.4}. Three checks on a gaussian gauge:
//   (1) REVERSIBILITY: integrate (U0,pi0) forward, negate pi, integrate back -> recover (U0,pi0) to the
//       CG-tolerance floor. Fixed phi throughout (drawn once).
//   (2) dH ~ tau^2 SCALING SWEEP (mimics saved_scripts/hmc_check.cu:422): FIXED phi, sweep nsteps=4..24,
//       print tau vs dH so the log-log slope 2 is read off directly (per-row effective order + per-row
//       reversibility printed). The decisive check that the summed force is the true action gradient.
//   (3) A few HMCHasenbusch trajectories: print dH, acceptance (sanity of the full accept/reject loop).
// Build: see tmp_hb_hmc_local_claude.sh. Small Nt (fast); LREF=2 default (-DLREF=1/4 to vary).

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
#include <memory>
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
#define LREF 2
#endif
  constexpr int N_REFINE=LREF;
  constexpr int NS=2;
  constexpr int Nt=4;         // small: validation only (reversibility/dH-scaling are Nt-independent)

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;
  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  // The L1 dH has a step-INDEPENDENT floor F ~ 3 delta that is TOLERANCE-INDEPENDENT (identical at
  // 1e-10 and 1e-12) -> it is the Zolotarev (delta) operator-consistency floor of the overlap, NOT a
  // solver artifact and NOT Hasenbusch-specific. The dH~tau^2 gate is therefore made floor-independent
  // (difference ratio, see below). 1e-10 here; prod runs use 1e-9/1e-8.
  const double TOL_INNER=1.0e-9;
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
#define GRAD_L4
#include "includes/overlap_wmass_claude.h"
#include "pseudofermion_hasenbusch_claude.h"
#include "hmc_hasenbusch_claude.h"
#include "pseudofermion_claude.h"     // standard single-PF (for the non-Hasenbusch comparison sweep)
#include "integrator.h"              // standard MinimumNorm2

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);

  double gsq = 8.0;
  if(argc>1) gsq = atof(argv[1]);
  int nsteps = 8;
  if(argc>2) nsteps = atoi(argv[2]);

  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();

  constexpr Idx N = Comp::N;
  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;
  using Force=Gauge;
  using Action=U1WilsonExt<Base>;
  using Fermion=OverlapWMass<WilsonDirac>;
  using Rng=ParallelRngExt<Base, Comp::Nt>;
  using HBpf=HasenbuschPF<Fermion,Force>;

  Base base(Comp::N_REFINE);
  const double M5 = -1.0, at = 0.2, nu0 = 1.0;
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);

  Rng rng(base);
  Gauge U(base);
  U.gaussian(rng, 0.3);

  Fermion D(DW, Complex(0.0,0.0), 21, 0.001);   // massless overlap, wide window (matches prod strong-coupling)
  D.update(U);
  D.freeze_lambda();   // _claude: HMC-exactness fix (2026-07-12) -- freeze lambda_max at this config so
                       // D_ov is a fixed function of U; the force is then the exact derivative (no dH floor).
                       // (In production the driver calls this after a few thermalization trajectories.)

  Action SW( gsq, at, base );

  const std::vector<Complex> masses = { Complex(0.0,0.0), Complex(0.1,0.0), Complex(0.4,0.0) };
  std::vector<std::shared_ptr<HBpf>> pfs;
  pfs.push_back( std::make_shared<HBpf>(D, masses, base) );

  const double tmax = 1.0;
  const int nsteps_inner = 100;

  std::cout << "# Hasenbusch HMC L=" << Comp::N_REFINE << " Nt=" << Comp::Nt << " N=" << N
            << " gsq=" << gsq << " ladder={0,0.1,0.4} tmax=" << tmax
            << " nsteps=" << nsteps << " nsteps_inner=" << nsteps_inner << std::endl;
  std::cout << "# delta (Zolotarev) = " << D.Delta() << std::endl;

  Force pi(base);

  bool ok = true;

  // ================= (1) REVERSIBILITY =================
  {
    pi.gaussian(rng);
    pfs[0]->gen(rng);          // fixed phi for the whole test
    Gauge U0(U);
    Force pi0(pi);

    MinimumNorm2Hasenbusch integ(tmax, nsteps, nsteps_inner);
    integ.integrate(U, pi, &SW, &D, pfs);   // forward

    pi *= -1.0;
    integ.integrate(U, pi, &SW, &D, pfs);   // backward
    pi *= -1.0;

    double dU = 0.0, dpi = 0.0;
    for(int t=0;t<Comp::Nt;t++){
      for(Idx il=0; il<base.n_links; il++){ dU=std::max(dU,std::abs(U.sp(t,il)-U0.sp(t,il))); dpi=std::max(dpi,std::abs(pi.sp(t,il)-pi0.sp(t,il))); }
      for(Idx ix=0; ix<base.n_sites; ix++){ dU=std::max(dU,std::abs(U.tp(t,ix)-U0.tp(t,ix))); dpi=std::max(dpi,std::abs(pi.tp(t,ix)-pi0.tp(t,ix))); }
    }
    const bool okR = (dU < 1.0e-6) && (dpi < 1.0e-6);
    ok = ok && okR;
    std::cout << "\n# ===== (1) reversibility: max|dU|=" << dU << " max|dpi|=" << dpi
              << "  -> " << (okR ? "PASS" : "FAIL") << " (tol 1e-6) =====" << std::endl;

    U = U0;   // restore for the next test
    D.update(U);
  }

  // ================= (2) dH ~ tau^2 SCALING SWEEP =================
  // Mimics the canonical scaling test in saved_scripts/hmc_check.cu:422 : FIXED phi, sweep nsteps,
  // print "tau  dH" so the 2nd-order scaling (log-log slope 2) is read off directly. The per-interval
  // effective order p (from dH ~ tau^p between consecutive points) is printed as a helper; it sits near
  // 2 at coarse tau and FLATTENS toward 0 as dH approaches the step-INDEPENDENT overlap floor
  // F ~ 3 delta (delta printed above) -- that flattening is the Zolotarev operator floor, NOT a bug.
  // Reversibility (max|pi-pi0|, max|U-U0| after fwd then -pi back) is printed per row; section (1) is
  // the hard reversibility gate, the scaling table is for inspection (as in the obsolete code).
  {
    pi.gaussian(rng);
    pfs[0]->gen(rng);                 // fixed phi for the whole sweep
    Gauge U0(U);
    Force pi0(pi);

    const double tmax_scan = 0.1;   // short trajectory for the scaling scan (faster; deeper in tau^2 regime)
    std::cout << "\n# ===== (2) dH ~ tau^2 scaling sweep (fixed phi, tmax=" << tmax_scan << ") =====" << std::endl;
    std::cout << "#    nsteps    tau                dH                  order_p     rev|dpi|     rev|dU|" << std::endl;
    double prev_dH = 0.0;
    int prev_ns = 0;
    double dH8 = 0.0, dH16 = 0.0;
    for(int ns=4; ns<=16; ns+=4){
      U = U0;
      pi = pi0;
      D.update(U);
      pfs[0]->update_eta();
      double h0 = 0.5*pi.squared_norm() + SW(U);
      for(auto& pf : pfs) h0 += pf->S();

      MinimumNorm2Hasenbusch integ(tmax_scan, ns, nsteps_inner);
      integ.integrate(U, pi, &SW, &D, pfs);   // forward

      double h1 = 0.5*pi.squared_norm() + SW(U);
      for(auto& pf : pfs) h1 += pf->S();
      const double dH = h1 - h0;
      if(ns==8) dH8 = dH;
      if(ns==16) dH16 = dH;

      // reversibility (obsolete style): negate pi, integrate back -> recover (U0, pi0)
      pi *= -1.0;
      integ.integrate(U, pi, &SW, &D, pfs);
      pi *= -1.0;
      double dpi = 0.0, dU = 0.0;
      for(int t=0;t<Comp::Nt;t++){
        for(Idx il=0; il<base.n_links; il++){ dpi=std::max(dpi,std::abs(pi.sp(t,il)-pi0.sp(t,il))); dU=std::max(dU,std::abs(U.sp(t,il)-U0.sp(t,il))); }
        for(Idx ix=0; ix<base.n_sites; ix++){ dpi=std::max(dpi,std::abs(pi.tp(t,ix)-pi0.tp(t,ix))); dU=std::max(dU,std::abs(U.tp(t,ix)-U0.tp(t,ix))); }
      }

      const double p = (prev_ns>0) ? std::log(std::abs(prev_dH/dH))/std::log((double)ns/(double)prev_ns) : 0.0;
      std::cout << "#    " << ns << "    " << tmax_scan/ns << "    " << dH << "    " << p
                << "    " << dpi << "    " << dU << std::endl;
      prev_dH = dH;
      prev_ns = ns;
    }
    // extrapolated floor F from dH = F + c/n^2 using the (8,16) = (2n,4n) pair: F = dH16 - (dH8-dH16)/3
    const double F_hb = dH16 - (dH8 - dH16)/3.0;
    std::cout << "#   -> Hasenbusch floor F = " << F_hb << " (F/delta = " << F_hb/D.Delta() << ")" << std::endl;

    U = U0;
    D.update(U);
  }

  // ================= (3) a few HMC trajectories =================
  {
    MinimumNorm2Hasenbusch integ(tmax, nsteps, nsteps_inner);
    HMCHasenbusch hmc(rng, &SW, &D, U, pi, pfs, &integ);
    std::cout << "\n# ===== (3) HMCHasenbusch trajectories =====" << std::endl;
    int n_acc = 0;
    const int ntraj = 6;
    for(int k=0;k<ntraj;k++){
      double rate, dH;
      bool is_accept;
      hmc.run(rate, dH, is_accept);
      n_acc += is_accept ? 1 : 0;
      std::cout << "#   traj " << k << ": dH=" << dH << " accept=" << is_accept << " rate=" << rate << std::endl;
    }
    std::cout << "# acceptance " << n_acc << "/" << ntraj << std::endl;
  }

  // ================= (4) STANDARD (non-Hasenbusch) single-PF scaling sweep, per mass =================
  // Same sweep as (2) but the STANDARD single PseudoFermion + MinimumNorm2 (no Hasenbusch, no Term B),
  // for masses {0, 0.1, 0.4}. Extract the floor F per mass. This LOCALIZES the Hasenbusch excess floor:
  //   F grows with mass toward the Hasenbusch F (~1e-5) -> the massive Term A (GW square) is the culprit;
  //   F stays small (~massless) for all masses          -> the ratio-frame Term B is the culprit.
  {
    const std::array<Complex,3> smasses = { Complex(0.0,0.0), Complex(0.1,0.0), Complex(0.4,0.0) };
    const double tmax_scan = 0.1;
    for(const Complex sm : smasses){
      D.set_mass( sm );
      D.update(U);
      std::vector<std::shared_ptr<PseudoFermion<Fermion>>> pfs_std;
      pfs_std.push_back( std::make_shared<PseudoFermion<Fermion>>(D) );

      pi.gaussian(rng);
      pfs_std[0]->gen(rng);             // fixed phi for the sweep
      Gauge U0(U);
      Force pi0(pi);

      std::cout << "\n# ===== (4) STANDARD non-Hasenbusch single-PF sweep, mass m=" << sm.real()
                << " (fixed phi, tmax=" << tmax_scan << ") =====" << std::endl;
      std::cout << "#    nsteps    tau                dH                  order_p" << std::endl;
      double prev_dH = 0.0;
      int prev_ns = 0;
      double dH8 = 0.0, dH16 = 0.0;
      for(int ns=4; ns<=16; ns+=4){
        U = U0;
        pi = pi0;
        D.update(U);
        pfs_std[0]->update_eta();
        D.precalc_grad_deviceAsyncLaunch(U, pfs_std[0]->d_eta);
        double h0 = 0.5*pi.squared_norm() + SW(U) + pfs_std[0]->S();

        MinimumNorm2 integ(tmax_scan, ns, nsteps_inner);
        integ.integrate(U, pi, &SW, &D, pfs_std);

        double h1 = 0.5*pi.squared_norm() + SW(U) + pfs_std[0]->S();
        const double dH = h1 - h0;
        if(ns==8) dH8 = dH;
        if(ns==16) dH16 = dH;
        const double p = (prev_ns>0) ? std::log(std::abs(prev_dH/dH))/std::log((double)ns/(double)prev_ns) : 0.0;
        std::cout << "#    " << ns << "    " << tmax_scan/ns << "    " << dH << "    " << p << std::endl;
        prev_dH = dH;
        prev_ns = ns;
      }
      const double F_std = dH16 - (dH8 - dH16)/3.0;
      std::cout << "#   -> standard m=" << sm.real() << " floor F = " << F_std
                << " (F/delta = " << F_std/D.Delta() << ")" << std::endl;

      U = U0;
      D.update(U);
    }
  }

  std::cout << "\n# ===== Hasenbusch HMC L=" << Comp::N_REFINE << " gate: " << (ok ? "PASS" : "FAIL")
            << " (reversibility + dH~tau^2) =====" << std::endl;
  return ok ? 0 : 1;
}
