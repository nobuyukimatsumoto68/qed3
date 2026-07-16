// test_hasenbusch_ml_claude.cu
// _claude: C6c validation of the multi-timescale MinimumNorm2ML (Grid-mimic nested Omelyan) for the
// Hasenbusch stack (hasenbusch_massless_impl_plan_claude.md). Three checks (massless, ladder {0,0.1,0.4},
// lambda_max FROZEN):
//   (1) 2-LEVEL REGRESSION: MinimumNorm2ML with levels [fermion(0..K), gauge] (mult 1, nsteps_inner)
//       MDsteps=nsteps MUST reproduce MinimumNorm2Hasenbusch(tmax,nsteps,nsteps_inner) -- same dH + final U.
//   (2) REVERSIBILITY of the 3-LEVEL split {light(0..K-1) supset heavy(K) supset gauge}: fwd then -pi
//       back -> recover (U0,pi0) to CG floor.
//   (3) dH ~ tau^2 for the 3-level ML (scale MDsteps -> dH/4 per doubling).
// Build: tmp_hb_ml_local_claude.sh.

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
  constexpr int Nt=4;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;
  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

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
#include "hmc_hasenbusch_claude.h"        // MinimumNorm2Hasenbusch (2-level reference)
#include "hmc_hasenbusch_ml_claude.h"     // MinimumNorm2ML (multi-timescale)

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);

  double gsq = 8.0;
  if(argc>1) gsq = atof(argv[1]);

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
  using PFptr=std::shared_ptr<HBpf>;

  Base base(Comp::N_REFINE);
  const double M5 = -1.0, at = 0.2, nu0 = 1.0;
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);

  Rng rng(base);
  Gauge U(base);
  U.gaussian(rng, 0.3);

  Fermion D(DW, Complex(0.0,0.0), 21, 0.001);
  D.update(U);
  D.freeze_lambda();                 // HMC-exactness fix: lambda_max fixed

  Action SW( gsq, at, base );

  const std::vector<Complex> masses = { Complex(0.0,0.0), Complex(0.1,0.0), Complex(0.4,0.0) };  // K=2
  const int K = (int)masses.size()-1;
  std::vector<PFptr> pfs;
  pfs.push_back( std::make_shared<HBpf>(D, masses, base) );

  const double tmax = 1.0;
  std::cout << "# ML integrator test L=" << Comp::N_REFINE << " Nt=" << Comp::Nt << " N=" << N
            << " gsq=" << gsq << " ladder{0,0.1,0.4} delta=" << D.Delta() << std::endl;

  Force pi(base);
  bool ok = true;

  // ===================== (1) 2-LEVEL REGRESSION =====================
  // ML levels [fermion(0..K) mult=1, gauge mult=nsteps_inner], MDsteps=nsteps  ==  MinimumNorm2Hasenbusch.
  {
    const int nsteps = 6, nsteps_inner = 20;
    pi.gaussian(rng);
    pfs[0]->gen(rng);                 // fixed phi
    Gauge U0(U);
    Force pi0(pi);

    // reference: MinimumNorm2Hasenbusch
    U = U0; pi = pi0; D.update(U); pfs[0]->update_eta();
    double h0r = 0.5*pi.squared_norm() + SW(U); for(auto& pf:pfs) h0r += pf->S();
    MinimumNorm2Hasenbusch integ_ref(tmax, nsteps, nsteps_inner);
    integ_ref.integrate(U, pi, &SW, &D, pfs);
    double h1r = 0.5*pi.squared_norm() + SW(U); for(auto& pf:pfs) h1r += pf->S();
    const double dH_ref = h1r - h0r;
    Gauge U_ref(U);

    // ML 2-level
    U = U0; pi = pi0; D.update(U); pfs[0]->update_eta();
    double h0m = 0.5*pi.squared_norm() + SW(U); for(auto& pf:pfs) h0m += pf->S();
    GaugeLevel<Action,Gauge,Force> gL(&SW, nsteps_inner);
    FermionGroupLevel<Fermion,PFptr,Gauge,Force> fAll(&D, pfs, 0, K, 1, base);
    std::vector<MDLevel<Gauge,Force>*> lv2 = { &fAll, &gL };
    MinimumNorm2ML<Gauge,Force> integ_ml(tmax, nsteps, lv2);
    integ_ml.integrate(U, pi);
    D.update(U); for(auto& pf:pfs) pf->update_eta();   // _claude (two-op split): action eta at final U
    double h1m = 0.5*pi.squared_norm() + SW(U); for(auto& pf:pfs) h1m += pf->S();
    const double dH_ml = h1m - h0m;

    double dU = 0.0;
    for(int t=0;t<Comp::Nt;t++){
      for(Idx il=0; il<base.n_links; il++) dU=std::max(dU,std::abs(U.sp(t,il)-U_ref.sp(t,il)));
      for(Idx ix=0; ix<base.n_sites; ix++) dU=std::max(dU,std::abs(U.tp(t,ix)-U_ref.tp(t,ix)));
    }
    const bool okR = (std::abs(dH_ml-dH_ref) < 1.0e-9) && (dU < 1.0e-9);
    ok = ok && okR;
    std::cout << "\n# ===== (1) 2-level regression: dH_ref=" << dH_ref << " dH_ml=" << dH_ml
              << " |ddH|=" << std::abs(dH_ml-dH_ref) << " max|dU|=" << dU
              << "  -> " << (okR?"PASS":"FAIL") << " (tol 1e-9) =====" << std::endl;
    U = U0; D.update(U);
  }

  // ===================== (2) 3-LEVEL REVERSIBILITY =====================
  // levels {light = ratio frames 0..K-1} supset {heavy = frame K} supset {gauge}
  {
    const int mult_light=1, mult_heavy=3, mult_gauge=8, MDsteps=4;
    pi.gaussian(rng);
    pfs[0]->gen(rng);
    Gauge U0(U);
    Force pi0(pi);

    GaugeLevel<Action,Gauge,Force> gL(&SW, mult_gauge);
    FermionGroupLevel<Fermion,PFptr,Gauge,Force> fLight(&D, pfs, 0, K-1, mult_light, base);
    FermionGroupLevel<Fermion,PFptr,Gauge,Force> fHeavy(&D, pfs, K, K,   mult_heavy, base);
    std::vector<MDLevel<Gauge,Force>*> lv3 = { &fLight, &fHeavy, &gL };
    MinimumNorm2ML<Gauge,Force> integ(tmax, MDsteps, lv3);

    integ.integrate(U, pi);       // forward
    pi *= -1.0;
    integ.integrate(U, pi);       // backward
    pi *= -1.0;
    double dU=0.0, dpi=0.0;
    for(int t=0;t<Comp::Nt;t++){
      for(Idx il=0; il<base.n_links; il++){ dU=std::max(dU,std::abs(U.sp(t,il)-U0.sp(t,il))); dpi=std::max(dpi,std::abs(pi.sp(t,il)-pi0.sp(t,il))); }
      for(Idx ix=0; ix<base.n_sites; ix++){ dU=std::max(dU,std::abs(U.tp(t,ix)-U0.tp(t,ix))); dpi=std::max(dpi,std::abs(pi.tp(t,ix)-pi0.tp(t,ix))); }
    }
    const bool okRev = (dU<1.0e-6) && (dpi<1.0e-6);
    ok = ok && okRev;
    std::cout << "\n# ===== (2) 3-level reversibility: max|dU|=" << dU << " max|dpi|=" << dpi
              << "  -> " << (okRev?"PASS":"FAIL") << " (tol 1e-6) =====" << std::endl;
    U = U0; D.update(U);
  }

  // ===================== (3) dH ~ tau^2 for the 3-level ML =====================
  {
    const int mult_light=1, mult_heavy=3, mult_gauge=8;
    const double tmax_scan=0.1;
    pi.gaussian(rng);
    pfs[0]->gen(rng);
    Gauge U0(U);
    Force pi0(pi);

    std::cout << "\n# ===== (3) dH ~ tau^2 (3-level ML, tmax=" << tmax_scan << ") =====" << std::endl;
    std::cout << "#   MDsteps     dH                  order_p" << std::endl;
    double prev=0.0; int prev_n=0;
    for(int MDsteps=2; MDsteps<=8; MDsteps*=2){
      U = U0; pi = pi0; D.update(U); pfs[0]->update_eta();
      double h0 = 0.5*pi.squared_norm() + SW(U); for(auto& pf:pfs) h0 += pf->S();
      GaugeLevel<Action,Gauge,Force> gL(&SW, mult_gauge);
      FermionGroupLevel<Fermion,PFptr,Gauge,Force> fLight(&D, pfs, 0, K-1, mult_light, base);
      FermionGroupLevel<Fermion,PFptr,Gauge,Force> fHeavy(&D, pfs, K, K,   mult_heavy, base);
      std::vector<MDLevel<Gauge,Force>*> lv3 = { &fLight, &fHeavy, &gL };
      MinimumNorm2ML<Gauge,Force> integ(tmax_scan, MDsteps, lv3);
      integ.integrate(U, pi);
      D.update(U); for(auto& pf:pfs) pf->update_eta();   // _claude (two-op split): action eta at final U
      double h1 = 0.5*pi.squared_norm() + SW(U); for(auto& pf:pfs) h1 += pf->S();
      const double dH = h1 - h0;
      const double p = (prev_n>0) ? std::log(std::abs(prev/dH))/std::log((double)MDsteps/(double)prev_n) : 0.0;
      std::cout << "#   " << MDsteps << "    " << dH << "    " << p << std::endl;
      prev = dH; prev_n = MDsteps;
    }
    U = U0; D.update(U);
  }

  std::cout << "\n# ===== ML integrator gate: " << (ok?"PASS":"FAIL")
            << " (2-level regression + 3-level reversibility) =====" << std::endl;
  return ok ? 0 : 1;
}
