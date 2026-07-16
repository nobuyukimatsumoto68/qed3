// test_hasenbusch_tune_claude.cu
// _claude: C6d nsteps TUNING for the FROZEN Hasenbusch ladder (hasenbusch_ladder(L)) at L=1/2/4
// (build -DLREF; hasenbusch_massless_impl_plan_claude.md). Reports the Osborn-style cost
//
//   Cost = C_S / ( tau^2 <P_acc> ),   C_S = N_CG/traj,   <P_acc> = min(1, exp(-dH)),
//
// so the cheapest per-stage step-count assignment at fixed trajectory length tau is the smallest Cost.
// Ref: Osborn et al., PoS LATTICE2024 052 "Automated tuning for HMC mass ratios"; ladder = M. Hasenbusch
// hep-lat/0107019.
//
// PROBE (SETUP COMPARISON): one config, shared momentum, FROZEN Zolotarev window (set_lambda). Compare a
// hand-picked list of (ladder, steps) setups. Per setup: fresh heatbath (phi mass-dependent), per-frame
// FORCE NORMS (L2/Linf), one MD integration -> dH / N_CG / Cost (+ per-stage CG). tau = hasenbusch_tau( Comp::N_REFINE ).
//
// Usage: ./bin [gsq=16] [Nf=2] [N_CFG=1] [N_DRAW=1]   (config dir Nf<Nf>_gsq<gsq>...nt128L1/)
// Build: tmp_hb_tune_l1_local_claude.sh (v2 -> test_hasenbusch_tune_l1_v2_claude.log).

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
#include <string>
#include <sstream>
#include <cmath>
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
  constexpr int Nt=128;         // production L1 temporal extent (matches the nt128L1 config files)

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;
  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;   // production tols
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
#include "frozen_window_claude.h"
#include "hasenbusch_ladder_claude.h"
#include "pseudofermion_hasenbusch_claude.h"
#include "hmc_hasenbusch_ml_claude.h"

// ---- one tuning setup: a mass ladder + a PER-STAGE step-count assignment ----
// Each frame i is its own timescale. fsteps[i] = absolute number of MD steps for frame i, ordered
// OUTER (i=0 = massless, most expensive solve, fewest steps) -> INNER (i=K = heaviest, cheapest, most
// steps). Nested-Omelyan mapping: MDsteps = fsteps[0], level-l multiplier = fsteps[l]/fsteps[l-1]
// (requires fsteps[l-1] | fsteps[l]). Gauge (innermost) sub-cycles the heaviest frame by mult_gauge.
struct Setup {
  std::string name;
  std::vector<Complex> masses;                 // {0, m_1, ..., m_K}
  std::vector<int> fsteps;                      // per-frame step counts, outer(massless)->inner(heavy); size K+1
  int mult_gauge;                               // gauge substeps per innermost-frame step
};

// L-infinity (max |component|) force norm over the whole lattice (all links + temporal sites, all t).
// Complements Force::squared_norm() (L2): L2 ~ typical force (sets timescale ratios n_i ~ sqrt||F||_2);
// Linf ~ worst link (bounds step-size stability, the max-force spikes). Templated on the driver types.
template<typename ForceT, typename BaseT>
static double force_inf_norm( const ForceT& f, const BaseT& base ){
  double mx = 0.0;
  for(int t=0; t<Comp::Nt; t++){
    for(Idx il=0; il<base.n_links; il++) mx = std::max( mx, std::abs(f.sp(t,il)) );
    for(Idx ix=0; ix<base.n_sites; ix++) mx = std::max( mx, std::abs(f.tp(t,ix)) );
  }
  return mx;
}

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(6);

  double gsq = 16.0;
  int    Nf  = 2;
  int    N_CFG  = 1;    // first try (per NM): one config, one momentum
  int    N_DRAW = 1;
  if(argc>1) gsq   = atof(argv[1]);
  if(argc>2) Nf    = atoi(argv[2]);
  if(argc>3) N_CFG = atoi(argv[3]);
  if(argc>4) N_DRAW= atoi(argv[4]);
  assert(Nf%2==0);

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

  Fermion D(DW, Complex(0.0,0.0), 21, 0.001);   // massless overlap; window overridden below
  double lmin_frozen, lmax_frozen;
  frozen_window( Comp::N_REFINE, lmin_frozen, lmax_frozen );
  D.set_lambda( lmin_frozen, lmax_frozen );      // FIXED Zolotarev window (frozen_window_claude.h), n=21
  std::cout << "# frozen window L=" << Comp::N_REFINE << ": lambda_min=" << lmin_frozen
            << " lambda_max=" << lmax_frozen << " k=" << (lmin_frozen/lmax_frozen)
            << " Delta=" << D.Delta() << std::endl;
  Action  SW( gsq, at, base );

  const double tmax = hasenbusch_tau( Comp::N_REFINE );   // = 1.2 (matches production/validation)

  // ----- locate the ensemble config dir and its last N_CFG configs -----
  // _claude: CFGDIR env overrides (point at any existing L4 config dir). Default = the current-naming dir
  // WITH mRe/mIm + the OLD-ladder _hb tag (that is where the existing L4 configs live; any thermalized L4
  // gauge is a valid probe for the per-frame force norms). Trailing '/' required.
  std::string cfgdir;
  const char* cfgdir_env = std::getenv("CFGDIR");
  if( cfgdir_env != nullptr ){
    cfgdir = std::string(cfgdir_env);
  } else {
    cfgdir = "Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)
           + "nu0"+std::to_string(nu0)+"mRe0.000000mIm0.000000nt"+std::to_string(Comp::Nt)
           + "L"+std::to_string(Comp::N_REFINE)+"_hb0.400000-1.000000/";
  }
  std::vector<int> ks;
  for(int k=1; k<=100000; k++){
    if(std::filesystem::exists(cfgdir+"ckpoint_lat."+std::to_string(k))) ks.push_back(k);
  }
  if((int)ks.size() < N_CFG){
    std::cout << "# ERROR: found only " << ks.size() << " configs in " << cfgdir
              << " (need " << N_CFG << ")" << std::endl;
    return 1;
  }
  std::vector<int> cfgk( ks.end()-N_CFG, ks.end() );   // the last N_CFG

  std::cout << "# ================ Hasenbusch L=" << Comp::N_REFINE << " massless tuning ================"
            << std::endl;
  std::cout << "# gsq=" << gsq << " Nf=" << Nf << " Nt=" << Comp::Nt << " N=" << N
            << " tau=" << tmax << " delta=" << D.Delta() << std::endl;
  std::cout << "# ensemble dir: " << cfgdir << std::endl;
  std::cout << "# probe: config " << cfgk.back() << " (last of {";
  for(int k : cfgk) std::cout << k << " ";
  std::cout << "}), one momentum; MDsteps swept per ladder  (Cost = N_CG/(tau^2 P_acc))" << std::endl;
  (void)N_DRAW;   // reserved (first try = 1 draw); kept as a CLI knob for later averaging

  // ----- SETUP COMPARISON (hand-picked ladder+steps): per setup, fresh heatbath (phi mass-dependent),
  // per-frame force norms (L2/Linf), one MD integration -> dH / N_CG / Cost. Config + momentum shared. -----
  // OLD (3-frame, gauge x20) vs NEW (5-frame, gauge x40) L4 ladder -- per-frame force norms show the split
  // balance. (Force norms depend only on the ladder masses, not steps; steps + MG affect the MD/dH/Cost part.)
  // BASELINE 3-frame MG20 vs the "safer" 3-frame + very fine gauge (MG100), motivated by Nf=6 (3 pf stacks
  // -> fermion forces add ~3x; gauge force is Nf-indep). Run at Nf=6 to see the balance + dH at 3 stacks.
  std::vector<Setup> setups = {
    { "base {0,0.4,1.0} MG20",  { Complex(0.0,0.0), Complex(0.4,0.0), Complex(1.0,0.0) }, { 4, 4, 4 }, 20  },
    { "safe {0,0.5,1.0} MG100", { Complex(0.0,0.0), Complex(0.5,0.0), Complex(1.0,0.0) }, { 4, 4, 4 }, 100 },
  };

  const int kcfg = cfgk.back();
  U.read( cfgdir+"ckpoint_lat."+std::to_string(kcfg) );
  D.update( U );
  Gauge Ucfg( U );
  Force pi0( base );
  pi0.gaussian( rng );

  std::cout << "\n# ===== L=" << Comp::N_REFINE << " setup comparison  config " << kcfg
            << "  gauge x(per-setup MG)  tau=" << tmax << " =====" << std::endl;
  { Force fg( base ); SW.get_force( fg, U );
    std::cout << "#   gauge force (mass-indep): L2=" << std::sqrt(fg.squared_norm())
              << " Linf=" << force_inf_norm(fg, base) << std::endl; }

  for(const Setup& st : setups){
    const int K = (int)st.masses.size()-1;
    assert( (int)st.fsteps.size() == K+1 );
    std::vector<int> mult(K+1);
    mult[0] = 1;
    for(int l=1; l<=K; l++){ assert(st.fsteps[l]%st.fsteps[l-1]==0); mult[l]=st.fsteps[l]/st.fsteps[l-1]; }

    std::vector<PFptr> pfs;
    for(int f=0; f<Nf/2; f++) pfs.push_back( std::make_shared<HBpf>(D, st.masses, base) );
    U = Ucfg; D.update( U );
    for( PFptr pf : pfs ) pf->gen( rng );          // heatbath (mass-dependent phi)
    // _claude: get_force_frames now reads the FORCE eta (d_eta_f); populate it (force op==action op here).
    for( PFptr pf : pfs ) pf->update_eta_force_frames( 0, (int)st.masses.size()-1 );

    std::cout << "\n#   --- " << st.name << "  ladder {";
    for(int i=0;i<=K;i++) std::cout << st.masses[i].real() << (i<K?",":"");
    std::cout << "}  steps {";
    for(int i=0;i<=K;i++) std::cout << st.fsteps[i] << (i<K?",":"");
    std::cout << "} ---" << std::endl;

    // per-frame force norms
    for(int i=0; i<=K; i++){
      Force f( base ), ft( base );
      bool first = true;
      for( PFptr pf : pfs ){
        if(first){ pf->get_force_frames( f, U, i, i ); first=false; }
        else     { pf->get_force_frames( ft, U, i, i ); f += ft; }
      }
      std::cout << "#     frame " << i << " (m=" << st.masses[i].real() << "): L2="
                << std::sqrt(f.squared_norm()) << " Linf=" << force_inf_norm(f, base) << std::endl;
    }

    // one MD integration
    U = Ucfg; D.update( U );
    Force pi( pi0 );
    for( PFptr pf : pfs ) pf->update_eta();
    std::vector<std::unique_ptr<MDLevel<Gauge,Force>>> owned;
    std::vector<MDLevel<Gauge,Force>*> levels;
    for(int l=0; l<=K; l++){
      owned.push_back( std::make_unique<FermionGroupLevel<Fermion,PFptr,Gauge,Force>>(&D,pfs,l,l,mult[l],base) );
      levels.push_back( owned.back().get() );
    }
    owned.push_back( std::make_unique<GaugeLevel<Action,Gauge,Force>>(&SW, st.mult_gauge) );   // per-setup MG
    levels.push_back( owned.back().get() );
    MinimumNorm2ML<Gauge,Force> integ( tmax, st.fsteps[0], levels );

    reset_cg_iters();
    double h0 = 0.5*pi.squared_norm() + SW(U);
    for( PFptr pf : pfs ) h0 += pf->S();
    integ.integrate( U, pi );
    // _claude (two-op split): the ML integrator refreshed only the FORCE eta; re-solve the ACTION eta at
    // the final U so S() (accept/reject) is accurate (mirrors HMCHasenbuschML::run).
    D.update( U );
    for( PFptr pf : pfs ) pf->update_eta();
    double h1 = 0.5*pi.squared_norm() + SW(U);
    for( PFptr pf : pfs ) h1 += pf->S();
    const unsigned long long ncg = get_cg_iters();
    const double dH   = h1 - h0;
    const double pacc = std::min( 1.0, std::exp(-dH) );
    const double cost = (double)ncg / ( tmax*tmax * pacc );

    std::cout << "#     result: N_CG=" << (double)ncg << " (per-stage";
    for(int l=0; l<=K; l++) std::cout << " f" << l << "=" << (double)levels[l]->cg_iters;
    std::cout << ")  dH=" << std::showpos << dH << std::noshowpos
              << "  P_acc=" << pacc << "  Cost=" << cost << std::endl;
  }

  std::cout << "\n# ================ setup comparison done ================" << std::endl;
  return 0;
}
