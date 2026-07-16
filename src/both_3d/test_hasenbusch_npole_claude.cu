// test_hasenbusch_npole_claude.cu
// _claude: FORCE cost test for the TWO-OPERATOR split-pole overlap HMC (two_operator_force_impl_plan_claude.md).
// ACTION op D (n_act(L)=25/25/31, full frozen window) does heatbath + accept/reject; FORCE op Df (n_f=hasenbusch_nforce=11)
// drives the MD force only, on a NARROWED window [w*lambda_min, lambda_max]. Exact by the Metropolis
// correction (Duane-Kennedy-Pendleton-Roweth, Phys. Lett. B 195 (1987) 216; cruder MD approx: Clark-Kennedy,
// PRL 98 (2007) 051601, hep-lat/0608015). The n_f-scan finding (2026-07-14): pole count is a weak speed knob
// (CG iters are npole-independent); the WINDOW factor w (larger k_force -> better inner conditioning) is the
// reversibility-clean lever -- so this test scans w.
//
// FINALIZED force (n_f=11 @ 2*lambda_min, action n=31). Hasenbusch SETUP comparison: for L4, compares the
// current {0,0.2,0.5}{2,2,2} vs the 2-stage {0,0.5}{2,4} (other L = header default). Per setup, ALWAYS the
// per-frame force comparison (CG iters / wall-time / L2 / Linf per frame), then:
//   * N_TRAJ>0 -> a REAL HMC chain: per-traj dH/P/sec + <|dH|>, acceptance, force CG/traj, sec/traj, Osborn Cost,
//   * N_TRAJ<=0 -> stop after the per-frame force (probe-only).
//
// Usage: ./bin [gsq=8] [Nf=2] [config_k=0(=last)] [N_TRAJ=16] [winfac_single=0]
//   N_TRAJ<=0 -> probe-only;  winfac_single is legacy (window factor is fixed at 2 in the setup comparison).
//   Build with -DIsVerbose2 (matpoly log macro) to get per-solve SOLVER/MULTISHIFT #iter + shift drift.
// Build: tmp_hb_npole_local_claude.sh.

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
  constexpr int Nt=128;

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
#include "frozen_window_claude.h"
#include "hasenbusch_ladder_claude.h"
#include "pseudofermion_hasenbusch_claude.h"
#include "hmc_hasenbusch_ml_claude.h"

// L-infinity (max |component|) force norm over the whole lattice (all links + temporal sites, all t).
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

  double gsq      = 8.0;
  int    Nf       = 2;
  int    config_k = 0;      // 0 -> last config of the massless ensemble
  int    N_TRAJ   = 16;
  int    winfac_single = 0; // 0 -> scan {1,2}; >0 -> just that window factor w (diagnostics)
  if(argc>1) gsq      = atof(argv[1]);
  if(argc>2) Nf       = atoi(argv[2]);
  if(argc>3) config_k = atoi(argv[3]);
  if(argc>4) N_TRAJ   = atoi(argv[4]);   // <=0 -> probe-only (per-frame CG/time diagnostic; skip the chain)
  if(argc>5) winfac_single = atoi(argv[5]);
  (void)winfac_single;   // legacy arg (window factor fixed at 2 in the setup comparison)
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

  // ACTION operator (n_act(L) = 25 for L1/L2, 31 for L4) + full frozen window (heatbath + accept/reject).
  Fermion D(DW, Complex(0.0,0.0), hasenbusch_naction(Comp::N_REFINE), 0.001);
  double lmin_frozen, lmax_frozen;
  frozen_window( Comp::N_REFINE, lmin_frozen, lmax_frozen );
  D.set_lambda( lmin_frozen, lmax_frozen );
  Action SW( gsq, at, base );

  // ----- massless config dir + a thermalized config (default = last) -----
  const std::string cfgdir = "Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)
                           + "nu0"+std::to_string(nu0)+"nt"+std::to_string(Comp::Nt)
                           + "L"+std::to_string(Comp::N_REFINE)+"/";
  if(config_k<=0){
    for(int k=1; k<=100000; k++)
      if(std::filesystem::exists(cfgdir+"ckpoint_lat."+std::to_string(k))) config_k = k;
  }
  const std::string str_lat = cfgdir+"ckpoint_lat."+std::to_string(config_k);
  if(!std::filesystem::exists(str_lat)){
    std::cout << "# ERROR: config not found: " << str_lat << std::endl;
    return 1;
  }
  U.read( str_lat );
  D.update( U );
  Gauge Ucfg( U );

  const int    MG   = 20;
  const double tmax = hasenbusch_tau( Comp::N_REFINE );
  const int    nf   = hasenbusch_nforce( Comp::N_REFINE );

  // FORCE op (built ONCE; ladder-independent): n_f=11 poles on the NARROWED window [2*lambda_min, lambda_max].
  const double lmin_f = 2.0 * lmin_frozen;
  Fermion Df(DW, Complex(0.0,0.0), nf, 0.001);
  Df.set_lambda( lmin_f, lmax_frozen );
  Df.update( Ucfg );

  // Hasenbusch SETUP comparison (ladder + steps). For L4 compare the current {0,0.2,0.5}{2,2,2} vs the
  // 2-stage {0,0.5}{2,4}; other L use the header default. Each setup runs an N_TRAJ chain with the FINALIZED
  // force (n_f=11 @ 2*lambda_min) + n=31 action, recording per-traj dH/P and <|dH|>/acceptance/Osborn Cost.
  struct Setup { std::string name; std::vector<Complex> masses; std::vector<int> steps; };
  std::vector<Setup> setups;
  if(Comp::N_REFINE==4){
#ifdef L4_STEPS224
    // -DL4_STEPS224 (e.g. the GPU1 parallel run): the single L4 setup {0,0.2,0.5} with steps {2,2,4}
    // (frame 2 sub-cycled x2). Kept separate so it does not duplicate the default two-setup comparison.
    setups.push_back( { "{0,0.2,0.5} steps{2,2,4}", { Complex(0,0), Complex(0.2,0), Complex(0.5,0) }, { 2, 2, 4 } } );
#else
    // setups.push_back( { "{0,0.2,0.5} steps{3,3,3}",     { Complex(0,0), Complex(0.2,0), Complex(0.5,0) },                { 3, 3, 3 } } );
    // setups.push_back( { "{0,0.1,0.3,0.5} steps{2,2,2,4}", { Complex(0,0), Complex(0.1,0), Complex(0.3,0), Complex(0.5,0) }, { 2, 2, 2, 4 } } );
    // setups.push_back( { "{0,0.2,0.5} steps{2,2,6}",     { Complex(0,0), Complex(0.2,0), Complex(0.5,0) },                { 2, 2, 6 } } );
    setups.push_back( { "{0,0.16,0.32,0.5} steps{2,2,2,4}", { Complex(0,0), Complex(0.16,0), Complex(0.32,0), Complex(0.5,0) }, { 2, 2, 2, 4 } } );
    setups.push_back( { "{0,0.2,0.4,0.5} steps{2,2,2,4}",   { Complex(0,0), Complex(0.2,0),  Complex(0.4,0),  Complex(0.5,0) }, { 2, 2, 2, 4 } } );
#endif
  } else if(Comp::N_REFINE==2){
#ifdef L2_COMPARE
    // -DL2_COMPARE: L2 setup comparison {0,0.5}{3,3} vs {0,0.5}{2,4} (opt-in; default L2 = single header).
    setups.push_back( { "{0,0.5} steps{3,3}", { Complex(0,0), Complex(0.5,0) }, { 3, 3 } } );
    setups.push_back( { "{0,0.5} steps{2,4}", { Complex(0,0), Complex(0.5,0) }, { 2, 4 } } );
#else
    setups.push_back( { "header", hasenbusch_ladder(2), hasenbusch_steps(2) } );
#endif
  } else {
    setups.push_back( { "header", hasenbusch_ladder(Comp::N_REFINE), hasenbusch_steps(Comp::N_REFINE) } );
  }

  std::cout << "# ================ Hasenbusch SETUP comparison  L=" << Comp::N_REFINE
            << " gsq=" << gsq << " Nf=" << Nf << " ================" << std::endl;
  std::cout << "# n_act=" << hasenbusch_naction(Comp::N_REFINE) << " action, n_f=" << nf
            << " force @ lambda_min=" << lmin_f << " (=2*" << lmin_frozen
            << "); Delta_action=" << D.Delta() << " Delta_force=" << Df.Delta() << std::endl;
  std::cout << "# tmax=" << tmax << " gauge x" << MG << "  config k=" << config_k
            << "  N_TRAJ=" << N_TRAJ << "  Cost = force_NCG/(tau^2 P_acc)." << std::endl;

  for(const Setup& st : setups){
    const std::vector<Complex>& ladder = st.masses;
    const std::vector<int>&     base_steps = st.steps;
    const int K = (int)ladder.size()-1;
    std::vector<int> base_mult(K+1);
    base_mult[0] = 1;
    for(int l=1; l<=K; l++){
      assert( base_steps[l] % base_steps[l-1] == 0 );
      base_mult[l] = base_steps[l] / base_steps[l-1];
    }
    const int MDsteps0 = base_steps[0];

    std::vector<PFptr> pfs;
    for(int f=0; f<Nf/2; f++) pfs.push_back( std::make_shared<HBpf>(D, Df, ladder, base) );

    std::cout << "\n# ===== setup " << st.name << "  MDsteps=" << MDsteps0 << " gauge x" << MG
              << " =====" << std::endl;

    // --- PER-FRAME force comparison (ALWAYS): isolates each frame's eta_f solve + grad -> CG iters +
    // wall-time + force norms, on the loaded config with a fixed heatbath. Then (if N_TRAJ>0) the chain. ---
    U = Ucfg; D.update( U ); Df.update( U );
    for( PFptr pf : pfs ) pf->gen( rng );
    std::cout << "#   per-frame force (mass coeff / CG iters / sec / L2 / Linf):" << std::endl;
    for(int i=0; i<=K; i++){
      reset_cg_iters();
      Timer ftime;
      for( PFptr pf : pfs ) pf->update_eta_force_frame( i );
      Force f( base ), ft( base );
      bool first = true;
      for( PFptr pf : pfs ){
        if(first){ pf->get_force_frames( f, U, i, i ); first=false; }
        else     { pf->get_force_frames( ft, U, i, i ); f += ft; }
      }
      const double sec = ftime.currentSeconds();
      const unsigned long long cg = get_cg_iters();
      std::cout << "#     frame " << i << " (c=" << ladder[i].real() << "): CG=" << (double)cg
                << "  sec=" << sec << "  L2=" << std::sqrt(f.squared_norm())
                << "  Linf=" << force_inf_norm(f, base) << std::endl;
    }

    if(N_TRAJ<=0){
      std::cout << "#   (probe-only: N_TRAJ<=0, chain skipped)" << std::endl;
      continue;
    }

    // --- REAL HMC chain (accept/reject with n=31, force n_f=11 @ 2*lambda_min) ---
    Force pi( base );
    std::vector<std::unique_ptr<MDLevel<Gauge,Force>>> owned;
    std::vector<MDLevel<Gauge,Force>*> levels;
    for(int l=0; l<=K; l++){
      owned.push_back( std::make_unique<FermionGroupLevel<Fermion,PFptr,Gauge,Force>>(
                         &Df, pfs, l, l, base_mult[l], base ) );   // level drives the FORCE op Df
      levels.push_back( owned.back().get() );
    }
    owned.push_back( std::make_unique<GaugeLevel<Action,Gauge,Force>>( &SW, MG ) );
    levels.push_back( owned.back().get() );
    MinimumNorm2ML<Gauge,Force> integrator( tmax, MDsteps0, levels );
    HMCHasenbuschML hmc( rng, &SW, &D, U, pi, pfs, &integrator );

    U = Ucfg; D.update( U ); Df.update( U );
    reset_cg_iters();
    double sum_absdH=0.0, sum_pacc=0.0, sum_sec=0.0;
    int    n_acc=0;
    std::cout << "#   chain (" << N_TRAJ << " trajectories, fresh momentum + phi each):" << std::endl;
    std::cout << "#     traj    dH                accept   P=min(1,e^-dH)   sec" << std::endl;
    for(int k=1; k<=N_TRAJ; k++){
      double rate, dH;
      bool is_accept;
      Timer traj_timer;
      hmc.run( rate, dH, is_accept );
      const double sec = traj_timer.currentSeconds();
      const double pacc = std::min( 1.0, std::exp(-dH) );
      sum_absdH += std::abs(dH);
      sum_pacc  += pacc;
      sum_sec   += sec;
      n_acc     += (is_accept?1:0);
      std::cout << "#     " << std::setw(4) << k << "    " << std::showpos << dH << std::noshowpos
                << "     " << (is_accept?"Y":"n") << "      " << pacc << "     " << sec << std::endl;
    }
    // FORCE CG = sum of the fermion levels' cg_iters (force eta-solve + grad, all at n_f). The action-eta
    // re-solve in run() (n=31, before h1) is OUTSIDE the levels; total - force = that action overhead.
    unsigned long long force_cg = 0;
    for(int l=0; l<=K; l++) force_cg += levels[l]->cg_iters;
    const unsigned long long total_cg = get_cg_iters();
    const double force_cg_traj = (double)force_cg / N_TRAJ;
    const double pacc_avg = sum_pacc / N_TRAJ;
    const double cost = force_cg_traj / ( tmax*tmax * pacc_avg );

    std::cout << "#   RESULT " << st.name << ":  <|dH|>=" << sum_absdH/N_TRAJ
              << "  acceptance=" << (double)n_acc/N_TRAJ
              << "  <P_acc>=" << pacc_avg << std::endl;
    std::cout << "#   cost:  force_NCG/traj=" << force_cg_traj
              << "  (total_NCG/traj=" << (double)total_cg/N_TRAJ
              << ", action re-solve=" << (double)(total_cg-force_cg)/N_TRAJ << ")"
              << "  <sec/traj>=" << sum_sec/N_TRAJ
              << "  Cost=" << cost << std::endl;
  }

  std::cout << "\n# ================ setup comparison done ================" << std::endl;
  return 0;
}
