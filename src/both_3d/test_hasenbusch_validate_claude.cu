// test_hasenbusch_validate_claude.cu
// _claude: validate the FROZEN production Hasenbusch overlap HMC (window + ladder + nsteps + tau, all from
// the single-source headers), L-generic (-DLREF). Loads a thermalized config and:
//   (1) runs a short REAL HMC chain (accept/reject) -> per-traj dH + acceptance + <P_acc>,
//   (2) a dH ~ tau^2 FLOOR-FREE check: same config + momentum + phi, integrate at nsteps {2,4,8} (no
//       accept) -> dH should drop ~4x per doubling (no step-independent floor).
// Uses the multi-timescale pipeline: MinimumNorm2ML (per-stage hasenbusch_steps) + HMCHasenbuschML -- needed
// for L4 {2,4} (heavy frame sub-cycled x2); L1/L2 {2,2} are the mult-1 special case.
//
// Usage: ./bin [gsq=8] [Nf=2] [config_k=0(=last)] [N_TRAJ=10]
// Build: tmp_hb_validate_local_claude.sh.

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
#include "hmc_hasenbusch_ml_claude.h"   // multi-timescale (handles L4 {2,4}; mult-1 for L1/L2 is a special case)

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(6);

  double gsq      = 8.0;
  int    Nf       = 2;
  int    config_k = 0;      // 0 -> last config of the massless ensemble
  int    N_TRAJ   = 10;
  if(argc>1) gsq      = atof(argv[1]);
  if(argc>2) Nf       = atoi(argv[2]);
  if(argc>3) config_k = atoi(argv[3]);
  if(argc>4) N_TRAJ   = atoi(argv[4]);
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

  Fermion D(DW, Complex(0.0,0.0), 21, 0.001);   // window overridden below
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

  const std::vector<Complex> ladder = hasenbusch_ladder( Comp::N_REFINE );
  // per-stage steps -> outer MDsteps + inner multipliers (for MinimumNorm2ML).
  const std::vector<int> base_steps = hasenbusch_steps( Comp::N_REFINE );
  const int K = (int)ladder.size()-1;
  std::vector<int> base_mult(K+1);
  base_mult[0] = 1;
  for(int l=1; l<=K; l++){
    assert( base_steps[l] % base_steps[l-1] == 0 );
    base_mult[l] = base_steps[l] / base_steps[l-1];
  }
  const int    MDsteps0 = base_steps[0];   // outer (massless) step count
  const int    MG       = 20;              // gauge substeps per heaviest-frame step (matches tuning)
  const double tmax     = hasenbusch_tau( Comp::N_REFINE );

  std::vector<PFptr> pfs;
  for(int f=0; f<Nf/2; f++) pfs.push_back( std::make_shared<HBpf>(D, ladder, base) );

  std::cout << "# ================ Hasenbusch VALIDATION  L=" << Comp::N_REFINE << " gsq=" << gsq
            << " Nf=" << Nf << " ================" << std::endl;
  std::cout << "# frozen window: lambda_min=" << lmin_frozen << " lambda_max=" << lmax_frozen
            << " k=" << (lmin_frozen/lmax_frozen) << " Delta=" << D.Delta() << std::endl;
  std::cout << "# ladder {";
  for(size_t i=0;i<ladder.size();i++) std::cout << ladder[i].real() << (i+1<ladder.size()?",":"");
  std::cout << "}  steps {";
  for(size_t i=0;i<base_steps.size();i++) std::cout << base_steps[i] << (i+1<base_steps.size()?",":"");
  std::cout << "}  tmax=" << tmax << " gauge x" << MG << "  config k=" << config_k << std::endl;

  // ============ (1) short REAL HMC chain (accept/reject) ============
  {
    Force pi( base );
    std::vector<std::unique_ptr<MDLevel<Gauge,Force>>> owned;
    std::vector<MDLevel<Gauge,Force>*> levels;
    for(int l=0; l<=K; l++){
      owned.push_back( std::make_unique<FermionGroupLevel<Fermion,PFptr,Gauge,Force>>(
                         &D, pfs, l, l, base_mult[l], base ) );
      levels.push_back( owned.back().get() );
    }
    owned.push_back( std::make_unique<GaugeLevel<Action,Gauge,Force>>( &SW, MG ) );
    levels.push_back( owned.back().get() );
    MinimumNorm2ML<Gauge,Force> integrator( tmax, MDsteps0, levels );
    HMCHasenbuschML hmc( rng, &SW, &D, U, pi, pfs, &integrator );
    D.update( U );

    std::cout << "\n# ---- (1) HMC chain (" << N_TRAJ << " trajectories) ----" << std::endl;
    std::cout << "#   traj    dH               accept   P_acc          sec" << std::endl;
    double sum_absdH=0.0, sum_pacc=0.0, sum_sec=0.0;
    int    n_acc=0;
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
      std::cout << "#   " << std::setw(4) << k << "    " << std::showpos << dH << std::noshowpos
                << "     " << (is_accept?"Y":"n") << "      " << pacc << "   " << sec << std::endl;
    }
    std::cout << "#   summary: <|dH|>=" << sum_absdH/N_TRAJ
              << "  acceptance=" << (double)n_acc/N_TRAJ
              << "  <P_acc>=" << sum_pacc/N_TRAJ
              << "  <sec/traj>=" << sum_sec/N_TRAJ << std::endl;
  }

  // ============ (2) dH ~ tau^2 FLOOR-FREE check ============
  // Same config + same momentum + same phi; integrate at MDsteps {1x,2x,4x} the tuned value (no accept).
  // dH should fall ~4x per doubling. A step-INDEPENDENT floor would signal a non-exact force (old lambda_max bug).
  {
    U = Ucfg; D.update( U );
    Force pi0( base );
    pi0.gaussian( rng );
    for( PFptr pf : pfs ) pf->gen( rng );          // fixed phi

    std::cout << "\n# ---- (2) dH ~ tau^2 (fixed config/momentum/phi; MDsteps scaled 1x,2x,4x) ----" << std::endl;
    std::cout << "#   MDsteps   dH                   ratio(prev/this)" << std::endl;
    double prev = 0.0;
    int    prev_n = 0;
    for(int fac=1; fac<=4; fac*=2){
      const int MD = MDsteps0*fac;
      U = Ucfg; D.update( U );
      Force pi( pi0 );
      for( PFptr pf : pfs ) pf->update_eta();
      double h0 = 0.5*pi.squared_norm() + SW(U);
      for( PFptr pf : pfs ) h0 += pf->S();

      std::vector<std::unique_ptr<MDLevel<Gauge,Force>>> owned;
      std::vector<MDLevel<Gauge,Force>*> levels;
      for(int l=0; l<=K; l++){
        owned.push_back( std::make_unique<FermionGroupLevel<Fermion,PFptr,Gauge,Force>>(
                           &D, pfs, l, l, base_mult[l], base ) );
        levels.push_back( owned.back().get() );
      }
      owned.push_back( std::make_unique<GaugeLevel<Action,Gauge,Force>>( &SW, MG ) );
      levels.push_back( owned.back().get() );
      MinimumNorm2ML<Gauge,Force> integ( tmax, MD, levels );
      integ.integrate( U, pi );

      // _claude (two-op split): re-solve the ACTION eta at the final U before S() (mirrors run()).
      D.update( U );
      for( PFptr pf : pfs ) pf->update_eta();
      double h1 = 0.5*pi.squared_norm() + SW(U);
      for( PFptr pf : pfs ) h1 += pf->S();
      const double dH = h1 - h0;
      const double ratio = (prev_n>0) ? (prev/dH) : 0.0;
      std::cout << "#   " << std::setw(3) << MD << "      " << std::showpos << dH << std::noshowpos
                << "     " << ratio << std::endl;
      prev = dH; prev_n = MD;
    }
    std::cout << "#   (ratio ~4 per doubling = dH~tau^2, floor-free; a flat ratio -> a dH floor)" << std::endl;
  }

  std::cout << "\n# ================ validation done ================" << std::endl;
  return 0;
}
