// test_hasenbusch_stot_claude.cu
// _claude: s_tot (trajectory length tau) scan for the massless Hasenbusch overlap HMC, L-generic
// (-DLREF; hasenbusch_massless_impl_plan_claude.md). On the STIFFEST config (min lambda_min per L) it
// measures how far the SMEARED PLAQUETTE moves per unit trajectory length -> pick the s_tot with the best
// decorrelation-per-cost, AFTER the window/ladder/nsteps are fixed.
//
// PROBE (per NM): load the min-lambda_min config; loop 4 REFRESHED momenta (redraw pi + pseudofermions);
// for each momentum use the SAME pi/phi across all s_tot (fair -- longer s_tot integrates the same start
// further). MDsteps are FIXED to the tuned hasenbusch_steps(L); ONLY trajL=s_tot varies -> eps=s_tot/nsteps
// grows with s_tot, so N_CG (cost) stays ~constant while dH grows. Pick the largest s_tot with acceptable dH.
// After each trajectory, flow the evolved config (spatial-only Wilson flow, flow_claude.h) by t_flow and
// take the smeared observable P = <spatial plaquette angle^2> (= the non-compact "smeared plaquette").
// Report per s_tot, averaged over the 4 momenta: <dH>, <P>, <|dP|> (dP = P - P_config), <|dP|>/s_tot and
// <|dP|>/<N_CG> (decorrelation per cost).
//
// Usage: ./bin [gsq=8] [Nf=2] [config_k=0(=min-lambda_min per L)] [N_MOM=4] [t_flow=1.0]
// Build: tmp_hb_stot_local_claude.sh.

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
#include "flow_claude.h"

// stiffest config (min lambda_min) per L, from the Wilson low-mode scan (Nf2 gsq8). CLI-overridable.
static int min_lambda_config( const int L ){
  if( L==1 ) return 2105;
  if( L==2 ) return 5;
  if( L==4 ) return 93;
  return 1;
}

// smeared "plaquette": mean SPATIAL plaquette angle^2 over the (flowed) lattice (non-compact -> theta^2 is
// the meaningful quantity; the mean angle averages to ~0). Spatial only, matching the spatial-only flow.
template<typename GaugeT, typename BaseT>
static double mean_sp_plaq2( const GaugeT& U, const BaseT& base ){
  double acc = 0.0;
  long   cnt = 0;
  for(int s=0; s<Comp::Nt; s++){
    for(Idx i=0; i<(Idx)base.faces.size(); i++){
      const double th = U.plaquette_angle( s, base.faces[i] );
      acc += th*th;
      cnt++;
    }
  }
  return (cnt>0) ? acc/(double)cnt : 0.0;
}

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(6);

  double gsq    = 8.0;
  int    Nf     = 2;
  int    config_k = 0;     // 0 -> min-lambda_min config for this L
  int    N_MOM  = 4;
  double t_flow = 1.0;
  if(argc>1) gsq      = atof(argv[1]);
  if(argc>2) Nf       = atoi(argv[2]);
  if(argc>3) config_k = atoi(argv[3]);
  if(argc>4) N_MOM    = atoi(argv[4]);
  if(argc>5) t_flow   = atof(argv[5]);
  if(config_k<=0) config_k = min_lambda_config( Comp::N_REFINE );
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
  D.set_lambda( lmin_frozen, lmax_frozen );      // FIXED Zolotarev window (frozen_window_claude.h)
  Action SW( gsq, at, base );

  // ----- massless config dir + the stiffest config -----
  const std::string cfgdir = "Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)
                           + "nu0"+std::to_string(nu0)+"nt"+std::to_string(Comp::Nt)
                           + "L"+std::to_string(Comp::N_REFINE)+"/";
  const std::string str_lat = cfgdir+"ckpoint_lat."+std::to_string(config_k);
  if(!std::filesystem::exists(str_lat)){
    std::cout << "# ERROR: config not found: " << str_lat << std::endl;
    return 1;
  }
  U.read( str_lat );
  D.update( U );
  Gauge Ucfg( U );

  // frozen ladder + base per-stage step counts -> fixed inner multipliers; outer MDsteps scales with s_tot.
  const std::vector<Complex> ladder     = hasenbusch_ladder( Comp::N_REFINE );
  const std::vector<int>     base_steps = hasenbusch_steps( Comp::N_REFINE );
  const int K = (int)ladder.size()-1;
  std::vector<int> base_mult(K+1);
  base_mult[0] = 1;
  for(int l=1; l<=K; l++){
    assert( base_steps[l] % base_steps[l-1] == 0 );
    base_mult[l] = base_steps[l] / base_steps[l-1];
  }
  const int MG = 20;    // gauge substeps per heaviest-frame step

  // Wilson flow (spatial-only) for the smeared plaquette; reference = the config itself.
  Flow<Action> flow( &SW, t_flow, 100 );
  Gauge Uref( Ucfg );
  flow( Uref );
  const double P0 = mean_sp_plaq2( Uref, base );

  std::vector<double> stots = { 1.0, 1.4, 1.8 };   // 2.2 omitted (per NM)
  const int nst = (int)stots.size();

  std::cout << "# ================ s_tot scan  L=" << Comp::N_REFINE << " gsq=" << gsq << " Nf=" << Nf
            << " ================" << std::endl;
  std::cout << "# stiffest config k=" << config_k << "  ladder {";
  for(int i=0;i<=K;i++) std::cout << ladder[i].real() << (i<K?",":"");
  std::cout << "}  base_steps {";
  for(int i=0;i<=K;i++) std::cout << base_steps[i] << (i<K?",":"");
  std::cout << "}  gauge x" << MG << "  N_MOM=" << N_MOM << "  t_flow=" << t_flow << std::endl;
  std::cout << "# smeared plaquette of the config: P0 = " << P0 << std::endl;

  std::vector<double> sdH(nst,0.0), sP(nst,0.0), sdP(nst,0.0), sNcg(nst,0.0);

  for(int m=0; m<N_MOM; m++){
    Force pi0( base );
    pi0.gaussian( rng );
    // one Hasenbusch stack per doublet pair; heatbath ONCE for this momentum (fixed phi across s_tot).
    std::vector<PFptr> pfs;
    for(int f=0; f<Nf/2; f++) pfs.push_back( std::make_shared<HBpf>(D, ladder, base) );
    U = Ucfg; D.update( U );
    for( PFptr pf : pfs ) pf->gen( rng );

    for(int si=0; si<nst; si++){
      const double s_tot = stots[si];
      U = Ucfg; D.update( U );
      Force pi( pi0 );                              // SAME momentum across all s_tot
      for( PFptr pf : pfs ) pf->update_eta();

      const int MDsteps = base_steps[0];   // FIXED to hasenbusch_steps(L); only trajL=s_tot varies (eps grows)

      std::vector<std::unique_ptr<MDLevel<Gauge,Force>>> owned;
      std::vector<MDLevel<Gauge,Force>*> levels;
      for(int l=0; l<=K; l++){
        owned.push_back( std::make_unique<FermionGroupLevel<Fermion,PFptr,Gauge,Force>>(
                           &D, pfs, l, l, base_mult[l], base ) );
        levels.push_back( owned.back().get() );
      }
      owned.push_back( std::make_unique<GaugeLevel<Action,Gauge,Force>>( &SW, MG ) );
      levels.push_back( owned.back().get() );
      MinimumNorm2ML<Gauge,Force> integ( s_tot, MDsteps, levels );

      reset_cg_iters();
      double h0 = 0.5*pi.squared_norm() + SW(U);
      for( PFptr pf : pfs ) h0 += pf->S();
      integ.integrate( U, pi );
      // _claude (two-op split): re-solve the ACTION eta at the final U before S() (mirrors run()).
      D.update( U );
      for( PFptr pf : pfs ) pf->update_eta();
      double h1 = 0.5*pi.squared_norm() + SW(U);
      for( PFptr pf : pfs ) h1 += pf->S();
      const unsigned long long ncg = get_cg_iters();
      const double dH = h1 - h0;

      Gauge Uflow( U );
      flow( Uflow );
      const double P  = mean_sp_plaq2( Uflow, base );
      const double dP = std::abs( P - P0 );

      sdH[si]  += dH;
      sP[si]   += P;
      sdP[si]  += dP;
      sNcg[si] += (double)ncg;

      std::cout << "#   m=" << m << " s_tot=" << s_tot << " MDsteps=" << MDsteps
                << " dH=" << std::showpos << dH << std::noshowpos
                << " P=" << P << " |dP|=" << dP << " N_CG=" << (double)ncg << std::endl;
    }
  }

  std::cout << "\n# s_tot    <dH>           <P>            <|dP|>         <|dP|>/s_tot     <|dP|>/<N_CG>" << std::endl;
  for(int si=0; si<nst; si++){
    const double mdH  = sdH[si]/N_MOM;
    const double mP   = sP[si]/N_MOM;
    const double mdP  = sdP[si]/N_MOM;
    const double mNcg = sNcg[si]/N_MOM;
    std::cout << "  " << stots[si]
              << "   " << std::showpos << mdH << std::noshowpos
              << "   " << mP
              << "   " << mdP
              << "   " << mdP/stots[si]
              << "   " << ( mNcg>0.0 ? mdP/mNcg : 0.0 ) << std::endl;
  }
  std::cout << "\n# ================ done (pick s_tot with best decorrelation-per-cost) ================" << std::endl;
  return 0;
}
