// hmc_hasenbusch_block_scc_claude.cu
// _scc: BU SCC copy of hmc_hasenbusch_block_claude.cu -- IDENTICAL physics/pipeline; the ONLY change is
// that the geometry paths are ABSOLUTE (SCC /projectnb/qfe/nmatsum/qed3/geometry) instead of the relative
// "../../geometry" of the base driver, so the binary can be launched from an SGE batch node with any cwd.
// Mirrors the _fermilab driver-copy convention. Keep this file in sync with the base driver except for the
// two geometry-path lines below. (The .dat geometry files are generated separately by NM into geometry/data/.)
//
// _claude: production overlap HMC with Hasenbusch mass preconditioning -- Nf-BLOCK (mrhs) version
// (nfblock_hasenbusch_impl_plan_claude.md). Same pipeline as the serial hmc_hasenbusch_claude.cu but the
// Nf/2 = NSTACK pseudofermion flavors are PACKED into one N*NSTACK block per frame and driven through
// BlockedMat / BlockedForce (inner Zolotarev multishift + outer block CG). UNIFIED: handles all Nf >= 2
// (Nf=2 => NSTACK=1, a no-op block). Nf is COMPILE-TIME (build -DNF=; NSTACK = Nf/2 is a template arg);
// build -DLREF=1/2/4 and mass-generic:
//   * frame 0 = the PHYSICAL target fermion mass (argv mass_re/mass_im; 0 = massless, !=0 = massive),
//   * frames 1..K = auxiliary M_mass coefficients from hasenbusch_ladder(L).
// Pipeline: HasenbuschPFBlock (two BlockedMats: action D + force Df; set_frame_mass) + MinimumNorm2ML +
// HMCHasenbuschML (multi-timescale, per-stage hasenbusch_steps(L); the ML templates are unchanged -- the
// block PF has the same interface as HasenbuschPF, driven via a size-1 pfs vector). FROZEN window.
// LOCAL geometry paths (../../geometry). Refs: M. Hasenbusch, hep-lat/0107019; B. Jegerlehner, hep-lat/9612014
// (block multishift). The serial hmc_hasenbusch_claude.cu is KEPT UNCHANGED (reference / current ensembles).

#include <typeinfo>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <filesystem>
#include <thread>
#include <chrono>
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
  constexpr int Nt=128;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;
  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-8;
}

// _scc: ABSOLUTE geometry paths for the BU SCC (was "../../geometry/..."). The .dat files under
// geometry/data/ are generated separately by NM (may not exist yet at build time; needed at run time).
const std::string dir = "/projectnb/qfe/nmatsum/qed3/geometry/data/";
#include "/projectnb/qfe/nmatsum/qed3/geometry/geodesic.h"

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
#include "blocked_mat_claude.h"                     // BlockMemPool + BlockedMat (mrhs; matpoly pulled the block kernels)
#include "overlap_force_Nfblock_claude.h"           // BlockedForce (Term A + block Term B); AFTER blocked_mat
#include "pseudofermion_hasenbusch_block_claude.h"  // HasenbuschPFBlock (Nf-BLOCK; replaces serial HasenbuschPF)
#include "hmc_hasenbusch_ml_claude.h"   // ALWAYS the multi-timescale integrator (per NM); mult-1 is a special case

// Nf fixed at COMPILE time for the block NSTACK = Nf/2 (BlockedMat is a template on NSTACK). Override -DNF=.
#ifndef NF
#define NF 2
#endif

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) == "-h") {
      printf("Usage: ./a.out [gsq] [Nf] [nu0] [mass_re] [mass_im] [max_sec]\n");
      printf("  gsq      Wilson coupling squared (default: 8.0)\n");
      printf("  Nf       number of fermion flavors (default: 2)\n");
      printf("  nu0      Wilson param (default: 1.0)\n");
      printf("  mass_re  PHYSICAL target fermion mass, real part (default: 0.0 = massless)\n");
      printf("  mass_im  PHYSICAL target fermion mass, imag part (default: 0.0)\n");
      printf("  max_sec  wall-time budget in seconds, 0 = unlimited (default: 0.0)\n");
      printf("  (Hasenbusch: frame 0 = physical m; auxiliary M_mass-coeff ladder per L, hasenbusch_ladder)\n");
      return 0;
    }
  }

  double gsq = 8.0;
  if(argc>1) gsq = atof(argv[1]);
  int Nf = 2;
  if(argc>2) Nf = atoi(argv[2]);
  double nu0 = 1.0;
  if(argc>3) nu0 = atof(argv[3]);
  double mass_re = 0.0;
  if(argc>4) mass_re = atof(argv[4]);
  double mass_im = 0.0;
  if(argc>5) mass_im = atof(argv[5]);
  double max_sec = 0.0;
  if(argc>6) max_sec = atof(argv[6]);
  const Complex mass = Complex(mass_re, mass_im);   // PHYSICAL target fermion mass (0 = massless)

  // Hasenbusch ladder per L: frame 0 = the PHYSICAL target mass, frames 1..K = auxiliary M_mass
  // COEFFICIENTS (hasenbusch_ladder). L-generic: build -DLREF=1/2/4 -> aux {0.4} / {0.1,0.4} / {0.02,0.4}.
  // (For a massive target, m*Abar/abar_s must stay below the first auxiliary coefficient for the split.)
  std::vector<Complex> masses = hasenbusch_ladder( Comp::N_REFINE );
  masses[0] = mass;   // frame 0: physical target (0 = massless); auxiliaries SHIFTED by resc below (after base)

  std::cout << "# gsq = " << gsq << " Nf = " << Nf << " nu0 = " << nu0 << " mass = " << mass
            << " (Hasenbusch L=" << Comp::N_REFINE << "; frame0=physical m, aux ladder = c_i + rescaled mass)" << std::endl;
  std::cout << "# max_sec = " << max_sec << " (wall-time budget; 0 = unlimited)" << std::endl;
  Timer wall_timer;

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();

  using BaseLink = std::array<Idx,2>;
  constexpr Idx N = Comp::N;
  constexpr int Nt = Comp::Nt;

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Force=GaugeExt<Base,Nt,Comp::is_compact>;
  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  using Action=U1WilsonExt<Base>;
  using Rng=ParallelRngExt<Base,Nt>;
  using Fermion=OverlapWMass<WilsonDirac>;
  static constexpr int NSTACK = NF/2;                     // Nf/2 pseudofermions packed into one block
  using HBpf=HasenbuschPFBlock<Fermion,Force,NSTACK>;     // Nf-BLOCK Hasenbusch PF (all NSTACK flavors)
  using PFptr=std::shared_ptr<HBpf>;

  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;

  // _claude (MASSIVE Hasenbusch, c += rescaled mass): shift the AUXILIARY coefficients UP by the physical
  // rescaled mass resc = m*Abar/abar_s = mass*mean_dual_area/mean_ell (the mass_coeff of the physical mass).
  // -> ladder in coeff space = {resc, c_1+resc, ...}, so the INTER-FRAME gaps (c_i, tuned massless) are
  // PRESERVED for any physical mass; no retuning + no degeneracy at heavy m. Frame 0 stays physical
  // (set_mass rescales internally). m=0 -> resc=0 -> massless ladder unchanged. Force-validated L1 (all
  // masses {0.1,0.5,1.0,1.5}, grad-vs-FD ~1e-8): test_hasenbusch_force_massive_l1_claude.cu.
  const Complex resc = mass * (base.mean_dual_area / base.mean_ell);
  for(size_t i=1; i<masses.size(); i++) masses[i] += resc;
  {
    std::string aux_str;
    for(size_t i=1; i<masses.size(); i++) aux_str += (i>1?",":"") + std::to_string(masses[i].real());
    std::cout << "# Hasenbusch ladder (resc = " << resc.real() << "): frame0 = phys m = " << mass.real()
              << ", aux coeffs {" << aux_str << "}" << std::endl;
  }

  const double r = 1.0;
  const double M5 = -1.0;
#ifndef AT_VAL
#define AT_VAL 0.2
#endif
  const double at = AT_VAL;   // temporal lattice spacing (a_t). -DAT_VAL=0.1 -> half-a_t ensembles (2026-07-22).
                              // dir3 embeds to_string(at) -> at0.100000 dirs, distinct from at0.200000 (no clash).
  assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);
  std::cout << "# DW set" << std::endl;

  Gauge U(base);
  Rng rng(base);

  // massless overlap with the FIXED FROZEN Zolotarev window per L (frozen_window_claude.h): set both edges
  // + rebuild the Zolotarev fit on k=lambda_min/lambda_max + freeze (no per-config recompute) -> exact HMC.
  Fermion D(DW, mass, hasenbusch_naction(Comp::N_REFINE), 0.001);  // action n_act(L); 0.001 overridden by set_lambda below
  double lmin_frozen, lmax_frozen;
  frozen_window( Comp::N_REFINE, lmin_frozen, lmax_frozen );
  D.set_lambda( lmin_frozen, lmax_frozen );
  std::cout << "# Dov set; M5 = " << M5 << std::endl;
  D.update(U);
  std::cout << "# frozen window: lambda_min " << D.lambda_min << " lambda_max " << D.lambda_max
            << " k " << D.lambda_min/D.lambda_max << " delta " << D.Delta() << std::endl;

  // _claude (two-operator split-pole, two_operator_force_impl_plan_claude.md): the FORCE operator Df uses
  // a cruder Zolotarev (n_f = hasenbusch_nforce(L) poles) on the SAME D_W + SAME frozen window. The MD
  // force is -dS_f/dU (Df); heatbath + accept/reject stay on the accurate n=21 D. Exact by Metropolis.
  // FORCE op: n_f poles on the NARROWED window [2*lambda_min, lambda_max] (better inner conditioning ->
  // cheaper, reversibility-clean; the sign approx below 2*lambda_min is cruder but Metropolis-corrected).
  const int nforce = hasenbusch_nforce( Comp::N_REFINE );
  const double lmin_force = 2.0 * lmin_frozen;
  Fermion Df(DW, mass, nforce, 0.001);
  Df.set_lambda( lmin_force, lmax_frozen );
  Df.update(U);
  std::cout << "# Dov FORCE op set; n_force = " << nforce << " force lambda_min " << lmin_force
            << " delta_force " << Df.Delta() << std::endl;

  Action SW( gsq, at, base );
  std::cout << "# alat = " << base.mean_ell << std::endl;

  // ONE Nf-BLOCK Hasenbusch stack packs all Nf/2 = NSTACK flavors (mrhs). Nf is compile-time (-DNF); the
  // runtime Nf must match so NSTACK = Nf/2 is consistent. (The size-1 pfs vector reuses the ML integrator's
  // per-pf loop unchanged -- HasenbuschPFBlock has the same interface as the serial HasenbuschPF.)
  assert(Nf%2==0);
  if(Nf != NF){
    std::cerr << "# ERROR: runtime Nf=" << Nf << " != compile-time NF=" << NF
              << " (rebuild with -DNF=" << Nf << ")" << std::endl;
    return 1;
  }
  std::vector<std::shared_ptr<HBpf>> pfs;
  pfs.push_back( std::make_shared<HBpf>(D, Df, masses, base) );

  Timer timer;

  // Hasenbusch-tagged output dir (tag = ladder coefficients) so it does NOT collide with the
  // non-Hasenbusch massless runs and is per-L-ladder distinct (e.g. L1 _hb0.4, L2 _hb0.1-0.4).
  std::string hbtag = "_hb";
  for(size_t i=1; i<masses.size(); i++) hbtag += (i>1?"-":"") + std::to_string(masses[i].real());
  std::string dir3;
  dir3="Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)
      +"mRe"+std::to_string(mass.real())+"mIm"+std::to_string(mass.imag())
      +"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+hbtag+"/";
  std::filesystem::create_directory(dir3);
  const int k_ckpoint=1;
#ifndef KRNG
#define KRNG 5
#endif
  const int k_ckpoint_rng=KRNG;   // -DKRNG=1 -> FULL checkpointing (keep every rng, no thinning)
#ifndef KMAX
#define KMAX 1200
#endif
  const int kmax=KMAX;   // -DKMAX bounds trajectory count (set >= target sample size per gsq)

  int k_tmp=0;
  {
    for(int k_scan=k_ckpoint; k_scan<=kmax; k_scan+=k_ckpoint ){
      const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k_scan);
      const std::string str_rng=dir3+"ckpoint_rng."+std::to_string(k_scan);
      const bool bool_lat = std::filesystem::exists(str_lat);
      const bool bool_rng = std::filesystem::exists(str_rng);
      if(!bool_lat) break;
      if(bool_lat && bool_rng) k_tmp = k_scan;
    }
    if(k_tmp>0){
      std::cout << "read from k_tmp = " << k_tmp << std::endl;
      const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k_tmp);
      const std::string str_rng=dir3+"ckpoint_rng."+std::to_string(k_tmp);
      U.read( str_lat );
      rng.read( str_rng );
    }
  }
  std::cout << "#starting from k_tmp = " << k_tmp << std::endl;

  Force pi( base );
  const double tmax = hasenbusch_tau( Comp::N_REFINE );   // trajectory length (s_tot scan)
  // Multi-timescale MinimumNorm2ML from the tuned per-stage hasenbusch_steps(L): outer MDsteps = steps[0],
  // level-l multiplier = steps[l]/steps[l-1], gauge sub-cycled MG per heaviest-frame step. Handles mult-1
  // (L1) and genuine multi-timescale (L2 {2,4}, L4 {2,2,4}) uniformly.
  const std::vector<int> base_steps = hasenbusch_steps( Comp::N_REFINE );
  const int K = (int)masses.size()-1;
  std::vector<int> base_mult(K+1);
  base_mult[0] = 1;
  for(int l=1; l<=K; l++){
    assert( base_steps[l] % base_steps[l-1] == 0 );
    base_mult[l] = base_steps[l] / base_steps[l-1];
  }
  const int MDsteps = base_steps[0];
  const int MG = hasenbusch_mg( Comp::N_REFINE );   // gauge substeps per heaviest-frame step (per-L; L4 -> 100)
  std::cout << "# tmax = " << tmax << "  steps {";
  for(size_t i=0;i<base_steps.size();i++) std::cout << base_steps[i] << (i+1<base_steps.size()?",":"");
  std::cout << "}  MDsteps = " << MDsteps << "  gauge x" << MG << std::endl;

  std::vector<std::unique_ptr<MDLevel<Gauge,Force>>> owned;
  std::vector<MDLevel<Gauge,Force>*> levels;
  for(int l=0; l<=K; l++){
    owned.push_back( std::make_unique<FermionGroupLevel<Fermion,PFptr,Gauge,Force>>(
                       &Df, pfs, l, l, base_mult[l], base ) );   // level drives the FORCE op Df
    levels.push_back( owned.back().get() );
  }
  owned.push_back( std::make_unique<GaugeLevel<Action,Gauge,Force>>( &SW, MG ) );
  levels.push_back( owned.back().get() );
  MinimumNorm2ML<Gauge,Force> integrator( tmax, MDsteps, levels );
  HMCHasenbuschML hmc( rng, &SW, &D, U, pi, pfs, &integrator );
  D.update( U );

  double rate, dH;
  bool is_accept;

  double r_mean = 0.0;
  double last_traj_sec = 0.0;
  for(int k=k_tmp+1; k<kmax; k++){
    if(max_sec > 0.0 && last_traj_sec > 0.0){
      const double elapsed = wall_timer.currentSeconds();
      const double est = 1.3*last_traj_sec;
      if(elapsed + est > max_sec){
        std::cout << "# wall budget reached: stopping before traj " << k
                  << " (elapsed " << elapsed << "s + est " << est << "s > budget " << max_sec << "s)" << std::endl;
        break;
      }
    }
    Timer timer;
    hmc.run( rate, dH, is_accept);
    std::cout << "# dH : " << dH
              << " is_accept : " << is_accept
              << " rate : " << rate << std::endl;
    r_mean += rate;
    last_traj_sec = timer.currentSeconds();
    std::cout << "# HMC : " << last_traj_sec << " sec" << std::endl;

    if(k%20==0) D.is_update = false;
    if(k%100==0){
      std::cout << "# k = " << k << std::endl;
    }

    if(k%k_ckpoint==0){
      const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k);
      const std::string str_rng=dir3+"ckpoint_rng."+std::to_string(k);
      U.ckpoint( str_lat );
      rng.ckpoint( str_rng );
      int k_prev = k - k_ckpoint;
      if(k_prev > 0 && k%k_ckpoint_rng != 0){
        std::error_code ec;
        std::filesystem::remove(dir3+"ckpoint_rng."+std::to_string(k_prev), ec);
        if(ec){
          std::cout << "# error removing ckpoint_rng." << k_prev << ": " << ec.message() << std::endl;
          assert(false);
        }
      }
    }
  }
  r_mean /= kmax;
  std::cout << "# r_mean = " << r_mean << std::endl;

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();
  return 0;
}
