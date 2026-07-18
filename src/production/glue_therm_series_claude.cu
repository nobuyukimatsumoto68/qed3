#include <typeinfo>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <filesystem>
#include <chrono>
#include <cstdint>
#include <complex>
#include <array>
#include <vector>
#include <map>
#include <set>
#include <Eigen/Dense>

// glue_therm_series_claude.cu
//
// GLUONIC THERMALIZATION MONITOR / FLOW-SCALE t0 (plan: redo_obs_measurement_impl_plan_claude.md Sec 4).
// Per-config FLOW TRAJECTORY of the SPATIAL average gauge-action density E_s(t_fl)
// (noncompact 0.5 beta_s theta^2 integrand of U1WilsonExt::operator(), magnetic/per-face
// channel ONLY, SAME beta_s weight as the glueball code; NM 2026-07-17). v6 recording:
// FIXED eps = 0.01, saved EVERY step, tmax = RUNTIME arg 8 (default 2.0, COMMON to all
// ensembles) -> t_fl = 0, 0.01, ..., tmax (tmax/eps + 1 points). Scale/coupling readout
// (Sec 4b of the plan): g_R^2 from the t^{3/2}E_s plateau; t0 from dimensionless t^2 E_s = c.
// Flow = per-timeslice SPATIAL gradient flow of the SIMULATED action (get_spatial, beta_s
// weight, NOT coupling-stripped). SWITCHED from full 3D (NM 2026-07-17): matches the glueball
// "on a timeslice" convention, and avoids the anisotropic two-rate kernel (beta_s ~ a_t vs
// beta_t ~ 1/a_t) of the full flow -- spatial-only has ONE uniform prefactor a_t/g^2, so the
// a_t bookkeeping drops out of the t/g^2 collapse. LO scaling unchanged: E_s ~ g^2/t^{3/2}.
// Wilson flow: M. Luscher, arXiv:1006.4518. Action density: qed3_v2-6.pdf Sec. IV,
// Eqs. (IV.24)-(IV.35) via action_ext_claude.h density_Ylm_spatial at l=m=0
// (Y_00 = 1/sqrt(4pi); multiply back by sqrt(4pi) -> plain average density).
//
// PURPOSE: (i) inspect typical trajectories to set c in the flow-scale condition
// t0^2 E_s(t0) = c (c = 0.3 first guess); (ii) the MC series of E_s at (later-chosen) t0
// is the thermalization monitor. Smooth E_s(t) is reconstructed from the saved points.
//
// Output (TEXT for this first check; production h5 layout = planned later):
// data_<ensemble-basename>/therm_series_claude.dat, one line per config:
//   k  E_s(0.0)  E_s(0.2)  ...  E_s(2.0)      (NREC+1 = 11 values)
// Resumable: existing k lines are skipped on re-run (append mode).
// Usage: ./a.out [gsq] [Nf] [nu0] [kmax_run] [kmin] [stride] [ens_dir]
//   ens_dir = explicit sea-config dir (REQUIRED for the redo `..._hb*` ensembles whose names
//   the legacy dir3 construction cannot form; omitted -> legacy both_3d-style dir3).
// The therm cut kmin per ensemble is then read off the series plateau (blackboard
// redo_ensembles_claude.txt).

using Double = double;
using Idx = std::int32_t;
using Complex = std::complex<double>;

using Link = std::array<Idx,2>; // <int,int>;
using Face = std::vector<Idx>;

using MS=Eigen::Matrix2cd;
using VD=Eigen::Vector2d;
using VE=Eigen::Vector3d;
using VC=Eigen::VectorXcd;

static constexpr int NS = 2;
static constexpr int DIM = 2;
static constexpr Complex I = Complex(0.0, 1.0);


// #define IS_DUAL
#define IS_OVERLAP

namespace Comp{
  constexpr bool is_compact=false;

  // d_DW.update() is always done independently
#ifdef IS_OVERLAP
  constexpr int NPARALLEL_DUPDATE=1;
  constexpr int NPARALLEL=1;
  constexpr int NSTREAMS=1;
#else
  constexpr int NPARALLEL_DUPDATE=12;
  constexpr int NPARALLEL=1;
  constexpr int NSTREAMS=12;
#endif
  constexpr int NPARALLEL_GAUGE=1;
  constexpr int NPARALLEL_SORT=1;

  // Refinement level L is a compile-time flag: -DN_REFINE_CLI=2 (or 4) for L=2,4.
  // Same #ifndef-guard convention as the L2/L4 stochastic Y_lm drivers.
#ifdef N_REFINE_CLI
  constexpr int N_REFINE=N_REFINE_CLI;
#else
  constexpr int N_REFINE=1;
#endif
  constexpr int NS=2;

  constexpr int Nt=128;

#ifdef IS_DUAL
  constexpr Idx N_SITES=20*N_REFINE*N_REFINE;
#else
  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
#endif
  constexpr int N_LINKS=30*N_REFINE*N_REFINE; // 30, 120, 480

  constexpr Idx Nx=NS*N_SITES; // matrix size of DW
  constexpr Idx N=Nx*Nt; // matrix size of DW

  const double TOL_INNER=1.0e-15;
  const double TOL_OUTER=1.0e-14;
}

const std::string dir = "../../geometry/data/";


#include "timer.h"

#include "s2n_simp.h"
#include "s2n_dual.h"
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
#include "action_ext_claude.h"   // density_Ylm_spatial / density_Ylm_temporal (l=m=0 -> avg density)

// SWITCHED to per-timeslice SPATIAL flow (NM 2026-07-17; matches the glueball convention).
// flow_claude.h default = get_spatial when FLOW_FULL is NOT defined. LO scaling unchanged:
// E_s(t) ~ g^2/t^{3/2} and t_phys ~ a_t t/g^2 (spatial force carries beta_s a^2 = a_t/g^2).
// // this driver ALWAYS flows in full 3D (no build-flag dependence)
// #define FLOW_FULL
#include "flow_claude.h"


// average spatial (magnetic) action density per face, timeslice-averaged:
// (1/Nt) sum_s sqrt(4pi) * density_Ylm_spatial(U, s, 0, 0)
template<typename Action, typename Gauge>
static double avg_density_spatial( const Action& SW, const Gauge& U ){
  const double sqrt4pi = std::sqrt( 4.0*M_PI );
  double res = 0.0;
  for(int s=0; s<U.Nt; s++){
    res += SW.density_Ylm_spatial( U, s, 0, 0 );
  }
  res *= sqrt4pi/U.Nt;
  return res;
}

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) == "-h") {
      printf("Usage: ./a.out [gsq] [Nf] [nu0] [kmax_run] [kmin] [stride] [ens_dir]\n");
      printf("  gsq      Wilson coupling squared (default: 8.0)\n");
      printf("  Nf       number of fermion flavors (default: 2)\n");
      printf("  nu0      mass parameter (default: 1.0)\n");
      printf("  kmax_run cap on configs processed (default: all)\n");
      printf("  kmin     first config index (default: 1)\n");
      printf("  stride   config stride (default: 1)\n");
      printf("  ens_dir  explicit config dir (redo _hb* names; omitted -> legacy dir3)\n");
      printf("  tmax     max flow time (default 2.0, common to all ensembles; eps=0.01, save every step)\n");
      return 0;
    }
  }

  double gsq = 8.0;
  if(argc>1) gsq = atof(argv[1]);
  int Nf = 2;
  if(argc>2) Nf = atoi(argv[2]);
  double nu0 = 1.0;
  if(argc>3) nu0 = atof(argv[3]);
  int kmax_run = 100000;
  if(argc>4) kmax_run = atoi(argv[4]);
  int kmin = 1;
  if(argc>5) kmin = atoi(argv[5]);
  int stride = 1;
  if(argc>6) stride = atoi(argv[6]);
  std::string ens_dir = "";
  if(argc>7) ens_dir = argv[7];
  double tmax = 2.0;
  if(argc>8) tmax = atof(argv[8]);
  std::cout << "# gsq = " << gsq << " Nf = " << Nf << " nu0 = " << nu0
            << " kmax_run = " << kmax_run << " kmin = " << kmin << " stride = " << stride
            << " ens_dir = " << (ens_dir.empty() ? std::string("<legacy dir3>") : ens_dir)
            << " tmax = " << tmax << std::endl;

  // ---------------------------------------

  constexpr Idx N = Comp::N;
  constexpr int Nt = Comp::Nt;

#ifdef IS_DUAL
  using Base=S2Trivalent;
#else
  using Base=S2Simp;
#endif

  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  using Action=U1WilsonExt<Base>;

  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;

  // ----------------------
  const double at = 0.2; // production a_t
  if(Nt!=1) assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);

  Action SW( gsq, at, base );
  std::cout << "# alat = " << base.mean_ell << std::endl;

  std::string dir3, dir4;
  if(!ens_dir.empty()){
    // explicit ensemble dir (redo production names Nf*_..._hb*); output shares the
    // data_<basename>/ dir with the fermionic obs (plan Sec 5).
    dir3 = ens_dir;
    if(dir3.back()!='/') dir3 += "/";
    std::string bn = dir3.substr(0, dir3.size()-1);
    const std::size_t pos = bn.find_last_of('/');
    if(pos!=std::string::npos) bn = bn.substr(pos+1);
    dir4 = "data_"+bn+"/";
    std::cout << "dir3 = " << dir3 << std::endl;
  }
  else if(Nf==0){
    dir3="gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
    dir4="data_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
    std::cout << "dir3 = " << dir3 << std::endl;
  }
  else{
    dir3="Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
    dir4="data_Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
  }
  std::filesystem::create_directory(dir4);

  Gauge U(base);

  const int kmax=1e5;

  int k_tmp=0;
  {
    // scan from kmin in steps of `stride` so non-contiguously-numbered ensembles are detected
    for(k_tmp=kmin; k_tmp<=kmax; k_tmp+=stride ){
      const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k_tmp);
      const bool bool_lat = std::filesystem::exists(str_lat);
      if(!bool_lat) break;
    }
    k_tmp -= stride;
  }
  if(k_tmp > kmax_run) k_tmp = kmax_run;
  std::cout << "# last config k = " << k_tmp << std::endl;

  // ---- trajectory recording (NM 2026-07-17 v6): FIXED eps = 0.01, save EVERY step
  //      (Dt_seg = 0.01); tmax = runtime arg (default 2.0), COMMON to all ensembles ----
  constexpr double FLOW_EPS = 0.01;
  constexpr int SAVE_EVERY = 1;
  const int nrec = (int)std::lround( tmax/(FLOW_EPS*SAVE_EVERY) );   // 200 for tmax = 2.0
  const Flow<Action> flow_seg( &SW, FLOW_EPS*SAVE_EVERY, SAVE_EVERY );

  // ---- resumable text output: collect k's already present ----
  // v3 = per-timeslice SPATIAL flow, r-scheme tmax (21 cols after k).
  // v2 (FULL 3D flow, same grid) and v1 (full flow, 0.2 grid, 11 cols) are superseded.
  // const std::string dat_path = dir4+"therm_series_claude.dat";
  // const std::string dat_path = dir4+"therm_series_v2_claude.dat";
  const std::string dat_path = dir4+"therm_series_v6_claude.dat";
  std::set<int> k_done;
  if(std::filesystem::exists(dat_path)){
    std::ifstream fin(dat_path);
    std::string line;
    while(std::getline(fin, line)){
      if(line.empty()) continue;
      if(line[0]=='#') continue;
      std::istringstream iss(line);
      int kk;
      if(iss >> kk) k_done.insert(kk);
    }
  }

  std::ofstream fout(dat_path, std::ios::app);
  fout << std::scientific << std::setprecision(15);
  if(k_done.empty()){
    fout << "# tmax = " << tmax << std::endl;   // parseable: analysis reconstructs tlist from this
    fout << "# k  E_s(t_fl) for t_fl = 0.0, " << (tmax/nrec) << ", ..., " << tmax << "  (" << (nrec+1) << " columns after k)" << std::endl;
    fout << "# E_s = SPATIAL avg action density (noncompact 0.5*beta_s*theta^2 integrand, per face," << std::endl;
    fout << "#       glueball-code weight); flow = per-timeslice SPATIAL gradient flow (get_spatial)" << std::endl;
    fout << "#       (Luscher arXiv:1006.4518), eps=" << FLOW_EPS << ", save every " << SAVE_EVERY << " steps." << std::endl;
    fout << "# flow-scale condition (analysis, 3D): t0^(3/2) E_s(t0) = c (c to be tuned)." << std::endl;
  }

  // serial over configs k; parallelism is ensemble-level (one process per ensemble).
  for(int k=kmin; k<=k_tmp; k+=stride ){
    if(k_done.count(k)!=0) continue;
    const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k);
    U.read( str_lat );

    std::vector<double> Es( nrec+1, 0.0 );
    Es[0] = avg_density_spatial( SW, U );

    Gauge Uflow = U;
    for(int i=1; i<=nrec; i++){
      flow_seg(Uflow);                     // t_fl = i * FLOW_EPS*SAVE_EVERY
      Es[i] = avg_density_spatial( SW, Uflow );
    }

    fout << k;
    for(int i=0; i<nrec+1; i++){
      fout << " " << Es[i];
    }
    fout << std::endl;   // endl = flush -> line-complete on disk per config
    std::cout << "# k = " << k << " done" << std::endl;
  }

  return 0;
}
