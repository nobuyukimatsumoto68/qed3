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
#include <highfive/H5File.hpp>

// glue_flow_weighted_claude.cu  (copy of glue_therm_flow_claude.cu; original untouched)
//
// LINK-WEIGHTED spatial Wilson flow, as defined in NM's paper note (Sec. IX.A.2, 2026-10-26):
//   d theta_l / d \hat\tau = - (\ell_l^2/(2 A_l)) d \hat E_t / d theta_l ,
//   \hat E_t = sum_{f in S^2(t)} (\bar a_s^2/A_f) (1/2) theta_f^2 ,   \bar a_s = base.mean_ell ,
// spatial links only, each timeslice independently, naive explicit Euler (includes/flow_weighted_claude.h).
// Continuum flow time tau = \bar a_s^2 \hat\tau. No coupling and no a_t enter the flow.
// Derivation: final/analysis_wilson/wilson_flow_notes_claude.md Secs. 6-7;
// plan: final/analysis_wilson/wilson_flow_impl_plan_claude.md.
// Wilson flow: M. Luscher, arXiv:1006.4518. Link weight = discrete Hodge star on 1-forms,
// M. Desbrun, A. N. Hirani, M. Leok, J. E. Marsden, arXiv:math/0508341.
//
// RECORDING (2026-10-26, NM): per-L tmax = R^2/(2 \bar a_s^2) (R = 1; flow radius sqrt(8 tau) = 2R at tmax),
//   200 Luscher-RK3 steps, eps = tmax/200, save EVERY step -> tlist = 0, eps, ..., tmax (201 points).
//   tmax: L1 0.408, L2 1.432, L3 3.158, L4 5.575; eps: L1 0.0020 ... L4 0.0279.
//   (Drafts: uniform tmax = 4.0; then 400 Euler steps.)
//
// H5 LAYOUT (data_<ensemble-basename>/flow_weighted_claude.h5; one file per ensemble):
//   tlist  [n_t]        lattice flow times \hat\tau (includes 0)
//   E      [n_t][n_cfg] <\hat E_t>_t (timeslice average) trajectories, config = FAST axis
//   klist  [n_cfg]      config indices, sorted ascending (columns of E)
//   eps, save_every, tmax, abar   scalars (grid metadata; abar = \bar a_s used in \hat E_t)
// RESUME: existing file is READ first; only configs missing from klist are flowed and merged
// (Wilson flow is the expensive part -- never recomputed). The file is REWRITTEN after EVERY
// new config via write-to-tmp + atomic rename (crash loses only the config in flight).
// GRID GUARD: if the existing file's tlist does not match the current (eps, tmax), the driver
// ABORTS rather than mix grids (remove the old file deliberately to change protocol).
//
// Usage: ./a.out [gsq] [Nf] [nu0] [kmax_run] [kmin] [stride] [ens_dir] [tmax]
//   ens_dir = explicit sea-config dir (redo `..._hb*` names; omitted -> legacy dir3).
//   Production: stride 1, kmin 1 (no thermalization cut at measurement; cut in analysis).

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
#include "action_ext_claude.h"   // density_Ylm_spatial (l=m=0 -> avg density)

// per-timeslice SPATIAL flow, OLD (beta-weighted action, no link weight):
// #include "flow_claude.h"
// per-timeslice SPATIAL flow, NEW (link-weighted, coupling-free \hat E_t):
#include "flow_weighted_claude.h"


// rewrite the per-ensemble h5 from the in-memory (k -> trajectory) map; atomic via tmp+rename
static void write_h5( const std::string& h5path,
                      const std::vector<double>& tlist,
                      const std::map<int, std::vector<double>>& data,
                      const double eps,
                      const int save_every,
                      const double tmax,
                      const double abar ){
  const std::string tmp_path = h5path + ".tmp";
  {
    HighFive::File h5( tmp_path, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate );
    std::vector<int> klist;
    for(const auto& kv : data){
      klist.push_back( kv.first );
    }
    const int n_t = (int)tlist.size();
    const int n_cfg = (int)klist.size();
    std::vector<std::vector<double>> E( n_t, std::vector<double>(n_cfg, 0.0) );
    for(int j=0; j<n_cfg; j++){
      const std::vector<double>& col = data.at( klist[j] );
      for(int it=0; it<n_t; it++){
        E[it][j] = col[it];
      }
    }
    h5.createDataSet( "tlist", tlist );
    h5.createDataSet( "klist", klist );
    h5.createDataSet( "E", E );
    h5.createDataSet( "eps", std::vector<double>{eps} );
    h5.createDataSet( "save_every", std::vector<int>{save_every} );
    h5.createDataSet( "tmax", std::vector<double>{tmax} );
    h5.createDataSet( "abar", std::vector<double>{abar} );
  }
  std::filesystem::rename( tmp_path, h5path );
}


// Derive a_t from the ensemble-dir string (parse "at<digits>"); returns -1 if absent.
// Same helper as glue2_msm_shapes_claude.cu (09-14 fix), so at0.1 configs are flowed with
// beta_s(0.1) instead of the old at=0.2 hardcode.
static double at_from_ensdir(const std::string& s){
  for(std::size_t p = 0; p + 2 < s.size(); ++p){
    const bool is_at = (s[p] == 'a' && s[p+1] == 't');
    const bool digit_after = (s[p+2] >= '0' && s[p+2] <= '9');
    if(is_at && digit_after) return std::stod(s.substr(p+2));
  }
  return -1.0;
}


int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) == "-h") {
      printf("Usage: ./a.out [gsq] [Nf] [nu0] [kmax_run] [kmin] [stride] [ens_dir] [tmax]\n");
      printf("  gsq      Wilson coupling squared (default: 8.0)\n");
      printf("  Nf       number of fermion flavors (default: 2)\n");
      printf("  nu0      mass parameter (default: 1.0)\n");
      printf("  kmax_run cap on configs processed (default: all)\n");
      printf("  kmin     first config index (default: 1)\n");
      printf("  stride   config stride (default: 1)\n");
      printf("  ens_dir  explicit config dir (redo _hb* names; omitted -> legacy dir3)\n");
      printf("  tmax     max lattice flow time \\hat\\tau (<=0 or omitted: R^2/(2 abar^2), R=1; RK3, eps=tmax/200, save every step)\n");
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
  // double tmax = 4.0;   // first draft: uniform tmax
  double tmax = -1.0;     // <= 0 -> per-L auto value R^2/(2 abar^2), set after the lattice is built
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

  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;

  // ----------------------
  // a_t does NOT enter the link-weighted flow; it is only used for the legacy dir names below.
  double at = 0.2;
  if(!ens_dir.empty()){
    const double at_ens = at_from_ensdir(ens_dir);
    if(at_ens > 0.0) at = at_ens;
  }
  std::cout << "# at = " << at << " (from ens_dir; not used by the flow)" << std::endl;
  std::cout << "# alat = abar = " << base.mean_ell << std::endl;
  if(tmax <= 0.0){
    const double R = 1.0;
    tmax = R*R / (2.0*base.mean_ell*base.mean_ell);   // flow radius sqrt(8 tau) = 2R at tmax
    std::cout << "# tmax (auto, R^2/(2 abar^2)) = " << tmax << std::endl;
  }

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

  // ---- recording grid (2026-10-26, NM): 200 Luscher-RK3 steps to the per-L tmax = R^2/(2 abar^2)
  //      (eps = tmax/200: L1 0.0020 ... L4 0.0279), save every step -> 201 points ----
  // constexpr int FLOW_NSTEP = 400;   // Euler draft
  constexpr int FLOW_NSTEP = 200;
  constexpr int SAVE_EVERY = 1;
  const double flow_eps = tmax/FLOW_NSTEP;
  const int nrec = FLOW_NSTEP/SAVE_EVERY;
  const FlowWeighted<Base> flow( base, flow_eps );
  std::cout << "# flow: link-weighted, Luscher RK3, eps = " << flow_eps << " nsteps = " << FLOW_NSTEP
            << " tmax = " << tmax << " abar = " << flow.abar << std::endl;

  std::vector<double> tlist( nrec+1, 0.0 );
  for(int i=0; i<=nrec; i++){
    tlist[i] = i*flow_eps*SAVE_EVERY;
  }

  // ---- read existing h5 (resume): flow ONLY configs missing from klist ----
  const std::string h5path = dir4+"flow_weighted_claude.h5";
  std::map<int, std::vector<double>> data;
  if(std::filesystem::exists(h5path)){
    HighFive::File fin( h5path, HighFive::File::ReadOnly );
    std::vector<double> tlist_in;
    fin.getDataSet("tlist").read(tlist_in);
    bool grid_ok = ( tlist_in.size() == tlist.size() );
    if(grid_ok){
      for(std::size_t i=0; i<tlist.size(); i++){
        if(std::abs(tlist_in[i]-tlist[i]) > 1.0e-12) grid_ok = false;
      }
    }
    if(!grid_ok){
      std::cerr << "### ERROR: existing " << h5path << " has a DIFFERENT flow grid (n_t="
                << tlist_in.size() << " vs " << tlist.size()
                << "). Refusing to mix; remove the old file deliberately to change protocol." << std::endl;
      return 1;
    }
    std::vector<int> klist_in;
    fin.getDataSet("klist").read(klist_in);
    std::vector<std::vector<double>> E_in;
    fin.getDataSet("E").read(E_in);
    for(std::size_t j=0; j<klist_in.size(); j++){
      std::vector<double> col( tlist.size(), 0.0 );
      for(std::size_t it=0; it<tlist.size(); it++){
        col[it] = E_in[it][j];
      }
      data[ klist_in[j] ] = col;
    }
    std::cout << "# resume: " << data.size() << " configs already in " << h5path << std::endl;
  }

  // serial over configs k; parallelism is ensemble-level (one process per ensemble).
  for(int k=kmin; k<=k_tmp; k+=stride ){
    if(data.count(k)!=0) continue;
    const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k);
    U.read( str_lat );

    std::vector<double> Es( nrec+1, 0.0 );
    Es[0] = flow.energy_avg( U );

    Gauge Uflow = U;
    for(int i=1; i<=nrec; i++){
      // flow.step_euler( Uflow );         // Euler draft
      flow.step_rk3( Uflow );              // \hat\tau = i * eps
      Es[i] = flow.energy_avg( Uflow );
    }

    data[k] = Es;
    write_h5( h5path, tlist, data, flow_eps, SAVE_EVERY, tmax, flow.abar );
    std::cout << "# k = " << k << " done (" << data.size() << " configs in h5)" << std::endl;
  }

  std::cout << "# finished: " << data.size() << " configs in " << h5path << std::endl;
  return 0;
}
