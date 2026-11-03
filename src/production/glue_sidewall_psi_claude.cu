// glue_sidewall_psi_claude.cu -- j_top Psi component (side-wall flux) at l=1,2, per config -> HDF5.
// Copy of glue2_msm_shapes_claude.cu (original untouched); the spatial-shape operators are replaced by
// the SIDE-WALL (temporal-plaquette) flux of j_top, paper-21 Eq. (V.14) with p = rectangle:
//   J_l(t) = U.plaquette_angle(t, link)   (noncompact -> linear; = E-flux along the spatial link l)
//   w_l    = link_volume[l] / ell[l]^2    (proportional to beta_t[l] of the gauge action = lattice Hodge metric)
//   A (exact gradient)  O_lm(t) = \sum_l w_l J_l(t) [ Y_lm(x1) - Y_lm(x0) ]  = \sum_x Y_lm(x) (dS_G/dA_t)(x)/c
//   B (face direction)  O_lm(t) = \sum_l w_l J_l(t) ell_l \Psi_lm(m_l) . N_l,  N_l = m_l \times (x1-x0)/|x1-x0|
//                        (\Psi.N = \nabla_{tan} Y . (x1-x0)/|x1-x0| at the link midpoint m_l; NM: 3D normal)
// A is exactly orthogonal (w metric) to the co-exact part of E fixed by the triangle circulations (Bianchi),
// i.e. no leakage from the radial F_12 tower; B = A + O(a^2).  Derivation + plan:
// final/analysis_gluonic/jtop_psi_ell1_note_claude.md, final/analysis_gluonic/jtop_psi_impl_plan_claude.md.
// Discrete Hodge decomposition: Desbrun, Hirani, Leok, Marsden, arXiv:math/0508341.
// Smearing levels (stored as the "shape" axis, ishape = ilev*N_VAR + ivar, flow-major like msm):
//   ilev 0 raw | 1 spatial flow tau=2.0 (get_spatial, = glue_msm_shapes) | 2..4 FULL 3D flow (get_force)
//   tau = 0.02, 0.05, 0.1 (cumulative).  Wilson/gradient flow: M. Luescher, arXiv:1006.4518.
// Output layout identical to glue_msm_shapes (F_corr_blk, F, n_lm, n_shapes, ell, em, complete) so
// glue_gevp_analysis_t0_claude.o reads it unchanged; prefix glue_sidewall_psi.<k>.h5.
#include <typeinfo>
#include <iostream>
#include <iomanip>
#include <fstream>
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
#include <random>
#include <Eigen/Dense>
#include <highfive/H5File.hpp>   // per-config HDF5 output (fast binary read for analysis)

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


// #define IS_FLOW

// #define IS_DUAL
#define IS_OVERLAP
// #define IS_DAGGER
// #undef _OPENMP

// #define GAUGE_TRSF


namespace Comp{
  constexpr bool is_compact=false;

  // d_DW.update() is always done independently
#ifdef IS_OVERLAP
  constexpr int NPARALLEL_DUPDATE=1;
  // constexpr int NPARALLEL=12; // 12
  // constexpr int NSTREAMS=4; // 4
  constexpr int NPARALLEL=1; // 12
  constexpr int NSTREAMS=1; // 4
#else
  constexpr int NPARALLEL_DUPDATE=12;
  constexpr int NPARALLEL=1; // 12
  constexpr int NSTREAMS=12; // for grad loop
#endif
  constexpr int NPARALLEL_GAUGE=1; // 12
  constexpr int NPARALLEL_SORT=1; // 12

  // Refinement level L is a compile-time flag: -DN_REFINE_CLI=2 (or 4) for L=2,4.
#ifdef N_REFINE_CLI
  constexpr int N_REFINE=N_REFINE_CLI;
#else
  constexpr int N_REFINE=1;
#endif
  constexpr int NS=2;

  // constexpr int Nt=24;
  // constexpr int Nt=48; // add 2
  // constexpr int Nt=64;
  // constexpr int Nt=96; // add 4
  // constexpr int Nt=120;
  // constexpr int Nt=144; // add 8
  // constexpr int Nt=168;

  // constexpr int Nt=24;
  // constexpr int Nt=192;
  constexpr int Nt=128;
  // constexpr int Nt=16;

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

// const std::string dir = "../../dats/";
const std::string dir = "../../geometry/data/";


// // #define IsVerbose
// #define IsVerbose2
// // #define InfoForce
// #define InfoDelta

#include "timer.h"

// #include "../../integrator/geodesic.h"

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
#include "action_ext.h"

// ======================================

#include "sparse_matrix.h"
#include "dirac_base.h"
#include "dirac_simp.h"
#include "dirac_dual.h"
#include "dirac_ext.h"
// // #include "pseudofermion.h"
// #include "dirac.h"

// #include "sparse_dirac.h"
#include "sparse_dirac_claude.h"   // O(len) bucketing CSR build (was O(N*len)); -DCSR_VERIFY to check
// #include "matpoly.h"
#include "matpoly_claude.h"

#include "dirac_pf.h"
#include "overlap.h"

// #include "flow.h"
#include "flow_claude.h"   // adds -DFLOW_FULL switch (full 3D flow via get_force vs spatial-only)

#include "icos_orbits_claude.h"    // native Ih orbit table on the lattice
#include "wilson_shapes_claude.h"  // generic shape (triangle/rectangle) orbit operators


// TODO: Cusparse for SparseMatrix::act_gpu, probably defining handle in matpoly.h
// make 2 streams in V Vdag in square in Overlap
// all the operation on GPU in Overlap::operator()
// pseudofermion
// 3d
// __m256 to vectorize with AVX2


// Derive a_t from the ensemble-dir string (parse "at<digits>"); returns -1 if absent.
// Mirrors the fermion valence fix at_from_ensdir (jj_local_ylm_scalar_conn_stoch_fnal_claude.cu),
// so the Wilson-flow smearing uses the ensemble's true a_t (0.1 at0.1) instead of a hardcode.
static double at_from_ensdir(const std::string& s){
  for(std::size_t p = 0; p + 2 < s.size(); ++p){
    const bool is_at = (s[p] == 'a' && s[p+1] == 't');
    const bool digit_after = (s[p+2] >= '0' && s[p+2] <= '9');
    if(is_at && digit_after) return std::stod(s.substr(p+2));
  }
  return -1.0;
}

// Explicit Euler gradient-flow steps, written here (not via flow_claude.h, whose spatial/full choice is a
// compile-time switch) so ONE binary produces both smearings.  Same update as Flow::operator().
template<typename Action, typename Gauge>
static void flow_spatial_steps(const Action& SW, Gauge& U, const double tau, const int nsteps){
  Gauge dSg(U.lattice);
  for(int i=0; i<nsteps; i++){
    SW.get_spatial( dSg, U );   // spatial-only (per-timeslice), temporal links untouched
    U += -tau * dSg;
  }
}

template<typename Action, typename Gauge>
static void flow_full_steps(const Action& SW, Gauge& U, const double tau, const int nsteps){
  Gauge dSg(U.lattice);
  for(int i=0; i<nsteps; i++){
    SW.get_force( dSg, U );     // full 3D flow (spatial + temporal links)
    U += -tau * dSg;
  }
}

// Side-wall projections at every smearing level.  obs[(ilev*n_var + ivar)*n_lm + ilm][t].
// Levels: kind 0 = raw U, 1 = spatial flow (cumulative on its own copy), 2 = full 3D flow (cumulative).
template<typename Action, typename Gauge>
static void measure_sidewall(const Action& SW,
                             const Gauge& U,
                             const std::vector<std::vector<std::vector<double>>>& coef,
                             const std::vector<Link>& link_of,
                             const double* lev_tau,
                             const int* lev_kind,
                             const int n_lev,
                             const double sp_dt,
                             const double full_dt,
                             std::vector<std::vector<double>>& obs){
  const int n_var = coef.size();
  const int n_lm = coef[0].size();
  const Idx n_links = link_of.size();
  const int Nt = Comp::Nt;
  obs.assign( n_lev*n_var*n_lm, std::vector<double>(Nt, 0.0) );

  Gauge Usp = U;     // spatial-flow copy
  Gauge Ufull = U;   // full-3D-flow copy
  double tau_sp = 0.0;
  double tau_full = 0.0;
  std::vector<double> J( n_links );
  for(int ilev=0; ilev<n_lev; ilev++){
    const Gauge* Um = &U;
    if(lev_kind[ilev]==1){
      const int nst = (int)((lev_tau[ilev] - tau_sp)/sp_dt + 0.5);
      flow_spatial_steps( SW, Usp, sp_dt, nst );
      tau_sp = lev_tau[ilev];
      Um = &Usp;
    }
    else if(lev_kind[ilev]==2){
      const int nst = (int)((lev_tau[ilev] - tau_full)/full_dt + 0.5);
      flow_full_steps( SW, Ufull, full_dt, nst );
      tau_full = lev_tau[ilev];
      Um = &Ufull;
    }

    for(int t=0; t<Nt; t++){
      for(Idx il=0; il<n_links; il++) J[il] = Um->plaquette_angle( t, link_of[il] );
      for(int ivar=0; ivar<n_var; ivar++){
        const int ishape = ilev*n_var + ivar;
        for(int ilm=0; ilm<n_lm; ilm++){
          double sum = 0.0;
          for(Idx il=0; il<n_links; il++) sum += coef[ivar][ilm][il] * J[il];
          obs[ishape*n_lm+ilm][t] = sum;
        }
      }
    }
  }
}

// -DSIDEWALL_TEST: checks on the first config (stdout only, no h5):
//  (a) Gauss law: c * A_raw == \sum_x Y_lm(x) (dS_G/dA_t)(x,t) from Action::get_force, c = 2/(a_t g^2)
//  (b) gauge invariance: random noncompact gauge transformation -> every level/variant unchanged
//  (c) A vs B per level: ratio and correlation over t (1 - corr = leakage measure, expect O(a^2))
//  (d) full 3D flow: action monotone in tau; dt convergence (dt vs dt/2) and dt = 10*dt (stability)
//  (e) temporal smearing radius sigma_t ~ sqrt(2 beta_t tau) [slices] per 3D level
template<typename Action, typename Gauge, typename BaseT>
static void run_sidewall_tests(const Action& SW,
                               const Gauge& U,
                               const BaseT& base,
                               const std::vector<std::vector<std::vector<double>>>& coef,
                               const std::vector<Link>& link_of,
                               const std::vector<std::array<int,2>>& lm_set,
                               const double* lev_tau,
                               const int* lev_kind,
                               const int n_lev,
                               const double sp_dt,
                               const double full_dt,
                               const std::vector<std::vector<double>>& obs0,
                               const double gsq,
                               const double at){
  const int n_var = coef.size();
  const int n_lm = coef[0].size();
  const Idx n_links = link_of.size();
  const int Nt = Comp::Nt;
  const double cnorm = 2.0/(at*gsq);   // beta_t[l] = cnorm * w_l
  std::cout << "# ===== SIDEWALL_TEST (gsq = " << gsq << ", at = " << at << ", n_links = " << n_links << ") =====" << std::endl;

  // (a) Gauss law
  {
    Gauge dSg(U.lattice);
    SW.get_force( dSg, U );
    double maxdiff = 0.0;
    double maxabs = 0.0;
    for(int t=0; t<Nt; t++){
      for(int ilm=0; ilm<n_lm; ilm++){
        double rhs = 0.0;
        for(Idx ix=0; ix<base.n_sites; ix++){
          const VE r = base.sites[ix] / base.sites[ix].norm();
          rhs += Ylm_real( lm_set[ilm][0], lm_set[ilm][1], r ) * dSg.tp( t, ix );
        }
        const double lhs = cnorm * obs0[(0*n_var + 0)*n_lm + ilm][t];
        maxdiff = std::max( maxdiff, std::abs(lhs - rhs) );
        maxabs = std::max( maxabs, std::abs(rhs) );
      }
    }
    std::cout << "# (a) Gauss law: max|c*A_raw - sum_x Y dS/dA_t| = " << maxdiff << "  (scale max|rhs| = " << maxabs << ")" << std::endl;
  }

  // (b) gauge invariance
  {
    std::mt19937 gen(12345);
    std::uniform_real_distribution<double> dist(-M_PI, M_PI);
    std::vector<std::vector<double>> lam( Nt, std::vector<double>(base.n_sites, 0.0) );
    for(int s=0; s<Nt; s++) for(Idx ix=0; ix<base.n_sites; ix++) lam[s][ix] = dist(gen);
    Gauge Ug = U;
    for(int s=0; s<Nt; s++){
      for(Idx il=0; il<n_links; il++){
        const Link& lk = link_of[il];
        const double sign = base.map2sign.at(lk);
        Ug.sp( s, il ) += sign * ( lam[s][lk[1]] - lam[s][lk[0]] );
      }
      for(Idx ix=0; ix<base.n_sites; ix++) Ug.tp( s, ix ) += lam[(s+1)%Nt][ix] - lam[s][ix];
    }
    std::vector<std::vector<double>> obsg;
    measure_sidewall( SW, Ug, coef, link_of, lev_tau, lev_kind, n_lev, sp_dt, full_dt, obsg );
    for(int ilev=0; ilev<n_lev; ilev++){
      double maxdiff = 0.0;
      double maxabs = 0.0;
      for(int ivar=0; ivar<n_var; ivar++){
        for(int ilm=0; ilm<n_lm; ilm++){
          const int op = (ilev*n_var + ivar)*n_lm + ilm;
          for(int t=0; t<Nt; t++){
            maxdiff = std::max( maxdiff, std::abs(obsg[op][t] - obs0[op][t]) );
            maxabs = std::max( maxabs, std::abs(obs0[op][t]) );
          }
        }
      }
      std::cout << "# (b) gauge inv. level " << ilev << " (kind " << lev_kind[ilev] << ", tau " << lev_tau[ilev]
                << "): max|diff| = " << maxdiff << "  scale " << maxabs << std::endl;
    }
  }

  // (c) A vs B
  for(int ilev=0; ilev<n_lev; ilev++){
    for(int ilm=0; ilm<n_lm; ilm++){
      const int opA = (ilev*n_var + 0)*n_lm + ilm;
      const int opB = (ilev*n_var + 1)*n_lm + ilm;
      double aa = 0.0;
      double bb = 0.0;
      double ab = 0.0;
      for(int t=0; t<Nt; t++){
        aa += obs0[opA][t] * obs0[opA][t];
        bb += obs0[opB][t] * obs0[opB][t];
        ab += obs0[opA][t] * obs0[opB][t];
      }
      std::cout << "# (c) A vs B level " << ilev << " (l,m)=(" << lm_set[ilm][0] << "," << lm_set[ilm][1]
                << "): <AB>/<AA> = " << ab/aa << "  1-corr = " << 1.0 - ab/std::sqrt(aa*bb) << std::endl;
    }
  }

  // (d) full 3D flow: action along tau, dt convergence, large-dt stability
  {
    Gauge Uf = U;
    double tau = 0.0;
    std::cout << "# (d) S_G(raw) = " << SW( U ) << std::endl;
    for(int ilev=0; ilev<n_lev; ilev++){
      if(lev_kind[ilev]!=2) continue;
      const int nst = (int)((lev_tau[ilev] - tau)/full_dt + 0.5);
      flow_full_steps( SW, Uf, full_dt, nst );
      tau = lev_tau[ilev];
      std::cout << "# (d) full flow tau = " << tau << " (dt " << full_dt << "): S_G = " << SW( Uf ) << std::endl;
    }
    Gauge Uh = U;
    const int nst_h = (int)(tau/(0.5*full_dt) + 0.5);
    flow_full_steps( SW, Uh, 0.5*full_dt, nst_h );
    double maxdiff = 0.0;
    double maxabs = 0.0;
    for(int t=0; t<Nt; t++){
      for(Idx il=0; il<n_links; il++){
        const double a = Uf.plaquette_angle( t, link_of[il] );
        const double b = Uh.plaquette_angle( t, link_of[il] );
        maxdiff = std::max( maxdiff, std::abs(a - b) );
        maxabs = std::max( maxabs, std::abs(a) );
      }
    }
    std::cout << "# (d) dt vs dt/2 at tau = " << tau << ": max|dJ| = " << maxdiff << "  scale " << maxabs << std::endl;
    Gauge Ub = U;
    const int nst_b = (int)(tau/(10.0*full_dt) + 0.5);
    flow_full_steps( SW, Ub, 10.0*full_dt, nst_b );
    std::cout << "# (d) dt x10 (" << 10.0*full_dt << ") at tau = " << tau << ": S_G = " << SW( Ub ) << "  (blow-up => dt too large)" << std::endl;
  }

  // (e) temporal smearing radius
  {
    double bt = 0.0;
    double bmax = 0.0;
    for(Idx il=0; il<n_links; il++){
      const double b = cnorm * base.link_volume[il] / std::pow(base.ell[il], 2);
      bt += b;
      bmax = std::max( bmax, b );
    }
    bt /= n_links;
    std::cout << "# (e) beta_t mean = " << bt << "  max = " << bmax << std::endl;
    for(int ilev=0; ilev<n_lev; ilev++){
      if(lev_kind[ilev]!=2) continue;
      std::cout << "# (e) tau = " << lev_tau[ilev] << ": sigma_t ~ sqrt(2 beta_t tau) = " << std::sqrt(2.0*bt*lev_tau[ilev]) << " slices" << std::endl;
    }
  }
  std::cout << "# ===== SIDEWALL_TEST done =====" << std::endl;
}

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) == "-h") {
      printf("Usage: ./a.out [gsq] [Nf] [nu0]\n");
      printf("  gsq  Wilson coupling squared (default: 8.0)\n");
      printf("  Nf   number of fermion flavors (default: 2)\n");
      printf("  nu0  mass parameter (default: 1.0)\n");
      return 0;
    }
  }

  double gsq = 8.0;
  if(argc>1) gsq = atof(argv[1]);
  int Nf = 2;
  if(argc>2) Nf = atoi(argv[2]);
  double nu0 = 1.0;
  if(argc>3) nu0 = atof(argv[3]);
  // config selection (added _claude): argv[4]=kmax_run cap, argv[5]=kmin, argv[6]=stride.
  // Disjoint (kmin,stride) across processes -> collision-free packing (per-k files).
  int kmax_run = 100000;
  if(argc>4) kmax_run = atoi(argv[4]);
  int kmin = 1;
  if(argc>5) kmin = atoi(argv[5]);
  int stride = 1;
  if(argc>6) stride = atoi(argv[6]);
  // arg7 = explicit ensemble dir (redo production `..._hb*` names, which the legacy dir3
  // construction below cannot form; plan redo_obs_measurement_impl_plan_claude.md Sec 5)
  std::string ens_dir = "";
  if(argc>7) ens_dir = argv[7];
  // arg8 = optional output-prefix suffix (e.g. "_v2.1") to tag a distinct data set (e.g. the
  // a_t-corrected at0.1 flow) WITHOUT clobbering existing h5.  Empty = production prefix.
  std::string out_suffix = "";
  if(argc>8) out_suffix = argv[8];
  std::cout << "# gsq = " << gsq << " Nf = " << Nf << " nu0 = " << nu0
            << " kmax_run = " << kmax_run << " kmin = " << kmin << " stride = " << stride
            << " ens_dir = " << (ens_dir.empty() ? std::string("<legacy dir3>") : ens_dir) << std::endl;


  // int igam=0;
  // if(argc>3) igam = atoi(argv[3]);
  // std::cout << "# igam = " << igam << " where 0=id., i=sigma_3" << std::endl;



  // int device;
  // CUDA_CHECK(cudaGetDeviceCount(&device));
  // cudaDeviceProp device_prop[device];
  // cudaGetDeviceProperties(&device_prop[0], 0);
  // std::cout << "# dev = " << device_prop[0].name << std::endl;
  // CUDA_CHECK(cudaSetDevice(0));// "TITAN V"
  // std::cout << "# (GPU device is set.)" << std::endl;

  // ---------------------------------------

  constexpr Idx N = Comp::N;
  constexpr int Nt = Comp::Nt;

#ifdef IS_DUAL
  using Base=S2Trivalent;
  // using WilsonDirac=DiracExt<Base, DiracS2Dual>;
#else
  using Base=S2Simp;
  // using WilsonDirac=DiracExt<Base, DiracS2Simp>;
#endif

  using Force=GaugeExt<Base,Nt,Comp::is_compact>;
  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  using Action=U1WilsonExt<Base>;

  using Rng=ParallelRngExt<Base,Nt>;
  // using FermionVector = FermionVector<Base>;


  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;

  // ----------------------
  // const double at = 0.5;
  // const double T = 0.2;
  // const double T = 24;
  // a_t enters ONLY the Wilson-flow smearing (Action SW -> beta_s = at/(vol*gsq), linear in at).
  // Derive it from the ensemble dir so at0.1 configs are flowed with beta_s(0.1); at0.2 dirs
  // resolve to 0.2 = the old hardcode (bit-identical).  Fallback 0.2 for the legacy/free path.
  // const double at = 0.2; // T/Comp::Nt; (OLD hardcode -- replaced by dir-derived at)
  double at = 0.2;
  if(!ens_dir.empty()){
    const double at_ens = at_from_ensdir(ens_dir);
    if(at_ens > 0.0) at = at_ens;
  }
  std::cout << "# at = " << at << " (from ens_dir)" << std::endl;
  if(Nt!=1) assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);


  // const double gsq = 0.05;
  // double at = 0.05; // base.mean_ell * 0.125 * ratio;
  // if(Comp::Nt==1) at=0.;
  Action SW( gsq, at, base );
  std::cout << "# alat = " << base.mean_ell << std::endl;

  std::string dir3, dir4;
  // #ifdef Nf2
  if(!ens_dir.empty()){
    // explicit ensemble dir; output h5 shares data_<basename>/ with the fermionic obs
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
    // dir4 gets a "_msm" suffix so the N_FLOW*n_lm operator output does not clobber
    // the 16-op data of glue2_claude.cu (dir3 input ckpoints are shared, unchanged).
    dir4="data_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";  // shared data_ dir (same as mesonic correlators)
    std::cout << "dir3 = " << dir3 << std::endl;
  }
  else{
    dir3="Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
    dir4="data_Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";  // shared data_ dir (same as mesonic correlators)
  }
  std::filesystem::create_directory(dir4);


  Gauge U(base);

  const int k_ckpoint=1;
  // const int kinit=10;
  const int kmax=1e5;

  int k_tmp=0;
  {
    // find the highest available config, scanning from kmin in steps of `stride` so ensembles that
    // are NOT numbered contiguously from 1 are detected (e.g. the free gauge set numbered 10,20,30,...).
    // (Previously scanned 1,2,3,... in steps of k_ckpoint=1 and stopped at the first gap -> k_tmp=0
    //  whenever ckpoint_lat.1 was absent, which silently measured nothing.)
    for(k_tmp=kmin; k_tmp<=kmax; k_tmp+=stride ){
      const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k_tmp);
      const bool bool_lat = std::filesystem::exists(str_lat);
      if(!bool_lat) break;
    }
    k_tmp -= stride;
  }

  // std::cout << "# kinit = " << kinit << " k_tmp = " << k_tmp << std::endl;

  // ---- smearing levels (stored on the "shape" axis, flow-major: ishape = ilev*N_VAR + ivar) ----
  // ilev 0: raw (exact Gauss-law operator); 1: spatial flow tau=2.0 (= glue_msm_shapes smearing, dt 0.01);
  // 2..4: SHORT-TIME FULL 3D flow, cumulative tau = 0.02, 0.05, 0.1 (NM 2026-10-12), dt 0.001 (temporal
  // stiffness ~ beta_t ~ 1/(a_t g^2) is larger than the spatial one; stability checked in the test).
  constexpr int N_LEV = 5;
  const double LEV_TAU[N_LEV] = { 0.0, 2.0, 0.02, 0.05, 0.1 };
  const int LEV_KIND[N_LEV] = { 0, 1, 2, 2, 2 };           // 0 raw, 1 spatial flow, 2 full 3D flow
  constexpr double SP_DT = 0.01;
  constexpr double FULL_DT = 0.001;
  constexpr int N_VAR = 2;                                 // 0 = A (exact gradient), 1 = B (face direction)
  const int n_shapes = N_LEV * N_VAR;
  const std::string H5PREFIX = "glue_sidewall_psi";

  // F-type (l,m) list, same as glue_msm_shapes (l=1 ground of interest, l=2 cheap).
  const std::vector<std::array<int,2>> lm_set = {
    {1,-1},{1,0},{1,1},
    {2,-2},{2,-1},{2,0},{2,1},{2,2},
  };
  const int n_lm = lm_set.size();
  const int nops = n_shapes * n_lm;        // op = ishape*n_lm + ilm
  std::cout << "# n_shapes = " << n_shapes << " (n_lev = " << N_LEV << " x n_var = " << N_VAR
            << ") n_lm = " << n_lm << " nops = " << nops << std::endl;

  // per-link projection coefficients coef[ivar][ilm][il] (config independent)
  const Idx n_links = base.n_links;
  std::vector<std::vector<std::vector<double>>> coef( N_VAR,
      std::vector<std::vector<double>>( n_lm, std::vector<double>( n_links, 0.0 ) ) );
  std::vector<Link> link_of( n_links );
  for(const Link& link : base.links) {
    const Idx il = base.map2il.at(link);
    link_of[il] = link;
    const VE r0 = base.sites[link[0]] / base.sites[link[0]].norm();
    const VE r1 = base.sites[link[1]] / base.sites[link[1]].norm();
    VE mid = r0 + r1;
    mid /= mid.norm();
    VE dhat = r1 - r0;
    dhat /= dhat.norm();
    const double w = base.link_volume[il] / std::pow(base.ell[il], 2);
    const double eps = 1.0e-5;
    VE rp = mid + eps * dhat;
    rp /= rp.norm();
    VE rm = mid - eps * dhat;
    rm /= rm.norm();
    for(int ilm=0; ilm<n_lm; ilm++){
      const int ell = lm_set[ilm][0];
      const int em  = lm_set[ilm][1];
      // A: exact discrete gradient along the link x0 -> x1
      coef[0][ilm][il] = w * ( Ylm_real(ell, em, r1) - Ylm_real(ell, em, r0) );
      // B: \Psi(mid).N |l| = \nabla_{tan}Y(mid).dhat * ell_geodesic (central difference on the sphere)
      const double dY = ( Ylm_real(ell, em, rp) - Ylm_real(ell, em, rm) ) / (2.0*eps);
      coef[1][ilm][il] = w * base.ell[il] * dY;
    }
  }

  if(k_tmp > kmax_run) k_tmp = kmax_run;
  // serial over configs k; parallelism is ensemble-level (one process per Nf).
  for(int k=kmin; k<=k_tmp; k+=stride ){
    // resume-safe: skip configs already fully measured (h5 with "complete" flag)
    const std::string h5path = dir4+H5PREFIX+out_suffix+"."+std::to_string(k)+".h5";
    {
      bool done=false;
      if(std::filesystem::exists(h5path)){
        try { HighFive::File f(h5path, HighFive::File::ReadOnly); done = f.exist("complete"); } catch(...) {}
      }
      if(done) continue;
    }
    std::cout << "# read from k = " << k << std::endl;
    const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k);
    U.read( str_lat );

    // obs[(ilev*N_VAR + ivar)*n_lm + ilm][t]
    std::vector<std::vector<double>> obs;
    measure_sidewall( SW, U, coef, link_of, LEV_TAU, LEV_KIND, N_LEV, SP_DT, FULL_DT, obs );

#ifdef SIDEWALL_TEST
    run_sidewall_tests( SW, U, base, coef, link_of, lm_set, LEV_TAU, LEV_KIND, N_LEV, SP_DT, FULL_DT, obs,
                        gsq, at );
    return 0;   // test mode: first config only, no h5 written
#endif

    // std::vector<double> plaq_avg(Comp::Nt);
    // for(int t=0; t<Comp::Nt; t++) plaq_avg[t] = U.plaquette_angle_avg(t);

    // for(int t=0; t<Comp::Nt; t++) flow_plaq_avg[t] = Uflow.plaquette_angle_avg(t);
    // std::vector<double> plaq_avg_00(Comp::Nt);
    // std::vector<double> plaq_avg_1m1(Comp::Nt);
    // std::vector<double> plaq_avg_10(Comp::Nt);
    // std::vector<double> plaq_avg_11(Comp::Nt);
    // std::vector<double> plaq_avg_2m2(Comp::Nt);
    // std::vector<double> plaq_avg_2m1(Comp::Nt);
    // std::vector<double> plaq_avg_20(Comp::Nt);
    // std::vector<double> plaq_avg_21(Comp::Nt);
    // std::vector<double> plaq_avg_22(Comp::Nt);
    // std::vector<double> plaq_avg_3m3(Comp::Nt);
    // std::vector<double> plaq_avg_3m2(Comp::Nt);
    // std::vector<double> plaq_avg_3m1(Comp::Nt);
    // std::vector<double> plaq_avg_30(Comp::Nt);
    // std::vector<double> plaq_avg_31(Comp::Nt);
    // std::vector<double> plaq_avg_32(Comp::Nt);
    // std::vector<double> plaq_avg_33(Comp::Nt);
    //
    // std::vector<double> plaq_avg_temporal_00(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_1m1(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_10(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_11(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_2m2(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_2m1(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_20(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_21(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_22(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_3m3(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_3m2(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_3m1(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_30(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_31(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_32(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_33(Comp::Nt);
    //
    // (named flow_plaq_avg_* vectors + obs_ptrs of the original now replaced by
    //  the obs[][] flow-checkpoint loop above; pristine version in glue2_claude.cu)

    // correlator matrix C(dt)[i][j] = (1/Nt) sum_t obs[i][t] obs[j][t+dt], and one-point <O_i>.
    // Only the small-separation window dt = 0..TMAX_CORR is stored: the GEVP uses dt up to tcut(~5),
    // and the backward fold slice is recovered LOSSLESSLY from the transpose of a stored forward
    // slice via the periodicity identity  C_ij(Nt-d) = C_ji(d)  (see glue_gevp_analysis_claude.cu
    // fold). This cuts both the correlator cost (Nt^2 -> Nt*TMAX) and the h5 size (~Nt/TMAX).
    constexpr int TMAX_CORR = 16;
    const int nsep = std::min(TMAX_CORR + 1, Comp::Nt); // stored separations dt = 0..nsep-1
    // Store ONLY the symmetry-allowed per-(l,m) blocks: for each ilm channel, the n_shapes x n_shapes
    // SHAPE correlator C(a*n_lm+ilm, b*n_lm+ilm). Cross-(l,m) entries vanish by rotational symmetry
    // (the analysis zero_noise) so they are neither computed nor stored -> ~n_lm x smaller + faster.
    // Layout: blk[dt][ilm*nob*nob + a*nob + b]. The analysis re-expands to the block-diagonal matrix.
    const int nob = n_shapes;
    std::vector<std::vector<double>> Fcorr( nsep, std::vector<double>(n_lm*nob*nob, 0.0) );
    for(int dt=0; dt<nsep; dt++){
      for(int t=0; t<Comp::Nt; t++) {
        const int tj = (t+dt)%Comp::Nt;
        for(int ilm=0; ilm<n_lm; ilm++){
          for(int a=0; a<nob; a++){
            const double oa = obs[a*n_lm+ilm][t];
            for(int b=0; b<nob; b++){
              Fcorr[dt][ilm*nob*nob + a*nob + b] += oa * obs[b*n_lm+ilm][tj];
            }
          }
        }
      }
      for(double& x : Fcorr[dt]) x /= Comp::Nt; // F_corr and F share the <..> scale for vacuum subtraction
    }

    std::vector<double> F1( nops, 0.0 );
    for(int t=0; t<Comp::Nt; t++) {
      for(int i=0; i<nops; i++) F1[i] += obs[i][t];
    }
    for(int i=0; i<nops; i++) F1[i] /= Comp::Nt;

    // write per-config HDF5 (F_corr_blk: nsep x n_lm*n_shapes^2 block-diagonal, F: nops); "complete" LAST
    HighFive::File h5( h5path, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate );
    h5.createDataSet( "F_corr_blk", Fcorr );  // per-(l,m) shape blocks (symmetry-allowed only)
    h5.createDataSet( "F", F1 );
    h5.createDataSet( "n_lm", std::vector<int>{n_lm} );  // (l,m) count so the analysis auto-adapts to the l-tower
    h5.createDataSet( "n_shapes", std::vector<int>{n_shapes} );  // shapes/orbit count -> nops = n_shapes*n_lm
    // smearing provenance: n_shapes = n_flow * n_shapes_geom, ishape = ilev*N_VAR + ivar (flow-major as msm);
    // n_shapes_geom = N_VAR (0 = A exact gradient, 1 = B face direction); flow_kind 0 raw / 1 spatial / 2 full 3D.
    h5.createDataSet( "n_flow", std::vector<int>{N_LEV} );
    h5.createDataSet( "n_shapes_geom", std::vector<int>{N_VAR} );
    h5.createDataSet( "flow_times", std::vector<double>( LEV_TAU, LEV_TAU + N_LEV ) );
    h5.createDataSet( "flow_kind", std::vector<int>( LEV_KIND, LEV_KIND + N_LEV ) );
    // explicit channel labels: ell[ilm], em[ilm] (op = iorbit*n_lm + ilm) so the analysis reads the
    // EXACT (l,m) list -- required now that l is non-contiguous (F saves l=1,2, no l=0).
    std::vector<int> ell_v( n_lm ), em_v( n_lm );
    for(int i=0; i<n_lm; i++){ ell_v[i] = lm_set[i][0]; em_v[i] = lm_set[i][1]; }
    h5.createDataSet( "ell", ell_v );
    h5.createDataSet( "em", em_v );
    h5.createDataSet( "complete", std::vector<int>{1} );
  } // end for k

  // ------------------


  return 0;

}

