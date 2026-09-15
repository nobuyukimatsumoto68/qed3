// jj_sigma_loops_stoch_claude.cu
// MASSLESS scalar zero-momentum loops for the sigma^2 - F^2 (0++) mixing program
// (Chester-Pufu, arXiv:1603.05582).  Measures, per config, the disconnected one-point loops
//   D_S   J_1   = tr[W_lm(t) D_ov^{-1}]              (PS vertex)
//   D_S   J_1mD = tr[W_lm(t) (1-D_ov) D_ov^{-1}]     (FS vertex, GW-furnished)
// AND the EXTENDED (one extra propagator hop) loops D'_S:
//   D'_1   (PS extended)  and  D'_1mD (FS extended).
//
// D'_S estimator: for each time+spin dilution class, take the class solution phi = D_ov^{-1} eta,
// project it back onto THIS dilution class (zero every timeslice outside t = t_s, t_s+interval, ...),
// then do ONE second solve D_ov^{-1}(phiP).  This reproduces the extended (double-hop) loop with the
// contamination from cross-class timeslices suppressed by the dilution interval (residual only at the
// t = Nt/2 wrap image).  W_lm = area-weighted real Y_lm identity vertex (mult_Ylm_real), NO mult_sigma.
//
// This driver uses ONLY the massless overlap D_ov (M5=-1); NO massive D_m, NO backward solve.  The RNG
// stream is folded with --seed-tag (default "sigmaloops") so its noise is INDEPENDENT of the disc-loop
// and FNAL-conn streams.  Boilerplate ported from src/both_3d/jj_local_ylm_scalar_disc_stoch_claude.cu;
// a_t resolution + massless operator construction ported from
// src/production/jj_local_ylm_scalar_conn_stoch_claude.cu.  Physics/design in the plan
// sigma_sigma_f2_mixing_impl_plan_claude.md ("Connected four-point" + "Measurement plan" sections).
#include <typeinfo>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <highfive/H5File.hpp>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <filesystem>
#include <memory>
#include <chrono>
#include <cstdint>
#include <complex>
#include <array>
#include <vector>
#include <map>
#include <string>
#include <random>
#include <Eigen/Dense>

using Double = double;
using Idx = std::int32_t;
using Complex = std::complex<double>;

using Link = std::array<Idx,2>;
using Face = std::vector<Idx>;

using MS=Eigen::Matrix2cd;
using VD=Eigen::Vector2d;
using VE=Eigen::Vector3d;
using VC=Eigen::VectorXcd;

static constexpr int NS = 2;
static constexpr int DIM = 2;
static constexpr Complex I = Complex(0.0, 1.0);

namespace Comp{
  constexpr bool is_compact=false;

  constexpr int NPARALLEL_DUPDATE=1;
  constexpr int NPARALLEL=NPARALLEL_DUPDATE;
  constexpr int NSTREAMS=4;
  constexpr int NPARALLEL_GAUGE=NPARALLEL_DUPDATE;
  constexpr int NPARALLEL_SORT=NPARALLEL_DUPDATE;

  // constexpr int N_REFINE=1;   // L1 default; now compile-time via -DN_REFINE_CLI
#ifndef N_REFINE_CLI
#define N_REFINE_CLI 1
#endif
  constexpr int N_REFINE=N_REFINE_CLI;
  constexpr int NS=2;
  constexpr int Nt=128;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;

  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-5;
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

#include "valence_claude.h"
#include "gauge_ext.h"
#include "action_ext.h"

#include "sparse_matrix.h"
#include "dirac_base.h"
#include "dirac_simp.h"
#include "dirac_dual.h"
#include "dirac_ext.h"

#include "sparse_dirac_claude.h"
#include "matpoly_claude.h"

#include "overlap_wmass_claude.h"        // complex-mass overlap (massless at mass=0)

//------------------------------------------
#include <getopt.h>

// Stable string -> RNG seed (std::seed_seq; reproducible from the stored rng_seed string).
static int seed_from_string(const std::string& s){
  std::seed_seq seq(s.begin(), s.end());
  std::uint32_t w;
  seq.generate(&w, &w + 1);
  return static_cast<int>(w);
}

// extract a_t from an ensemble dir name.  The a_t token is the unique 'at' immediately
// followed by a digit, e.g. "...gsq2.000000at0.100000nu0..." -> 0.1 .  Returns -1.0 if not found.
// Used so the valence temporal spacing ALWAYS tracks the sea ensemble (no silent 0.2 mismatch).
static double at_from_ensdir(const std::string& s){
  for(std::size_t p = 0; p + 2 < s.size(); ++p){
    const bool is_at = (s[p] == 'a' && s[p+1] == 't');
    const bool digit_after = (s[p+2] >= '0' && s[p+2] <= '9');
    if(is_at && digit_after) return std::stod(s.substr(p+2));
  }
  return -1.0;
}

// Project a full-lattice spinor onto ONE time-dilution class: copy `in` only on the class timeslices
// t = t_s, t_s+interval, ... (BOTH spin components, all spatial sites), zero everywhere else.  Mirrors
// the timeslice stride of FermionVector::time_spin_dilution (t=t_s; t<Nt; t+=interval).  This is the
// middle "projection to the class" that turns the single solve phi = D_ov^{-1} eta into the extended
// (double-hop) loop D'_S when followed by a second solve.  File-scope static helper (no lambda).
static void project_to_class(FermionVector& out, const FermionVector& in,
                             const int t_s, const int interval){
  memset(out.field, 0, Comp::N*CD);
  for(int t=t_s; t<Comp::Nt; t+=interval){
    for(Idx ix=0; ix<Comp::N_SITES; ix++){
      for(int i=0; i<NS; i++){
        out(t, ix, i) = in(t, ix, i);
      }
    }
  }
}

void PrintHelp(){
  printf("jj_sigma_loops_stoch: MASSLESS scalar zero-momentum loops D_S (J_1, J_1mD) and extended D'_S\n");
  printf("  (D'_1, D'_1mD) for the sigma^2 - F^2 (0++) mixing program (arXiv:1603.05582).\n");
  printf("  --gsq <x>            Wilson coupling squared (ensemble id; default 8.0)\n");
  printf("  --Nf <n>             number of fermion flavors (ensemble id; default 2)\n");
  printf("  --nu0 <x>            sea quark asymmetry (ensemble id; default 1.0)\n");
  printf("  --nu1 <x>            valence Wilson-Dirac asymmetry (operator; default nu0)\n");
  printf("  --at <x>             valence temporal spacing a_t (default: derive from ens-dir; free field 0.2)\n");
  printf("  --mass-re <x>        IGNORED (massless driver; forced to 0)\n");
  printf("  --mass-im <y>        IGNORED (massless driver; forced to 0)\n");
  printf("  --ens-dir <path>     sea config dir; OMIT => free field (U=1)\n");
  printf("  --nhits <n>          stochastic hits (default 2)\n");
  printf("  --disc-tblock <tb>   timeslices per time-dilution class (default 2; interval=Nt/tb classes)\n");
  printf("  --stride <s>         ensemble config stride (default 10)\n");
  printf("  --kmin <a> --kmax <b> config range [a,b) (default 0..1e6)\n");
  printf("  --seed-tag <str>     RNG stream tag (default 'sigmaloops'; keeps noise independent of other streams)\n");
  printf("  -h, --help\n");
  exit(0);
}

void ParseArgs(int argc, char* argv[],
               double& gsq, int& Nf, double& nu0, double& nu1,
               double& mass_re, double& mass_im,
               std::string& ens_dir, int& nhits, int& disc_tblock, int& stride,
               int& kmin, int& kmax, double& at_cli, std::string& seed_tag){
  static struct option long_opts[] = {
    {"gsq",     required_argument, nullptr, 'g'},
    {"Nf",      required_argument, nullptr, 'N'},
    {"nu0",     required_argument, nullptr, 'n'},
    {"nu1",     required_argument, nullptr, 'm'},
    {"at",      required_argument, nullptr, 'A'},
    {"mass-re", required_argument, nullptr, 'r'},
    {"mass-im", required_argument, nullptr, 'i'},
    {"ens-dir", required_argument, nullptr, 'e'},
    {"nhits",   required_argument, nullptr, 'H'},
    {"disc-tblock", required_argument, nullptr, 'T'},
    {"stride",  required_argument, nullptr, 'I'},
    {"kmin",    required_argument, nullptr, 'a'},
    {"kmax",    required_argument, nullptr, 'b'},
    {"seed-tag", required_argument, nullptr, 'G'},
    {"help",    no_argument,       nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
  };
  int opt, idx;
  while((opt = getopt_long(argc, argv, "g:N:n:m:A:r:i:e:H:T:I:a:b:G:h", long_opts, &idx)) != -1){
    switch(opt){
    case 'g': gsq     = std::stod(optarg); break;
    case 'N': Nf      = std::stoi(optarg); break;
    case 'n': nu0     = std::stod(optarg); break;
    case 'm': nu1     = std::stod(optarg); break;
    case 'A': at_cli  = std::stod(optarg); break;
    case 'r': mass_re = std::stod(optarg); break;
    case 'i': mass_im = std::stod(optarg); break;
    case 'e': ens_dir = optarg; break;
    case 'H': nhits   = std::stoi(optarg); break;
    case 'T': disc_tblock = std::stoi(optarg); break;
    case 'I': stride  = std::stoi(optarg); break;
    case 'a': kmin    = std::stoi(optarg); break;
    case 'b': kmax    = std::stoi(optarg); break;
    case 'G': seed_tag = optarg; break;
    case 'h':
    case '?':
    default:  PrintHelp(); break;
    }
  }
}
//------------------------------------------

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  double gsq=8.0;  int Nf=2;  double nu0=1.0;  double nu1=-1.0;
  double mass_re=0.0, mass_im=0.0;
  std::string ens_dir="";
  int nhits=2;
  int disc_tblock=2;         // interval = Nt/disc_tblock time-dilution classes
  int stride=10;
  int kmin=0;
  int kmax=1000000;
  double at_cli=-1.0;        // -1 => derive a_t from ens-dir (free field: 0.2)
  std::string seed_tag="sigmaloops";

  ParseArgs(argc, argv, gsq, Nf, nu0, nu1, mass_re, mass_im, ens_dir, nhits, disc_tblock, stride,
            kmin, kmax, at_cli, seed_tag);
  if(nu1 < 0.0) nu1 = nu0;

  // MASSLESS driver: mass CLI is ignored (forced to 0); D_ov is the only operator used.
  mass_re = 0.0;
  mass_im = 0.0;
  const bool free_field = ens_dir.empty();

  // resolve valence temporal spacing a_t: --at overrides; else derive from the ensemble dir name
  // (so a_t ALWAYS tracks the sea ensemble); free field defaults to 0.2.
  double at = at_cli;
  if(at < 0.0) at = free_field ? 0.2 : at_from_ensdir(ens_dir);
  assert(at > 0.0 && "could not resolve a_t from ens-dir; pass --at explicitly");
  assert(disc_tblock >= 1 && Comp::Nt % disc_tblock == 0 && "Nt % disc_tblock must be 0");

  std::cout << "# gsq="<<gsq<<" Nf="<<Nf<<" nu0="<<nu0<<" nu1="<<nu1<<" at="<<at
            << " (massless D_ov)"
            << " ens_dir="<<(free_field?std::string("<free-field U=1>"):ens_dir)
            << " nhits="<<nhits<<" disc_tblock="<<disc_tblock
            << " seed_tag='"<<seed_tag<<"'" << std::endl;

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();

  int device;
  CUDA_CHECK(cudaGetDeviceCount(&device));
  cudaDeviceProp device_prop[device];
  cudaGetDeviceProperties(&device_prop[0], 0);
  std::cout << "# dev = " << device_prop[0].name << std::endl;
  CUDA_CHECK(cudaSetDevice(0));

  constexpr Idx N  = Comp::N;
  constexpr int Nt = Comp::Nt;

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  using Rng=ParallelRngExt<Base,Nt>;
  using Fermion=OverlapWMass<WilsonDirac>;

  Base base(Comp::N_REFINE);
  std::cout << "# lattice set." << std::endl;

  const double M5 = -1.0;
  if(Nt!=1) assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu1);
  std::cout << "# DW set." << std::endl;

  Gauge U(base);
  Rng rng(base, 1234);

  // D = D_ov (massless overlap).  This driver uses ONLY the massless operator: forward solve
  // D_ov^{-1} b = op_Dsq^{-1}(op_DH b) (apply D_ov^dag, then CG on D_ov^dag D_ov; overlap is normal so exact),
  // plus the GW factor (1-D_ov) via op_oneMinusD.
  Fermion D(DW, Complex(0.0), 11);
  std::cout << "# overlap operator set: D_ov (massless, M5="<<M5<<")." << std::endl;

  // D_ov^{-1} via op_Dsq (CG on D_ov^dag D_ov = D_ov D_ov^dag, D_ov normal) + op_DH (RHS):
  auto f_D   = std::bind(&Fermion::mult_deviceAsyncLaunch_ms, &D, std::placeholders::_1, std::placeholders::_2);
  auto f_DH  = std::bind(&Fermion::adj_deviceAsyncLaunch_ms,  &D, std::placeholders::_1, std::placeholders::_2);
  auto f_Dsq = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms,  &D, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_D(f_D), M_DH(f_DH), M_Dsq(f_Dsq);
  MatPoly op_D;   op_D.push_back(cplx(1.0), {&M_D});
  MatPoly op_DH;  op_DH.push_back(cplx(1.0), {&M_DH});
  MatPoly op_Dsq; op_Dsq.push_back(cplx(1.0), {&M_Dsq});

  // massless D_ov apply: op_oneMinusD : v -> (1 - D_ov) v  (identity term + a -D_ov term).
  MatPoly op_oneMinusD;
  op_oneMinusD.push_back(cplx( 1.0), {});        // identity term (empty product)
  op_oneMinusD.push_back(cplx(-1.0), {&M_D});    // - D_ov
  // op_oneMinusDdag : v -> (1 - D_ov^dag) v.  The FS scalar vertex is \tilde S = -(1 - D_ov^dag) (Eq. 5.3;
  // sigma_FS = eta^dag xi - xi^dag (1 - D_ov^dag) eta).  D_S folds this via the (1-D_ov)+conj one-point trick,
  // but the TWO-vertex D'_S factor needs the literal (1 - D_ov^dag) at each vertex (conj does not recover it).
  MatPoly op_oneMinusDdag;
  op_oneMinusDdag.push_back(cplx( 1.0), {});     // identity term
  op_oneMinusDdag.push_back(cplx(-1.0), {&M_DH}); // - D_ov^dag

  constexpr int L_MAX_YLM = 3;
  const int interval = Nt / disc_tblock;   // number of time-dilution classes

  // ---- output: data_<ESNID>/corr_sigma_loops_tb<tb>_nhits<H>/corr.<k>.h<h>.h5
  std::string ens_base = ens_dir;
  if(!ens_base.empty() && ens_base.back()=='/') ens_base.pop_back();
  { const auto slash = ens_base.find_last_of('/'); if(slash!=std::string::npos) ens_base = ens_base.substr(slash+1); }
  const std::string esnid = (free_field ? std::string("free") : ens_base)
                          + "_vmRe"+std::to_string(mass_re)+"vmIm"+std::to_string(mass_im);
  const std::string dir_out = "data_"+esnid+"/corr_sigma_loops_tb"+std::to_string(disc_tblock)
                            + "_nhits"+std::to_string(nhits)+"/";
  std::filesystem::create_directories(dir_out);
  std::cout << "# dir_out = " << dir_out << std::endl;

  // per-(l,m) complex length-Nt loop vector, written split real/imag (RAW, no fold -- matches disc driver).
  auto write_vec = [&](HighFive::File& h5, const std::string& key, const std::vector<Complex>& C){
    std::vector<double> re(C.size()), im(C.size());
    for(size_t t=0;t<C.size();t++){ re[t]=C[t].real(); im[t]=C[t].imag(); }
    h5.createDataSet(key+"/real", re);  h5.createDataSet(key+"/imag", im);
  };

  // work vectors: eta = diluted Z2 source; phi = D_ov^{-1} eta; dphi = (1-D_ov) phi (D_S J_1mD, no dagger);
  // phiP/dphiP = class-projected middles; chi1 = D_ov^{-1} phiP; chimD = D_ov^{-1} dphiP; dchi = (1-D_ov^dag) chimD;
  // dphidag = (1-D_ov^dag) phi (FS D'_S middle vertex).
  FermionVector eta, tmp, phi, dphi, dphidag, phiP, dphiP, chi1, chimD, dchi, Gamma;

  const int k_ckpoint = free_field ? 1 : stride;
  const int k_lo      = free_field ? 0 : kmin;
  const int k_hi      = free_field ? 1 : kmax;

  for(int k = k_lo; k < k_hi; k += k_ckpoint){
    std::string str_lat;
    if(!free_field){
      str_lat = ens_dir + "ckpoint_lat." + std::to_string(k);
      if(!std::filesystem::exists(str_lat)){ if(k==0) continue; else break; }
    }
    // CHEAP pre-skip (BEFORE any construction): skip the whole config if every hit's DS group is already
    // present -- only the output .h5 is needed, NOT U.read or D.update (the expensive lambda_min/max).
    {
      bool all_done = true;
      for(int h=0; h<nhits; h++){
        const std::string h5p = dir_out + "corr." + std::to_string(k) + ".h" + std::to_string(h) + ".h5";
        bool done_h = false;
        if(std::filesystem::exists(h5p)){
          try {
            HighFive::File f(h5p, HighFive::File::ReadOnly);
            // non-throwing DS check (navigate groups; getDataSet throws + HDF5 spams stderr when absent).
            if(f.exist("h0")){
              auto g0 = f.getGroup("h0");
              if(g0.exist("sigma_loops")) done_h = g0.getGroup("sigma_loops").exist("DS");
            }
          } catch(...) {}
        }
        if(!done_h){ all_done = false; break; }
      }
      if(all_done){ std::cout<<"# skip k="<<k<<" (all "<<nhits<<" hits done; no U.read/update)"<<std::endl; continue; }
    }
    if(!free_field) U.read(str_lat);
    D.update(U);
    std::cout << "# k="<<k<<(free_field?" (free field)":"")
              << "  lambda_min/max="<<D.lambda_min<<"/"<<D.lambda_max<<std::endl;

    for(int h=0; h<nhits; h++){
      const std::string h5path_h = dir_out + "corr." + std::to_string(k) + ".h" + std::to_string(h) + ".h5";
      // Per-hit skip: skip if this hit's file already has the DS group.
      if(std::filesystem::exists(h5path_h)){
        bool has_DS=false;
        try {
          HighFive::File f(h5path_h, HighFive::File::ReadOnly);
          if(f.exist("h0")){
            auto g0 = f.getGroup("h0");
            if(g0.exist("sigma_loops")) has_DS = g0.getGroup("sigma_loops").exist("DS");
          }
        } catch(...) {}
        if(has_DS){ std::cout<<"# skip k="<<k<<" hit "<<h<<" (DS present)"<<std::endl; continue; }
      }
      // seed folds in seed_tag so this stream's noise is INDEPENDENT of disc-loop / FNAL-conn streams.
      const std::string seed_str = esnid + "_" + seed_tag + "_k" + std::to_string(k) + "_h" + std::to_string(h);
      rng.reseed(seed_from_string(seed_str));
      const auto t_hit0 = std::chrono::steady_clock::now();
      std::cout << "# k="<<k<<" hit "<<(h+1)<<"/"<<nhits<<"  (disc_tblock="<<disc_tblock
                <<", "<<interval<<" time x "<<NS<<" spin classes, seed='"<<seed_str<<"')" << std::endl;

      // accumulators (per l,m; length-Nt complex), summed over the (t_s, spin) dilution classes.
      //   JS    = D_S  J_1   = tr[W_lm D_ov^{-1}]              (PS one-point loop)
      //   JS1mD = D_S  J_1mD = tr[W_lm (1-D_ov) D_ov^{-1}]     (FS one-point loop)
      //   Dp1   = D'_S PS extended;   Dp1mD = D'_S FS extended.
      std::vector<std::vector<std::vector<Complex>>> JS(L_MAX_YLM+1), JS1mD(L_MAX_YLM+1),
                                                     Dp1(L_MAX_YLM+1), Dp1mD(L_MAX_YLM+1);
      for(int l=0;l<=L_MAX_YLM;l++){
        JS[l]   .assign(2*l+1, std::vector<Complex>(Nt, Complex(0,0)));
        JS1mD[l].assign(2*l+1, std::vector<Complex>(Nt, Complex(0,0)));
        Dp1[l]  .assign(2*l+1, std::vector<Complex>(Nt, Complex(0,0)));
        Dp1mD[l].assign(2*l+1, std::vector<Complex>(Nt, Complex(0,0)));
      }

      // ===== TIME+SPIN dilution sweep (spin loop MANDATORY, as in the disc driver) =====
      for(int t_s=0; t_s<interval; t_s++){
        for(int spin=0; spin<NS; spin++){
          eta.time_spin_dilution(rng, t_s, disc_tblock, spin);   // volume Z2, this (t-class, spin)

          // (A) first solve: phi = D_ov^{-1} eta = op_Dsq^{-1}(op_DH eta).
          op_DH.from_cpu<N>(tmp.field, eta.field);
          op_Dsq.solve<N>(phi.field, tmp.field, Comp::TOL_OUTER);

          // (B) dphi = (1 - D_ov) phi  (mat-vec; base for the FS one-point loop and the FS extended loop).
          op_oneMinusD.from_cpu<N>(dphi.field, phi.field);

          // (C) D_S one-point loops: JS from phi, JS1mD from dphi.
          for(int l=0; l<=L_MAX_YLM; l++){
            for(int m=-l; m<=l; m++){
              Gamma = phi;
              Gamma.mult_Ylm_real(l, m, base);
              eta.accumulate_loop_raw(JS[l][m+l], Gamma, t_s, disc_tblock, spin);
              Gamma = dphi;
              Gamma.mult_Ylm_real(l, m, base);
              eta.accumulate_loop_raw(JS1mD[l][m+l], Gamma, t_s, disc_tblock, spin);
            }
          }

          // (D) PS extended loop D'_1: project phi to the class, second solve, then Y_lm + accumulate.
          project_to_class(phiP, phi, t_s, interval);
          op_DH.from_cpu<N>(tmp.field, phiP.field);
          op_Dsq.solve<N>(chi1.field, tmp.field, Comp::TOL_OUTER);    // chi1 = D_ov^{-1} phiP
          for(int l=0; l<=L_MAX_YLM; l++){
            for(int m=-l; m<=l; m++){
              Gamma = chi1;
              Gamma.mult_Ylm_real(l, m, base);
              eta.accumulate_loop_raw(Dp1[l][m+l], Gamma, t_s, disc_tblock, spin);
            }
          }

          // (E) FS extended loop D'_1mD: \tilde S = -(1-D_ov^dag) at BOTH vertices (Eq. 5.3).  The two vertex
          // signs cancel ((-1)^2=+1), so apply (1-D_ov^dag) twice with no explicit minus.  Middle vertex:
          // (1-D_ov^dag) phi; project to the class; third solve; sink vertex: (1-D_ov^dag); accumulate.
          op_oneMinusDdag.from_cpu<N>(dphidag.field, phi.field);      // (1-D_ov^dag) phi   (middle vertex)
          project_to_class(dphiP, dphidag, t_s, interval);
          op_DH.from_cpu<N>(tmp.field, dphiP.field);
          op_Dsq.solve<N>(chimD.field, tmp.field, Comp::TOL_OUTER);   // chimD = D_ov^{-1} (P_class (1-D_ov^dag) phi)
          op_oneMinusDdag.from_cpu<N>(dchi.field, chimD.field);       // dchi = (1-D_ov^dag) chimD   (sink vertex)
          for(int l=0; l<=L_MAX_YLM; l++){
            for(int m=-l; m<=l; m++){
              Gamma = dchi;
              Gamma.mult_Ylm_real(l, m, base);
              eta.accumulate_loop_raw(Dp1mD[l][m+l], Gamma, t_s, disc_tblock, spin);
            }
          }
        }
      }

      // ---- write fresh file per (k,h): .tmp then atomic rename ----
      const std::string hp = "h0/sigma_loops/";
      const std::string h5tmp = h5path_h + ".tmp";
      std::unique_ptr<HighFive::File> h5p = std::make_unique<HighFive::File>(h5tmp,
              HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
      HighFive::File& h5 = *h5p;
      h5.createDataSet("nhits",       std::vector<int>{nhits});
      h5.createDataSet("hit",         std::vector<int>{h});
      h5.createDataSet("disc_tblock", std::vector<int>{disc_tblock});
      h5.createDataSet("L_MAX_YLM",   std::vector<int>{L_MAX_YLM});
      h5.createDataSet("at",          std::vector<double>{at});
      h5.createDataSet("M5",          std::vector<double>{M5});
      h5.createDataSet("rng_seed",    seed_str);
      h5.createDataSet("seed_tag",    seed_tag);
      for(int l=0; l<=L_MAX_YLM; l++){
        for(int m=-l; m<=l; m++){
          const std::string lm = "l"+std::to_string(l)+"/m"+std::to_string(m)+"/J";
          write_vec(h5, hp+"DS/"    +lm, JS[l][m+l]);
          write_vec(h5, hp+"DS_1mD/"+lm, JS1mD[l][m+l]);
          write_vec(h5, hp+"Dp/"    +lm, Dp1[l][m+l]);
          write_vec(h5, hp+"Dp_1mD/"+lm, Dp1mD[l][m+l]);
        }
      }
      h5p.reset();
      std::filesystem::rename(h5tmp, h5path_h);
      const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now()-t_hit0).count();
      std::cout << "#   hit "<<(h+1)<<" done (D_S J_1/J_1mD + D'_S PS/FS, l<=3 per-m) ["<<secs<<" s] -> "
                << h5path_h << std::endl;
    } // hits
  } // k

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();
  return 0;
}
