// jj_conn_tpproj_axial_claude.cu
// Connected AXIAL current-current correlator C_{A+-}, TEMPORAL site-diagonal projection G_t
// (Eq. 4.32 of qed3int_v2-10.pdf).  Plan: conserved_current_correlators_impl_plan_v3_claude.md
// Sec. 3.3 (axial code realization) -- the axial sibling of jj_conn_tpproj_claude.cu (vector).
//
// Estimator (Sec. 3.3, fix-0/loop-t, shared buffers):
//   C_{A+-}(0->t) = sum_n w_n psi_n^dag phi_n(t),
//     phi'    = D_m^{-1} eta                                  (shared forward leg)
//     chi     = (1 - D_ov) phi'                               (shared; GW factor, ALWAYS massless)
//     psi_n   = [D_m or tilde D_{m_P}]^{-1} (1-D_ov) K^dag(n,0) eta  (fixed-0, one solve / site)
//     phi_n(t)= K^dag(n,t) chi                                (looped kernel, no inversion)
// Both GW factors collapse to the non-dagger (1-D_ov); parity (m purely imaginary) switches ONLY
// the psi_n forward solve D_m -> tilde D_{m_P} (Eq. 3.68).  Only C_{A+-} is written (Apm); the
// other ordering C_{A-+}(dt)=C_{A+-}(Nt-dt) is recovered by reflection downstream (Eq. 3.57).
//
// Conventions (v3): valence mass via --mass-re/--mass-im (D_m = D_ov + m); ensemble selected by
// --ens-dir (omit => free field, U = 1); kernel K is the massless-form Noether kernel.

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
  constexpr int NSTREAMS=4; // NPARALLEL_DUPDATE;
  constexpr int NPARALLEL_GAUGE=NPARALLEL_DUPDATE;
  constexpr int NPARALLEL_SORT=NPARALLEL_DUPDATE;

  constexpr int N_REFINE=1;
  constexpr int NS=2;
  constexpr int Nt=128;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;

  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;   // inner Zolotarev pole solves
  const double TOL_OUTER=1.0e-5;   // outer CG; above the ~1e-15 machine-precision residual floor
                                   // for the small-norm sink RHS (plenty for the correlator)
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

#include "sparse_dirac.h"
// #include "matpoly.h"
#include "matpoly_claude.h"

#include "overlap_wmass_claude.h"        // complex-mass overlap (massless at mass=0)
#include "conserved_current_claude.h"   // ConservedCurrent: apply_k / apply_k_dag

//------------------------------------------
#include <getopt.h>

void PrintHelp(){
  printf("Usage: ./a.out [options]\n");
  printf("  --gsq <gsq>          Wilson coupling squared (ensemble id; default: 8.0)\n");
  printf("  --Nf <Nf>            number of fermion flavors (ensemble id; default: 2)\n");
  printf("  --nu0 <nu0>          sea quark asymmetry (ensemble id; default: 1.0)\n");
  printf("  --nu1 <nu1>          valence Wilson-Dirac asymmetry (operator; default: nu0)\n");
  printf("  --mass-re <x>        valence mass Re (default: 0.0)\n");
  printf("  --mass-im <y>        valence mass Im (default: 0.0)\n");
  printf("  --current <c>        axial (this binary; default: axial)\n");
  printf("  --ens-dir <path>     sea config directory; OMIT => free field (U=1) check\n");
  printf("  --nhits <n>          stochastic hits (default: 1)\n");
  printf("  --tmax <T>           t-loop computes |dt| in [0,T] and [Nt-T,Nt-1] only (default: Nt/4)\n");
  printf("  -h, --help           show this help\n");
  exit(0);
}

void ParseArgs(int argc, char* argv[],
               double& gsq, int& Nf, double& nu0, double& nu1,
               double& mass_re, double& mass_im,
               std::string& current, std::string& ens_dir, int& nhits, int& tmax){
  static struct option long_opts[] = {
    {"gsq",     required_argument, nullptr, 'g'},
    {"Nf",      required_argument, nullptr, 'N'},
    {"nu0",     required_argument, nullptr, 'n'},
    {"nu1",     required_argument, nullptr, 'm'},
    {"mass-re", required_argument, nullptr, 'r'},
    {"mass-im", required_argument, nullptr, 'i'},
    {"current", required_argument, nullptr, 'c'},
    {"ens-dir", required_argument, nullptr, 'e'},
    {"nhits",   required_argument, nullptr, 'H'},
    {"tmax",    required_argument, nullptr, 'T'},
    {"help",    no_argument,       nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
  };
  int opt, idx;
  while((opt = getopt_long(argc, argv, "g:N:n:m:r:i:c:e:H:T:h", long_opts, &idx)) != -1){
    switch(opt){
    case 'g': gsq     = std::stod(optarg); break;
    case 'N': Nf      = std::stoi(optarg); break;
    case 'n': nu0     = std::stod(optarg); break;
    case 'm': nu1     = std::stod(optarg); break;
    case 'r': mass_re = std::stod(optarg); break;
    case 'i': mass_im = std::stod(optarg); break;
    case 'c': current = optarg; break;
    case 'e': ens_dir = optarg; break;
    case 'H': nhits   = std::stoi(optarg); break;
    case 'T': tmax    = std::stoi(optarg); break;
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

  double gsq=8.0;  int Nf=2;  double nu0=1.0;  double nu1=-1.0;   // nu1<0 => default to nu0
  double mass_re=0.0, mass_im=0.0;
  std::string current="axial";
  std::string ens_dir="";     // empty => free-field mode
  int nhits=1;
  int tmax=-1;   // <0 sentinel => default Nt/4 (set after Nt is known)

  ParseArgs(argc, argv, gsq, Nf, nu0, nu1, mass_re, mass_im, current, ens_dir, nhits, tmax);
  if(nu1 < 0.0) nu1 = nu0;    // valence asymmetry defaults to the sea value nu0 (knob retained)

  const Complex valence_mass(mass_re, mass_im);
  const bool free_field = ens_dir.empty();
  assert(current=="axial");

  // parity case: purely imaginary valence mass -> dagger leg uses \tilde D_{m_P} (Eq. 3.68).
  const bool parity = (std::abs(mass_im) > 1.0e-15) && (std::abs(mass_re) <= 1.0e-15);
  // flavor case: purely real valence mass.  At m_F != 0 the axial current is NOT conserved, so the
  // only well-defined correlator is the massless-limit expression -> BOTH propagator legs use the
  // massless D_ov (D instead of D_m).  Axial handles massless, pure m_F, or pure m_P -- not mixed.
  const bool flavor = (std::abs(mass_re) > 1.0e-15) && (std::abs(mass_im) <= 1.0e-15);
  assert( std::abs(mass_re) <= 1.0e-15 || std::abs(mass_im) <= 1.0e-15 );

  std::cout << "# gsq="<<gsq<<" Nf="<<Nf<<" nu0="<<nu0<<" nu1="<<nu1
            << " mass=("<<mass_re<<","<<mass_im<<")"
            << " current="<<current
            << " ens_dir="<<(free_field?std::string("<free-field U=1>"):ens_dir)
            << " nhits="<<nhits << std::endl;

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();

  int device;
  CUDA_CHECK(cudaGetDeviceCount(&device));
  cudaDeviceProp device_prop[device];
  cudaGetDeviceProperties(&device_prop[0], 0);
  std::cout << "# dev = " << device_prop[0].name << std::endl;
  CUDA_CHECK(cudaSetDevice(0));

  constexpr Idx N  = Comp::N;
  constexpr int Nt = Comp::Nt;

  // tmax cap: t-loop computes |dt| in [0,tmax] and [Nt-tmax,Nt-1] only (skip the noise middle).
  if(tmax < 0)     tmax = Nt/4;     // default
  if(tmax > Nt-1)  tmax = Nt-1;     // clamp
  std::cout << "# tmax=" << tmax << " (t-loop = [0,"<<tmax<<"] U ["<<(Nt-tmax)<<","<<(Nt-1)<<"])" << std::endl;

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  using Rng=ParallelRngExt<Base,Nt>;
  using Fermion=OverlapWMass<WilsonDirac>;

  // ---- operators (grouped as in hmc_w_mass_claude.cu) -----------------------------
  Base base(Comp::N_REFINE);
  std::cout << "# lattice set." << std::endl;

  const double r  = 1.0;
  const double M5 = -1.0;
  const double at = 0.2;
  if(Nt!=1) assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu1);
  std::cout << "# DW set." << std::endl;

  Gauge U(base);                       // free field (U=1) unless a config is read below
  Rng rng(base, 1234);

  // Overlap operators (all three built by value, sharing DW; updated each config):
  //   D    = D_ov              (massless; the GW factors (1-D_ov) are ALWAYS massless, Eq. 3.55)
  //   Dm   = D_ov + m          (forward leg D_m^{-1}; Eqs. 3.55, 3.68)
  //   Dtil = D_ov + m/(1-m)    (tilde D_{m_P}; parity dagger leg, Eq. 3.68)
  Fermion D   (DW, Complex(0.0), 21);
  Fermion Dm  (DW, valence_mass, 21);
  Fermion Dtil(DW, valence_mass / (Complex(1.0) - valence_mass), 21);
  std::cout << "# overlap operators set: D_ov, D_m, tilde D_{m_P} (M5="<<M5<<")." << std::endl;

  ConservedCurrent<Fermion,Gauge> kop(Dm);   // K is mass-independent; LinOp configured per (link,dag)
  MatPoly op_K; op_K.push_back(cplx(1.0), {&kop});   // apply K via op_K.from_cpu (no raw device buffers)

  // Operator set (Sec. 2.1/2.2).  Forward solve  X^{-1} b  via op_XH (RHS X^dag b) + op_Xsq (CG).
  // massless D_ov series (op_D / op_DH / op_Dsq).  Used by:
  //   - the GW factor (1 - D_ov), built once as op_oneMinusD (identity + (-D_ov), via M_D); and
  //   - the FLAVOR case, where the axial current is not conserved and the only well-defined
  //     correlator is the massless-limit expression -> BOTH legs use the massless forward solve
  //     op_DH + op_Dsq (Sec. 3.3).
  auto f_D   = std::bind(&Fermion::mult_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);
  auto f_DH  = std::bind(&Fermion::adj_deviceAsyncLaunch,  &D, std::placeholders::_1, std::placeholders::_2);
  auto f_Dsq = std::bind(&Fermion::DDH_deviceAsyncLaunch,  &D, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_D(f_D), M_DH(f_DH), M_Dsq(f_Dsq);
  MatPoly op_oneMinusD;
  op_oneMinusD.push_back(cplx( 1.0), {});       // identity term (empty product)
  op_oneMinusD.push_back(cplx(-1.0), {&M_D});   // - D_ov   => op_oneMinusD : v -> (1 - D_ov) v
  MatPoly op_DH;  op_DH.push_back(cplx(1.0), {&M_DH});
  MatPoly op_Dsq; op_Dsq.push_back(cplx(1.0), {&M_Dsq});

  // D_m forward solve  (phi' and the non-parity sink):
  auto f_DmH  = std::bind(&Fermion::adj_deviceAsyncLaunch, &Dm, std::placeholders::_1, std::placeholders::_2);
  auto f_Dmsq = std::bind(&Fermion::DDH_deviceAsyncLaunch, &Dm, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_DmH(f_DmH), M_Dmsq(f_Dmsq);
  MatPoly op_DmH;  op_DmH.push_back(cplx(1.0), {&M_DmH});
  MatPoly op_Dmsq; op_Dmsq.push_back(cplx(1.0), {&M_Dmsq});

  // tilde D_{m_P} forward solve (parity sink leg only):
  auto f_tilDmH  = std::bind(&Fermion::adj_deviceAsyncLaunch, &Dtil, std::placeholders::_1, std::placeholders::_2);
  auto f_tilDmsq = std::bind(&Fermion::DDH_deviceAsyncLaunch, &Dtil, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_tilDmH(f_tilDmH), M_tilDmsq(f_tilDmsq);
  MatPoly op_tilDmH;  op_tilDmH.push_back(cplx(1.0), {&M_tilDmH});
  MatPoly op_tilDmsq; op_tilDmsq.push_back(cplx(1.0), {&M_tilDmsq});

  // geometry weights for the temporal projection (Eq. 4.32): w_n = A_n / kappa_t(n)^2
  std::vector<double> w_tp(base.n_sites);
  for(Idx n=0; n<base.n_sites; n++){ const double kt=DW.kappa_t[n]; w_tp[n]=base.dual_areas[n]/(kt*kt); }
  const double inv4pi = 1.0/(4.0*std::acos(-1.0));
  const int n_sites = static_cast<int>(base.n_sites);

  // output: data_<ESNID>/tp_axial/conn.<config>.<hit>.h5   (datasets: "Apm/real", "Apm/imag")
  std::string ens_base = ens_dir;
  if(!ens_base.empty() && ens_base.back()=='/') ens_base.pop_back();
  { const auto slash = ens_base.find_last_of('/'); if(slash!=std::string::npos) ens_base = ens_base.substr(slash+1); }
  const std::string esnid = (free_field ? std::string("free") : ens_base)
                          + "_vmRe"+std::to_string(mass_re)+"vmIm"+std::to_string(mass_im);
  const std::string dir_out = "data_"+esnid+"/tp_"+current+"_tmax"+std::to_string(tmax)+"/";
  std::filesystem::create_directories(dir_out);
  std::cout << "# dir_out = " << dir_out << std::endl;

  // ---- host buffers (Sec. 2.2 symbols + axial GW intermediate chi) -------------------------
  //   phi  := phi'   (shared forward leg D_m^{-1} eta),   chi := (1 - D_ov) phi'  (shared GW factor),
  //   psi  := psi_n  (fixed-0 sink leg),   phit := phi_n(t) (looped current),
  //   rho  := the K-applied (and GW-applied) source,   tmp := the preconditioned CG RHS.
  FermionVector eta, phi, chi, tmp, rho, psi, phit;

  // free field: single deterministic config (k=0), U=1.  ensemble: loop ckpoint_lat.k in ens_dir.
  const int k_ckpoint = free_field ? 1 : 10;
  const int kmax      = free_field ? 0 : 1000;

  for(int k = 0; k <= kmax; k += k_ckpoint){
    if(!free_field){
      const std::string str_lat = ens_dir + "ckpoint_lat." + std::to_string(k);
      if(!std::filesystem::exists(str_lat)){ if(k==0) continue; else break; }
      U.read(str_lat);
    }
    D.update(U);  Dm.update(U);  Dtil.update(U);
    std::cout << "# k="<<k<<(free_field?" (free field)":"")
              << "  lambda_min/max="<<Dm.lambda_min<<"/"<<Dm.lambda_max<<std::endl;

    for(int h=0; h<nhits; h++){
      // one stochastic hit -> one file conn.<k>.<h>.h5 (single-realization estimate; average in post)
      const std::string h5path = dir_out + "conn." + std::to_string(k) + "." + std::to_string(h) + ".h5";
      if(std::filesystem::exists(h5path)){ std::cout<<"# skip k="<<k<<" h="<<h<<" (done)"<<std::endl; continue; }
      std::cout << "# k="<<k<<" hit "<<(h+1)<<"/"<<nhits<<" : start ("<<n_sites<<" sites x "<<Nt<<" t)" << std::endl;
      const auto t_hit0 = std::chrono::steady_clock::now();

      eta.fill_z2_source(rng);
      std::vector<double> pmRe(Nt), pmIm(Nt);

      // shared legs:  phi' = X^{-1} eta  (forward solve; massless D_ov in the flavor case, else D_m);
      //               chi  = (1 - D_ov) phi'  (massless GW factor, reused for all sites and t)
      if(flavor){
        op_DH.from_cpu<N>(tmp.field, eta.field);                 // tmp = D_ov^dag eta
        op_Dsq.solve<N>(phi.field, tmp.field, Comp::TOL_OUTER);  // phi := phi' = D_ov^{-1} eta (massless)
      } else {
        op_DmH.from_cpu<N>(tmp.field, eta.field);                // tmp = D_m^dag eta
        op_Dmsq.solve<N>(phi.field, tmp.field, Comp::TOL_OUTER); // phi := phi' = D_m^{-1} eta
      }
      op_oneMinusD.from_cpu<N>(chi.field, phi.field);            // chi := (1 - D_ov) phi'
      std::cout << "#   phi' and chi = (1-D_ov) phi' solved (reused for all sites)" << std::endl;

      std::vector<Complex> Apm(Nt, Complex(0.0,0.0));
      // explicit loop over spatial sites (site-diagonal temporal projection, Eq. 4.32)
      for(int n=0; n<n_sites; n++){
        std::cout << "#   site "<<(n+1)<<"/"<<n_sites<<" : sink solve + "<<Nt<<" kernel applies" << std::endl;
        // fixed at (n,0): psi_n = [D_ov (flavor) | tilde (parity) | D_m]^{-1} (1-D_ov) K^dag(n,0) eta
        kop.set_temporal(U, 0, (Idx)n, /*dag=*/true);
        op_K.from_cpu<N>(rho.field, eta.field);                  // rho = K^dag(n,0) eta
        op_oneMinusD.from_cpu<N>(rho.field, rho.field);          // rho = (1 - D_ov) rho   (in place)
        if(flavor){
          op_DH.from_cpu<N>(tmp.field, rho.field);               // tmp = D_ov^dag rho
          op_Dsq.solve<N>(psi.field, tmp.field, Comp::TOL_OUTER);       // psi := D_ov^{-1} (...) (massless)
        } else if(parity){
          op_tilDmH.from_cpu<N>(tmp.field, rho.field);           // tmp = tilde^dag rho
          op_tilDmsq.solve<N>(psi.field, tmp.field, Comp::TOL_OUTER);   // psi := tilde^{-1} (...)
        } else {
          op_DmH.from_cpu<N>(tmp.field, rho.field);              // tmp = D_m^dag rho
          op_Dmsq.solve<N>(psi.field, tmp.field, Comp::TOL_OUTER);      // psi := D_m^{-1} (...)
        }
        // looped kernel (no inversion): phi_n(t) = K^dag(n,t) chi ;  Apm[t] += w_n psi_n^dag phi_n(t)
        for(int t=0; t<Nt; t++){
          if(t>tmax && t<Nt-tmax) continue;   // tmax cap: skip the noise middle (both signal ends kept)
          kop.set_temporal(U, t, (Idx)n, /*dag=*/true);
          op_K.from_cpu<N>(phit.field, chi.field);               // phit = K^dag(n,t) chi
          Apm[t] += w_tp[n] * psi.dag(phit);
        }
        const double se = std::chrono::duration<double>(std::chrono::steady_clock::now()-t_hit0).count();
        std::cout << "#   site "<<(n+1)<<"/"<<n_sites<<" done  ["<<se<<" s]" << std::endl;
      }
      for(int t=0;t<Nt;t++){ const Complex g = inv4pi*Apm[t]; pmRe[t]=g.real(); pmIm[t]=g.imag(); }

      HighFive::File h5(h5path, HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
      h5.createDataSet("Apm/real", pmRe);   // C_{A+-} connected (mixed ordering)
      h5.createDataSet("Apm/imag", pmIm);
      h5.createDataSet("tmax", std::vector<int>{tmax});   // computed region = [0,tmax] U [Nt-tmax,Nt-1]
      const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now()-t_hit0).count();
      std::cout << "# wrote "<<h5path<<"  Apm(dt=0)=("<<pmRe[0]<<","<<pmIm[0]<<")  ["<<secs<<" s]"<<std::endl;
    } // hits
  } // k

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();
  return 0;
}
