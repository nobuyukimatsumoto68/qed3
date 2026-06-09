// jj_conn_spproj_claude.cu
// Connected current-current correlator, SPATIAL equal-link projection G_s (Eq. 4.29 of
// qed3int_v2-10.pdf).  Plan: conserved_current_correlators_impl_plan_v3_claude.md Sec. 3.5 --
// the spatial sibling of jj_conn_tpproj_claude.cu (temporal).  Differences from tp (Sec. 3.5):
//   - insertion points are SPATIAL LINKS (loop base.links; kernel via kop.set_spatial, not set_temporal);
//   - weight w_{nn'} = A_{nn'}/kappa^{(0)2}_{nn'} = base.link_volume[il]/(DW.bd.kappa[il])^2 (Eq. 4.29);
//   - output data_<ESNID>/sp_<current>/.
// Everything else (operator set, fix-0/loop-t estimator, (--)=adjoint-mirror, buffers, CLI) is identical.
//
// Solver: binds the overlap multi-shift device entry points (*_deviceAsyncLaunch_ms; ~4x over the
// per-pole loop) -- the default for new programs.  Multi-RHS batching of the 1+n_links sink solves
// is a later pass; the insertion-point loop is the batching lever.
//
// Conventions (v3): valence mass via --mass-re/--mass-im (D_m = D_ov + m); ensemble selected by
// --ens-dir (omit => free field, U = 1); kernel K is always the massless-form Noether kernel.
// Parity ( m purely imaginary ) dagger leg uses \tilde D_{m_P} = D_ov + m_P/(1-m_P).

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
#include "conserved_current_claude.h"   // ConservedCurrent: apply_k / apply_k_dag; BaseLink

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
  printf("  --current <c>        vector | axial (default: vector)\n");
  printf("  --ens-dir <path>     sea config directory; OMIT => free field (U=1) check\n");
  printf("  --nhits <n>          stochastic hits (default: 1)\n");
  printf("  --n-t0 <N>           number of source-time origins t0=b*(Nt/N), b=0..N-1 (default: 2)\n");
  printf("  -h, --help           show this help\n");
  exit(0);
}

void ParseArgs(int argc, char* argv[],
               double& gsq, int& Nf, double& nu0, double& nu1,
               double& mass_re, double& mass_im,
               std::string& current, std::string& ens_dir, int& nhits, int& n_t0){
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
    {"n-t0",    required_argument, nullptr, 'T'},
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
    case 'T': n_t0    = std::stoi(optarg); break;
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
  std::string current="vector";
  std::string ens_dir="";     // empty => free-field mode
  int nhits=1;
  int n_t0=2;   // number of source-time origins (Sec. 3.7)

  ParseArgs(argc, argv, gsq, Nf, nu0, nu1, mass_re, mass_im, current, ens_dir, nhits, n_t0);
  if(nu1 < 0.0) nu1 = nu0;    // valence asymmetry defaults to the sea value nu0 (knob retained)

  const Complex valence_mass(mass_re, mass_im);
  const bool free_field = ens_dir.empty();
  assert(current=="vector" || current=="axial");

  // parity case: purely imaginary valence mass -> dagger leg uses \tilde D_{m_P}
  const bool parity = (std::abs(mass_im) > 1.0e-15) && (std::abs(mass_re) <= 1.0e-15);

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

  // source-time origins (Sec. 3.7): n_t0 evenly-spaced t0 = b*(Nt/n_t0); FULL dt-loop dt=0..Nt-1;
  // store every (hit, t0) raw (no in-program averaging); average + jackknife downstream.
  assert(n_t0 >= 1 && Nt % n_t0 == 0 && "Nt must be divisible by n_t0");
  const int t0_spacing = Nt / n_t0;
  std::vector<int> t0s(n_t0);
  for(int b=0; b<n_t0; b++) t0s[b] = b*t0_spacing;
  std::cout << "# n_t0=" << n_t0 << " source origins (t0=b*"<<t0_spacing<<"), full dt; one file/config, keys per (h,t0)" << std::endl;

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

  // Overlap operators (Dtil is built unconditionally; unused unless parity, but harmless):
  //   Dm   = D_ov + m           (Eq. 3.60; (++) uses it on both legs, Eq. 3.64)
  //   Dtil = D_ov + m/(1-m)     (\tilde D_{m_P}, Eq. 3.63 first form; parity (-) dagger leg)
  // Fermion D   (DW, Complex(0.0), 21);   // massless -- axial-only prototype
  Fermion Dm  (DW, valence_mass, 21);
  Fermion Dtil(DW, valence_mass / (Complex(1.0) - valence_mass), 21);
  std::cout << "# overlap operators set: D_m, tilde D_{m_P} (M5="<<M5<<")." << std::endl;

  ConservedCurrent<Fermion,Gauge> kop(Dm);   // K is mass-independent; LinOp configured per (link,dag)
  MatPoly op_K; op_K.push_back(cplx(1.0), {&kop});   // apply K via op_K.from_cpu (no raw device buffers)

  // Uniform operator set: for each overlap (massless D, D_m, tilde D) the three pieces
  //   mult = X,  H = X^dag,  sq = X^dag X = X X^dag  (equal since D_ov+m is normal; fused DDH=DHD).
  // Solves:  X^{-1} b  via op_XH (RHS X^dag b) + op_Xsq (CG);  X^{-dag} b via op_X (RHS X b) + op_Xsq.
  // Bound to the MULTI-SHIFT entry points (*_deviceAsyncLaunch_ms; ~4x), Sec. 3.5.
  // massless D_ov -- UNUSED by the vector currents (axial-case prototype).

  // D_m = D_ov + m
  auto f_Dm   = std::bind(&Fermion::mult_deviceAsyncLaunch_ms, &Dm, std::placeholders::_1, std::placeholders::_2);
  auto f_DmH  = std::bind(&Fermion::adj_deviceAsyncLaunch_ms,  &Dm, std::placeholders::_1, std::placeholders::_2);
  auto f_Dmsq = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms,  &Dm, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_Dm(f_Dm), M_DmH(f_DmH), M_Dmsq(f_Dmsq);
  MatPoly op_Dm;   op_Dm.push_back(cplx(1.0), {&M_Dm});
  MatPoly op_DmH;  op_DmH.push_back(cplx(1.0), {&M_DmH});
  MatPoly op_Dmsq; op_Dmsq.push_back(cplx(1.0), {&M_Dmsq});

  // tilde D_{m_P} = D_ov + m/(1-m)
  auto f_tilDm   = std::bind(&Fermion::mult_deviceAsyncLaunch_ms, &Dtil, std::placeholders::_1, std::placeholders::_2);
  auto f_tilDmH  = std::bind(&Fermion::adj_deviceAsyncLaunch_ms,  &Dtil, std::placeholders::_1, std::placeholders::_2);
  auto f_tilDmsq = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms,  &Dtil, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_tilDm(f_tilDm), M_tilDmH(f_tilDmH), M_tilDmsq(f_tilDmsq);
  MatPoly op_tilDm;   op_tilDm.push_back(cplx(1.0), {&M_tilDm});
  MatPoly op_tilDmH;  op_tilDmH.push_back(cplx(1.0), {&M_tilDmH});
  MatPoly op_tilDmsq; op_tilDmsq.push_back(cplx(1.0), {&M_tilDmsq});

  // geometry weights for the spatial projection (Eq. 4.29): w_{nn'} = A_{nn'} / kappa^{(0)}_{nn'}^2
  // indexed by link il = base.map2il.at(lk);  A_{nn'} = link_volume[il], kappa = DW.bd.kappa[il].
  std::vector<double> w_sp(base.n_links);
  for(Idx il=0; il<base.n_links; il++){ const double ks=DW.bd.kappa[il]; w_sp[il]=base.link_volume[il]/(ks*ks); }
  const double inv4pi = 1.0/(4.0*std::acos(-1.0));
  const int n_links = static_cast<int>(base.links.size());

  // output: data_<ESNID>/sp_<current>/conn.<config>.<hit>.h5
  std::string ens_base = ens_dir;
  if(!ens_base.empty() && ens_base.back()=='/') ens_base.pop_back();
  { const auto slash = ens_base.find_last_of('/'); if(slash!=std::string::npos) ens_base = ens_base.substr(slash+1); }
  const std::string esnid = (free_field ? std::string("free") : ens_base)
                          + "_vmRe"+std::to_string(mass_re)+"vmIm"+std::to_string(mass_im);
  const std::string dir_out = "data_"+esnid+"/sp_"+current+"_nt0"+std::to_string(n_t0)+"_nhits"+std::to_string(nhits)+"/";
  std::filesystem::create_directories(dir_out);
  std::cout << "# dir_out = " << dir_out << std::endl;

  if(current!="vector"){
    std::cout << "# NOTE: only --current vector is implemented in this build; axial is "
              << "jj_conn_spproj_axial_claude.cu. Exiting cleanly." << std::endl;
    for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
    return 0;
  }

  // ---- host buffers (shared by BOTH channels; Sec. 3.2 symbols, no +-/subscript) ----------
  //   phi  := phi'   (forward leg, reused over links and t),   psi := psi_a (fixed-0 leg),
  //   phit := phi(t) (looped current),   rho := the K-applied source,
  //   tmp  := the preconditioned CG RHS (an operator X^{(dag)} applied to eta or rho, either channel).
  FermionVector eta, phi, tmp, rho, psi, phit;

  // free field: single deterministic config (k=0), U=1.  ensemble: loop ckpoint_lat.k in ens_dir.
  const int k_ckpoint = free_field ? 1 : 10;
  const int kmax      = free_field ? 0 : 1000;

  for(int k = 0; k <= kmax; k += k_ckpoint){
    if(!free_field){
      const std::string str_lat = ens_dir + "ckpoint_lat." + std::to_string(k);
      if(!std::filesystem::exists(str_lat)){ if(k==0) continue; else break; }
      U.read(str_lat);
    }
    // ONE file per config (Sec. 3.7).  Resume: skip ONLY if the file is COMPLETE -- HighFive writes
    // datasets sequentially, so the LAST key present implies the whole config is there; an interrupted
    // write leaves it absent -> recompute.  Read-only open (never truncates a maybe-good file).
    const std::string h5path = dir_out + "conn." + std::to_string(k) + ".h5";
    if(std::filesystem::exists(h5path)){
      bool complete = false;
      try {
        HighFive::File h5check(h5path, HighFive::File::ReadOnly);
        const std::string last_ds = "h"+std::to_string(nhits-1)+"/t0_"+std::to_string(n_t0-1)+"/Vmm/imag";
        complete = h5check.exist(last_ds);
      } catch(...) {}
      if(complete){ std::cout<<"# skip k="<<k<<" (complete)"<<std::endl; continue; }
      std::cout<<"# k="<<k<<" exists but INCOMPLETE -> recompute"<<std::endl;
    }
    Dm.update(U);  Dtil.update(U);
    std::cout << "# k="<<k<<(free_field?" (free field)":"")
              << "  lambda_min/max="<<Dm.lambda_min<<"/"<<Dm.lambda_max<<std::endl;

    HighFive::File h5(h5path, HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
    h5.createDataSet("t0s",   t0s);
    h5.createDataSet("n_t0",  std::vector<int>{n_t0});
    h5.createDataSet("nhits", std::vector<int>{nhits});

    for(int h=0; h<nhits; h++){
      const auto t_hit0 = std::chrono::steady_clock::now();
      std::cout << "# k="<<k<<" hit "<<(h+1)<<"/"<<nhits<<" : "<<n_t0<<" origins x "<<n_links<<" links x "<<Nt<<" dt" << std::endl;
      eta.fill_z2_source(rng);
      const std::string hp = "h" + std::to_string(h) + "/";   // key prefix /h{h}/

      // (++) channel: phi'_++ = D_m^{-1} eta (shared across origins/links)
      op_DmH.from_cpu<N>(tmp.field, eta.field);
      op_Dmsq.solve<N>(phi.field, tmp.field, Comp::TOL_OUTER);
      for(int b=0; b<n_t0; b++){
        const int t0 = t0s[b];
        std::vector<Complex> Cpp(Nt, Complex(0.0,0.0));
        for(int a=0; a<n_links; a++){
          const BaseLink lk = base.links[a];
          const Idx il = base.map2il.at(lk);
          // psi_++ = D_m^{-dag} K^dag(lk,t0) eta
          kop.set_spatial(U, t0, lk, /*dag=*/true);
          op_K.from_cpu<N>(rho.field, eta.field);
          op_Dm.from_cpu<N>(tmp.field, rho.field);
          op_Dmsq.solve<N>(psi.field, tmp.field, Comp::TOL_OUTER);
          for(int dt=0; dt<Nt; dt++){
            kop.set_spatial(U, (t0+dt)%Nt, lk, /*dag=*/false);
            op_K.from_cpu<N>(phit.field, phi.field);
            Cpp[dt] += w_sp[il] * psi.dag(phit);
          }
        }
        std::vector<double> ppRe(Nt), ppIm(Nt);
        for(int dt=0;dt<Nt;dt++){ const Complex g=inv4pi*Cpp[dt]; ppRe[dt]=g.real(); ppIm[dt]=g.imag(); }
        const std::string kp = hp + "t0_" + std::to_string(b) + "/";   // /h{h}/t0_{b}/
        h5.createDataSet(kp+"Vpp/real", ppRe);
        h5.createDataSet(kp+"Vpp/imag", ppIm);
        if(!parity){   // massless / m_F: Vmm = conj(Vpp)
          std::vector<double> mmRe(Nt), mmIm(Nt);
          for(int dt=0;dt<Nt;dt++){ mmRe[dt]=ppRe[dt]; mmIm[dt]=-ppIm[dt]; }
          h5.createDataSet(kp+"Vmm/real", mmRe);
          h5.createDataSet(kp+"Vmm/imag", mmIm);
        }
      }

      // (--) parity channel: phi'_-- = tilde^{-dag} eta; independent tilde computation per origin
      if(parity){
        op_tilDm.from_cpu<N>(tmp.field, eta.field);
        op_tilDmsq.solve<N>(phi.field, tmp.field, Comp::TOL_OUTER);
        for(int b=0; b<n_t0; b++){
          const int t0 = t0s[b];
          std::vector<Complex> Cmm(Nt, Complex(0.0,0.0));
          for(int a=0; a<n_links; a++){
            const BaseLink lk = base.links[a];
            const Idx il = base.map2il.at(lk);
            // psi_-- = tilde^{-1} K(lk,t0) eta  (forward solve; K, not K^dag)
            kop.set_spatial(U, t0, lk, /*dag=*/false);
            op_K.from_cpu<N>(rho.field, eta.field);
            op_tilDmH.from_cpu<N>(tmp.field, rho.field);
            op_tilDmsq.solve<N>(psi.field, tmp.field, Comp::TOL_OUTER);
            for(int dt=0; dt<Nt; dt++){
              kop.set_spatial(U, (t0+dt)%Nt, lk, /*dag=*/true);
              op_K.from_cpu<N>(phit.field, phi.field);
              Cmm[dt] += w_sp[il] * psi.dag(phit);
            }
          }
          std::vector<double> mmRe(Nt), mmIm(Nt);
          for(int dt=0;dt<Nt;dt++){ const Complex g=inv4pi*Cmm[dt]; mmRe[dt]=g.real(); mmIm[dt]=g.imag(); }
          const std::string kp = hp + "t0_" + std::to_string(b) + "/";
          h5.createDataSet(kp+"Vmm/real", mmRe);
          h5.createDataSet(kp+"Vmm/imag", mmIm);
        }
      }
      const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now()-t_hit0).count();
      std::cout << "#   hit "<<(h+1)<<" done ["<<secs<<" s]" << std::endl;
    } // hits
    std::cout << "# wrote " << h5path << std::endl;
  } // k

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();
  return 0;
}
