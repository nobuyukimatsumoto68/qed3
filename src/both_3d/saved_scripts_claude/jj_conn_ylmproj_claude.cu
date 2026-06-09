// jj_conn_ylmproj_claude.cu
// Connected current-current correlator, SPHERICAL-HARMONIC m-summed TOWER G_l(t)
// (Eq. 4.36 of qed3int_v2-10.pdf).  Plan: conserved_current_correlators_impl_plan_v3_claude.md
// Sec. 3.6.  VECTOR only this pass ((++) and parity (--); axial ylm later).
//
// ylm folds the real Y_{lm} over ALL sites at BOTH currents.  The correlator is diagonal in (l,m), so
// (Sec. 3.6) we keep l_1=l_2=l and SUM both m_1,m_2 over -l..l => a TOWER G_l(t), l=0,1,2 (n_ell=3).
// Because the inner product is bilinear and the inverse is linear, the double m-sum collapses BEFORE the
// solve, giving ONE inversion per l (not per (l,m)):
//   G_l(t) = sum_{m1,m2} G_{(l m1)(l m2)}(t) = ( sum_{m1} psi_{(l m1)} )^dag ( sum_{m2} Phi_{(l m2)} )
//          = Psi_l^dag PhiL_l ,   with the COMBINED real weight
//   W_ell[l][n] = (A_n / kappa_t(n)) * sum_{m=-l}^{l} Y_lm(n^)   (kappa^1, Eq. 4.36; NOT kappa^2).
//   (Real Y_{lm} via the free Ylm_real(l,m,VE); NOT FermionVector::mult_Ylm_real -- the overlap kernel
//    spreads, so each per-site kernel vector must be weighted by its CURRENT-site scalar.)
// The pre-summed estimator includes off-diagonal m1!=m2 (vanish in expectation) -> unbiased estimate of
// the diagonal trace sum_m G_{(lm)(lm)}, carrying only extra off-diagonal noise (l=0 ~ 0, a check).
// Estimator (Sec. 2.2 op-set + op_K):
//   phi' = D_m^{-1} eta  (shared)
//   Psi_l = D_m^{-dag} sum_n W_ell[l][n] K^dag(n,t0) eta   (source; ONE inversion / l)
//   PhiL_l(t) = sum_n W_ell[l][n] K(n,t) phi'              (looped kernel; no inversion)
//   G^{++}_l(t) = (1/4pi) Psi_l^dag PhiL_l(t)
// (--) = operator-adjoint mirror: massless/m_F -> G^{--}_l=conj(G^{++}_l); parity -> tilde.
//
// Solver: overlap multishift (*_deviceAsyncLaunch_ms) + the kernel is the multishift apply_k_ms via
// ConservedCurrent::operator() (inherited through op_K.from_cpu).  --n-t0 sets the source origins; the
// full dt=0..Nt-1 is computed and every (h,t0) tower is stored raw (averaged downstream).

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
  printf("  --current <c>        vector (this binary; default: vector)\n");
  printf("  --ens-dir <path>     sea config directory; OMIT => free field (U=1) check\n");
  printf("  --nhits <n>          stochastic hits (default: 1)\n");
  printf("  --n-t0 <N>           number of source-time origins t0=b*(Nt/N), full dt-loop (default: 2)\n");
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
  int n_t0=2;    // number of source-time origins t0=b*(Nt/n_t0); full dt-loop, no in-program averaging

  ParseArgs(argc, argv, gsq, Nf, nu0, nu1, mass_re, mass_im, current, ens_dir, nhits, n_t0);
  if(nu1 < 0.0) nu1 = nu0;    // valence asymmetry defaults to the sea value nu0 (knob retained)

  const Complex valence_mass(mass_re, mass_im);
  const bool free_field = ens_dir.empty();
  assert(current=="vector");   // axial ylm is a later pass

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

  // source-time origins (Sec. 3.7): n_t0 evenly-spaced t0=b*(Nt/n_t0); FULL dt=0..Nt-1; every (h,t0)
  // stored raw (no in-program averaging) -- the notebook averages over h and t0 per config, then jackknifes.
  assert(n_t0 >= 1 && Nt % n_t0 == 0 && "Nt must be divisible by n_t0");
  const int t0_spacing = Nt / n_t0;
  std::vector<int> t0s(n_t0);
  for(int b=0; b<n_t0; b++) t0s[b] = b*t0_spacing;
  std::cout << "# n_t0=" << n_t0 << " source origins (t0=b*"<<t0_spacing<<"), full dt-loop, raw per (h,t0)" << std::endl;

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

  // Overlap operators (Dtil built unconditionally; unused unless parity, but harmless):
  //   Dm   = D_ov + m           ((++) uses it on both legs)
  //   Dtil = D_ov + m/(1-m)     (\tilde D_{m_P}; parity (-) dagger leg)
  Fermion Dm  (DW, valence_mass, 21);
  Fermion Dtil(DW, valence_mass / (Complex(1.0) - valence_mass), 21);
  std::cout << "# overlap operators set: D_m, tilde D_{m_P} (M5="<<M5<<")." << std::endl;

  ConservedCurrent<Fermion,Gauge> kop(Dm);   // K is mass-independent; multishift apply_k_ms via operator()
  MatPoly op_K; op_K.push_back(cplx(1.0), {&kop});

  // Uniform operator set bound to the multishift entry points (*_deviceAsyncLaunch_ms; ~4x):
  //   X^{-1} b = op_XH (RHS X^dag b) + op_Xsq (CG);  X^{-dag} b = op_X (RHS X b) + op_Xsq.
  auto f_Dm   = std::bind(&Fermion::mult_deviceAsyncLaunch_ms, &Dm, std::placeholders::_1, std::placeholders::_2);
  auto f_DmH  = std::bind(&Fermion::adj_deviceAsyncLaunch_ms,  &Dm, std::placeholders::_1, std::placeholders::_2);
  auto f_Dmsq = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms,  &Dm, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_Dm(f_Dm), M_DmH(f_DmH), M_Dmsq(f_Dmsq);
  MatPoly op_Dm;   op_Dm.push_back(cplx(1.0), {&M_Dm});
  MatPoly op_DmH;  op_DmH.push_back(cplx(1.0), {&M_DmH});
  MatPoly op_Dmsq; op_Dmsq.push_back(cplx(1.0), {&M_Dmsq});

  auto f_tilDm   = std::bind(&Fermion::mult_deviceAsyncLaunch_ms, &Dtil, std::placeholders::_1, std::placeholders::_2);
  auto f_tilDmH  = std::bind(&Fermion::adj_deviceAsyncLaunch_ms,  &Dtil, std::placeholders::_1, std::placeholders::_2);
  auto f_tilDmsq = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms,  &Dtil, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_tilDm(f_tilDm), M_tilDmH(f_tilDmH), M_tilDmsq(f_tilDmsq);
  MatPoly op_tilDm;   op_tilDm.push_back(cplx(1.0), {&M_tilDm});
  MatPoly op_tilDmH;  op_tilDmH.push_back(cplx(1.0), {&M_tilDmH});
  MatPoly op_tilDmsq; op_tilDmsq.push_back(cplx(1.0), {&M_tilDmsq});

  // ---- spherical-harmonic m-summed tower (l_max=2, keep l=0): n_ell=3 towers l=0,1,2 -------------
  constexpr int L_MAX = 2;
  constexpr int N_ELL = L_MAX + 1;            // l = 0,1,2  (l=0 ~ 0 by charge conservation, a check)
  std::vector<int> ls(N_ELL); for(int l=0;l<N_ELL;l++) ls[l]=l;

  const int n_sites = static_cast<int>(base.n_sites);
  const double inv4pi = 1.0/(4.0*std::acos(-1.0));

  // combined real weight W_ell[l][n] = (A_n / kappa_t(n)) * sum_{m=-l}^{l} Y_lm(n^)  (Sec. 3.6; kappa^1).
  // Summing m here is what collapses the 2l+1 (l,m)-channel solves into ONE solve per l.
  std::vector<std::vector<double>> W_ell(N_ELL, std::vector<double>(n_sites, 0.0));
  for(int l=0; l<N_ELL; l++)
    for(int n=0; n<n_sites; n++){
      const double kt = DW.kappa_t[n];
      double s = 0.0;
      for(int m=-l; m<=l; m++) s += Ylm_real(l, m, base.sites[n]);
      W_ell[l][n] = base.dual_areas[n] * s / kt;
    }

  // output: data_<ESNID>/ylm_<current>_nt0<N>_nhits<H>/conn.<config>.h5
  std::string ens_base = ens_dir;
  if(!ens_base.empty() && ens_base.back()=='/') ens_base.pop_back();
  { const auto slash = ens_base.find_last_of('/'); if(slash!=std::string::npos) ens_base = ens_base.substr(slash+1); }
  const std::string esnid = (free_field ? std::string("free") : ens_base)
                          + "_vmRe"+std::to_string(mass_re)+"vmIm"+std::to_string(mass_im);
  const std::string dir_out = "data_"+esnid+"/ylm_"+current+"_nt0"+std::to_string(n_t0)+"_nhits"+std::to_string(nhits)+"/";
  std::filesystem::create_directories(dir_out);
  std::cout << "# dir_out = " << dir_out << "  (n_ell="<<N_ELL<<" towers l=0,1,2)" << std::endl;

  // ---- host buffers ----------------------------------------------------------------------------
  //   eta = source; phi = phi' (shared forward leg); rho = K-applied source vector (per site);
  //   kphi = K-applied looped vector (per site); tmp = preconditioned CG RHS.
  //   srcL[l] = sum_n W_ell[l][n] K(n,t0) eta (accumulated); psiL[l] = inverted source tower (held over t);
  //   PhiL[l] = sum_n W_ell[l][n] K(n,t) phi' (per t).  Tower arrays use std::array (FermionVector has no
  //   safe copy ctor, so never std::vector-push_back them).
  FermionVector eta, phi, rho, kphi, tmp;
  std::array<FermionVector, N_ELL> srcL, psiL, PhiL;

  // free field: single deterministic config (k=0), U=1.  ensemble: loop ckpoint_lat.k in ens_dir.
  const int k_ckpoint = free_field ? 1 : 10;
  const int kmax      = free_field ? 0 : 1000;

  // helper: write the per-l tower vectors (Re/Im) under /h{h}/t0_{b}/<chan>/l{l}/{real,imag}
  auto write_tower = [&](HighFive::File& h5, const std::string& kp, const std::string& chan,
                         const std::vector<std::vector<double>>& Re,
                         const std::vector<std::vector<double>>& Im){
    for(int l=0;l<N_ELL;l++){
      const std::string key = kp + chan + "/l" + std::to_string(ls[l]) + "/";
      h5.createDataSet(key+"real", Re[l]);
      h5.createDataSet(key+"imag", Im[l]);
    }
  };

  for(int k = 0; k <= kmax; k += k_ckpoint){
    if(!free_field){
      const std::string str_lat = ens_dir + "ckpoint_lat." + std::to_string(k);
      if(!std::filesystem::exists(str_lat)){ if(k==0) continue; else break; }
      U.read(str_lat);
    }
    // ONE file per config (Sec. 3.7).  Resume: skip ONLY if the file is COMPLETE -- HighFive writes
    // datasets sequentially, so the LAST key present implies the whole config is there; an interrupted
    // write leaves it absent -> recompute.  Read-only open (never truncates a maybe-good file).
    // Last write = Vmm of the last (h,t0), last tower l=ls.back().
    const std::string h5path = dir_out + "conn." + std::to_string(k) + ".h5";
    if(std::filesystem::exists(h5path)){
      bool complete = false;
      try {
        HighFive::File h5check(h5path, HighFive::File::ReadOnly);
        const std::string last_ds = "h"+std::to_string(nhits-1)+"/t0_"+std::to_string(n_t0-1)
                                   + "/Vmm/l" + std::to_string(ls.back()) + "/imag";
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
    h5.createDataSet("ls", ls);           // tower l list (0,1,2)

    for(int h=0; h<nhits; h++){
      const auto t_hit0 = std::chrono::steady_clock::now();
      std::cout << "# k="<<k<<" hit "<<(h+1)<<"/"<<nhits<<" : "<<n_t0<<" origins x "<<N_ELL
                << " tower-inversions x "<<n_sites<<" sites x "<<Nt<<" dt" << std::endl;
      eta.fill_z2_source(rng);
      const std::string hp = "h" + std::to_string(h) + "/";   // key prefix /h{h}/

      // (++) tower: phi'_++ = D_m^{-1} eta (shared across origins)
      op_DmH.from_cpu<N>(tmp.field, eta.field);
      op_Dmsq.solve<N>(phi.field, tmp.field, Comp::TOL_OUTER);

      for(int b=0;b<n_t0;b++){
        const int t0 = t0s[b];
        const std::string kp = hp + "t0_" + std::to_string(b) + "/";   // /h{h}/t0_{b}/
        std::vector<std::vector<Complex>> Gpp(N_ELL, std::vector<Complex>(Nt, Complex(0.0,0.0)));
        // source towers at t0: srcL[l] = sum_n W_ell[l][n] K^dag(n,t0) eta ;  psiL[l] = D_m^{-dag} srcL[l]
        for(int l=0;l<N_ELL;l++) memset(srcL[l].field, 0, Comp::N*CD);
        for(int n=0;n<n_sites;n++){
          kop.set_temporal(U, t0, (Idx)n, /*dag=*/true);
          op_K.from_cpu<N>(rho.field, eta.field);                 // rho = K^dag(n,t0) eta
          for(int l=0;l<N_ELL;l++){ const double w=W_ell[l][n]; for(Idx i=0;i<N;i++) srcL[l].field[i]+=w*rho.field[i]; }
        }
        for(int l=0;l<N_ELL;l++){
          op_Dm.from_cpu<N>(tmp.field, srcL[l].field);            // tmp = D_m srcL[l]
          op_Dmsq.solve<N>(psiL[l].field, tmp.field, Comp::TOL_OUTER);  // psiL[l] = D_m^{-dag} srcL[l]
        }
        // full dt: PhiL[l] = sum_n W_ell[l][n] K(n,(t0+dt)%Nt) phi' ;  Gpp_l[dt] += psiL[l]^dag PhiL[l]
        for(int dt=0; dt<Nt; dt++){
          const int tsink = (t0+dt)%Nt;
          for(int l=0;l<N_ELL;l++) memset(PhiL[l].field, 0, Comp::N*CD);
          for(int n=0;n<n_sites;n++){
            kop.set_temporal(U, tsink, (Idx)n, /*dag=*/false);
            op_K.from_cpu<N>(kphi.field, phi.field);              // kphi = K(n,tsink) phi'
            for(int l=0;l<N_ELL;l++){ const double w=W_ell[l][n]; for(Idx i=0;i<N;i++) PhiL[l].field[i]+=w*kphi.field[i]; }
          }
          for(int l=0;l<N_ELL;l++) Gpp[l][dt] += psiL[l].dag(PhiL[l]);
        }
        std::vector<std::vector<double>> ppRe(N_ELL, std::vector<double>(Nt,0.0)), ppIm(N_ELL, std::vector<double>(Nt,0.0));
        for(int l=0;l<N_ELL;l++) for(int t=0;t<Nt;t++){ const Complex g=inv4pi*Gpp[l][t]; ppRe[l][t]=g.real(); ppIm[l][t]=g.imag(); }
        write_tower(h5, kp, "Vpp", ppRe, ppIm);
        // (--) for massless / m_F: G^{--}_l=conj(G^{++}_l) (parity handled in its own loop below)
        if(!parity){
          std::vector<std::vector<double>> mmRe(N_ELL, std::vector<double>(Nt,0.0)), mmIm(N_ELL, std::vector<double>(Nt,0.0));
          for(int l=0;l<N_ELL;l++) for(int t=0;t<Nt;t++){ mmRe[l][t]=ppRe[l][t]; mmIm[l][t]=-ppIm[l][t]; }
          write_tower(h5, kp, "Vmm", mmRe, mmIm);
        }
      }

      // (--) parity tower: independent via tilde.  phi'_-- = tilde^{-dag} eta (shared across origins)
      if(parity){
        op_tilDm.from_cpu<N>(tmp.field, eta.field);
        op_tilDmsq.solve<N>(phi.field, tmp.field, Comp::TOL_OUTER);
        for(int b=0;b<n_t0;b++){
          const int t0 = t0s[b];
          const std::string kp = hp + "t0_" + std::to_string(b) + "/";
          std::vector<std::vector<Complex>> Gmm(N_ELL, std::vector<Complex>(Nt, Complex(0.0,0.0)));
          // source at t0: srcL[l] = sum_n W_ell[l][n] K(n,t0) eta ;  psiL[l] = tilde^{-1} srcL[l]  (forward solve)
          for(int l=0;l<N_ELL;l++) memset(srcL[l].field, 0, Comp::N*CD);
          for(int n=0;n<n_sites;n++){
            kop.set_temporal(U, t0, (Idx)n, /*dag=*/false);
            op_K.from_cpu<N>(rho.field, eta.field);                 // rho = K(n,t0) eta
            for(int l=0;l<N_ELL;l++){ const double w=W_ell[l][n]; for(Idx i=0;i<N;i++) srcL[l].field[i]+=w*rho.field[i]; }
          }
          for(int l=0;l<N_ELL;l++){
            op_tilDmH.from_cpu<N>(tmp.field, srcL[l].field);         // tmp = tilde^dag srcL[l]
            op_tilDmsq.solve<N>(psiL[l].field, tmp.field, Comp::TOL_OUTER);  // psiL[l] = tilde^{-1} srcL[l]
          }
          for(int dt=0; dt<Nt; dt++){
            const int tsink = (t0+dt)%Nt;
            for(int l=0;l<N_ELL;l++) memset(PhiL[l].field, 0, Comp::N*CD);
            for(int n=0;n<n_sites;n++){
              kop.set_temporal(U, tsink, (Idx)n, /*dag=*/true);
              op_K.from_cpu<N>(kphi.field, phi.field);              // kphi = K^dag(n,tsink) phi'_--
              for(int l=0;l<N_ELL;l++){ const double w=W_ell[l][n]; for(Idx i=0;i<N;i++) PhiL[l].field[i]+=w*kphi.field[i]; }
            }
            for(int l=0;l<N_ELL;l++) Gmm[l][dt] += psiL[l].dag(PhiL[l]);
          }
          std::vector<std::vector<double>> mmRe(N_ELL, std::vector<double>(Nt,0.0)), mmIm(N_ELL, std::vector<double>(Nt,0.0));
          for(int l=0;l<N_ELL;l++) for(int t=0;t<Nt;t++){ const Complex g=inv4pi*Gmm[l][t]; mmRe[l][t]=g.real(); mmIm[l][t]=g.imag(); }
          write_tower(h5, kp, "Vmm", mmRe, mmIm);
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
