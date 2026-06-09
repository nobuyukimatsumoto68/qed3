// jj_disc_claude.cu
// DISCONNECTED current-current correlator -- standalone (plan: conserved_current_correlators_impl_plan_v3
// _claude.md Sec. 3.1).  Disc is treated SEPARATELY from the connected unified file (jj_correlators_
// claude.cu): the disc factorizes into single-time traces and is apply-only (one solve per hit), while the
// connected ties two propagators and must solve the current-dressed source -- nothing to share but phi'.
//
// Vector only (disc has no axial piece, Sec. 1).  Single-time trace, Z_2 noise, shared phi' = D_m^{-1} eta:
//   T(a,t) = tr[D_m^{-1} K(a,t)] ~ eta^dag K(a,t) phi'   (Eqs. 3.53, 3.65).
// PROJECTED single-current vectors are dumped (one Nt-vector per hit; RAW, no 1/4pi):
//   J_tp(t)  = sum_n w_tp[n]  T(n,t)        (temporal site-diagonal, Eq. 4.32)
//   J_ell(t) = sum_n W_ell[l][n] T(n,t)     (ylm m-summed tower, Sec. 3.6; same temporal pass)
//   J_sp(t)  = sum_l w_sp[il] T(l,t)        (spatial link, Eq. 4.29)
// Disconnected correlator downstream = (1/4pi) <J(t0) J(t1)>_conn; with nhits>1 pair DIFFERENT hits for the
// two factors to remove the same-hit connected bias.  Parity ( m purely imaginary ): dagger-leg tilde trace
//   \tilde T(a,t) = (K(a,t) tilphi)^dag eta,  tilphi = tilde D_{m_P}^{-1} eta  (massless/m_F: \tilde T = T^*).
//
// Solver: overlap multishift (*_deviceAsyncLaunch_ms); kernel via the multishift apply_k_ms (op_K.from_cpu).

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
  const double TOL_OUTER=1.0e-5;   // outer CG (only phi' here; plenty for the single-time traces)
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
#include "conserved_current_claude.h"   // ConservedCurrent: apply_k / apply_k_dag (multishift apply_k_ms)

//------------------------------------------
#include <getopt.h>

void PrintHelp(){
  printf("Usage: ./a.out [options]   (DISCONNECTED: projected single-current traces, vector only)\n");
  printf("  --gsq <gsq>          Wilson coupling squared (ensemble id; default: 8.0)\n");
  printf("  --Nf <Nf>            number of fermion flavors (ensemble id; default: 2)\n");
  printf("  --nu0 <nu0>          sea quark asymmetry (ensemble id; default: 1.0)\n");
  printf("  --nu1 <nu1>          valence Wilson-Dirac asymmetry (operator; default: nu0)\n");
  printf("  --mass-re <x>        valence mass Re (default: 0.0)\n");
  printf("  --mass-im <y>        valence mass Im (default: 0.0)\n");
  printf("  --ens-dir <path>     sea config directory; OMIT => free field (U=1) check\n");
  printf("  --nhits <n>          stochastic hits (default: 1)\n");
  printf("  --ninter <N>         ensemble config stride: measure ckpoint_lat.k for k=0,N,2N,... (default: 10)\n");
  printf("  -h, --help           show this help\n");
  exit(0);
}

void ParseArgs(int argc, char* argv[],
               double& gsq, int& Nf, double& nu0, double& nu1,
               double& mass_re, double& mass_im,
               std::string& ens_dir, int& nhits, int& ninter){
  static struct option long_opts[] = {
    {"gsq",     required_argument, nullptr, 'g'},
    {"Nf",      required_argument, nullptr, 'N'},
    {"nu0",     required_argument, nullptr, 'n'},
    {"nu1",     required_argument, nullptr, 'm'},
    {"mass-re", required_argument, nullptr, 'r'},
    {"mass-im", required_argument, nullptr, 'i'},
    {"ens-dir", required_argument, nullptr, 'e'},
    {"nhits",   required_argument, nullptr, 'H'},
    {"ninter",  required_argument, nullptr, 'I'},
    {"help",    no_argument,       nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
  };
  int opt, idx;
  while((opt = getopt_long(argc, argv, "g:N:n:m:r:i:e:H:I:h", long_opts, &idx)) != -1){
    switch(opt){
    case 'g': gsq     = std::stod(optarg); break;
    case 'N': Nf      = std::stoi(optarg); break;
    case 'n': nu0     = std::stod(optarg); break;
    case 'm': nu1     = std::stod(optarg); break;
    case 'r': mass_re = std::stod(optarg); break;
    case 'i': mass_im = std::stod(optarg); break;
    case 'e': ens_dir = optarg; break;
    case 'H': nhits   = std::stoi(optarg); break;
    case 'I': ninter  = std::stoi(optarg); break;
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
  std::string ens_dir="";     // empty => free-field mode
  int nhits=1;
  int ninter=10;   // ensemble config stride (k = 0, ninter, 2*ninter, ...)

  ParseArgs(argc, argv, gsq, Nf, nu0, nu1, mass_re, mass_im, ens_dir, nhits, ninter);
  if(nu1 < 0.0) nu1 = nu0;    // valence asymmetry defaults to the sea value nu0 (knob retained)

  const Complex valence_mass(mass_re, mass_im);
  const bool free_field = ens_dir.empty();

  // parity case: purely imaginary valence mass -> dagger-leg tilde trace uses \tilde D_{m_P}
  const bool parity = (std::abs(mass_im) > 1.0e-15) && (std::abs(mass_re) <= 1.0e-15);

  std::cout << "# gsq="<<gsq<<" Nf="<<Nf<<" nu0="<<nu0<<" nu1="<<nu1
            << " mass=("<<mass_re<<","<<mass_im<<")"
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

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  using Rng=ParallelRngExt<Base,Nt>;
  using Fermion=OverlapWMass<WilsonDirac>;

  // ---- operators (vector only: D_m, and tilde D_{m_P} for parity) ------------------------------
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

  Fermion Dm  (DW, valence_mass, 21);
  Fermion Dtil(DW, valence_mass / (Complex(1.0) - valence_mass), 21);
  std::cout << "# overlap operators set: D_m, tilde D_{m_P} (M5="<<M5<<")." << std::endl;

  ConservedCurrent<Fermion,Gauge> kop(Dm);   // K is mass-independent; multishift apply_k_ms via operator()
  MatPoly op_K; op_K.push_back(cplx(1.0), {&kop});

  // multishift (_ms) op-set for the single forward solves phi' / tilphi:
  //   X^{-1} b = op_XH (RHS X^dag b) + op_Xsq (CG).
  auto f_DmH  = std::bind(&Fermion::adj_deviceAsyncLaunch_ms,  &Dm, std::placeholders::_1, std::placeholders::_2);
  auto f_Dmsq = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms,  &Dm, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_DmH(f_DmH), M_Dmsq(f_Dmsq);
  MatPoly op_DmH;  op_DmH.push_back(cplx(1.0), {&M_DmH});
  MatPoly op_Dmsq; op_Dmsq.push_back(cplx(1.0), {&M_Dmsq});

  auto f_tilDmH  = std::bind(&Fermion::adj_deviceAsyncLaunch_ms,  &Dtil, std::placeholders::_1, std::placeholders::_2);
  auto f_tilDmsq = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms,  &Dtil, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_tilDmH(f_tilDmH), M_tilDmsq(f_tilDmsq);
  MatPoly op_tilDmH;  op_tilDmH.push_back(cplx(1.0), {&M_tilDmH});
  MatPoly op_tilDmsq; op_tilDmsq.push_back(cplx(1.0), {&M_tilDmsq});

  // ---- geometry weights ------------------------------------------------------------------------
  const int n_sites = static_cast<int>(base.n_sites);
  const int n_links = static_cast<int>(base.links.size());

  std::vector<double> w_tp(base.n_sites);    // temporal (Eq. 4.32): A_n / kappa_t(n)^2
  for(Idx n=0; n<base.n_sites; n++){ const double kt=DW.kappa_t[n]; w_tp[n]=base.dual_areas[n]/(kt*kt); }

  std::vector<double> w_sp(base.n_links);    // spatial (Eq. 4.29): link_volume[il] / kappa[il]^2
  for(Idx il=0; il<base.n_links; il++){ const double ks=DW.bd.kappa[il]; w_sp[il]=base.link_volume[il]/(ks*ks); }

  // ylm m-summed tower (Sec. 3.6): W_ell[l][n] = (A_n / kappa_t(n)) * sum_{m=-l}^{l} Y_lm(n^)  (kappa^1)
  constexpr int L_MAX = 2;
  constexpr int N_ELL = L_MAX + 1;
  std::vector<int> ls(N_ELL); for(int l=0;l<N_ELL;l++) ls[l]=l;
  std::vector<std::vector<double>> W_ell(N_ELL, std::vector<double>(n_sites, 0.0));
  for(int l=0; l<N_ELL; l++)
    for(int n=0; n<n_sites; n++){
      const double kt = DW.kappa_t[n];
      double s = 0.0;
      for(int m=-l; m<=l; m++) s += Ylm_real(l, m, base.sites[n]);
      W_ell[l][n] = base.dual_areas[n] * s / kt;
    }

  // ---- output: data_<ESNID>/disc_nhits<H>/disc.<config>.h5 -------------------------------------
  std::string ens_base = ens_dir;
  if(!ens_base.empty() && ens_base.back()=='/') ens_base.pop_back();
  { const auto slash = ens_base.find_last_of('/'); if(slash!=std::string::npos) ens_base = ens_base.substr(slash+1); }
  const std::string esnid = (free_field ? std::string("free") : ens_base)
                          + "_vmRe"+std::to_string(mass_re)+"vmIm"+std::to_string(mass_im);
  const std::string dir_out = "data_"+esnid+"/disc_nhits"+std::to_string(nhits)+"/";
  std::filesystem::create_directories(dir_out);
  std::cout << "# dir_out = " << dir_out
            << "  (n_sites="<<n_sites<<", n_links="<<n_links<<", n_ell="<<N_ELL<<")" << std::endl;

  // ---- host buffers ----------------------------------------------------------------------------
  //   eta = Z_2 source; phi = phi' = D_m^{-1} eta (forward leg); tilphi = tilde D_{m_P}^{-1} eta (parity);
  //   kphi = K(.,t) phi' (looped apply); tmp = preconditioned CG RHS.
  FermionVector eta, phi, tilphi, kphi, tmp;

  // helper: write a length-Nt complex vector (RAW -- no 1/4pi) under <key>/{real,imag}.
  auto write_vec = [&](HighFive::File& h5, const std::string& key, const std::vector<Complex>& v){
    std::vector<double> re(Nt), im(Nt);
    for(int t=0;t<Nt;t++){ re[t]=v[t].real(); im[t]=v[t].imag(); }
    h5.createDataSet(key+"/real", re);
    h5.createDataSet(key+"/imag", im);
  };

  // free field: single deterministic config (k=0), U=1.  ensemble: loop ckpoint_lat.k in ens_dir
  // with stride ninter (--ninter; default 10).
  const int k_ckpoint = free_field ? 1 : ninter;
  const int kmax      = free_field ? 0 : 1000;

  for(int k = 0; k <= kmax; k += k_ckpoint){
    if(!free_field){
      const std::string str_lat = ens_dir + "ckpoint_lat." + std::to_string(k);
      if(!std::filesystem::exists(str_lat)){ if(k==0) continue; else break; }
      U.read(str_lat);
    }
    // ONE file per config.  Resume: skip ONLY if the "complete" sentinel (written LAST) is present.
    const std::string h5path = dir_out + "disc." + std::to_string(k) + ".h5";
    if(std::filesystem::exists(h5path)){
      bool complete = false;
      try { HighFive::File h5c(h5path, HighFive::File::ReadOnly); complete = h5c.exist("complete"); }
      catch(...) {}
      if(complete){ std::cout<<"# skip k="<<k<<" (complete)"<<std::endl; continue; }
      std::cout<<"# k="<<k<<" exists but INCOMPLETE -> recompute"<<std::endl;
    }
    Dm.update(U);  Dtil.update(U);
    std::cout << "# k="<<k<<(free_field?" (free field)":"")
              << "  lambda_min/max="<<Dm.lambda_min<<"/"<<Dm.lambda_max<<std::endl;

    HighFive::File h5(h5path, HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
    h5.createDataSet("nhits", std::vector<int>{nhits});
    h5.createDataSet("ls",    ls);

    for(int h=0; h<nhits; h++){
      const auto t_hit0 = std::chrono::steady_clock::now();
      std::cout << "# k="<<k<<" hit "<<(h+1)<<"/"<<nhits<<std::endl;
      eta.fill_z2_source(rng);
      const std::string hp = "h" + std::to_string(h) + "/";   // key prefix /h{h}/

      // shared forward leg phi' = D_m^{-1} eta
      op_DmH.from_cpu<N>(tmp.field, eta.field);
      op_Dmsq.solve<N>(phi.field, tmp.field, Comp::TOL_OUTER);

      // single-time trace T(a,t) = eta^dag K(a,t) phi', folded with the projection weights in-program:
      //   J_tp(t)=sum_n w_tp[n] T(n,t),  J_ell(t)=sum_n W_ell[l][n] T(n,t)  (temporal pass);
      //   J_sp(t)=sum_l w_sp[il] T(l,t)  (spatial pass).  RAW (no 1/4pi).
      std::vector<Complex> Jtp(Nt, Complex(0,0)), Jsp(Nt, Complex(0,0));
      std::vector<std::vector<Complex>> Jyl(N_ELL, std::vector<Complex>(Nt, Complex(0,0)));
      for(int n=0; n<n_sites; n++)
        for(int t=0; t<Nt; t++){
          kop.set_temporal(U, t, (Idx)n, /*dag=*/false);
          op_K.from_cpu<N>(kphi.field, phi.field);            // kphi = K(n,t) phi'
          const Complex Tnt = eta.dag(kphi);                  // T(n,t) = eta^dag kphi
          Jtp[t] += w_tp[n] * Tnt;
          for(int l=0; l<N_ELL; l++) Jyl[l][t] += W_ell[l][n] * Tnt;
        }
      for(int a=0; a<n_links; a++){
        const BaseLink lk = base.links[a];
        const Idx il = base.map2il.at(lk);
        for(int t=0; t<Nt; t++){
          kop.set_spatial(U, t, lk, /*dag=*/false);
          op_K.from_cpu<N>(kphi.field, phi.field);            // kphi = K(l,t) phi'
          Jsp[t] += w_sp[il] * eta.dag(kphi);
        }
      }
      write_vec(h5, hp+"tp/J", Jtp);
      for(int l=0;l<N_ELL;l++) write_vec(h5, hp+"ylm/l"+std::to_string(l)+"/J", Jyl[l]);
      write_vec(h5, hp+"sp/J", Jsp);

      // parity: dagger-leg tilde single-trace \tilde T(a,t) = (K(a,t) tilphi)^dag eta, tilphi=tilde D^{-1} eta
      // (massless / m_F: \tilde T = T^* -> reconstructed downstream, no dump).
      if(parity){
        op_tilDmH.from_cpu<N>(tmp.field, eta.field);
        op_tilDmsq.solve<N>(tilphi.field, tmp.field, Comp::TOL_OUTER);   // tilphi = tilde D_{m_P}^{-1} eta
        std::vector<Complex> JtpT(Nt, Complex(0,0)), JspT(Nt, Complex(0,0));
        std::vector<std::vector<Complex>> JylT(N_ELL, std::vector<Complex>(Nt, Complex(0,0)));
        for(int n=0; n<n_sites; n++)
          for(int t=0; t<Nt; t++){
            kop.set_temporal(U, t, (Idx)n, /*dag=*/false);
            op_K.from_cpu<N>(kphi.field, tilphi.field);       // kphi = K(n,t) tilphi
            const Complex TtilNT = kphi.dag(eta);             // \tilde T(n,t) = (K tilphi)^dag eta
            JtpT[t] += w_tp[n] * TtilNT;
            for(int l=0; l<N_ELL; l++) JylT[l][t] += W_ell[l][n] * TtilNT;
          }
        for(int a=0; a<n_links; a++){
          const BaseLink lk = base.links[a];
          const Idx il = base.map2il.at(lk);
          for(int t=0; t<Nt; t++){
            kop.set_spatial(U, t, lk, /*dag=*/false);
            op_K.from_cpu<N>(kphi.field, tilphi.field);
            JspT[t] += w_sp[il] * kphi.dag(eta);
          }
        }
        write_vec(h5, hp+"tp/Jtil", JtpT);
        for(int l=0;l<N_ELL;l++) write_vec(h5, hp+"ylm/l"+std::to_string(l)+"/Jtil", JylT[l]);
        write_vec(h5, hp+"sp/Jtil", JspT);
      }

      const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now()-t_hit0).count();
      std::cout << "#   hit "<<(h+1)<<" done ["<<secs<<" s]  Jtp(0)=("<<Jtp[0].real()<<","<<Jtp[0].imag()<<")" << std::endl;
    } // hits

    h5.createDataSet("complete", std::vector<int>{1});   // sentinel: ALL datasets present (written LAST)
    std::cout << "# wrote " << h5path << std::endl;
  } // k

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();
  return 0;
}
