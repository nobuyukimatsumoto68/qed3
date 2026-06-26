// jj_local_ylm_disc_stoch_claude.cu
// Chunk B2: STOCHASTIC DISCONNECTED one-point loops J^a_{lm}(t) (per-m, vector a=1,2,3, l<=3) of the
// LOCAL current, PLUS the sigma_PS condensate -- both share one phi = D_m^{-1} eta solve.  Volume
// TIME+SPIN dilution sweep (model: saved_scripts_claude/disc_claude.cu).  Massless / m_F; m_P OUT OF SCOPE.
//
//   J^a_{lm}(t) = tr[ Sigma^a_{lm}(t) D_m^{-1} ],  Sigma^a_{lm}(t)=sum_n A_n Y_lm(n^) sigma_a(n,t).
// One-point estimator with the volume source eta and phi = D_m^{-1} eta (per dilution pattern):
//   J^a_{lm}(t) += sum_n A_n Y_lm(n^) conj(eta(t,n,spin)) (sigma_a phi)(t,n,spin)   (RAW; summed over patterns).
// The spin loop reconstructs the FULL spin trace (sigma_a mixes spin; each pattern reads one spin row),
// so spin dilution is MANDATORY here (as in disc_claude.cu).  Analysis forms
//   C_disc,l(dt) = (1/(2l+1)) sum_m two_point(J^a_{lm}),  physical = -C_conn + C_disc (Eq. 3.39).
// sigma_PS condensate (a=0 scalar, area-weighted, spacetime-summed): etadag_xi = -sum_d eta_A^dag phi
//   (eta_A = A_n eta); sigma_PS = etadag_xi + conj(etadag_xi) for massless/m_F.  NO extra solve.
//
// Disc is intrinsically noisy (one volume hit/config; config-average-dominated); = 0 in free.
// --disc-tblock (default 8) = timeslices per time-dilution class (interval = Nt/t_block classes).
// Plan: jj_local_ylm_impl_plan_claude.md (Chunk B2).  Per-hit atomic + "complete"-gated.
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

static int seed_from_string(const std::string& s){
  std::seed_seq seq(s.begin(), s.end());
  std::uint32_t w;
  seq.generate(&w, &w + 1);
  return static_cast<int>(w);
}

void PrintHelp(){
  printf("jj_local_ylm_disc_stoch: stochastic DISCONNECTED per-m Y_lm one-point loops + sigma_PS.\n");
  printf("  --gsq <x>            Wilson coupling squared (ensemble id; default 8.0)\n");
  printf("  --Nf <n>             number of fermion flavors (ensemble id; default 2)\n");
  printf("  --nu0 <x>            sea quark asymmetry (ensemble id; default 1.0)\n");
  printf("  --nu1 <x>            valence Wilson-Dirac asymmetry (operator; default nu0)\n");
  printf("  --mass-re <x>        valence mass Re (default 0.0; m_F is real)\n");
  printf("  --mass-im <y>        valence mass Im (default 0.0; m_P parity NOT supported here)\n");
  printf("  --ens-dir <path>     sea config dir; OMIT => free field (U=1)\n");
  printf("  --nhits <n>          stochastic hits (default 1)\n");
  printf("  --disc-tblock <tb>   timeslices per time-dilution class (default 8; interval=Nt/tb classes)\n");
  printf("  --stride <s>         ensemble config stride (default 10)\n");
  printf("  --kmin <a> --kmax <b> config range [a,b) (default 0..1e6)\n");
  printf("  -h, --help\n");
  exit(0);
}

void ParseArgs(int argc, char* argv[],
               double& gsq, int& Nf, double& nu0, double& nu1,
               double& mass_re, double& mass_im,
               std::string& ens_dir, int& nhits, int& disc_tblock, int& stride,
               int& kmin, int& kmax){
  static struct option long_opts[] = {
    {"gsq",     required_argument, nullptr, 'g'},
    {"Nf",      required_argument, nullptr, 'N'},
    {"nu0",     required_argument, nullptr, 'n'},
    {"nu1",     required_argument, nullptr, 'm'},
    {"mass-re", required_argument, nullptr, 'r'},
    {"mass-im", required_argument, nullptr, 'i'},
    {"ens-dir", required_argument, nullptr, 'e'},
    {"nhits",   required_argument, nullptr, 'H'},
    {"disc-tblock", required_argument, nullptr, 'T'},
    {"stride",  required_argument, nullptr, 'I'},
    {"kmin",    required_argument, nullptr, 'a'},
    {"kmax",    required_argument, nullptr, 'b'},
    {"help",    no_argument,       nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
  };
  int opt, idx;
  while((opt = getopt_long(argc, argv, "g:N:n:m:r:i:e:H:T:I:a:b:h", long_opts, &idx)) != -1){
    switch(opt){
    case 'g': gsq     = std::stod(optarg); break;
    case 'N': Nf      = std::stoi(optarg); break;
    case 'n': nu0     = std::stod(optarg); break;
    case 'm': nu1     = std::stod(optarg); break;
    case 'r': mass_re = std::stod(optarg); break;
    case 'i': mass_im = std::stod(optarg); break;
    case 'e': ens_dir = optarg; break;
    case 'H': nhits   = std::stoi(optarg); break;
    case 'T': disc_tblock = std::stoi(optarg); break;
    case 'I': stride  = std::stoi(optarg); break;
    case 'a': kmin    = std::stoi(optarg); break;
    case 'b': kmax    = std::stoi(optarg); break;
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
  int nhits=1;
  int disc_tblock=8;
  int stride=10;
  int kmin=0;
  int kmax=1000000;

  ParseArgs(argc, argv, gsq, Nf, nu0, nu1, mass_re, mass_im, ens_dir, nhits, disc_tblock, stride, kmin, kmax);
  if(nu1 < 0.0) nu1 = nu0;

  const Complex valence_mass(mass_re, mass_im);
  const bool free_field = ens_dir.empty();
  const bool parity = (std::abs(mass_im) > 1.0e-15) && (std::abs(mass_re) <= 1.0e-15);
  assert(!parity && "m_P needs the tilde-D backward leg -- OUT OF SCOPE for B2 (vector disc + sigma_PS, massless/m_F)");
  assert(disc_tblock >= 1 && Comp::Nt % disc_tblock == 0 && "Nt % disc_tblock must be 0");

  std::cout << "# gsq="<<gsq<<" Nf="<<Nf<<" nu0="<<nu0<<" nu1="<<nu1
            << " mass=("<<mass_re<<","<<mass_im<<")"
            << " ens_dir="<<(free_field?std::string("<free-field U=1>"):ens_dir)
            << " nhits="<<nhits<<" disc_tblock="<<disc_tblock << std::endl;

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
  const double at = 0.2;
  if(Nt!=1) assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu1);
  std::cout << "# DW set." << std::endl;

  Gauge U(base);
  Rng rng(base, 1234);

  Fermion Dm(DW, valence_mass, 11);
  std::cout << "# overlap operator set: D_m (M5="<<M5<<")." << std::endl;

  // phi = D_m^{-1} eta = op_Dmsq^{-1}(op_DmH eta).  (No D_m^{-dag} needed for disc/condensate.)
  auto f_DmH  = std::bind(&Fermion::adj_deviceAsyncLaunch_ms,  &Dm, std::placeholders::_1, std::placeholders::_2);
  auto f_Dmsq = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms,  &Dm, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_DmH(f_DmH), M_Dmsq(f_Dmsq);
  MatPoly op_DmH;  op_DmH.push_back(cplx(1.0), {&M_DmH});
  MatPoly op_Dmsq; op_Dmsq.push_back(cplx(1.0), {&M_Dmsq});

  constexpr int L_MAX_YLM = 3;
  const int n_sites = static_cast<int>(base.n_sites);
  const int interval = Nt / disc_tblock;   // number of time-dilution classes

  // LOCAL site current area weight A_n = dual_areas (bare sigma, no kappa).
  std::vector<double> w_site(base.n_sites);
  for(Idx n=0; n<base.n_sites; n++) w_site[n]=base.dual_areas[n];

  // ---- output: data_<ESNID>/corr_ylm_disc_tb<tb>_nhits<H>/corr.<k>.h<h>.h5
  std::string ens_base = ens_dir;
  if(!ens_base.empty() && ens_base.back()=='/') ens_base.pop_back();
  { const auto slash = ens_base.find_last_of('/'); if(slash!=std::string::npos) ens_base = ens_base.substr(slash+1); }
  const std::string esnid = (free_field ? std::string("free") : ens_base)
                          + "_vmRe"+std::to_string(mass_re)+"vmIm"+std::to_string(mass_im);
  const std::string dir_out = "data_"+esnid+"/corr_ylm_disc_tb"+std::to_string(disc_tblock)
                            + "_nhits"+std::to_string(nhits)+"/";
  std::filesystem::create_directories(dir_out);
  std::cout << "# dir_out = " << dir_out << std::endl;

  // disc J(t) and condensate are RAW (no inv4pi fold), matching disc_claude.cu / the dilute disc.
  auto write_vec = [&](HighFive::File& h5, const std::string& key, const std::vector<Complex>& C){
    std::vector<double> re(C.size()), im(C.size());
    for(size_t t=0;t<C.size();t++){ re[t]=C[t].real(); im[t]=C[t].imag(); }
    h5.createDataSet(key+"/real", re);  h5.createDataSet(key+"/imag", im);
  };

  FermionVector eta, eta_A, tmp, phi, Gamma_phi;

  const int k_ckpoint = free_field ? 1 : stride;
  const int k_lo      = free_field ? 0 : kmin;
  const int k_hi      = free_field ? 1 : kmax;

  for(int k = k_lo; k < k_hi; k += k_ckpoint){
    if(!free_field){
      const std::string str_lat = ens_dir + "ckpoint_lat." + std::to_string(k);
      if(!std::filesystem::exists(str_lat)){ if(k==0) continue; else break; }
      U.read(str_lat);
    }
    Dm.update(U);
    std::cout << "# k="<<k<<(free_field?" (free field)":"")
              << "  lambda_min/max="<<Dm.lambda_min<<"/"<<Dm.lambda_max<<std::endl;

    for(int h=0; h<nhits; h++){
      const std::string h5path_h = dir_out + "corr." + std::to_string(k) + ".h" + std::to_string(h) + ".h5";
      if(std::filesystem::exists(h5path_h)){
        bool c=false; try { HighFive::File f(h5path_h,HighFive::File::ReadOnly); c=f.exist("complete"); } catch(...) {}
        if(c){ std::cout<<"# skip k="<<k<<" hit "<<h<<" (complete)"<<std::endl; continue; }
      }
      const std::string seed_str = esnid + "_k" + std::to_string(k) + "_h" + std::to_string(h);
      rng.reseed(seed_from_string(seed_str));
      const auto t_hit0 = std::chrono::steady_clock::now();
      std::cout << "# k="<<k<<" hit "<<(h+1)<<"/"<<nhits<<"  (disc_tblock="<<disc_tblock
                <<", "<<interval<<" time x "<<NS<<" spin classes, seed='"<<seed_str<<"')" << std::endl;

      // accumulators J[a-1][l][m+l][t] (per-m one-point), and the condensate scalar (SUMMED over patterns).
      std::vector<std::vector<std::vector<std::vector<Complex>>>> J(3);
      for(int a=0;a<3;a++){ J[a].resize(L_MAX_YLM+1);
        for(int l=0;l<=L_MAX_YLM;l++) J[a][l].assign(2*l+1, std::vector<Complex>(Nt, Complex(0,0))); }
      Complex acc_etadag_xi(0,0);

      // ===== TIME+SPIN dilution sweep (spin loop MANDATORY: reconstructs the sigma_a spin trace) =====
      for(int t_s=0; t_s<interval; t_s++){
        for(int spin=0; spin<NS; spin++){
          eta.time_spin_dilution(rng, t_s, disc_tblock, spin);   // volume Z2, this (t-class, spin)
          // phi = D_m^{-1} eta
          op_DmH.from_cpu<N>(tmp.field, eta.field);
          op_Dmsq.solve<N>(phi.field, tmp.field, Comp::TOL_OUTER);

          // CONDENSATE sigma_PS: eta_A = A_n eta; acc += eta_A^dag phi  -> tr[A D_m^{-1}] (both spins via patterns)
          for(int tt=0; tt<Nt; tt++)
            for(int n=0; n<n_sites; n++){
              const double aw = w_site[n];
              for(int i=0;i<NS;i++) eta_A(tt,(Idx)n,i) = aw*eta(tt,(Idx)n,i);
            }
          acc_etadag_xi += eta_A.dag(phi);

          // DISC currents J^a_{lm}(t): per (a,l,m), W = A_n Y_lm sigma_a phi; J += sum_n conj(eta) W.
          for(int a=1; a<=3; a++){
            for(int l=0; l<=L_MAX_YLM; l++){
              for(int m=-l; m<=l; m++){
                Gamma_phi = phi;
                Gamma_phi.mult_sigma(a);
                Gamma_phi.mult_Ylm_real(l, m, base);   // folds A_n Y_lm
                eta.accumulate_loop_raw(J[a-1][l][m+l], Gamma_phi, t_s, disc_tblock, spin);
              }
            }
          }
        }
      }

      // sigma_PS = etadag_xi + conj(etadag_xi); etadag_xi = -sum_d eta_A^dag phi (loop sign), matches the dilute.
      const Complex etadag_xi = -acc_etadag_xi;

      // ---- write ----
      const std::string hp = "h0/";
      const std::string h5tmp = h5path_h + ".tmp";
      auto h5p = std::make_unique<HighFive::File>(h5tmp,
                   HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
      HighFive::File& h5 = *h5p;
      h5.createDataSet("nhits", std::vector<int>{nhits});
      h5.createDataSet("hit",   std::vector<int>{h});
      h5.createDataSet("rng_seed", seed_str);
      h5.createDataSet("disc_tblock", std::vector<int>{disc_tblock});
      h5.createDataSet("L_MAX_YLM",   std::vector<int>{L_MAX_YLM});
      for(int a=1; a<=3; a++){
        const std::string chan="s"+std::to_string(a);
        for(int l=0; l<=L_MAX_YLM; l++){
          for(int m=-l; m<=l; m++){
            const std::string kp=hp+"disc/ylm/"+chan+"/l"+std::to_string(l)+"/m"+std::to_string(m)+"/J";
            write_vec(h5, kp, J[a-1][l][m+l]);
          }
        }
      }
      write_vec(h5, hp+"condensate/etadag_xi", std::vector<Complex>{ etadag_xi });
      h5.createDataSet("complete", std::vector<int>{1});
      h5p.reset();
      std::filesystem::rename(h5tmp, h5path_h);
      const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now()-t_hit0).count();
      std::cout << "#   hit "<<(h+1)<<" done (disc ylm s1/s2/s3 l<=3 per-m + sigma_PS) ["<<secs<<" s] -> "<<h5path_h << std::endl;
    } // hits
  } // k

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();
  return 0;
}
