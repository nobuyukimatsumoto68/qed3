// jj_local_ylm_conn_stoch_claude.cu
// Chunk B1: STOCHASTIC CONNECTED Y_lm tower of the LOCAL current correlator, per-m, VECTOR + AXIAL.
// One-end trick with a SINGLE-TIMESLICE wall source at one origin t0 (variance-optimal for the
// connected piece); NO conserved-K apply (the local sigma_a sink is element-wise).  Vector only,
// massless / m_F; m_P (parity) OUT OF SCOPE.  Per-m output (analysis does (1/(2l+1)) sum_m).
//
// Estimator (per Pauli a in {1,2,3}, real Y_lm, A_n = dual_areas; Sigma^a_{lm}(t)=sum_n A_n Y_lm(n^) sigma_a(n,t)):
//   phi            = D_m^{-1} eta                         (eta = Z2 wall at t0)
//   psi^a_{lm}     = D_m^{-dag} Sigma^a_{lm}(t0) eta
//   g^a_{lm}(t-t0) = psi^a_{lm}^dag [ Sigma^a_{lm}(t) phi ]   (estimates tr[Sigma_0 D_m^{-1} Sigma_t D_m^{-1}])
// tp tower = s3 (G_t), sp tower = (s1+s2)/2 (G_s); analysis forms g^a_l = (1/(2l+1)) sum_m g^a_{lm}.
// D_m^{-1} b   = op_Dmsq^{-1}(op_DmH b)  (apply D_m^dag, then CG on D_m^dag D_m).
// D_m^{-dag} b = op_Dm(op_Dmsq^{-1} b)   (CG on D_m^dag D_m, then apply D_m; D_m normal so DDH=DHD).
//
// Spin dilution (--spin-dilution): two single-spin wall classes summed (unbiased).  SINGLE origin
// (--t0, default 0); NO multi-t0, NO master-field superposition.  Per-hit atomic + "complete"-gated.
// Deterministic ground truth: jj_local_deter_claude.cu (ylm/s{a}/l{l}/m{m}).  Plan:
// jj_local_ylm_impl_plan_claude.md (Chunk B1).  Wall/Y_lm/sigma mechanics: meson_pq_wall_v2_claude.cu.
// Solve/output scaffolding: jj_corr_dilute_claude.cu.
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

void PrintHelp(){
  printf("jj_local_ylm_conn_stoch: stochastic CONNECTED per-m Y_lm tower of the LOCAL vector current.\n");
  printf("  --gsq <x>            Wilson coupling squared (ensemble id; default 8.0)\n");
  printf("  --Nf <n>             number of fermion flavors (ensemble id; default 2)\n");
  printf("  --nu0 <x>            sea quark asymmetry (ensemble id; default 1.0)\n");
  printf("  --nu1 <x>            valence Wilson-Dirac asymmetry (operator; default nu0)\n");
  printf("  --mass-re <x>        valence mass Re (default 0.0; m_F is real)\n");
  printf("  --mass-im <y>        valence mass Im (default 0.0; m_P parity NOT supported here)\n");
  printf("  --ens-dir <path>     sea config dir; OMIT => free field (U=1)\n");
  printf("  --nhits <n>          stochastic hits (default 1)\n");
  printf("  --t0 <t>             SINGLE source-time origin (default 0)\n");
  printf("  --stride <s>         ensemble config stride (default 10)\n");
  printf("  --kmin <a> --kmax <b> config range [a,b) (default 0..1e6)\n");
  printf("  --spin-dilution      two single-spin wall classes (default OFF = both spins, 1 class)\n");
  printf("  -h, --help\n");
  exit(0);
}

void ParseArgs(int argc, char* argv[],
               double& gsq, int& Nf, double& nu0, double& nu1,
               double& mass_re, double& mass_im,
               std::string& ens_dir, int& nhits, int& t0, int& stride,
               int& kmin, int& kmax, bool& spin_dilution){
  static struct option long_opts[] = {
    {"gsq",     required_argument, nullptr, 'g'},
    {"Nf",      required_argument, nullptr, 'N'},
    {"nu0",     required_argument, nullptr, 'n'},
    {"nu1",     required_argument, nullptr, 'm'},
    {"mass-re", required_argument, nullptr, 'r'},
    {"mass-im", required_argument, nullptr, 'i'},
    {"ens-dir", required_argument, nullptr, 'e'},
    {"nhits",   required_argument, nullptr, 'H'},
    {"t0",      required_argument, nullptr, 'T'},
    {"stride",  required_argument, nullptr, 'I'},
    {"kmin",    required_argument, nullptr, 'a'},
    {"kmax",    required_argument, nullptr, 'b'},
    {"spin-dilution", no_argument, nullptr, 's'},
    {"help",    no_argument,       nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
  };
  int opt, idx;
  while((opt = getopt_long(argc, argv, "g:N:n:m:r:i:e:H:T:I:a:b:sh", long_opts, &idx)) != -1){
    switch(opt){
    case 'g': gsq     = std::stod(optarg); break;
    case 'N': Nf      = std::stoi(optarg); break;
    case 'n': nu0     = std::stod(optarg); break;
    case 'm': nu1     = std::stod(optarg); break;
    case 'r': mass_re = std::stod(optarg); break;
    case 'i': mass_im = std::stod(optarg); break;
    case 'e': ens_dir = optarg; break;
    case 'H': nhits   = std::stoi(optarg); break;
    case 'T': t0      = std::stoi(optarg); break;
    case 'I': stride  = std::stoi(optarg); break;
    case 'a': kmin    = std::stoi(optarg); break;
    case 'b': kmax    = std::stoi(optarg); break;
    case 's': spin_dilution = true; break;
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
  int t0=0;            // SINGLE source-time origin
  int stride=10;
  int kmin=0;
  int kmax=1000000;
  bool spin_dilution=false;

  ParseArgs(argc, argv, gsq, Nf, nu0, nu1, mass_re, mass_im, ens_dir, nhits, t0, stride, kmin, kmax,
            spin_dilution);
  if(nu1 < 0.0) nu1 = nu0;

  const Complex valence_mass(mass_re, mass_im);
  const bool free_field = ens_dir.empty();
  const bool parity = (std::abs(mass_im) > 1.0e-15) && (std::abs(mass_re) <= 1.0e-15);
  assert(!parity && "m_P (purely imaginary mass) needs the tilde-D leg -- OUT OF SCOPE for B1 (vector, massless/m_F)");
  assert(t0>=0 && t0<Comp::Nt && "t0 out of range");

  std::cout << "# gsq="<<gsq<<" Nf="<<Nf<<" nu0="<<nu0<<" nu1="<<nu1
            << " mass=("<<mass_re<<","<<mass_im<<")"
            << " ens_dir="<<(free_field?std::string("<free-field U=1>"):ens_dir)
            << " nhits="<<nhits<<" t0="<<t0<<" spin_dilution="<<(spin_dilution?1:0) << std::endl;

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

  // D_m = D_ov + m (Eq. 3.59).  Only operator needed: vector local current, massless/m_F.
  Fermion Dm(DW, valence_mass, 11);
  std::cout << "# overlap operator set: D_m (M5="<<M5<<")." << std::endl;

  // X^{-1} via op_Xsq (CG on D_m^dag D_m = D_m D_m^dag, D_m normal) + op_DmH/op_Dm (RHS):
  auto f_Dm   = std::bind(&Fermion::mult_deviceAsyncLaunch_ms, &Dm, std::placeholders::_1, std::placeholders::_2);
  auto f_DmH  = std::bind(&Fermion::adj_deviceAsyncLaunch_ms,  &Dm, std::placeholders::_1, std::placeholders::_2);
  auto f_Dmsq = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms,  &Dm, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_Dm(f_Dm), M_DmH(f_DmH), M_Dmsq(f_Dmsq);
  MatPoly op_Dm;   op_Dm.push_back(cplx(1.0), {&M_Dm});
  MatPoly op_DmH;  op_DmH.push_back(cplx(1.0), {&M_DmH});
  MatPoly op_Dmsq; op_Dmsq.push_back(cplx(1.0), {&M_Dmsq});

  const double inv4pi = 1.0/(4.0*std::acos(-1.0));
  constexpr int L_MAX_YLM = 3;
  const int n_spin = spin_dilution ? NS : 1;

  // ---- output: data_<ESNID>/corr_ylm_conn_t0<t0>_nhits<H>_s<0/1>/corr.<k>.h<h>.h5
  std::string ens_base = ens_dir;
  if(!ens_base.empty() && ens_base.back()=='/') ens_base.pop_back();
  { const auto slash = ens_base.find_last_of('/'); if(slash!=std::string::npos) ens_base = ens_base.substr(slash+1); }
  const std::string esnid = (free_field ? std::string("free") : ens_base)
                          + "_vmRe"+std::to_string(mass_re)+"vmIm"+std::to_string(mass_im);
  const std::string dir_out = "data_"+esnid+"/corr_ylm_conn_t0"+std::to_string(t0)
                            + "_nhits"+std::to_string(nhits)+"_s"+std::to_string(spin_dilution?1:0)+"/";
  std::filesystem::create_directories(dir_out);
  std::cout << "# dir_out = " << dir_out << std::endl;

  // write per-m correlator with the inv4pi fold (matches jj_local_deter_claude.cu::write_corr).
  auto write_corr = [&](HighFive::File& h5, const std::string& key, const std::vector<Complex>& C, bool conj){
    std::vector<double> re(Nt), im(Nt);
    for(int t=0;t<Nt;t++){ const Complex g=inv4pi*C[t]; re[t]=g.real(); im[t]=conj?-g.imag():g.imag(); }
    h5.createDataSet(key+"/real", re);  h5.createDataSet(key+"/imag", im);
  };

  // work vectors  (psi = vector source leg D_m^{-dag} src; psiA = axial source leg D_m^{-1} src)
  FermionVector eta, tmp, phi, src, psi, psiA, Phi;

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
      std::cout << "# k="<<k<<" hit "<<(h+1)<<"/"<<nhits<<"  (t0="<<t0<<", seed='"<<seed_str<<"')" << std::endl;

      // accumulators g (vector) and gA (axial), [a-1][l][m+l][dt], SUMMED over the n_spin dilution patterns.
      std::vector<std::vector<std::vector<std::vector<Complex>>>> g(3), gA(3);
      for(int a=0;a<3;a++){ g[a].resize(L_MAX_YLM+1); gA[a].resize(L_MAX_YLM+1);
        for(int l=0;l<=L_MAX_YLM;l++){ g[a][l].assign(2*l+1, std::vector<Complex>(Nt, Complex(0,0)));
                                        gA[a][l].assign(2*l+1, std::vector<Complex>(Nt, Complex(0,0))); } }

      for(int sc=0; sc<n_spin; sc++){
        if(spin_dilution) eta.fill_z2_wall_spin(rng, t0, sc);   // single-spin wall at t0
        else              eta.fill_z2_wall_source(rng, t0);     // both-spin wall at t0
        // phi = D_m^{-1} eta
        op_DmH.from_cpu<N>(tmp.field, eta.field);
        op_Dmsq.solve<N>(phi.field, tmp.field, Comp::TOL_OUTER);

        for(int a=1; a<=3; a++){
          for(int l=0; l<=L_MAX_YLM; l++){
            for(int m=-l; m<=l; m++){
              // src = Sigma^a_{lm}(t0) eta  (shared by both source legs; mult_Ylm_real folds A_n)
              src = eta;
              src.mult_sigma(a);
              src.mult_Ylm_real(l, m, base);
              // VECTOR source: psi   = D_m^{-dag} src  (op_Dmsq solve, then op_Dm)
              op_Dmsq.solve<N>(tmp.field, src.field, Comp::TOL_OUTER);
              op_Dm.from_cpu<N>(psi.field, tmp.field);
              // AXIAL  source: psi_A = D_m^{-1}  src  (op_DmH, then op_Dmsq solve) -- forward leg, NO (1-D_ov)
              op_DmH.from_cpu<N>(tmp.field, src.field);
              op_Dmsq.solve<N>(psiA.field, tmp.field, Comp::TOL_OUTER);
              // shared sink leg: Phi = Sigma^a_{lm}(t) phi  (all t, element-wise)
              Phi = phi;
              Phi.mult_sigma(a);
              Phi.mult_Ylm_real(l, m, base);
              // g[a,l,m](t-t0) += psi^dag [Sigma_t phi]; gA same with the axial source leg
              auto& gv  = g [a-1][l][m+l];
              auto& gva = gA[a-1][l][m+l];
              for(int t=0;t<Nt;t++){
                const int dt=((t-t0)%Nt+Nt)%Nt;
                const VC Phit = Phi.slice(t);                 // VC::dot = conj(.).Phit over slice t
                gv [dt] += psi .slice(t).dot(Phit);
                gva[dt] += psiA.slice(t).dot(Phit);
              }
            }
          }
        }
      }

      // ---- write per-m (Vpp + Vmm=conj for non-parity), under h0/ylm/s{a}/l{l}/m{m}/ ----
      const std::string hp = "h0/";
      const std::string h5tmp = h5path_h + ".tmp";
      auto h5p = std::make_unique<HighFive::File>(h5tmp,
                   HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
      HighFive::File& h5 = *h5p;
      h5.createDataSet("t0",    std::vector<int>{t0});
      h5.createDataSet("nhits", std::vector<int>{nhits});
      h5.createDataSet("hit",   std::vector<int>{h});
      h5.createDataSet("rng_seed", seed_str);
      h5.createDataSet("spin_dilution", std::vector<int>{spin_dilution?1:0});
      h5.createDataSet("L_MAX_YLM", std::vector<int>{L_MAX_YLM});
      for(int a=1; a<=3; a++){
        const std::string chan="s"+std::to_string(a);
        for(int l=0; l<=L_MAX_YLM; l++){
          for(int m=-l; m<=l; m++){
            const std::string kp =hp+"ylm/"      +chan+"/l"+std::to_string(l)+"/m"+std::to_string(m)+"/";
            const std::string kpA=hp+"ylm_axial/"+chan+"/l"+std::to_string(l)+"/m"+std::to_string(m)+"/";
            write_corr(h5, kp +"Vpp", g [a-1][l][m+l], false);
            write_corr(h5, kp +"Vmm", g [a-1][l][m+l], true);
            write_corr(h5, kpA+"Vpp", gA[a-1][l][m+l], false);
            write_corr(h5, kpA+"Vmm", gA[a-1][l][m+l], true);
          }
        }
      }
      h5.createDataSet("complete", std::vector<int>{1});
      h5p.reset();
      std::filesystem::rename(h5tmp, h5path_h);
      const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now()-t_hit0).count();
      std::cout << "#   hit "<<(h+1)<<" done (conn ylm vector+axial, s1/s2/s3, l<=3, per-m) ["<<secs<<" s] -> "<<h5path_h << std::endl;
    } // hits
  } // k

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();
  return 0;
}
