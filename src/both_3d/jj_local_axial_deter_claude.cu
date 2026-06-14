// jj_local_axial_deter_claude.cu
// -----------------------------------------------------------------------------
// LOCAL (ultralocal) AXIAL current correlators with the EXACT all-to-all propagator.  AXIAL sibling of
// jj_local_deter_claude.cu (vector); same bare on-site Pauli sigma^a (a=1,2,3), same SITE weights, same
// channels (s1,s2,s3) + ylm tower -- only the contraction differs.
//
// The axial reduction "replace K -> sigma, rip off (1-D_ov)" of Eq.(3.55) gives, per site n:
//   C^loc_{A+-}(t0,t) = sum_n w_site[n] tr[ sigma^a(n,t0) P^dag sigma^a(n,t) P ]
// vs the vector  C^loc_{V++}(t0,t) = sum_n w_site[n] tr[ sigma^a(n,t0) P sigma^a(n,t) P ].
// The ONLY change is the t0-leg propagator P -> P^dag (= D_ov^{-dag}); the sink kernel keeps sigma^a
// (Hermitian, sigma^{a dag} = sigma^a, so no sink change for the local current).  P^dag_{ab}=conj(P_{ba}),
// a conjugate-transpose lookup of the loaded dense P -- no extra solve.
//   s3 = <sigma_3 sigma_3> (temporal, G_t); s1,s2 = spatial (G_s = s1 + s2).  Ref: qed3int_v2-13.pdf
//   Sec. 3.3.2; plan jj_axial_trio_deter_impl_plan_claude.md.
//
// C_{A+-} is NOT self-reflection-even (reflects to C_{A-+}, Eq. 3.57) -> SINGLE complex channel "Apm"
// (no Vmm).
//
// L COMPILE-TIME (-DN_REFINE_CLI).  CLI: --ens-dir(omit=free) --mass-re/--mass-im (P + esnid) --n-t0
//   --ninter --gpu --ins(single) --prop-file --out-tag.
//   Output: data_<ESNID>/corr_deter_local_axial[1][_<tag>]_L<L>/corr.<k>.h5, keys h0/t0_b/{s1,s2,s3}/Apm
//   + h0/t0_b/ylm/l{l}/Apm + h0/disc/{s1,s2,s3}/J.
// -----------------------------------------------------------------------------

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
#include <getopt.h>
#include <Eigen/Dense>

using Double = double;
using Idx = std::int32_t;
using Complex = std::complex<double>;
using Link = std::array<Idx,2>;
using Face = std::vector<Idx>;
using MS=Eigen::Matrix2cd; using VD=Eigen::Vector2d; using VE=Eigen::Vector3d; using VC=Eigen::VectorXcd;
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
// #include "sparse_dirac.h"
#include "sparse_dirac_claude.h"   // O(len) bucketing CSR build (was O(N*len)); -DCSR_VERIFY to check
#include "matpoly_claude.h"
#include "dirac_pf.h"
#include "overlap_wmass_claude.h"
#include "conserved_current_claude.h"

// one sparse matrix entry of a current insertion W: contributes v to W[i,j].
struct Ent{ Idx i,j; Complex v; };

// LOCAL (site) current insertion: bare Pauli sigma^a (a=1,2,3) at dual site ix, timeslice t.  (Identical to
// the vector code: the axial difference is in the contraction, not the insertion -- sigma is Hermitian.)
template<typename WilsonDirac>
static void local_W_sigma(std::vector<Ent>& en, const WilsonDirac& DW, int a, int t, Idx ix){
  en.clear();
  const Idx off = (Idx)Comp::Nx*t + NS*ix;
  const MS s = DW.sigma[a];
  en.push_back({off,   off,   s(0,0)});
  en.push_back({off,   off+1, s(0,1)});
  en.push_back({off+1, off,   s(1,0)});
  en.push_back({off+1, off+1, s(1,1)});
}

// Ylm-weighted TEMPORAL current Sigma_{l,m}(t) = sum_n A_n Y_lm(n^) sigma_3(n,t) (diagonal real entries).
template<typename Base>
static void build_Sigma_ylm(std::vector<Ent>& en, const Base& base, const std::vector<double>& w_site,
                            int l, int m, int t){
  en.clear();
  const Idx Nx=Comp::Nx;
  for(int n=0; n<(int)base.n_sites; n++){
    const double wy = w_site[n] * Ylm_real(l, m, base.sites[n]);
    const Idx off = Nx*t + NS*(Idx)n;
    en.push_back({off,   off,   Complex(+wy, 0.0)});
    en.push_back({off+1, off+1, Complex(-wy, 0.0)});
  }
}

// AXIAL insertion-DIAGONAL trace for ONE insertion:
//   tr[ W(t0) P0^dag W(t) P ] = sum_{(i,j,v0)} sum_{(k,l,vt)} v0 P0^dag_{jk} vt P_{li},
//   P0^dag_{jk} = conj(P0_{kj}).  The FIRST (t0-leg = backward propagator eta xi^dag) is P0 daggered; the
//   sink-leg (forward xi eta^dag) is P.  massless/m_F: P0=P=D_m^{-1}.  m_P (Eqs. 3.60/3.61): P0=Dtil_inv
//   -> P0^dag = tilde D_{m_P}^{-dag}, P=D_m^{-1}, and the caller multiplies by (1+m_P)^{-1} (bare current).
static Complex tr_WPWP_axial(const std::vector<Ent>& E0, const std::vector<Ent>& Et,
                             const std::vector<Complex>& P0, const std::vector<Complex>& P){
  const Idx N=Comp::N;
  Complex s(0,0);
  for(const auto& e0 : E0)
    for(const auto& et : Et)
      s += e0.v * std::conj(P0[(size_t)et.i*N + e0.j]) * et.v * P[(size_t)et.j*N + e0.i];
  return s;
}

static void load_mat(HighFive::File& f, const std::string& key, std::vector<Complex>& M){
  // CHUNKED lean load: read real/imag in hyperslab blocks (~0.5 GB) straight into the complex M, so the
  // transient double buffer is ONE chunk, not the whole N^2 (13.75 GB at L=4).  Peak RAM = M (N^2 complex)
  // + one chunk ~ 28 GB at L=4 (vs ~41 GB for a full real/imag read).  Byte-identical to the full read.
  HighFive::DataSet dre = f.getDataSet(key+"/real");
  HighFive::DataSet dim = f.getDataSet(key+"/imag");
  size_t n = 1; for(auto d : dre.getDimensions()) n *= d;   // flat element count (1D N^2 here)
  M.resize(n);
  const size_t CH = (size_t)1 << 26;                         // 64 Mi doubles per block = 0.5 GB
  std::vector<double> buf;
  for(size_t off=0; off<n; off+=CH){
    const size_t cnt = std::min(CH, n-off);
    dre.select(std::vector<size_t>{off}, std::vector<size_t>{cnt}).read(buf);
    for(size_t i=0;i<cnt;i++) M[off+i].real(buf[i]);
    dim.select(std::vector<size_t>{off}, std::vector<size_t>{cnt}).read(buf);
    for(size_t i=0;i<cnt;i++) M[off+i].imag(buf[i]);
  }
}

// AXIAL single complex channel "Apm" (no Vmm); 1/(4pi) folded.
static void write_corr_axial(HighFive::File& h5, const std::string& key, const std::vector<Complex>& C){
  const int Nt=Comp::Nt;
  const double inv4pi=1.0/(4.0*M_PI);
  std::vector<double> re(Nt),im(Nt);
  for(int t=0;t<Nt;t++){ Complex g=inv4pi*C[t]; re[t]=g.real(); im[t]=g.imag(); }
  h5.createDataSet(key+"/real",re);
  h5.createDataSet(key+"/imag",im);
}

static void write_vec(HighFive::File& h5, const std::string& key, const std::vector<Complex>& C){
  const int Nt=Comp::Nt;
  std::vector<double> re(Nt),im(Nt);
  for(int t=0;t<Nt;t++){ re[t]=C[t].real(); im[t]=C[t].imag(); }
  h5.createDataSet(key+"/real",re);
  h5.createDataSet(key+"/imag",im);
}

struct Args{ double nu0=1.0,nu1=-1.0,mass_re=0.0,mass_im=0.0; std::string ens_dir,prop_file,out_tag; int n_t0=2,ninter=10,gpu=0,ins=-1; };
void PrintHelp(){ printf("jj_local_axial_deter: AXIAL C_{A+-} (bare sigma, P^dag t0-leg).  --mass-re --mass-im\n"
                         "  --ens-dir(omit=free) --n-t0 --ninter --nu0 --nu1 --gpu\n"
                         "  --prop-file <path>  read P from this exact file (e.g. cont_prop_L<L>/Dinv.0.h5)\n"
                         "  --out-tag <tag>     corr_deter_local_axial_<tag>_L<L>\n"
                         "  --ins <i>           SINGLE-insertion at site i (no sum/weight; -> corr_deter_local_axial1[_<tag>]_L<L>;\n"
                         "                      Ylm tower skipped).  default -1 = full sum over sites.\n"); }
void ParseArgs(int argc,char**argv,Args&a){
  static struct option lo[]={{"nu0",required_argument,0,'n'},{"nu1",required_argument,0,'m'},
    {"mass-re",required_argument,0,'r'},{"mass-im",required_argument,0,'i'},{"ens-dir",required_argument,0,'e'},
    {"n-t0",required_argument,0,'T'},{"ninter",required_argument,0,'I'},{"gpu",required_argument,0,'G'},
    {"prop-file",required_argument,0,'P'},{"out-tag",required_argument,0,'O'},{"ins",required_argument,0,'A'},
    {"help",no_argument,0,'h'},{0,0,0,0}};
  int opt,idx;
  while((opt=getopt_long(argc,argv,"n:m:r:i:e:T:I:G:P:O:A:h",lo,&idx))!=-1){ switch(opt){
    case 'n':a.nu0=std::stod(optarg);break; case 'm':a.nu1=std::stod(optarg);break;
    case 'r':a.mass_re=std::stod(optarg);break; case 'i':a.mass_im=std::stod(optarg);break;
    case 'e':a.ens_dir=optarg;break; case 'T':a.n_t0=std::stoi(optarg);break;
    case 'I':a.ninter=std::stoi(optarg);break; case 'G':a.gpu=std::stoi(optarg);break;
    case 'P':a.prop_file=optarg;break; case 'O':a.out_tag=optarg;break; case 'A':a.ins=std::stoi(optarg);break;
    case 'h':default:PrintHelp();std::exit(0);} }
}

int main(int argc,char* argv[]){
  std::cout<<std::scientific<<std::setprecision(15);
  Args a; ParseArgs(argc,argv,a); if(a.nu1<0.0) a.nu1=a.nu0;
  (void)a.gpu; CUDA_CHECK(cudaSetDevice(0));
  cudaDeviceProp prop; cudaGetDeviceProperties(&prop,0); std::cout<<"# dev = "<<prop.name<<"\n";
  constexpr Idx N=Comp::N; constexpr int Nt=Comp::Nt;
  const bool free_field=a.ens_dir.empty();
  std::cout<<"# LOCAL AXIAL current (bare sigma, P^dag t0-leg): mass=("<<a.mass_re<<","<<a.mass_im<<")  N="<<N
           <<(free_field?"  [free]":"  [interacting]")<<"  per-t build\n";

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;

  Base base(Comp::N_REFINE);
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();
  const double r=1.0,M5=-1.0,at=0.2;
  WilsonDirac DW(base,0.0,r,M5,at,a.nu1);                 // supplies the Pauli sigma matrices
  Gauge U(base);

  const int n_sites=(int)base.n_sites;
  std::vector<double> w_site(n_sites);
  for(int n=0;n<n_sites;n++){ w_site[n]=base.dual_areas[n]; }

  const bool single=(a.ins>=0);
  if(single && (a.ins>=n_sites)){
    std::cout<<"# ERROR: --ins="<<a.ins<<" out of range (n_sites="<<n_sites<<")\n"; return 1; }
  const int n_lo = single?a.ins   : 0;
  const int n_hi = single?a.ins+1 : n_sites;

  int n_t0=a.n_t0; std::vector<int> t0s(n_t0); for(int b=0;b<n_t0;b++) t0s[b]=b*(Nt/n_t0);

  std::string ens_base=a.ens_dir;
  if(!ens_base.empty()&&ens_base.back()=='/') ens_base.pop_back();
  { auto s=ens_base.find_last_of('/'); if(s!=std::string::npos) ens_base=ens_base.substr(s+1); }
  const std::string tag=(free_field?std::string("free"):ens_base);
  const std::string esnid=tag+"_vmRe"+std::to_string(a.mass_re)+"vmIm"+std::to_string(a.mass_im);
  const std::string propdir="data_"+esnid+"/prop_deter_L"+std::to_string(Comp::N_REFINE)+"/";
  const std::string locbase=single?std::string("corr_deter_local_axial1"):std::string("corr_deter_local_axial");
  const std::string locname=a.out_tag.empty()?locbase:locbase+"_"+a.out_tag;
  const std::string outdir ="data_"+esnid+"/"+locname+"_L"+std::to_string(Comp::N_REFINE)+"/";
  std::filesystem::create_directories(outdir);

  Timer timer;
  const int k_lo=0, k_step=free_field?1:a.ninter, k_hi=free_field?1:1000000;
  for(int k=k_lo;k<k_hi;k+=k_step){
    if(!free_field){
      const std::string lat=a.ens_dir+"ckpoint_lat."+std::to_string(k);
      if(!std::filesystem::exists(lat)){ if(k==0) continue; else break; }
      U.read(lat);
    }
    const std::string pfile=a.prop_file.empty()?propdir+"Dinv."+std::to_string(k)+".h5":a.prop_file;
    if(!std::filesystem::exists(pfile)){ std::cout<<"# k="<<k<<" no propagator "<<pfile<<"\n"; if(free_field) break; else continue; }
    const std::string h5path=outdir+"corr."+std::to_string(k)+".h5";
    if(std::filesystem::exists(h5path)){ bool c=false; try{HighFive::File f(h5path,HighFive::File::ReadOnly);c=f.exist("complete");}catch(...){}
      if(c){ std::cout<<"# skip k="<<k<<" (complete)\n"; if(free_field) break; else continue; } }

    std::vector<Complex> P; { HighFive::File f(pfile,HighFive::File::ReadOnly); load_mat(f,"Dm_inv",P); }
    std::cout<<"# k="<<k<<"  loaded P\n";
    // m_P: the t0-leg (backward propagator, Eq. 3.61) is (1+m_P)^{-1} tilde D_{m_P}^{-dag}.  Load Dtil_inv
    // -> P0 = Dtil_inv (tr_WPWP_axial daggers it to tilde D^{-dag}); the (1+m_P)^{-1} factor is applied to
    // the conn below.  massless/m_F: P0 = P = D_m^{-1} (t0-leg = D_m^{-dag}).
    const bool parity = (a.mass_im != 0.0);
    std::vector<Complex> Dtil;
    if(parity){
      HighFive::File f(pfile,HighFive::File::ReadOnly);
      if(!f.exist("Dtil_inv")){ std::cout<<"# parity but no Dtil_inv in "<<pfile<<" (rerun jj_propagator_deter)\n"; return 1; }
      load_mat(f,"Dtil_inv",Dtil);
      std::cout<<"# m_P: t0-leg uses tilde D^{-dag} (Dtil_inv) + (1+m_P)^{-1}\n";
    }
    const std::vector<Complex>& P0 = parity ? Dtil : P;     // t0-leg propagator (daggered inside tr_WPWP_axial)
    const Complex inv1pmP = parity ? Complex(1.0,0.0)/(Complex(1.0,0.0)+Complex(a.mass_re,a.mass_im))
                                   : Complex(1.0,0.0);

    const std::string h5tmp=h5path+".tmp";
    auto h5p=std::make_unique<HighFive::File>(h5tmp,HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
    HighFive::File& h5=*h5p;
    h5.createDataSet("t0s",t0s); h5.createDataSet("n_t0",std::vector<int>{n_t0});
    h5.createDataSet("nhits",std::vector<int>{1}); h5.createDataSet("ls",std::vector<int>{0,1,2});
    h5.createDataSet("ins",std::vector<int>{a.ins});

    // Three diagonal Pauli channels.  AXIAL: conn(t0,t) = sum_n w_site[n] tr[ sigma^a(n,t0) P^dag sigma^a(n,t) P ].
    // disc(t) = sum_n w_site[n] tr[ sigma^a(n,t) P ] (single-current vev; diagnostic only, same as vector).
    for(int a=1; a<=3; a++){
      const std::string chan = "s"+std::to_string(a);
      std::vector<std::vector<Complex>> Cpp(n_t0,std::vector<Complex>(Nt,Complex(0,0)));
      std::vector<Complex> discvec(Nt,Complex(0,0));
      std::vector<std::vector<Ent>> E0(n_t0);
      std::vector<Ent> Et;
      for(int n=n_lo; n<n_hi; n++){
        const double w = single?1.0:w_site[n];
        for(int b=0;b<n_t0;b++) local_W_sigma(E0[b], DW, a, t0s[b], (Idx)n);
        for(int t=0;t<Nt;t++){
          local_W_sigma(Et, DW, a, t, (Idx)n);
          Complex d(0,0);
          for(const auto& e : Et) d += e.v * P[(size_t)e.j*Comp::N + e.i];   // tr[sigma^a(n,t) P]
          discvec[t] += w * d;
          for(int b=0;b<n_t0;b++){ const int dt=((t-t0s[b])%Nt+Nt)%Nt; Cpp[b][dt] += w * tr_WPWP_axial(E0[b], Et, P0, P); }
        }
      }
      if(parity) for(int b=0;b<n_t0;b++) for(int t=0;t<Nt;t++) Cpp[b][t] *= inv1pmP;   // m_P (1+m_P)^{-1}
      for(int b=0;b<n_t0;b++){ const std::string kp="h0/t0_"+std::to_string(b)+"/";
        write_corr_axial(h5,kp+chan+"/Apm",Cpp[b]); }
      write_vec(h5,"h0/disc/"+chan+"/J",discvec);
      std::cout<<"#   "<<chan<<": disc(0)=("<<discvec[0].real()<<","<<discvec[0].imag()
               <<")  conn(dt=4)="<<Cpp[0][4].real()<<"  ["<<timer.currentSeconds()<<" s]\n";
    }

    // ---- axial ylm tower (Eq. 4.36): g_l(t) = (1/(2l+1)) sum_m tr[ Sigma_lm(t0) P^dag Sigma_lm(t) P ].
    constexpr int L_MAX_YLM = 2;
    for(int l=0; !single && l<=L_MAX_YLM; l++){
      std::vector<std::vector<Complex>> Cyl(n_t0,std::vector<Complex>(Nt,Complex(0,0)));
      std::vector<std::vector<Ent>> Sig0(n_t0);
      std::vector<Ent> Sigt;
      for(int m=-l; m<=l; m++){
        for(int b=0;b<n_t0;b++) build_Sigma_ylm(Sig0[b], base, w_site, l, m, t0s[b]);
        for(int t=0;t<Nt;t++){
          build_Sigma_ylm(Sigt, base, w_site, l, m, t);
          for(int b=0;b<n_t0;b++){ const int dt=((t-t0s[b])%Nt+Nt)%Nt; Cyl[b][dt] += tr_WPWP_axial(Sig0[b], Sigt, P0, P); }
        }
      }
      const double inv2lp1 = 1.0/(2.0*l+1.0);
      for(int b=0;b<n_t0;b++) for(int t=0;t<Nt;t++) Cyl[b][t] *= inv2lp1*inv1pmP;   // m_P: also (1+m_P)^{-1}
      for(int b=0;b<n_t0;b++){ const std::string kp="h0/t0_"+std::to_string(b)+"/ylm/l"+std::to_string(l)+"/";
        write_corr_axial(h5,kp+"Apm",Cyl[b]); }
      std::cout<<"#   ylm l="<<l<<": conn(dt=4)="<<Cyl[0][4].real()<<"  ["<<timer.currentSeconds()<<" s]\n";
    }
    h5.createDataSet("complete",std::vector<int>{1});
    h5p.reset(); std::filesystem::rename(h5tmp,h5path);
    std::cout<<"# wrote "<<h5path<<"\n";
    if(free_field) break;
  }
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
  return 0;
}
