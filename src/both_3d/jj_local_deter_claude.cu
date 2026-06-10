// jj_local_deter_claude.cu
// -----------------------------------------------------------------------------
// LOCAL (ultralocal) current correlators with the EXACT all-to-all propagator.
// (jj_exact_freefield_impl_plan_claude.md, "LOCAL (ultralocal) current path".)
//
// Replaces the full overlap conserved current K (Eq. 3.27 = ultralocal W^{wz} + non-local resolvents)
// by the BARE vector current \bar\psi \sigma^a \psi, with NO overlap and NO multishift.
// The LOCAL current is a pure SITE (ultralocal) object: at dual site n the insertion is the bare
// Pauli \sigma^a (a=1,2,3) -- NO link hop, NO \kappa, NO i, NO gauge phase, NO \Omega, NO Wilson
// -r\sigma_0.  Both propagator legs attach at the SAME site n (on-site G((n,t0),(n,t))).  We measure
// the three DIAGONAL channels SEPARATELY (no tp/sp projection):
//   s3 = <\sigma_3 \sigma_3>             (the old "tp"/temporal projection; G_t, Eq. 4.31)
//   s1 = <\sigma_1 \sigma_1>, s2 = <\sigma_2 \sigma_2>   (spatial; G_s = s1 + s2, Eq. 4.28)
// All channels use the SITE weight w_site[n] = dual_areas[n] (no 1/\kappa^2; the bare \sigma carries
// no \kappa to cancel).  Overall scale is irrelevant for the 4.28/4.31 shape and the G_s/G_t ratio.
//
// Pipeline: load P=D_m^{-1} (jj_propagator_deter) -> for each channel a (1,2,3) and dual site n:
//   conn(t0,t) = sum_n w_site[n] tr[ \sigma^a(n,t0) P \sigma^a(n,t) P ]
//   disc(t)    = sum_n w_site[n] tr[ \sigma^a(n,t) P ]
// W(n,t) is built at EACH t explicitly (NOT a cyclic time-shift): the shift breaks at the
// ANTIPERIODIC time boundary.  Cheap here because the bare-\sigma build is O(N).
// AXIAL (local): the GW factor (1-D_ov) -> 1 (D_ov = O(a), dropped), so the axial local current is
//   the same structure with no overlap dressing -- TODO (next; vector s1/s2/s3 here).
//
// L COMPILE-TIME (-DN_REFINE_CLI).  CLI: --ens-dir(omit=free) --mass-re/--mass-im (P + esnid) --n-t0
//   --ninter --gpu.   Output: data_<ESNID>/corr_deter_local_L<L>/corr.<k>.h5  (jj-style; atomic).
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
#include "sparse_dirac.h"
#include "matpoly_claude.h"
#include "dirac_pf.h"
#include "overlap_wmass_claude.h"
#include "conserved_current_claude.h"

// one sparse matrix entry of a current insertion W: contributes v to W[i,j].
struct Ent{ Idx i,j; Complex v; };

// LOCAL (site) current insertion: bare Pauli \sigma^a (a=1,2,3) at dual site ix, timeslice t.
// A pure SITE object for ALL channels -- both propagator legs attach at site ix (no link hop, no
// \kappa, no i, no gauge phase, no \Omega, no Wilson -r\sigma_0).  \sigma_1,\sigma_2 are the spatial
// channels (G_s); \sigma_3 the temporal (G_t).  DW.sigma[a]: 1=\sigma_x, 2=\sigma_y, 3=\sigma_z.
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

// Ylm-weighted TEMPORAL current Sigma_{l,m}(t) = sum_n A_n Y_lm(n^) sigma_3(n,t), as Ent entries:
// the diagonal sigma_3 = diag(+1,-1) site block scaled by A_n Y_lm(n^).  Feeds the Eq.(4.35) tower
// g_l(t) = (1/(2l+1)) sum_m tr[Sigma_{l,m}(t0) P Sigma_{l,m}(t) P].  Ylm_real (valence_claude.h) is the
// real spherical harmonic (no Condon-Shortley), App. C convention; base.sites[n] is the unit VE.
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

// Insertion-DIAGONAL current-current trace for ONE site insertion n:
//   tr[ W(n,t0) P W(n,t) P ] = sum_{(i,j,v0) in E0} sum_{(k,l,vt) in Et} v0 P_{jk} vt P_{li}.
// E0 = W(n,t0) entries, Et = W(n,t) entries, P row-major (N^2).  W is the 2x2 spin block at
// site n (|E0|,|Et| <= 4), so this is O(16) dense-P lookups -- no N x N dense build.
static Complex tr_WPWP(const std::vector<Ent>& E0, const std::vector<Ent>& Et,
                       const std::vector<Complex>& P){
  const Idx N=Comp::N;
  Complex s(0,0);
  for(const auto& e0 : E0)
    for(const auto& et : Et)
      s += e0.v * P[(size_t)e0.j*N + et.i] * et.v * P[(size_t)et.j*N + e0.i];
  return s;
}

static void load_mat(HighFive::File& f, const std::string& key, std::vector<Complex>& M){
  // Lean load: scope each real/imag buffer so it frees before the next read.  Peak RAM = M (N^2
  // complex) + ONE N^2 double buffer, not both real+imag at once -- lets the dense P fit at L=4
  // (~41 GB instead of ~55 GB; cont_prop_L4 is 26 GB on disk).
  size_t n=0;
  { std::vector<double> re; f.getDataSet(key+"/real").read(re); n=re.size();
    M.resize(n); for(size_t i=0;i<n;i++) M[i].real(re[i]); }
  { std::vector<double> im; f.getDataSet(key+"/imag").read(im);
    for(size_t i=0;i<n;i++) M[i].imag(im[i]); }
}

static void write_corr(HighFive::File& h5, const std::string& key, const std::vector<Complex>& C, bool conj){
  const int Nt=Comp::Nt;
  const double inv4pi=1.0/(4.0*M_PI);
  std::vector<double> re(Nt),im(Nt);
  for(int t=0;t<Nt;t++){ Complex g=inv4pi*C[t]; re[t]=g.real(); im[t]=conj?-g.imag():g.imag(); }
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

struct Args{ double nu0=1.0,nu1=-1.0,mass_re=0.0,mass_im=0.0; std::string ens_dir,prop_file,out_tag; int n_t0=2,ninter=10,gpu=0; };
void PrintHelp(){ printf("jj_local_deter: --mass-re --mass-im --ens-dir(omit=free) --n-t0 --ninter --nu0 --nu1 --gpu\n"
                         "  --prop-file <path>  read P from this exact file (e.g. continuum cont_prop_L<L>/Dinv.0.h5)\n"
                         "  --out-tag <tag>     output to corr_deter_local_<tag>_L<L> (e.g. cont) instead of corr_deter_local_L<L>\n"); }
void ParseArgs(int argc,char**argv,Args&a){
  static struct option lo[]={{"nu0",required_argument,0,'n'},{"nu1",required_argument,0,'m'},
    {"mass-re",required_argument,0,'r'},{"mass-im",required_argument,0,'i'},{"ens-dir",required_argument,0,'e'},
    {"n-t0",required_argument,0,'T'},{"ninter",required_argument,0,'I'},{"gpu",required_argument,0,'G'},
    {"prop-file",required_argument,0,'P'},{"out-tag",required_argument,0,'O'},
    {"help",no_argument,0,'h'},{0,0,0,0}};
  int opt,idx;
  while((opt=getopt_long(argc,argv,"n:m:r:i:e:T:I:G:P:O:h",lo,&idx))!=-1){ switch(opt){
    case 'n':a.nu0=std::stod(optarg);break; case 'm':a.nu1=std::stod(optarg);break;
    case 'r':a.mass_re=std::stod(optarg);break; case 'i':a.mass_im=std::stod(optarg);break;
    case 'e':a.ens_dir=optarg;break; case 'T':a.n_t0=std::stoi(optarg);break;
    case 'I':a.ninter=std::stoi(optarg);break; case 'G':a.gpu=std::stoi(optarg);break;
    case 'P':a.prop_file=optarg;break; case 'O':a.out_tag=optarg;break;
    case 'h':default:PrintHelp();std::exit(0);} }
}

int main(int argc,char* argv[]){
  std::cout<<std::scientific<<std::setprecision(15);
  Args a; ParseArgs(argc,argv,a); if(a.nu1<0.0) a.nu1=a.nu0;
  (void)a.gpu; CUDA_CHECK(cudaSetDevice(0));
  cudaDeviceProp prop; cudaGetDeviceProperties(&prop,0); std::cout<<"# dev = "<<prop.name<<"\n";
  constexpr Idx N=Comp::N; constexpr int Nt=Comp::Nt;
  const bool free_field=a.ens_dir.empty();
  const bool parity=(a.mass_im!=0.0);
  std::cout<<"# LOCAL current (bare e.sigma, no Omega): mass=("<<a.mass_re<<","<<a.mass_im<<")  N="<<N
           <<(free_field?"  [free]":"  [interacting]")<<"  per-t build\n";

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;

  Base base(Comp::N_REFINE);
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();
  const double r=1.0,M5=-1.0,at=0.2;
  WilsonDirac DW(base,0.0,r,M5,at,a.nu1);                 // supplies gamma/kappa for the bare current
  Gauge U(base);

  const int n_sites=(int)base.n_sites;
  // LOCAL current is a SITE object for ALL channels => single site weight, no link weight.
  // Bare geometric measure dual_areas[n] (no 1/\kappa^2: the bare \sigma carries no \kappa).
  std::vector<double> w_site(n_sites);
  for(int n=0;n<n_sites;n++){ w_site[n]=base.dual_areas[n]; }

  int n_t0=a.n_t0; std::vector<int> t0s(n_t0); for(int b=0;b<n_t0;b++) t0s[b]=b*(Nt/n_t0);

  std::string ens_base=a.ens_dir;
  if(!ens_base.empty()&&ens_base.back()=='/') ens_base.pop_back();
  { auto s=ens_base.find_last_of('/'); if(s!=std::string::npos) ens_base=ens_base.substr(s+1); }
  const std::string tag=(free_field?std::string("free"):ens_base);
  const std::string esnid=tag+"_vmRe"+std::to_string(a.mass_re)+"vmIm"+std::to_string(a.mass_im);
  const std::string propdir="data_"+esnid+"/prop_deter_L"+std::to_string(Comp::N_REFINE)+"/";
  // --out-tag: corr_deter_local_<tag>_L<L> (e.g. cont) instead of corr_deter_local_L<L>
  const std::string locname=a.out_tag.empty()?std::string("corr_deter_local")
                                             :std::string("corr_deter_local_")+a.out_tag;
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
    // --prop-file overrides the propagator source (e.g. continuum cont_prop_L<L>/Dinv.0.h5)
    const std::string pfile=a.prop_file.empty()?propdir+"Dinv."+std::to_string(k)+".h5":a.prop_file;
    if(!std::filesystem::exists(pfile)){ std::cout<<"# k="<<k<<" no propagator "<<pfile<<"\n"; if(free_field) break; else continue; }
    const std::string h5path=outdir+"corr."+std::to_string(k)+".h5";
    if(std::filesystem::exists(h5path)){ bool c=false; try{HighFive::File f(h5path,HighFive::File::ReadOnly);c=f.exist("complete");}catch(...){}
      if(c){ std::cout<<"# skip k="<<k<<" (complete)\n"; if(free_field) break; else continue; } }

    std::vector<Complex> P; { HighFive::File f(pfile,HighFive::File::ReadOnly); load_mat(f,"Dm_inv",P); }
    std::cout<<"# k="<<k<<"  loaded P\n";

    const std::string h5tmp=h5path+".tmp";
    auto h5p=std::make_unique<HighFive::File>(h5tmp,HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
    HighFive::File& h5=*h5p;
    h5.createDataSet("t0s",t0s); h5.createDataSet("n_t0",std::vector<int>{n_t0});
    h5.createDataSet("nhits",std::vector<int>{1}); h5.createDataSet("ls",std::vector<int>{0,1,2});

    // Three diagonal Pauli channels, each a SITE object summed over dual sites n with w_site[n]:
    //   conn(t0,t) = sum_n w_site[n] tr[ \sigma^a(n,t0) P \sigma^a(n,t) P ]   (same site n both times)
    //   disc(t)    = sum_n w_site[n] tr[ \sigma^a(n,t) P ]
    // s3 = old "tp" (G_t); s1,s2 = spatial (G_s = s1 + s2).
    for(int a=1; a<=3; a++){
      const std::string chan = "s"+std::to_string(a);
      std::vector<std::vector<Complex>> Cpp(n_t0,std::vector<Complex>(Nt,Complex(0,0)));
      std::vector<Complex> discvec(Nt,Complex(0,0));
      std::vector<std::vector<Ent>> E0(n_t0);   // \sigma^a(n,t0_b) for the current site n
      std::vector<Ent> Et;
      for(int n=0; n<n_sites; n++){
        const double w = w_site[n];
        for(int b=0;b<n_t0;b++) local_W_sigma(E0[b], DW, a, t0s[b], (Idx)n);
        for(int t=0;t<Nt;t++){
          local_W_sigma(Et, DW, a, t, (Idx)n);
          Complex d(0,0);
          for(const auto& e : Et) d += e.v * P[(size_t)e.j*Comp::N + e.i];   // tr[\sigma^a(n,t) P]
          discvec[t] += w * d;
          for(int b=0;b<n_t0;b++){ const int dt=((t-t0s[b])%Nt+Nt)%Nt; Cpp[b][dt] += w * tr_WPWP(E0[b], Et, P); }
        }
      }
      for(int b=0;b<n_t0;b++){ const std::string kp="h0/t0_"+std::to_string(b)+"/";
        write_corr(h5,kp+chan+"/Vpp",Cpp[b],false);
        if(!parity) write_corr(h5,kp+chan+"/Vmm",Cpp[b],true); }
      write_vec(h5,"h0/disc/"+chan+"/J",discvec);
      std::cout<<"#   "<<chan<<": disc(0)=("<<discvec[0].real()<<","<<discvec[0].imag()
               <<")  conn(dt=4)="<<Cpp[0][4].real()<<"  ["<<timer.currentSeconds()<<" s]\n";
    }

    // ---- Ylm tower (Eq. 4.35): spherical-harmonic descendants of the temporal sigma_3 correlator.
    // Diagonal-m Legendre coefficient (connected only):
    //   g_l(t) = (1/(2l+1)) sum_{m=-l}^{l} tr[ Sigma_{l,m}(t0) P Sigma_{l,m}(t) P ]
    // with Sigma_{l,m}(t) = sum_n A_n Y_lm(n^) sigma_3(n,t).  write_corr folds 1/(4pi) so g_l matches
    // jj_cft_ylm_check_claude.cc / Eq.(4.35): rates (l=0->0, l=1->e^{-2t}, l=2->e^{-3t}),
    // G22 e^{3t}/G11 e^{2t} -> 12/5 = 2.4.  Off-site pairs enter via tr_WPWP (Sigma touches all sites).
    constexpr int L_MAX_YLM = 2;
    for(int l=0; l<=L_MAX_YLM; l++){
      std::vector<std::vector<Complex>> Cyl(n_t0,std::vector<Complex>(Nt,Complex(0,0)));
      std::vector<std::vector<Ent>> Sig0(n_t0);
      std::vector<Ent> Sigt;
      for(int m=-l; m<=l; m++){
        for(int b=0;b<n_t0;b++) build_Sigma_ylm(Sig0[b], base, w_site, l, m, t0s[b]);
        for(int t=0;t<Nt;t++){
          build_Sigma_ylm(Sigt, base, w_site, l, m, t);
          for(int b=0;b<n_t0;b++){ const int dt=((t-t0s[b])%Nt+Nt)%Nt; Cyl[b][dt] += tr_WPWP(Sig0[b], Sigt, P); }
        }
      }
      const double inv2lp1 = 1.0/(2.0*l+1.0);
      for(int b=0;b<n_t0;b++) for(int t=0;t<Nt;t++) Cyl[b][t] *= inv2lp1;
      for(int b=0;b<n_t0;b++){ const std::string kp="h0/t0_"+std::to_string(b)+"/ylm/l"+std::to_string(l)+"/";
        write_corr(h5,kp+"Vpp",Cyl[b],false);
        if(!parity) write_corr(h5,kp+"Vmm",Cyl[b],true); }
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
