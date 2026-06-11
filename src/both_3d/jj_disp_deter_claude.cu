// jj_disp_deter_claude.cu
// -----------------------------------------------------------------------------
// DISPLACED (link) current correlator with the EXACT all-to-all propagator.  BOTH channels:
//   sp = point-split SPATIAL link current ; tp = point-split TEMPORAL (time-link) current.
// Standalone, parallel to jj_local_deter_claude.cu (loc = on-site sigma^a) and the exact-K path.
//
// The displaced current is the genuine Wilson POINT-SPLIT link current (qed3int_v2-12, Eqs. 3.33-3.35):
//   W^{wz}_{xy} = (1/lambda_M)(d_xw d_yz C_wz - d_xz d_yw C_zw),  C_xy = -kappa_xy P_xy U_xy,
//   P_xy = (1/2)(1 - e^a_xy sigma_a) Omega_xy.
// We use the current with kappa DIVIDED OUT (Eq. 3.46 / classical_JV: J^{wz}_V ~ kappa^(0)_wz e^a_wz j^a;
// removing the spatially-dependent kappa^(0) coupling):
//   j^disp(x,y) = psibar(x) [ W^{wz}_xy / kappa_xy ] psi(y) = psibar(x) [ -P_xy U_xy ] psi(y).
// Legs sit on TWO distinct sites x != y (DISPLACED) -- this is the un-reduced (older) form of the loc
// current, which keeps the full Dirac structure (-r sigma_0 + e^a sigma_a), the spin connection Omega,
// and the gauge link U on the link, NOT collapsed to the on-site sigma of jj_local_deter.
//
// The bare kernel C/kappa = -P U is a routine of K: ConservedCurrent::build_W_ov_kappa (d_coo_format = i*C
// divided by i*kappa).  disp constructs a ConservedCurrent and exercises that routine -- the W
// DEFINITION lives in K; kappa is the relative lattice->continuum weight, factored out.
//
// Spatially-projected estimator (qed3int_v2-12, eq:spatially_projected_JJ):
//   sum_{(n,n') in tri} (ell ell* / kappa^(0)2) J_V(t0) J_V(t) ~ C_J (delta^ab - e3 e3) j^a j^b A_tri.
// With j^disp = J_V/kappa the 1/kappa^2 is already in the operator, so the sp per-link weight is
//   w_disp[il] = link_volume[il] = ell ell* = A_{nn'}  (the dual LINK area; sums to 4pi).
// conn_sp(t0,t) = sum_il w_disp[il] tr[ W_d(il,t0) P W_d(il,t) P ];  write_corr folds 1/(4pi).
//
// tp = TEMPORAL displaced current: W_d^t = C^t/kappa_t on the time-link (t,t+1) at site ix (build_W_disp_t,
// DiracExt::d_coo_format on {t,ix} / (i*kappa_t[ix])).  Sum over SITES with the dual SITE weight
// w_disp_t[ix] = dual_areas[ix].  disp thus supplies BOTH channels => its own G_s/G_t -> -(D-1) = -2.
//
// L COMPILE-TIME (-DN_REFINE_CLI).  CLI: --ens-dir(omit=free) --mass-re/--mass-im --n-t0 --ninter --gpu
//   --prop-file <path> (e.g. continuum cont_prop_L<L>/Dinv.0.h5) --out-tag <tag> --ins <i> (single).
// Output: data_<ESNID>/corr_deter_disp[1][_<tag>]_L<L>/corr.<k>.h5, keys h0/t0_b/{sp,tp}/{Vpp,Vmm}.
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

using BaseLink = std::array<Idx,2>;

// one sparse matrix entry of a current insertion W: contributes v to W[i,j].
struct Ent{ Idx i,j; Complex v; };

// DISPLACED (link) sp current W_d = C/kappa = -P U (Eqs. 3.33-3.35 / kappa_xy): the genuine Wilson
// point-split LINK current -(1/2)(1 - e^a sigma_a) Omega U on link (ix,iy), legs on TWO distinct sites.
// The kernel DEFINITION now lives in K (ConservedCurrent::build_W_ov_kappa) -- disp just exercises that
// routine of K, with kappa factored OUT as a weight.  These wrappers only repack the host COOEntry
// (CuC) into the local Ent (Complex) for the dense-propagator trace.
template<typename Kop, typename Gauge>
static void build_W_disp(std::vector<Ent>& en, const Kop& kop, const Gauge& U, int t, BaseLink lk){
  std::vector<COOEntry> coo;
  kop.build_W_ov_kappa(coo, U, std::pair<int,BaseLink>{t, lk});   // = C/kappa (normalized kernel W_ov_kappa)
  en.clear();
  for(const auto& c : coo) en.push_back({ c.i, c.j, Complex(c.v.x, c.v.y) });
}

// DISPLACED (link) TEMPORAL current W_d^t = C^t/kappa_t: the genuine Wilson point-split current along
// the TIME link, legs on the two timeslices t,t+1 at the SAME spatial site ix.  Same K routine
// (ConservedCurrent::build_W_ov_kappa, temporal overload {t,ix}); this wrapper only repacks COOEntry->Ent.
template<typename Kop, typename Gauge>
static void build_W_disp_t(std::vector<Ent>& en, const Kop& kop, const Gauge& U, int t, Idx ix){
  std::vector<COOEntry> coo;
  kop.build_W_ov_kappa(coo, U, std::pair<int,Idx>{t, ix});        // = C^t/kappa_t (normalized kernel W_ov_kappa)
  en.clear();
  for(const auto& c : coo) en.push_back({ c.i, c.j, Complex(c.v.x, c.v.y) });
}

// Insertion-DIAGONAL current-current trace for ONE link insertion:
//   tr[ W_d(t0) P W_d(t) P ] = sum_{(i,j,v0) in E0} sum_{(k,l,vt) in Et} v0 P_{jk} vt P_{li}.
// W_d touches the 2 sites of the link (|E0|,|Et| = 8), so this is O(64) dense-P lookups -- no dense build.
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
  // Lean load: scope each real/imag buffer so it frees before the next read (peak ~ M + one N^2
  // double, not M+real+imag) -- lets the dense P fit at L=4 (~41 GB).
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

struct Args{ double nu0=1.0,nu1=-1.0,mass_re=0.0,mass_im=0.0; std::string ens_dir,prop_file,out_tag; int n_t0=2,ninter=10,gpu=0,ins=-1; };
void PrintHelp(){ printf("jj_disp_deter: --mass-re --mass-im --ens-dir(omit=free) --n-t0 --ninter --nu0 --nu1 --gpu\n"
                         "  --prop-file <path>  read P from this exact file (e.g. continuum cont_prop_L<L>/Dinv.0.h5)\n"
                         "  --out-tag <tag>     output to corr_deter_disp_<tag>_L<L> (e.g. cont) instead of corr_deter_disp_L<L>\n"
                         "  --ins <i>           SINGLE-insertion mode at link i (no sum, no w_disp; RAW trace,\n"
                         "                      matches jj_exact_diag_deter_free sp).  -> corr_deter_disp1[_<tag>]_L<L>;\n"
                         "                      default -1 = full sum over links (current behaviour).\n"); }
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
  const bool parity=(a.mass_im!=0.0);
  std::cout<<"# DISPLACED (link) current W^{wz}/kappa (sp) + C^t/kappa_t (tp): mass=("<<a.mass_re<<","<<a.mass_im<<")  N="<<N
           <<(free_field?"  [free]":"  [interacting]")<<"  per-t build\n";

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;
  using Fermion=OverlapWMass<WilsonDirac>;

  Base base(Comp::N_REFINE);
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();
  const double r=1.0,M5=-1.0,at=0.2;
  WilsonDirac DW(base,0.0,r,M5,at,a.nu1);                 // supplies the Wilson W operator (kappa/gamma/Omega)
  Gauge U(base);
  // The bare current kernel is a ConservedCurrent routine (build_W_ov_kappa): disp exercises it.  Only DW
  // (kappa/gamma/Omega) is used by build_W_ov_kappa -- no overlap solve / D.update() is needed here; the
  // ConservedCurrent ctor's multishift scratch (~29 MB at L=4) is the accepted minor cost of homing W in K.
  Fermion D(DW, Complex(0.0), 21);
  ConservedCurrent<Fermion,Gauge> kop(D);

  const int n_links=(int)base.n_links;
  const int n_sites=(int)base.n_sites;
  // sp per-link weight = link_volume[il] = ell ell* = A_{nn'} (dual LINK area; the 1/kappa^2 of the
  // spatially-projected estimator is folded into W_d = W/kappa).
  std::vector<double> w_disp(n_links);
  for(int il=0;il<n_links;il++){ w_disp[il]=base.link_volume[il]; }
  // tp per-site weight = dual_areas[ix] (dual SITE area; the temporal current sums over sites, same
  // convention as loc's s3 tp weight w_site).
  std::vector<double> w_disp_t(n_sites);
  for(int ix=0;ix<n_sites;ix++){ w_disp_t[ix]=base.dual_areas[ix]; }

  // SINGLE-insertion mode (--ins i >= 0): one link (sp) / one site (tp), RAW trace (no weight).  Mirrors
  // jj_exact_diag_deter_free so disp1 at link/site i lines up with exact1 sp/tp at the same index i.
  // (n_sites < n_links always, so ins<n_sites is the binding constraint for both channels.)
  const bool single=(a.ins>=0);
  if(single && (a.ins>=n_sites)){
    std::cout<<"# ERROR: --ins="<<a.ins<<" out of range (n_sites="<<n_sites<<", n_links="<<n_links<<")\n"; return 1; }
  const int il_lo = single?a.ins   : 0;
  const int il_hi = single?a.ins+1 : n_links;
  const int it_lo = single?a.ins   : 0;
  const int it_hi = single?a.ins+1 : n_sites;

  int n_t0=a.n_t0; std::vector<int> t0s(n_t0); for(int b=0;b<n_t0;b++) t0s[b]=b*(Nt/n_t0);

  std::string ens_base=a.ens_dir;
  if(!ens_base.empty()&&ens_base.back()=='/') ens_base.pop_back();
  { auto s=ens_base.find_last_of('/'); if(s!=std::string::npos) ens_base=ens_base.substr(s+1); }
  const std::string tag=(free_field?std::string("free"):ens_base);
  const std::string esnid=tag+"_vmRe"+std::to_string(a.mass_re)+"vmIm"+std::to_string(a.mass_im);
  const std::string propdir="data_"+esnid+"/prop_deter_L"+std::to_string(Comp::N_REFINE)+"/";
  // --out-tag: corr_deter_disp_<tag>_L<L> (e.g. cont) instead of corr_deter_disp_L<L>.
  // single-insertion mode adds a "1": corr_deter_disp1[_<tag>]_L<L> (parallel to corr_deter_exact1).
  const std::string dispbase=single?std::string("corr_deter_disp1"):std::string("corr_deter_disp");
  const std::string dispname=a.out_tag.empty()?dispbase:dispbase+"_"+a.out_tag;
  const std::string outdir ="data_"+esnid+"/"+dispname+"_L"+std::to_string(Comp::N_REFINE)+"/";
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
    h5.createDataSet("nhits",std::vector<int>{1}); h5.createDataSet("ins",std::vector<int>{a.ins});

    // DISPLACED sp current, insertion-diagonal over spatial links a:
    //   conn(t0,t) = sum_a w_disp[a] tr[ W_d(a,t0) P W_d(a,t) P ],  disc(t) = sum_a w_disp[a] tr[ W_d(a,t) P ].
    // single-insertion mode (--ins i): one link i, w=1 (RAW), to match jj_exact_diag_deter_free sp.
    std::vector<std::vector<Complex>> Cpp(n_t0,std::vector<Complex>(Nt,Complex(0,0)));
    std::vector<Complex> discvec(Nt,Complex(0,0));
    std::vector<std::vector<Ent>> E0(n_t0);   // W_d(a,t0_b) for the current link a
    std::vector<Ent> Et;
    for(int il=il_lo; il<il_hi; il++){
      const double w = single?1.0:w_disp[il];
      for(int b=0;b<n_t0;b++) build_W_disp(E0[b], kop, U, t0s[b], base.links[il]);
      for(int t=0;t<Nt;t++){
        build_W_disp(Et, kop, U, t, base.links[il]);
        Complex d(0,0);
        for(const auto& e : Et) d += e.v * P[(size_t)e.j*Comp::N + e.i];   // tr[W_d(a,t) P]
        discvec[t] += w * d;
        for(int b=0;b<n_t0;b++){ const int dt=((t-t0s[b])%Nt+Nt)%Nt; Cpp[b][dt] += w * tr_WPWP(E0[b], Et, P); }
      }
    }
    for(int b=0;b<n_t0;b++){ const std::string kp="h0/t0_"+std::to_string(b)+"/sp/";
      write_corr(h5,kp+"Vpp",Cpp[b],false);
      if(!parity) write_corr(h5,kp+"Vmm",Cpp[b],true); }
    write_vec(h5,"h0/disc/sp/J",discvec);
    std::cout<<"#   disp(sp): disc(0)=("<<discvec[0].real()<<","<<discvec[0].imag()
             <<")  conn(dt=4)="<<Cpp[0][4].real()<<"  ["<<timer.currentSeconds()<<" s]\n";

    // DISPLACED tp current, insertion-diagonal over SITES ix (temporal time-link at ix):
    //   conn(t0,t) = sum_ix w_disp_t[ix] tr[ W_d^t(ix,t0) P W_d^t(ix,t) P ],  disc(t) likewise.
    // single-insertion mode (--ins i): one site i, w=1 (RAW), to match jj_exact_diag_deter_free tp.
    // disp now supplies BOTH channels, so disp's own G_s/G_t ratio -> -(D-1) = -2 as L grows.
    std::vector<std::vector<Complex>> Cpp_t(n_t0,std::vector<Complex>(Nt,Complex(0,0)));
    std::vector<Complex> discvec_t(Nt,Complex(0,0));
    std::vector<std::vector<Ent>> E0t(n_t0);   // W_d^t(ix,t0_b) for the current site ix
    std::vector<Ent> Ett;
    for(int ix=it_lo; ix<it_hi; ix++){
      const double w = single?1.0:w_disp_t[ix];
      for(int b=0;b<n_t0;b++) build_W_disp_t(E0t[b], kop, U, t0s[b], (Idx)ix);
      for(int t=0;t<Nt;t++){
        build_W_disp_t(Ett, kop, U, t, (Idx)ix);
        Complex d(0,0);
        for(const auto& e : Ett) d += e.v * P[(size_t)e.j*Comp::N + e.i];   // tr[W_d^t(ix,t) P]
        discvec_t[t] += w * d;
        for(int b=0;b<n_t0;b++){ const int dt=((t-t0s[b])%Nt+Nt)%Nt; Cpp_t[b][dt] += w * tr_WPWP(E0t[b], Ett, P); }
      }
    }
    for(int b=0;b<n_t0;b++){ const std::string kp="h0/t0_"+std::to_string(b)+"/tp/";
      write_corr(h5,kp+"Vpp",Cpp_t[b],false);
      if(!parity) write_corr(h5,kp+"Vmm",Cpp_t[b],true); }
    write_vec(h5,"h0/disc/tp/J",discvec_t);
    std::cout<<"#   disp(tp): disc(0)=("<<discvec_t[0].real()<<","<<discvec_t[0].imag()
             <<")  conn(dt=4)="<<Cpp_t[0][4].real()<<"  ["<<timer.currentSeconds()<<" s]\n";

    h5.createDataSet("complete",std::vector<int>{1});
    h5p.reset(); std::filesystem::rename(h5tmp,h5path);
    std::cout<<"# wrote "<<h5path<<"\n";
    if(free_field) break;
  }
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
  return 0;
}
