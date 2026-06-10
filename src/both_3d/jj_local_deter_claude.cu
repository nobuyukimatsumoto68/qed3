// jj_local_deter_claude.cu
// -----------------------------------------------------------------------------
// LOCAL (ultralocal) current correlators with the EXACT all-to-all propagator.
// (jj_exact_freefield_impl_plan_claude.md, "LOCAL (ultralocal) current path".)
//
// Replaces the full overlap conserved current K (Eq. 3.27 = ultralocal W^{wz} + non-local resolvents)
// by the BARE point-split vector current \hat e^a \bar\psi \sigma^a \psi -- i.e. ONLY the \hat e.\sigma
// hopping structure, NO overlap, NO multishift, and (per the note) NO spin connection \Omega and NO
// Wilson -r\sigma_0 term.
//   spatial link {ix,iy}: 0.5 kappa gamma(ix,iy) i exp(i u_sp)   (pos at ix, neg at iy)
//   temporal link n,t   : the sigma_3 part of the temporal hopping (i exp(i u_tp))
// This is DiracExt::d_coo_format with \Omega -> 1 and the -r\sigma_0 term dropped (vector current =
// the gamma part only).  We do NOT use ConservedCurrent::build_W (which dresses with \Omega + the
// Wilson term).  The projection weights w_tp/w_sp are kept (geometric projector, \Omega-independent);
// the overall scale is irrelevant for the 4.28/4.31 shape and the G_s/G_t ratio.  Sparse => feasible
// to L=4, free AND interacting (no N x n_ins x Nt multishift).
//
// Pipeline: load P=D_m^{-1} (jj_propagator_deter) -> build the SPARSE projected current
//   W^P(t) = sum_a w^P_a W(a,t)   (tp: sites/w_tp ; sp: links/w_sp)
// -> A^P(t) = W^P(t).P  (sparse row-scatter: A[i,:] += v P[j,:]) -> disc=tr(A), conn=tr(A(t0)A(t)).
// We build A(t) at EACH t explicitly (NOT a free-field cyclic time-shift of A(0)): the cyclic shift
// breaks at the ANTIPERIODIC time boundary.  Cheap here because the bare-gamma build is O(N).
// AXIAL (local): the GW factor (1-D_ov) -> 1 (D_ov = O(a), dropped), so the axial local current is
//   the same structure with no overlap dressing -- TODO (next; vector tp/sp here).
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

using BaseLink = std::array<Idx,2>;
using BaseFace = std::vector<Idx>;

// one entry of the sparse projected current W^P: A[i,:] += v * P[j,:].
struct Ent{ Idx i,j; Complex v; };

// --- bare-gamma local current entries (NO \Omega, NO -r\sigma_0) ------------------------------------
// Spatial link {ix,iy} at timeslice t.  Mirrors DiracExt::d_coo_format(pair<int,BaseLink>) with
// \Omega -> 1 and the Wilson -r\sigma_0 term dropped: the vector current is the gamma part only.
template<typename WilsonDirac, typename Gauge>
static void local_W_sp(std::vector<Ent>& en, const WilsonDirac& DW, const Gauge& U, int t, BaseLink lk){
  en.clear();
  const Idx ix=lk[0], iy=lk[1];
  const Idx il=DW.lattice.map2il.at(BaseLink{ix,iy});
  const Idx off=(Idx)Comp::Nx*t;
  // pos: row ix, col iy
  const MS p = 0.5 * DW.bd.kappa[il] * DW.bd.gamma(ix,iy) * I*std::exp( I*U.sp(t,BaseLink{ix,iy}) );
  en.push_back({off+NS*ix,   off+NS*iy,   p(0,0)});
  en.push_back({off+NS*ix,   off+NS*iy+1, p(0,1)});
  en.push_back({off+NS*ix+1, off+NS*iy,   p(1,0)});
  en.push_back({off+NS*ix+1, off+NS*iy+1, p(1,1)});
  // neg: row iy, col ix
  const MS q = -0.5 * DW.bd.kappa[il] * DW.bd.gamma(iy,ix) * I*std::exp( I*U.sp(t,BaseLink{iy,ix}) );
  en.push_back({off+NS*iy,   off+NS*ix,   q(0,0)});
  en.push_back({off+NS*iy,   off+NS*ix+1, q(0,1)});
  en.push_back({off+NS*iy+1, off+NS*ix,   q(1,0)});
  en.push_back({off+NS*iy+1, off+NS*ix+1, q(1,1)});
}

// Temporal link at dual site ix, timeslice t.  Mirrors DiracExt::d_coo_format(pair<int,Idx>) with the
// -r\sigma_0 term dropped: the vector current is the sigma_3 part only.  (No \Omega in the temporal
// hopping anyway.)  The wrap %N carries the t -> t+1 hop across the time boundary.
template<typename WilsonDirac, typename Gauge>
static void local_W_tp(std::vector<Ent>& en, const WilsonDirac& DW, const Gauge& U, int t, Idx ix){
  en.clear();
  const Idx Nx=Comp::Nx, N=Comp::N; const int Nt=Comp::Nt;
  int sign=1;
  if(t==Nt-1) sign=-1;
  const MS tmpP = -0.5 * sign * DW.kappa_t[ix] * ( -DW.sigma[3] ) * I*std::exp( -I*U.tp(t,ix) );
  const MS tmpM =  0.5 * sign * DW.kappa_t[ix] * (  DW.sigma[3] ) * I*std::exp(  I*U.tp(t,ix) );
  en.push_back({(Idx)((Nx*(t+1)+NS*ix)%N),   (Idx)(Nx*t+NS*ix),               tmpP(0,0)});
  en.push_back({(Idx)((Nx*(t+1)+NS*ix)%N),   (Idx)(Nx*t+NS*ix+1),             tmpP(0,1)});
  en.push_back({(Idx)((Nx*t+NS*ix + N)%N),   (Idx)((Nx*(t+1)+NS*ix)%N),       tmpM(0,0)});
  en.push_back({(Idx)((Nx*t+NS*ix + N)%N),   (Idx)((Nx*(t+1)+NS*ix+1)%N),     tmpM(0,1)});
  en.push_back({(Idx)((Nx*(t+1)+NS*ix+1)%N), (Idx)(Nx*t+NS*ix),               tmpP(1,0)});
  en.push_back({(Idx)((Nx*(t+1)+NS*ix+1)%N), (Idx)(Nx*t+NS*ix+1),             tmpP(1,1)});
  en.push_back({(Idx)((Nx*t+NS*ix+1 + N)%N), (Idx)((Nx*(t+1)+NS*ix)%N),       tmpM(1,0)});
  en.push_back({(Idx)((Nx*t+NS*ix+1 + N)%N), (Idx)((Nx*(t+1)+NS*ix+1)%N),     tmpM(1,1)});
}

// projected current W^P(t) = sum_a w^P_a W(a,t).  tp: sum over dual sites (w_tp); sp: over links (w_sp).
template<typename WilsonDirac, typename Gauge, typename Base>
static void build_Wp_local(const std::string& proj, int t, const WilsonDirac& DW, const Gauge& U,
                           const Base& base, const std::vector<double>& w_tp,
                           const std::vector<double>& w_sp, std::vector<Ent>& Wp){
  Wp.clear();
  std::vector<Ent> en;
  if(proj=="tp"){
    for(int n=0;n<(int)base.n_sites;n++){
      const double w=w_tp[n];
      local_W_tp(en, DW, U, t, (Idx)n);
      for(const auto& e : en) Wp.push_back({e.i, e.j, w*e.v});
    }
  } else {
    for(int il=0;il<(int)base.n_links;il++){
      const double w=w_sp[il];
      local_W_sp(en, DW, U, t, base.links[il]);
      for(const auto& e : en) Wp.push_back({e.i, e.j, w*e.v});
    }
  }
}

// A = W^P . P  (P row-major N^2):  A[i,c] += v * P[j,c].  A returned row-major (N^2), zero-init.
static void apply_Wp(const std::vector<Ent>& Wp, const std::vector<Complex>& P, std::vector<Complex>& A){
  const Idx N=Comp::N;
  A.assign((size_t)N*N, Complex(0,0));
  for(const auto& e : Wp){
    const size_t ai=(size_t)e.i*N, pj=(size_t)e.j*N;
    for(Idx c=0;c<N;c++) A[ai+c] += e.v * P[pj+c];
  }
}

static Complex trace(const std::vector<Complex>& A){
  const Idx N=Comp::N;
  Complex s(0,0);
  for(Idx i=0;i<N;i++) s+=A[(size_t)i*N+i];
  return s;
}

static Complex tr_AB(const std::vector<Complex>& A, const std::vector<Complex>& B){
  const Idx N=Comp::N;
  Complex s(0,0);
  for(Idx i=0;i<N;i++) for(Idx j=0;j<N;j++) s+=A[(size_t)i*N+j]*B[(size_t)j*N+i];
  return s;
}

static void load_mat(HighFive::File& f, const std::string& key, std::vector<Complex>& M){
  std::vector<double> re,im;
  f.getDataSet(key+"/real").read(re);
  f.getDataSet(key+"/imag").read(im);
  M.resize(re.size());
  for(size_t i=0;i<re.size();i++) M[i]=Complex(re[i],im[i]);
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

  const int n_sites=(int)base.n_sites, n_links=(int)base.n_links;
  std::vector<double> w_tp(n_sites), w_sp(n_links);
  for(int n=0;n<n_sites;n++){ const double kt=DW.kappa_t[n]; w_tp[n]=base.dual_areas[n]/(kt*kt); }
  for(int il=0;il<n_links;il++){ const double ks=DW.bd.kappa[il]; w_sp[il]=base.link_volume[il]/(ks*ks); }

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

    for(const std::string proj : {std::string("tp"), std::string("sp")}){
      std::vector<std::vector<Complex>> Cpp(n_t0,std::vector<Complex>(Nt,Complex(0,0)));
      std::vector<Complex> discvec(Nt,Complex(0,0));
      // Build the current at EACH t explicitly and form conn(t0,t)=tr(A(t0)A(t)).  No free-field
      // cyclic time-shift: it breaks at the ANTIPERIODIC time boundary.  Cheap (bare-gamma build).
      std::vector<Ent> Wp;
      std::vector<std::vector<Complex>> Ab(n_t0);
      for(int b=0;b<n_t0;b++){ build_Wp_local(proj,t0s[b],DW,U,base,w_tp,w_sp,Wp); apply_Wp(Wp,P,Ab[b]); }   // A(t0_b)
      for(int t=0;t<Nt;t++){ std::vector<Complex> At; build_Wp_local(proj,t,DW,U,base,w_tp,w_sp,Wp); apply_Wp(Wp,P,At);   // A(t)
        discvec[t]=trace(At);
        for(int b=0;b<n_t0;b++){ const int dt=((t-t0s[b])%Nt+Nt)%Nt; Cpp[b][dt]=tr_AB(Ab[b],At); } }
      for(int b=0;b<n_t0;b++){ const std::string kp="h0/t0_"+std::to_string(b)+"/";
        write_corr(h5,kp+proj+"/Vpp",Cpp[b],false);
        if(!parity) write_corr(h5,kp+proj+"/Vmm",Cpp[b],true); }
      write_vec(h5,"h0/disc/"+proj+"/J",discvec);
      std::cout<<"#   "<<proj<<": disc(0)=("<<discvec[0].real()<<","<<discvec[0].imag()
               <<")  conn(dt=4)="<<Cpp[0][4].real()<<"  ["<<timer.currentSeconds()<<" s]\n";
    }
    h5.createDataSet("complete",std::vector<int>{1});
    h5p.reset(); std::filesystem::rename(h5tmp,h5path);
    std::cout<<"# wrote "<<h5path<<"\n";
    if(free_field) break;
  }
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
  return 0;
}
