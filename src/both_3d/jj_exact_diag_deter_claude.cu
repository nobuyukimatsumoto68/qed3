// jj_exact_diag_deter_claude.cu
// -----------------------------------------------------------------------------
// INSERTION-DIAGONAL deterministic EXACT conserved-current (K) correlator, tp + sp.
// Aligns with jj_local_deter (loc) and jj_disp_deter (disp) and the production stochastic
// jj_corr_block_t -- and reproduces Eqs. (4.31)/(4.28) -- UNLIKE jj_kbuild_exact+jj_contract_deter
// which contract the SUMMED K^P (a q=0 DOUBLE sum over insertions).
//
//   C^P(t0->t) = sum_a w^P_a tr[ K(a,t0) P K(a,t) P ]            (P=Dm^{-1}=D_ov^{-1} at m=0)
//   tp: a over sites, w_tp = dual_areas/kappa_t^2  (Eq. 4.32)
//   sp: a over links, w_sp = link_volume/kappa^2   (Eq. 4.29)
// K(a,t) = full overlap conserved current (ConservedCurrent kop + op_K, multishift resolvent) -- the
// SAME kernel/weights as jj_kbuild_exact; only the contraction is diagonal (same insertion a at t0,t).
//
// Per insertion a: B_a = K(a,0).P (apply op_K to each COLUMN of P; N applies => total N*n_ins, the
// SAME op_K count as jj_kbuild_exact).  Free field (time-translation invariant):
//   C^P(dt) = sum_a w_a conn_shift(B_a, dt)   (ANTIPERIODIC time-shift; from jj_contract_deter)
// disc^P(t) = sum_a w_a tr(B_a) (time-independent for free).  One B_a (N x N) live at a time.
//
// L COMPILE-TIME (-DN_REFINE_CLI).  CLI: --ens-dir(omit=free; interacting NOT supported here, free-
//   field check only) --mass-re/--mass-im (selects P dir + esnid) --prop-file --out-tag --n-t0 --gpu.
// Output: data_<ESNID>/corr_deter_exactdiag[_<tag>]_L<L>/corr.<k>.h5, keys h0/t0_b/{tp,sp}/{Vpp,Vmm}
//         + h0/disc/{tp,sp}/J  (jj-style; DISTINCT from the double-sum corr_deter_exact_L*).
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

// ---- time-shift helpers (ANTIPERIODIC; from jj_contract_deter): time is the OUTER block of size Nx ----
static inline Idx tshift(Idx i, int dt){
  const Idx Nx=Comp::Nx; const int Nt=Comp::Nt;
  const int t=(int)(i/Nx); const Idx x=i%Nx;
  const int tt=((t+dt)%Nt+Nt)%Nt;
  return (Idx)tt*Nx + x;
}
static inline int twrap_sign(int t, int dt){
  const int Nt=Comp::Nt; int s=t+dt, w=0;
  while(s<0){ s+=Nt; w++; } while(s>=Nt){ s-=Nt; w++; }
  return (w%2==0)?1:-1;
}

static void load_mat(HighFive::File& f, const std::string& key, std::vector<Complex>& M){
  size_t n=0;
  { std::vector<double> re; f.getDataSet(key+"/real").read(re); n=re.size();
    M.resize(n); for(size_t i=0;i<n;i++) M[i].real(re[i]); }
  { std::vector<double> im; f.getDataSet(key+"/imag").read(im);
    for(size_t i=0;i<n;i++) M[i].imag(im[i]); }
}

static void write_corr(HighFive::File& h5, const std::string& key, const std::vector<Complex>& C, bool conj){
  const int Nt=Comp::Nt; const double inv4pi=1.0/(4.0*M_PI);
  std::vector<double> re(Nt),im(Nt);
  for(int t=0;t<Nt;t++){ Complex g=inv4pi*C[t]; re[t]=g.real(); im[t]=conj?-g.imag():g.imag(); }
  h5.createDataSet(key+"/real",re); h5.createDataSet(key+"/imag",im);
}
static void write_vec(HighFive::File& h5, const std::string& key, const std::vector<Complex>& C){
  const int Nt=Comp::Nt; std::vector<double> re(Nt),im(Nt);
  for(int t=0;t<Nt;t++){ re[t]=C[t].real(); im[t]=C[t].imag(); }
  h5.createDataSet(key+"/real",re); h5.createDataSet(key+"/imag",im);
}

static Complex trace(const std::vector<Complex>& A){
  const Idx N=Comp::N; Complex s(0,0);
  for(Idx i=0;i<N;i++) s+=A[(size_t)i*N+i];
  return s;
}

// conn(dt) = tr(B . S_dt B S_-dt) for the free field (antiperiodic), B = K(a,0) P (from jj_contract_deter).
static Complex conn_shift(const std::vector<Complex>& A, int dt){
  const Idx N=Comp::N, Nx=Comp::Nx; const int Nt=Comp::Nt;
  std::vector<int> sgn(Nt); for(int t=0;t<Nt;t++) sgn[t]=twrap_sign(t,-dt);
  Complex s(0,0);
  for(Idx p=0;p<N;p++){
    const int sp=sgn[(int)(p/Nx)];
    for(Idx q=0;q<N;q++){
      const int sq=sgn[(int)(q/Nx)];
      s += A[(size_t)p*N+q] * (double)(sp*sq) * A[(size_t)tshift(q,-dt)*N + tshift(p,-dt)];
    }
  }
  return s;
}

struct Args{ double nu0=1.0,nu1=-1.0,mass_re=0.0,mass_im=0.0; std::string ens_dir,prop_file,out_tag; int n_t0=2,ninter=10,gpu=0; };
void PrintHelp(){ printf("jj_exact_diag_deter: --mass-re --mass-im --ens-dir(omit=free) --n-t0 --ninter --nu0 --nu1 --gpu\n"
                         "  --prop-file <path>  read P from this exact file (e.g. continuum cont_prop_L<L>/Dinv.0.h5)\n"
                         "  --out-tag <tag>     output to corr_deter_exactdiag_<tag>_L<L> instead of corr_deter_exactdiag_L<L>\n"); }
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
  if(!free_field){ std::cout<<"# ERROR: interacting not supported (free-field check only; needs per-t K).\n"; return 1; }
  std::cout<<"# EXACT-DIAG conserved current (insertion-diagonal): mass=("<<a.mass_re<<","<<a.mass_im<<")  N="<<N<<"  [free]\n";

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;
  using Fermion=OverlapWMass<WilsonDirac>;

  Base base(Comp::N_REFINE);
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();
  const double r=1.0,M5=-1.0,at=0.2;
  WilsonDirac DW(base,0.0,r,M5,at,a.nu1);
  Gauge U(base);
  Fermion D(DW, Complex(0.0), 21);                 // massless D_ov; K is mass-independent
  ConservedCurrent<Fermion,Gauge> kop(D);
  MatPoly op_K; op_K.push_back(cplx(1.0), {&kop});

  const int n_sites=(int)base.n_sites, n_links=(int)base.n_links;
  std::vector<double> w_tp(n_sites), w_sp(n_links);                          // Eq. 4.32 tp, Eq. 4.29 sp
  for(int n=0;n<n_sites;n++){ const double kt=DW.kappa_t[n];  w_tp[n]=base.dual_areas[n]/(kt*kt); }
  for(int il=0;il<n_links;il++){ const double ks=DW.bd.kappa[il]; w_sp[il]=base.link_volume[il]/(ks*ks); }

  int n_t0=a.n_t0; std::vector<int> t0s(n_t0); for(int b=0;b<n_t0;b++) t0s[b]=b*(Nt/n_t0);

  std::string ens_base=a.ens_dir;
  if(!ens_base.empty()&&ens_base.back()=='/') ens_base.pop_back();
  { auto s=ens_base.find_last_of('/'); if(s!=std::string::npos) ens_base=ens_base.substr(s+1); }
  const std::string tag=(free_field?std::string("free"):ens_base);
  const std::string esnid=tag+"_vmRe"+std::to_string(a.mass_re)+"vmIm"+std::to_string(a.mass_im);
  const std::string propdir="data_"+esnid+"/prop_deter_L"+std::to_string(Comp::N_REFINE)+"/";
  const std::string dname=a.out_tag.empty()?std::string("corr_deter_exactdiag")
                                           :std::string("corr_deter_exactdiag_")+a.out_tag;
  const std::string outdir ="data_"+esnid+"/"+dname+"_L"+std::to_string(Comp::N_REFINE)+"/";
  std::filesystem::create_directories(outdir);

  Timer timer;
  const int k_lo=0, k_step=1, k_hi=1;     // free: single config
  for(int k=k_lo;k<k_hi;k+=k_step){
    const std::string pfile=a.prop_file.empty()?propdir+"Dinv."+std::to_string(k)+".h5":a.prop_file;
    if(!std::filesystem::exists(pfile)){ std::cout<<"# no propagator "<<pfile<<"\n"; break; }
    const std::string h5path=outdir+"corr."+std::to_string(k)+".h5";
    if(std::filesystem::exists(h5path)){ bool c=false; try{HighFive::File f(h5path,HighFive::File::ReadOnly);c=f.exist("complete");}catch(...){}
      if(c){ std::cout<<"# skip k="<<k<<" (complete)\n"; break; } }

    std::vector<Complex> P; { HighFive::File f(pfile,HighFive::File::ReadOnly); load_mat(f,"Dm_inv",P); }
    std::cout<<"# k="<<k<<"  loaded P\n";
    D.update(U);

    const std::string h5tmp=h5path+".tmp";
    auto h5p=std::make_unique<HighFive::File>(h5tmp,HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
    HighFive::File& h5=*h5p;
    h5.createDataSet("t0s",t0s); h5.createDataSet("n_t0",std::vector<int>{n_t0});
    h5.createDataSet("nhits",std::vector<int>{1}); h5.createDataSet("ls",std::vector<int>{0,1,2});

    std::vector<Complex> colj(N), out(N), B((size_t)N*N);
    for(const std::string proj : {std::string("tp"), std::string("sp")}){
      const int nins=(proj=="tp")?n_sites:n_links;
      std::vector<Complex> cdt(Nt, Complex(0,0));   // free: conn depends only on dt
      Complex disc0(0,0);
      for(int ins=0; ins<nins; ins++){
        if(proj=="tp") kop.set_temporal(U, 0, (Idx)ins, /*dag=*/false);
        else           kop.set_spatial (U, 0, base.links[ins], /*dag=*/false);
        // B = K(ins,0) . P : apply op_K to each column of P
        for(Idx j=0;j<N;j++){
          for(Idx i=0;i<N;i++) colj[i]=P[(size_t)i*N+j];
          op_K.from_cpu<N>(out.data(), colj.data());
          for(Idx i=0;i<N;i++) B[(size_t)i*N+j]=out[i];
        }
        const double w=(proj=="tp")?w_tp[ins]:w_sp[ins];
        disc0 += w*trace(B);
        // each dt independent (cdt[dt] unique) -> parallel over the Nt time-shifts (CPU)
        #pragma omp parallel for schedule(dynamic)
        for(int dt=0;dt<Nt;dt++) cdt[dt] += w*conn_shift(B, dt);
        std::cout<<"#   "<<proj<<" ins="<<ins<<"/"<<nins<<"  ["<<timer.currentSeconds()<<" s]\n";
      }
      std::vector<std::vector<Complex>> Cpp(n_t0,std::vector<Complex>(Nt,Complex(0,0)));
      std::vector<Complex> discvec(Nt, disc0);
      for(int b=0;b<n_t0;b++) for(int t=0;t<Nt;t++){ const int dt=((t-t0s[b])%Nt+Nt)%Nt; Cpp[b][dt]=cdt[dt]; }
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
  }
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
  return 0;
}
