// jj_exact_diag_deter_free_claude.cu
// -----------------------------------------------------------------------------
// FREE-FIELD exact conserved-current correlator <J(a,0) J(a,t)>, tp + sp.  TWO modes:
//   default       : SINGLE insertion (--ins i), dense K CACHING.  (The single-link/site object; the
//                   q=0 DOUBLE sum of jj_kbuild_exact/jj_contract_deter was wrong.)
//   --sum         : DIAGONALLY-n-SUMMED (Eq.4.29) over ALL sites(tp)/links(sp), area-weighted
//                   (dual_areas/link_volume) + 1/4pi, K_ov_kappa per insertion, build-use-discard
//                   (no cache).  -> corr_deter_exactsum_L<L>.  EXPENSIVE: (n_sites+n_links)*N solves.
//                   CHECK: G^sp/G^tp -> -(D-1) = -2 (geometric sum restores the D-1 directions).
//
// For ONE insertion a (tp: site `ins`; sp: link base.links[ins]; default ins=0):
//   1. Build the DENSE K(a,0) (N x N) col-by-col: K[:,j] = op_K(a,0) e_j   (N op_K applies).
//      CACHE it to data_free_Kcache_L<L>/K_ins<ins>.h5 (datasets tp,sp).  K is mass-INDEPENDENT and
//      time-translatable, so it is built ONCE and reused: if the cache exists, READ it and SKIP the
//      (expensive) op_K applies entirely.
//   2. Divide K by the insertion kappa -> BARE conserved current K/kappa (kappa_sp[link] for sp,
//      kappa_t[site] for tp).  K is LINEAR in the single insertion W (W = -C/lambda_M, C = -kappa P U),
//      so this single scalar division removes the (kappa_sp/kappa_t)^2 ~ (a_t nu0)^2 anisotropy that
//      otherwise swamps the sp/tp ratio; K/kappa is the exact-current analog of disp's W_d = C/kappa.
//      (The K cache stays RAW/kappa-in -> reusable; the division is applied on every run.)
//   3. Load P = D_m^{-1} (= D_ov^{-1} at m=0).  A0 = (K/kappa) . P  (dense cuBLAS Zgemm).
//   4. <J(a,0)J(a,t)> = tr[(K/kappa)(0) P (K/kappa)(t) P] = conn_shift(A0, dt)   (free-field
//      antiperiodic time-shift: K(a,t) = S_t K(a,0) S_-t).  No projection weight, no sum.
//
// L COMPILE-TIME (-DN_REFINE_CLI).  FREE ONLY.  CLI: --ins <i> (default 0) --mass-re/--mass-im (P dir)
//   --prop-file --out-tag --n-t0 --gpu.  Output: data_<ESNID>/corr_deter_exact1[_<tag>]_L<L>/corr.0.h5,
//   keys h0/t0_b/{tp,sp}/{Vpp,Vmm} + h0/disc/{tp,sp}/J (jj-style; "1" = single insertion) + kappa_sp,
//   kappa_tp (multiply the correlator back by kappa^2 to recover the raw-K result).
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
static void save_mat(HighFive::File& f, const std::string& key, const std::vector<Complex>& M){
  std::vector<double> re(M.size()), im(M.size());
  for(size_t i=0;i<M.size();i++){ re[i]=M[i].real(); im[i]=M[i].imag(); }
  f.createDataSet(key+"/real",re); f.createDataSet(key+"/imag",im);
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

// conjugate-transpose of a row-major NxN: (M^dag)_{ij} = conj(M_{ji}).  Used for the m_P "-" channel:
// the dagger kernel K^dag = conjT(K) and tilde D_{m_P}^{-dag} = conjT(Dtil_inv) (no extra solves).
static std::vector<Complex> conj_transpose(const std::vector<Complex>& M){
  const Idx N=Comp::N;
  std::vector<Complex> Md((size_t)N*N);
  for(Idx i=0;i<N;i++) for(Idx j=0;j<N;j++) Md[(size_t)i*N+j]=std::conj(M[(size_t)j*N+i]);
  return Md;
}

// A = K . P  (both row-major) via cuBLAS Zgemm.  (P,K) to a column-major gemm gives row-major K P.
static void matmul_A(cublasHandle_t cub, CuC* d_K, CuC* d_P, CuC* d_A,
                     const std::vector<Complex>& K, const std::vector<Complex>& P, std::vector<Complex>& A){
  const int n=(int)Comp::N;
  const CuC one=make_cuDoubleComplex(1.0,0.0), zero=make_cuDoubleComplex(0.0,0.0);
  CUDA_CHECK(cudaMemcpy(d_K, reinterpret_cast<const CuC*>(K.data()), (size_t)n*n*sizeof(CuC), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_P, reinterpret_cast<const CuC*>(P.data()), (size_t)n*n*sizeof(CuC), cudaMemcpyHostToDevice));
  cublasZgemm(cub, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &one, d_P, n, d_K, n, &zero, d_A, n);
  A.resize((size_t)n*n);
  CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(A.data()), d_A, (size_t)n*n*sizeof(CuC), cudaMemcpyDeviceToHost));
}

// conn(dt) = tr(A . S_dt A S_-dt) for the free field (antiperiodic), A = K(a,0) P.
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

struct Args{ double nu0=1.0,nu1=-1.0,mass_re=0.0,mass_im=0.0; std::string ens_dir,prop_file,out_tag; int n_t0=2,gpu=0,ins=0; bool sum=false; };
void PrintHelp(){ printf("jj_exact_diag_deter_free: FREE exact <J(a,0)J(a,t)>, tp+sp, with K cache.\n"
                         "  --ins <i>          single-insertion index (tp: site i; sp: link i).  default 0\n"
                         "  --sum              DIAGONALLY-n-SUMMED (Eq.4.29): sum over ALL sites(tp)/links(sp),\n"
                         "                     area-weighted (dual_areas/link_volume) + 1/4pi, K_ov_kappa per\n"
                         "                     insertion.  build-use-discard (no cache).  -> corr_deter_exactsum_L<L>.\n"
                         "                     EXPENSIVE: (n_sites+n_links)*N op_K solves (L=1 min, L=2 hours).\n"
                         "  --mass-re/--mass-im  selects the P dir + esnid\n"
                         "  --prop-file <path> read P from this exact file (e.g. cont_prop_L<L>/Dinv.0.h5)\n"
                         "  --out-tag <tag>    corr_deter_exact{1,sum}_<tag>_L<L>;  --n-t0 --gpu\n"); }
void ParseArgs(int argc,char**argv,Args&a){
  static struct option lo[]={{"nu0",required_argument,0,'n'},{"nu1",required_argument,0,'m'},
    {"mass-re",required_argument,0,'r'},{"mass-im",required_argument,0,'i'},{"ens-dir",required_argument,0,'e'},
    {"n-t0",required_argument,0,'T'},{"gpu",required_argument,0,'G'},{"ins",required_argument,0,'A'},
    {"prop-file",required_argument,0,'P'},{"out-tag",required_argument,0,'O'},{"sum",no_argument,0,'S'},
    {"help",no_argument,0,'h'},{0,0,0,0}};
  int opt,idx;
  while((opt=getopt_long(argc,argv,"n:m:r:i:e:T:G:A:P:O:Sh",lo,&idx))!=-1){ switch(opt){
    case 'n':a.nu0=std::stod(optarg);break; case 'm':a.nu1=std::stod(optarg);break;
    case 'r':a.mass_re=std::stod(optarg);break; case 'i':a.mass_im=std::stod(optarg);break;
    case 'e':a.ens_dir=optarg;break; case 'T':a.n_t0=std::stoi(optarg);break;
    case 'G':a.gpu=std::stoi(optarg);break; case 'A':a.ins=std::stoi(optarg);break;
    case 'P':a.prop_file=optarg;break; case 'O':a.out_tag=optarg;break; case 'S':a.sum=true;break;
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
  if(!free_field){ std::cout<<"# ERROR: free-field only.\n"; return 1; }
  if(!a.sum) std::cout<<"# EXACT single-insertion (ins="<<a.ins<<") <J(a,0)J(a,t)>, tp+sp:  N="<<N<<"  [free]\n";

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;
  using Fermion=OverlapWMass<WilsonDirac>;

  Base base(Comp::N_REFINE);
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();
  const double r=1.0,M5=-1.0,at=0.2;
  WilsonDirac DW(base,0.0,r,M5,at,a.nu1);
  Gauge U(base);

  const int n_sites=(int)base.n_sites, n_links=(int)base.n_links;
  if(!a.sum && (a.ins<0 || a.ins>=n_sites || a.ins>=n_links)){
    std::cout<<"# ERROR: ins="<<a.ins<<" out of range (n_sites="<<n_sites<<", n_links="<<n_links<<")\n"; return 1; }

  int n_t0=a.n_t0; std::vector<int> t0s(n_t0); for(int b=0;b<n_t0;b++) t0s[b]=b*(Nt/n_t0);

  // K operator.  Built unconditionally so kop.insertion_kappa (the single kappa lookup, used for the
  // K_ov_kappa normalization below) is available on the cache-HIT path too; the D.update()/solves only
  // run on a cache miss.  (massless D_ov; K is mass-independent.)
  Fermion D(DW, Complex(0.0), 21);
  ConservedCurrent<Fermion,Gauge> kop(D);

  // ===================== DIAGONALLY-n-SUMMED mode (--sum), Eq.(4.29) =====================
  // G^{tp/sp}(t) = (1/4pi) sum_{n in sites/links} A_n tr[ K_ov_kappa,n(0) P K_ov_kappa,n(t) P ],
  //   A_n = dual_areas (tp) / link_volume (sp);  1/(4pi) folded by write_corr.
  // Per insertion: build dense K (N op_K applies) / kappa = K_ov_kappa, A = K.P, accumulate
  //   A_n * conn_shift(A, dt).  BUILD-USE-DISCARD (no cache: all-n caching ~300 GB at L=2).
  // EXPENSIVE: (n_sites+n_links)*N solves.  CHECK: G^sp/G^tp -> -(D-1) = -2 (geometric sum restores the
  //   D-1 transverse directions; single-insertion was -1).  Single-insertion path below is untouched.
  if(a.sum){
    std::cout<<"# EXACT DIAGONAL-SUM (--sum): area-weighted Eq.(4.29), build-use-discard.  N="<<N
             <<"  n_sites="<<n_sites<<" n_links="<<n_links<<"\n";
    MatPoly op_K; op_K.push_back(cplx(1.0), {&kop});
    D.update(U);

    const std::string esnid="free_vmRe"+std::to_string(a.mass_re)+"vmIm"+std::to_string(a.mass_im);
    const std::string propdir="data_"+esnid+"/prop_deter_L"+std::to_string(Comp::N_REFINE)+"/";
    const std::string pfile=a.prop_file.empty()?propdir+"Dinv.0.h5":a.prop_file;
    if(!std::filesystem::exists(pfile)){ std::cout<<"# no propagator "<<pfile<<"\n"; return 1; }
    std::vector<Complex> P; { HighFive::File f(pfile,HighFive::File::ReadOnly); load_mat(f,"Dm_inv",P); }
    std::cout<<"# loaded P "<<pfile<<"\n";
    const std::string dname=a.out_tag.empty()?std::string("corr_deter_exactsum")
                                             :std::string("corr_deter_exactsum_")+a.out_tag;
    const std::string outdir="data_"+esnid+"/"+dname+"_L"+std::to_string(Comp::N_REFINE)+"/";
    std::filesystem::create_directories(outdir);
    const std::string h5path=outdir+"corr.0.h5";

    cublasHandle_t cub; cublasCreate(&cub);
    CuC *d_K,*d_P,*d_A;
    CUDA_CHECK(cudaMalloc(&d_K,(size_t)N*N*sizeof(CuC)));
    CUDA_CHECK(cudaMalloc(&d_P,(size_t)N*N*sizeof(CuC)));
    CUDA_CHECK(cudaMalloc(&d_A,(size_t)N*N*sizeof(CuC)));

    Timer timer;
    const std::string h5tmp=h5path+".tmp";
    auto h5p=std::make_unique<HighFive::File>(h5tmp,HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
    HighFive::File& h5=*h5p;
    h5.createDataSet("t0s",t0s); h5.createDataSet("n_t0",std::vector<int>{n_t0});
    h5.createDataSet("nhits",std::vector<int>{1}); h5.createDataSet("summed",std::vector<int>{1});

    std::vector<Complex> K((size_t)N*N), ej(N), out(N), A0;
    const long total_cols=(long)(n_sites+n_links)*(long)N;   // total op_K solves (1 per K column)
    long done_cols=0;
    std::cout<<"#   total op_K solves = (n_sites+n_links)*N = "<<total_cols<<"  (build-use-discard)\n";

    // Ylm tower (Eq.4.36): accumulate Sigma_{l,m} = sum_n A_n Y_lm(n) K^t_ov_kappa,n DURING the tp site
    // loop (reuses the per-site K -- no extra solves).  l=0,1,2 -> n_lm=9 dense (l,m) channels.
    // Memory: n_lm * N^2 complex ~ 17 GB at L=2 (L=1 trivial).
    constexpr int L_MAX_YLM = 2;
    std::vector<std::pair<int,int>> lm;
    for(int l=0;l<=L_MAX_YLM;l++) for(int m=-l;m<=l;m++) lm.push_back({l,m});
    const int n_lm=(int)lm.size();
    std::vector<std::vector<Complex>> Sigma(n_lm);
    for(int c=0;c<n_lm;c++) Sigma[c].assign((size_t)N*N, Complex(0,0));

    for(int which=0; which<2; which++){
      const std::string proj=(which==0)?"tp":"sp";
      const int n_ins=(which==0)?n_sites:n_links;
      std::vector<Complex> cdt(Nt, Complex(0,0));     // summed conn(dt) = sum_n A_n conn_n(dt)
      Complex disc(0,0);
      for(int ins=0; ins<n_ins; ins++){
        double w, kappa;
        if(which==0){ kop.set_temporal(U, 0, (Idx)ins, false);
                      w=base.dual_areas[ins];  kappa=kop.insertion_kappa(std::pair<int,Idx>{0,(Idx)ins}); }
        else        { kop.set_spatial (U, 0, base.links[ins], false);
                      w=base.link_volume[ins]; kappa=kop.insertion_kappa(std::pair<int,BaseLink>{0,base.links[ins]}); }
        // dense K_ov_kappa(ins,0) col-by-col  (each column = one op_K multishift solve)
        for(Idx j=0;j<N;j++){
          std::fill(ej.begin(), ej.end(), Complex(0,0)); ej[j]=Complex(1,0);
          op_K.from_cpu<N>(out.data(), ej.data());
          for(Idx i=0;i<N;i++) K[(size_t)i*N+j]=out[i]/kappa;
          ++done_cols;
          if(j%256==0) std::cout<<"#   "<<proj<<" ins "<<(ins+1)<<"/"<<n_ins<<" col "<<j<<"/"<<N
                                <<"  ("<<(int)(100.0*done_cols/total_cols)<<"% of solves)  ["
                                <<timer.currentSeconds()<<" s]\n";
        }
        if(which==0){   // Ylm: accumulate Sigma_{l,m} += A_n Y_lm(n) K^t_ov_kappa,n  (n = site ins)
          const VE site=base.sites[ins];
          for(int c=0;c<n_lm;c++){
            const double wy = w * Ylm_real(lm[c].first, lm[c].second, site);
            Complex* S=Sigma[c].data();
            #pragma omp parallel for schedule(static)
            for(long idx=0; idx<(long)N*N; idx++) S[idx]+=wy*K[(size_t)idx];
          }
        }
        matmul_A(cub, d_K, d_P, d_A, K, P, A0);       // A0 = K_ov_kappa . P
        disc += w*trace(A0);
        std::vector<Complex> cn(Nt);
        #pragma omp parallel for schedule(dynamic)
        for(int dt=0;dt<Nt;dt++) cn[dt]=conn_shift(A0,dt);
        for(int dt=0;dt<Nt;dt++) cdt[dt]+=w*cn[dt];
        std::cout<<"#   "<<proj<<" ins "<<(ins+1)<<"/"<<n_ins<<" DONE  conn(dt=4)+="<<(w*cn[4]).real()
                 <<"  ("<<(int)(100.0*done_cols/total_cols)<<"% of solves)  ["<<timer.currentSeconds()<<" s]\n";
      }
      // free field: time-translation invariant -> map the single conn(dt) onto each t0 origin.
      std::vector<std::vector<Complex>> Cpp(n_t0,std::vector<Complex>(Nt,Complex(0,0)));
      for(int b=0;b<n_t0;b++) for(int t=0;t<Nt;t++){ const int dt=((t-t0s[b])%Nt+Nt)%Nt; Cpp[b][dt]=cdt[dt]; }
      std::vector<Complex> discvec(Nt, disc);
      for(int b=0;b<n_t0;b++){ const std::string kp="h0/t0_"+std::to_string(b)+"/";
        write_corr(h5,kp+proj+"/Vpp",Cpp[b],false);
        if(!parity) write_corr(h5,kp+proj+"/Vmm",Cpp[b],true); }
      write_vec(h5,"h0/disc/"+proj+"/J",discvec);
      std::cout<<"#   "<<proj<<" SUM done: conn(dt=4)="<<cdt[4].real()<<"  ["<<timer.currentSeconds()<<" s]\n";
    }

    // ---- Ylm tower output (Eq.4.36) from the accumulated Sigma_{l,m}:
    //   g_l(t) = (1/(2l+1)) sum_m tr[ Sigma_lm(0) P Sigma_lm(t) P ] = (1/(2l+1)) sum_m conn_shift(Sigma_lm . P).
    for(int l=0; l<=L_MAX_YLM; l++){
      std::vector<Complex> gl(Nt, Complex(0,0));
      for(int c=0;c<n_lm;c++){
        if(lm[c].first!=l) continue;
        std::vector<Complex> A0lm; matmul_A(cub, d_K, d_P, d_A, Sigma[c], P, A0lm);
        #pragma omp parallel for schedule(dynamic)
        for(int dt=0;dt<Nt;dt++) gl[dt]+=conn_shift(A0lm,dt);
      }
      const double inv2lp1=1.0/(2.0*l+1.0);
      for(int t=0;t<Nt;t++) gl[t]*=inv2lp1;
      std::vector<std::vector<Complex>> Cyl(n_t0,std::vector<Complex>(Nt,Complex(0,0)));
      for(int b=0;b<n_t0;b++) for(int t=0;t<Nt;t++){ const int dt=((t-t0s[b])%Nt+Nt)%Nt; Cyl[b][dt]=gl[dt]; }
      for(int b=0;b<n_t0;b++){ const std::string kp="h0/t0_"+std::to_string(b)+"/ylm/l"+std::to_string(l)+"/";
        write_corr(h5,kp+"Vpp",Cyl[b],false);
        if(!parity) write_corr(h5,kp+"Vmm",Cyl[b],true); }
      std::cout<<"#   exact ylm l="<<l<<": conn(dt=4)="<<Cyl[0][4].real()<<"  ["<<timer.currentSeconds()<<" s]\n";
    }

    h5.createDataSet("complete",std::vector<int>{1});
    h5p.reset(); std::filesystem::rename(h5tmp,h5path);
    std::cout<<"# wrote "<<h5path<<"\n";
    cudaFree(d_K); cudaFree(d_P); cudaFree(d_A); cublasDestroy(cub);
    for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
    return 0;
  }
  // ======================================================================================

  // ---- K cache (mass-INDEPENDENT, so a mass-free dir; per insertion) ----
  const std::string kcdir = "data_free_Kcache_L"+std::to_string(Comp::N_REFINE)+"/";
  std::filesystem::create_directories(kcdir);
  const std::string kfile = kcdir + "K_ins"+std::to_string(a.ins)+".h5";

  std::vector<Complex> K_tp, K_sp;
  bool have_cache=false;
  if(std::filesystem::exists(kfile)){
    try{ HighFive::File kf(kfile,HighFive::File::ReadOnly);
      if(kf.exist("complete")){ load_mat(kf,"tp",K_tp); load_mat(kf,"sp",K_sp); have_cache=true; } }catch(...){}
  }
  if(have_cache){
    std::cout<<"# K cache HIT  "<<kfile<<"  (skipping op_K applies)\n";
  } else {
    std::cout<<"# K cache MISS "<<kfile<<"  -> building dense K(a,0) for tp+sp (N op_K applies each) ...\n";
    MatPoly op_K; op_K.push_back(cplx(1.0), {&kop});
    D.update(U);
    Timer kt;
    std::vector<Complex> ej(N), out(N);
    for(int which=0; which<2; which++){
      const std::string proj = (which==0)?"tp":"sp";
      std::vector<Complex>& K = (which==0)?K_tp:K_sp;
      K.assign((size_t)N*N, Complex(0,0));
      if(proj=="tp") kop.set_temporal(U, 0, (Idx)a.ins, /*dag=*/false);
      else           kop.set_spatial (U, 0, base.links[a.ins], /*dag=*/false);
      for(Idx j=0;j<N;j++){
        std::fill(ej.begin(), ej.end(), Complex(0,0)); ej[j]=Complex(1,0);
        op_K.from_cpu<N>(out.data(), ej.data());       // out = K(a,0) e_j  = column j of K
        for(Idx i=0;i<N;i++) K[(size_t)i*N+j]=out[i];
        if(j%512==0) std::cout<<"#   "<<proj<<" col "<<j<<"/"<<N<<"  ["<<kt.currentSeconds()<<" s]\n";
      }
    }
    const std::string ktmp=kfile+".tmp";
    { HighFive::File kf(ktmp,HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
      kf.createDataSet("ins",std::vector<int>{a.ins}); kf.createDataSet("N",std::vector<int>{(int)N});
      save_mat(kf,"tp",K_tp); save_mat(kf,"sp",K_sp);
      kf.createDataSet("complete",std::vector<int>{1}); }
    std::filesystem::rename(ktmp,kfile);
    std::cout<<"# wrote K cache "<<kfile<<"  ["<<kt.currentSeconds()<<" s]\n";
  }

  // ---- normalize to K_ov_kappa = K / kappa (the bare conserved current) ----
  // K is LINEAR in the single current insertion W (every apply_k term applies W or W^dag once),
  // and W = -C/lambda_M with C = -kappa P U, so K_sp \propto kappa[link] and K_tp \propto kappa_t[site].
  // The JJ correlator is then \propto kappa^2, and the raw sp/tp ratio is swamped by (kappa_sp/kappa_t)^2
  // ~ (a_t nu0)^2.  Dividing the dense K by its insertion kappa (a single scalar, by linearity) gives
  // K_ov_kappa = K/kappa, directly comparable to the displaced normalized kernel W_ov_kappa = C/kappa
  // (disp).  kappa comes from kop.insertion_kappa -- the SAME single lookup K uses for build_W_ov_kappa.
  // The K cache is kept RAW (kappa-in) so it stays reusable; the division is applied on every run.
  const double kappa_sp = kop.insertion_kappa(std::pair<int,BaseLink>{0, base.links[a.ins]});
  const double kappa_tp = kop.insertion_kappa(std::pair<int,Idx>{0, (Idx)a.ins});
  for(auto& z : K_sp) z /= kappa_sp;
  for(auto& z : K_tp) z /= kappa_tp;
  std::cout<<"# K_ov_kappa = K/kappa:  kappa_sp(link "<<a.ins<<")="
           <<kappa_sp<<"  kappa_t(site "<<a.ins<<")="<<kappa_tp<<"\n";

  // ---- propagator P ----
  std::string ens_base=a.ens_dir;
  if(!ens_base.empty()&&ens_base.back()=='/') ens_base.pop_back();
  { auto s=ens_base.find_last_of('/'); if(s!=std::string::npos) ens_base=ens_base.substr(s+1); }
  const std::string tag=std::string("free");
  const std::string esnid=tag+"_vmRe"+std::to_string(a.mass_re)+"vmIm"+std::to_string(a.mass_im);
  const std::string propdir="data_"+esnid+"/prop_deter_L"+std::to_string(Comp::N_REFINE)+"/";
  const std::string pfile=a.prop_file.empty()?propdir+"Dinv.0.h5":a.prop_file;
  if(!std::filesystem::exists(pfile)){ std::cout<<"# no propagator "<<pfile<<"\n"; return 1; }
  const std::string dname=a.out_tag.empty()?std::string("corr_deter_exact1")
                                           :std::string("corr_deter_exact1_")+a.out_tag;
  const std::string outdir="data_"+esnid+"/"+dname+"_L"+std::to_string(Comp::N_REFINE)+"/";
  std::filesystem::create_directories(outdir);
  const std::string h5path=outdir+"corr.0.h5";

  std::vector<Complex> P; { HighFive::File f(pfile,HighFive::File::ReadOnly); load_mat(f,"Dm_inv",P); }
  std::cout<<"# loaded P "<<pfile<<"\n";
  // m_P: the "-" channel (Eq. 3.65) uses tilde D_{m_P}^{-dag} = conjT(Dtil_inv) and K^dag = conjT(K).
  std::vector<Complex> Ptil;
  if(parity){
    std::vector<Complex> Dtil;
    { HighFive::File f(pfile,HighFive::File::ReadOnly);
      if(!f.exist("Dtil_inv")){ std::cout<<"# parity but no Dtil_inv in "<<pfile<<" (rerun jj_propagator_deter)\n"; return 1; }
      load_mat(f,"Dtil_inv",Dtil); }
    Ptil = conj_transpose(Dtil);                 // tilde D_{m_P}^{-dag}
    std::cout<<"# loaded Dtil_inv (parity, m_P) -> Ptil = tilde D^{-dag}\n";
  }

  // ---- A0 = K . P  (dense) ; conn(dt) = conn_shift(A0,dt) = tr[K(a,0) P K(a,dt) P] (single insertion) ----
  cublasHandle_t cub; cublasCreate(&cub);
  CuC *d_K,*d_P,*d_A; CUDA_CHECK(cudaMalloc(&d_K,(size_t)N*N*sizeof(CuC)));
  CUDA_CHECK(cudaMalloc(&d_P,(size_t)N*N*sizeof(CuC))); CUDA_CHECK(cudaMalloc(&d_A,(size_t)N*N*sizeof(CuC)));

  Timer timer;
  const std::string h5tmp=h5path+".tmp";
  auto h5p=std::make_unique<HighFive::File>(h5tmp,HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
  HighFive::File& h5=*h5p;
  h5.createDataSet("t0s",t0s); h5.createDataSet("n_t0",std::vector<int>{n_t0});
  h5.createDataSet("nhits",std::vector<int>{1}); h5.createDataSet("ins",std::vector<int>{a.ins});
  // kappa used to make K bare (K/kappa); multiply back by kappa^2 to recover the raw-K correlator.
  h5.createDataSet("kappa_sp",std::vector<double>{kappa_sp});
  h5.createDataSet("kappa_tp",std::vector<double>{kappa_tp});

  for(int which=0; which<2; which++){
    const std::string proj=(which==0)?"tp":"sp";
    const std::vector<Complex>& K=(which==0)?K_tp:K_sp;
    std::vector<Complex> A0; matmul_A(cub, d_K, d_P, d_A, K, P, A0);
    const Complex disc=trace(A0);
    std::vector<Complex> cdt(Nt);
    #pragma omp parallel for schedule(dynamic)
    for(int dt=0;dt<Nt;dt++) cdt[dt]=conn_shift(A0,dt);
    std::vector<std::vector<Complex>> Cpp(n_t0,std::vector<Complex>(Nt,Complex(0,0)));
    std::vector<Complex> discvec(Nt, disc);
    for(int b=0;b<n_t0;b++) for(int t=0;t<Nt;t++){ const int dt=((t-t0s[b])%Nt+Nt)%Nt; Cpp[b][dt]=cdt[dt]; }
    for(int b=0;b<n_t0;b++){ const std::string kp="h0/t0_"+std::to_string(b)+"/";
      write_corr(h5,kp+proj+"/Vpp",Cpp[b],false); }
    write_vec(h5,"h0/disc/"+proj+"/J",discvec);
    if(!parity){   // massless / m_F: the "-" channel is the elementwise conjugate of "+"
      for(int b=0;b<n_t0;b++) write_corr(h5,"h0/t0_"+std::to_string(b)+"/"+proj+"/Vmm",Cpp[b],true);
    } else {       // m_P: Vmm (Eq. 3.65) = conn_shift( K^dag . tilde D^{-dag} ) ; disc_mm = tr(.)
      std::vector<Complex> Kdag = conj_transpose(K);
      std::vector<Complex> A_mm; matmul_A(cub, d_K, d_P, d_A, Kdag, Ptil, A_mm);
      const Complex disc_mm = trace(A_mm);
      std::vector<Complex> cdt_mm(Nt);
      #pragma omp parallel for schedule(dynamic)
      for(int dt=0;dt<Nt;dt++) cdt_mm[dt]=conn_shift(A_mm,dt);
      std::vector<std::vector<Complex>> Cmm(n_t0,std::vector<Complex>(Nt,Complex(0,0)));
      for(int b=0;b<n_t0;b++) for(int t=0;t<Nt;t++){ const int dt=((t-t0s[b])%Nt+Nt)%Nt; Cmm[b][dt]=cdt_mm[dt]; }
      for(int b=0;b<n_t0;b++) write_corr(h5,"h0/t0_"+std::to_string(b)+"/"+proj+"/Vmm",Cmm[b],false);
      write_vec(h5,"h0/disc/"+proj+"/Jmm",std::vector<Complex>(Nt,disc_mm));
    }
    std::cout<<"#   "<<proj<<": disc(0)=("<<discvec[0].real()<<","<<discvec[0].imag()
             <<")  conn(dt=4)="<<Cpp[0][4].real()<<"  ["<<timer.currentSeconds()<<" s]\n";
  }
  h5.createDataSet("complete",std::vector<int>{1});
  h5p.reset(); std::filesystem::rename(h5tmp,h5path);
  std::cout<<"# wrote "<<h5path<<"\n";

  cudaFree(d_K); cudaFree(d_P); cudaFree(d_A); cublasDestroy(cub);
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
  return 0;
}
