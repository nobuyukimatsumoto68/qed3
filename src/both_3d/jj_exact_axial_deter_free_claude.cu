// jj_exact_axial_deter_free_claude.cu
// -----------------------------------------------------------------------------
// FREE-FIELD exact AXIAL conserved-current correlator  C_{A+-}(t0 -> t), tp + sp.  AXIAL sibling of
// jj_exact_diag_deter_free_claude.cu (vector C_{V++}); same dense / time-shift machinery, same K cache.
//
//   C_{A+-,c}(t0 -> t) = tr[ K(t0) (1 - D_ov^dag) D_ov^{-dag} K^dag(t) (1 - D_ov) D_ov^{-1} ]   (Eq. 3.55)
//
// vs the vector C_{V++,c}(t0 -> t) = tr[ K(t0) D_ov^{-1} K(t) D_ov^{-1} ] (Eq. 3.52).  The axial differs in
// three places: GW factors (1 - D_ov^dag) on t0 and (1 - D_ov) on t; the t0-leg propagator is D_ov^{-dag}
// (= P^dag) not P; the SINK kernel is K^dag(t) not K(t).
//   Ref: qed3int_v2-13.pdf Sec. 3.3.2; audit jj_corr_pdf_audit_claude.md; plan
//        jj_axial_trio_deter_impl_plan_claude.md.  Dense free-field contraction = our own conn_shift
//        (antiperiodic time-shift) from jj_exact_diag_deter_free_claude.cu.
//
// Dense build (per insertion a):
//   1. K(a,0) dense col-by-col (N op_K applies) -> divide by kappa = K_ov_kappa (normalized current,
//      exactly as the vector exactsum, so the axial trio compares normalized currents).  In --sum:
//      build-use-discard (no cache).  In single --ins: read/write data_free_Kcache_L<L>/K_ins<i>.h5.
//   2. G = (1 - D_ov) dense col-by-col (N op_oneMinusD applies), built ONCE (mass- and insertion-
//      INDEPENDENT) -> cache data_free_Gcache_L<L>/G.h5, reused across insertions and tp/sp.
//   3. K^dag = conj-transpose(K) (no extra applies).  Source  B_src = K G^dag P^dag (two matmuls);
//      sink B_snk = K^dag G P (two matmuls).  G^dag, P^dag are host conj-transposes.
//   4. C_{A+-}(dt) = conn_shift2(B_src, B_snk, dt) = tr[B_src S_dt B_snk S_-dt]
//                  = tr[K(0) G^dag P^dag K^dag(dt) G P]  (Eq. 3.55 with t0=0; valid since G,P are
//                    time-translation invariant -> B_snk(t)=S_t B_snk(0) S_-t).
//
// C_{A+-} is NOT self-reflection-even (reflects to C_{A-+}, Eq. 3.57) -> SINGLE complex channel "Apm"
// (no Vmm).  --sum also builds the axial ylm tower (Eq. 4.36) riding the temporal Sigma_{l,m}.
//
// L COMPILE-TIME (-DN_REFINE_CLI).  FREE ONLY.  CLI: --ins <i> --sum --mass-re/--mass-im (P dir)
//   --prop-file --out-tag --n-t0 --gpu.  Output: data_<ESNID>/corr_deter_exactsum_axial[_<tag>]_L<L>/corr.0.h5
//   (--sum) or corr_deter_exact1_axial[_<tag>]_L<L>/corr.0.h5 (single), keys h0/t0_b/{tp,sp}/Apm
//   + h0/t0_b/ylm/l{l}/Apm + h0/disc/{tp,sp}/J + kappa_sp,kappa_tp.
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

// AXIAL is a single complex channel "Apm" (C_{A+-} reflects to C_{A-+}, Eq. 3.57; no Vmm).  1/(4pi) folded.
static void write_corr_axial(HighFive::File& h5, const std::string& key, const std::vector<Complex>& C){
  const int Nt=Comp::Nt; const double inv4pi=1.0/(4.0*M_PI);
  std::vector<double> re(Nt),im(Nt);
  for(int t=0;t<Nt;t++){ Complex g=inv4pi*C[t]; re[t]=g.real(); im[t]=g.imag(); }
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

// host conjugate-transpose:  At = A^dag,  At[i,j] = conj(A[j,i])  (row-major N x N).
static void conj_transpose(const std::vector<Complex>& A, std::vector<Complex>& At){
  const Idx N=Comp::N; At.resize((size_t)N*N);
  #pragma omp parallel for schedule(static)
  for(Idx i=0;i<N;i++)
    for(Idx j=0;j<N;j++) At[(size_t)i*N+j]=std::conj(A[(size_t)j*N+i]);
}

// C = L . R  (all row-major) via cuBLAS Zgemm.  Passing (R,L) to a column-major gemm gives row-major L R.
static void matmul(cublasHandle_t cub, CuC* d_L, CuC* d_R, CuC* d_C,
                   const std::vector<Complex>& L, const std::vector<Complex>& R, std::vector<Complex>& C){
  const int n=(int)Comp::N;
  const CuC one=make_cuDoubleComplex(1.0,0.0), zero=make_cuDoubleComplex(0.0,0.0);
  CUDA_CHECK(cudaMemcpy(d_L, reinterpret_cast<const CuC*>(L.data()), (size_t)n*n*sizeof(CuC), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_R, reinterpret_cast<const CuC*>(R.data()), (size_t)n*n*sizeof(CuC), cudaMemcpyHostToDevice));
  cublasZgemm(cub, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &one, d_R, n, d_L, n, &zero, d_C, n);
  C.resize((size_t)n*n);
  CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(C.data()), d_C, (size_t)n*n*sizeof(CuC), cudaMemcpyDeviceToHost));
}

// conn_shift2(X,Y,dt) = tr[ X . S_dt Y S_-dt ] for the free field (antiperiodic time-shift).  Two-matrix
// generalization of jj_exact_diag_deter_free's conn_shift (which uses X==Y for the vector).
static Complex conn_shift2(const std::vector<Complex>& X, const std::vector<Complex>& Y, int dt){
  const Idx N=Comp::N, Nx=Comp::Nx; const int Nt=Comp::Nt;
  std::vector<int> sgn(Nt); for(int t=0;t<Nt;t++) sgn[t]=twrap_sign(t,-dt);
  Complex s(0,0);
  for(Idx p=0;p<N;p++){
    const int sp=sgn[(int)(p/Nx)];
    for(Idx q=0;q<N;q++){
      const int sq=sgn[(int)(q/Nx)];
      s += X[(size_t)p*N+q] * (double)(sp*sq) * Y[(size_t)tshift(q,-dt)*N + tshift(p,-dt)];
    }
  }
  return s;
}

// Build the AXIAL dressed source/sink and contract: given the dense kernel KK (= K_ov_kappa or Sigma),
// the GW factor G=(1-D_ov), and the propagator P (+ its dagger), accumulate the n_t0-mapped C_{A+-}(dt).
//   B_src = KK . G^dag . P^dag ;  B_snk = KK^dag . G . P ;  C(dt) = conn_shift2(B_src,B_snk,dt).
// (KK^dag is formed by conj_transpose.)  Returns the per-dt connected vector in cdt, and tr[B_snk] in disc.
static void axial_contract(cublasHandle_t cub, CuC* d_L, CuC* d_R, CuC* d_C,
                           const std::vector<Complex>& KK, const std::vector<Complex>& G,
                           const std::vector<Complex>& Gdag, const std::vector<Complex>& P,
                           const std::vector<Complex>& Pdag,
                           std::vector<Complex>& cdt, Complex& disc){
  const int Nt=Comp::Nt;
  std::vector<Complex> tmp, B_src, KKdag, B_snk;
  matmul(cub, d_L, d_R, d_C, KK,  Gdag, tmp);      // tmp   = KK . G^dag
  matmul(cub, d_L, d_R, d_C, tmp, Pdag, B_src);    // B_src = (KK G^dag) . P^dag
  conj_transpose(KK, KKdag);
  matmul(cub, d_L, d_R, d_C, KKdag, G, tmp);       // tmp   = KK^dag . G
  matmul(cub, d_L, d_R, d_C, tmp,   P, B_snk);     // B_snk = (KK^dag G) . P
  disc = trace(B_snk);
  cdt.assign(Nt, Complex(0,0));
  #pragma omp parallel for schedule(dynamic)
  for(int dt=0;dt<Nt;dt++) cdt[dt]=conn_shift2(B_src,B_snk,dt);
}

struct Args{ double nu0=1.0,nu1=-1.0,mass_re=0.0,mass_im=0.0; std::string ens_dir,prop_file,out_tag; int n_t0=2,gpu=0,ins=0; bool sum=false; };
void PrintHelp(){ printf("jj_exact_axial_deter_free: FREE exact AXIAL C_{A+-}(t0->t), tp+sp (Eq.3.55).\n"
                         "  --ins <i>          single-insertion index (tp: site i; sp: link i).  default 0\n"
                         "  --sum              DIAGONALLY-n-SUMMED (Eq.4.29) over ALL sites(tp)/links(sp),\n"
                         "                     area-weighted + 1/4pi, K_ov_kappa per insertion, + ylm tower.\n"
                         "                     build-use-discard.  -> corr_deter_exactsum_axial_L<L>.\n"
                         "  --mass-re/--mass-im  selects the P dir + esnid\n"
                         "  --prop-file <path> read P from this exact file (e.g. cont_prop_L<L>/Dinv.0.h5)\n"
                         "  --out-tag <tag>    corr_deter_exact{1,sum}_axial_<tag>_L<L>;  --n-t0 --gpu\n"); }
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

// Build (or load from cache) the dense G = (1 - D_ov) (N op_oneMinusD applies).  Mass- and insertion-
// INDEPENDENT, so a single shared cache data_free_Gcache_L<L>/G.h5.
static void build_or_load_G(MatPoly& op_oneMinusD, std::vector<Complex>& G){
  const Idx N=Comp::N;
  const std::string gdir="data_free_Gcache_L"+std::to_string(Comp::N_REFINE)+"/";
  std::filesystem::create_directories(gdir);
  const std::string gfile=gdir+"G.h5";
  if(std::filesystem::exists(gfile)){
    try{ HighFive::File gf(gfile,HighFive::File::ReadOnly);
      if(gf.exist("complete")){ load_mat(gf,"G",G); std::cout<<"# G=(1-D_ov) cache HIT "<<gfile<<"\n"; return; } }catch(...){}
  }
  std::cout<<"# G=(1-D_ov) cache MISS "<<gfile<<" -> building dense (N op_oneMinusD applies) ...\n";
  Timer gt; std::vector<Complex> ej(N), out(N); G.assign((size_t)N*N, Complex(0,0));
  for(Idx j=0;j<N;j++){
    std::fill(ej.begin(), ej.end(), Complex(0,0)); ej[j]=Complex(1,0);
    op_oneMinusD.from_cpu<N>(out.data(), ej.data());       // out = (1 - D_ov) e_j = column j of G
    for(Idx i=0;i<N;i++) G[(size_t)i*N+j]=out[i];
    if(j%512==0) std::cout<<"#   G col "<<j<<"/"<<N<<"  ["<<gt.currentSeconds()<<" s]\n";
  }
  const std::string gtmp=gfile+".tmp";
  { HighFive::File gf(gtmp,HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
    gf.createDataSet("N",std::vector<int>{(int)N}); save_mat(gf,"G",G);
    gf.createDataSet("complete",std::vector<int>{1}); }
  std::filesystem::rename(gtmp,gfile);
  std::cout<<"# wrote G cache "<<gfile<<"  ["<<gt.currentSeconds()<<" s]\n";
}

int main(int argc,char* argv[]){
  std::cout<<std::scientific<<std::setprecision(15);
  Args a; ParseArgs(argc,argv,a); if(a.nu1<0.0) a.nu1=a.nu0;
  (void)a.gpu; CUDA_CHECK(cudaSetDevice(0));
  cudaDeviceProp prop; cudaGetDeviceProperties(&prop,0); std::cout<<"# dev = "<<prop.name<<"\n";
  constexpr Idx N=Comp::N; constexpr int Nt=Comp::Nt;
  const bool free_field=a.ens_dir.empty();
  if(!free_field){ std::cout<<"# ERROR: free-field only.\n"; return 1; }
  if(!a.sum) std::cout<<"# EXACT AXIAL single-insertion (ins="<<a.ins<<") C_{A+-}, tp+sp:  N="<<N<<"  [free]\n";

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

  // K operator + the (1 - D_ov) GW factor.  K is mass-independent; massless D_ov for the GW factor.
  Fermion D(DW, Complex(0.0), 21);
  ConservedCurrent<Fermion,Gauge> kop(D);
  auto f_D=std::bind(&Fermion::mult_deviceAsyncLaunch_ms, &D, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_D(f_D);
  MatPoly op_oneMinusD;
  op_oneMinusD.push_back(cplx( 1.0), {});       // identity
  op_oneMinusD.push_back(cplx(-1.0), {&M_D});   // - D_ov   => (1 - D_ov)
  D.update(U);                                  // set the overlap (free U=1) once: needed for op_K AND
                                                // op_oneMinusD (the G build) on every code path below.

  cublasHandle_t cub; cublasCreate(&cub);
  CuC *d_L,*d_R,*d_C;
  CUDA_CHECK(cudaMalloc(&d_L,(size_t)N*N*sizeof(CuC)));
  CUDA_CHECK(cudaMalloc(&d_R,(size_t)N*N*sizeof(CuC)));
  CUDA_CHECK(cudaMalloc(&d_C,(size_t)N*N*sizeof(CuC)));

  // ===================== DIAGONALLY-n-SUMMED mode (--sum), Eq.(4.29) =====================
  if(a.sum){
    std::cout<<"# EXACT AXIAL DIAGONAL-SUM (--sum): area-weighted Eq.(4.29), build-use-discard.  N="<<N
             <<"  n_sites="<<n_sites<<" n_links="<<n_links<<"\n";
    MatPoly op_K; op_K.push_back(cplx(1.0), {&kop});   // D already updated above

    std::vector<Complex> G, Gdag; build_or_load_G(op_oneMinusD, G); conj_transpose(G, Gdag);

    const std::string esnid="free_vmRe"+std::to_string(a.mass_re)+"vmIm"+std::to_string(a.mass_im);
    const std::string propdir="data_"+esnid+"/prop_deter_L"+std::to_string(Comp::N_REFINE)+"/";
    const std::string pfile=a.prop_file.empty()?propdir+"Dinv.0.h5":a.prop_file;
    if(!std::filesystem::exists(pfile)){ std::cout<<"# no propagator "<<pfile<<"\n"; return 1; }
    std::vector<Complex> P, Pdag;
    { HighFive::File f(pfile,HighFive::File::ReadOnly); load_mat(f,"Dm_inv",P); }
    conj_transpose(P, Pdag);
    std::cout<<"# loaded P "<<pfile<<"  (+ P^dag)\n";

    const std::string dname=a.out_tag.empty()?std::string("corr_deter_exactsum_axial")
                                             :std::string("corr_deter_exactsum_axial_")+a.out_tag;
    const std::string outdir="data_"+esnid+"/"+dname+"_L"+std::to_string(Comp::N_REFINE)+"/";
    std::filesystem::create_directories(outdir);
    const std::string h5path=outdir+"corr.0.h5";

    Timer timer;
    const std::string h5tmp=h5path+".tmp";
    auto h5p=std::make_unique<HighFive::File>(h5tmp,HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
    HighFive::File& h5=*h5p;
    h5.createDataSet("t0s",t0s); h5.createDataSet("n_t0",std::vector<int>{n_t0});
    h5.createDataSet("nhits",std::vector<int>{1}); h5.createDataSet("summed",std::vector<int>{1});

    std::vector<Complex> K((size_t)N*N), ej(N), out(N);
    const long total_cols=(long)(n_sites+n_links)*(long)N;
    long done_cols=0;
    std::cout<<"#   total op_K solves = (n_sites+n_links)*N = "<<total_cols<<"  (build-use-discard)\n";

    // axial ylm tower: accumulate Sigma_{l,m} = sum_n A_n Y_lm(n) K^t_ov_kappa,n during the tp site loop.
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
        std::vector<Complex> cn; Complex dn;
        axial_contract(cub, d_L, d_R, d_C, K, G, Gdag, P, Pdag, cn, dn);
        disc += w*dn;
        for(int dt=0;dt<Nt;dt++) cdt[dt]+=w*cn[dt];
        std::cout<<"#   "<<proj<<" ins "<<(ins+1)<<"/"<<n_ins<<" DONE  conn(dt=4)+="<<(w*cn[4]).real()
                 <<"  ("<<(int)(100.0*done_cols/total_cols)<<"% of solves)  ["<<timer.currentSeconds()<<" s]\n";
      }
      std::vector<std::vector<Complex>> Cpp(n_t0,std::vector<Complex>(Nt,Complex(0,0)));
      for(int b=0;b<n_t0;b++) for(int t=0;t<Nt;t++){ const int dt=((t-t0s[b])%Nt+Nt)%Nt; Cpp[b][dt]=cdt[dt]; }
      std::vector<Complex> discvec(Nt, disc);
      for(int b=0;b<n_t0;b++){ const std::string kp="h0/t0_"+std::to_string(b)+"/";
        write_corr_axial(h5,kp+proj+"/Apm",Cpp[b]); }
      write_vec(h5,"h0/disc/"+proj+"/J",discvec);
      std::cout<<"#   "<<proj<<" SUM done: conn(dt=4)="<<cdt[4].real()<<"  ["<<timer.currentSeconds()<<" s]\n";
    }

    // axial ylm tower (Eq.4.36): g_l(t) = (1/(2l+1)) sum_m conn_shift2 of the dressed Sigma_{l,m}.
    for(int l=0; l<=L_MAX_YLM; l++){
      std::vector<Complex> gl(Nt, Complex(0,0));
      for(int c=0;c<n_lm;c++){
        if(lm[c].first!=l) continue;
        std::vector<Complex> cn; Complex dn;
        axial_contract(cub, d_L, d_R, d_C, Sigma[c], G, Gdag, P, Pdag, cn, dn);
        for(int dt=0;dt<Nt;dt++) gl[dt]+=cn[dt];
      }
      const double inv2lp1=1.0/(2.0*l+1.0);
      for(int t=0;t<Nt;t++) gl[t]*=inv2lp1;
      std::vector<std::vector<Complex>> Cyl(n_t0,std::vector<Complex>(Nt,Complex(0,0)));
      for(int b=0;b<n_t0;b++) for(int t=0;t<Nt;t++){ const int dt=((t-t0s[b])%Nt+Nt)%Nt; Cyl[b][dt]=gl[dt]; }
      for(int b=0;b<n_t0;b++){ const std::string kp="h0/t0_"+std::to_string(b)+"/ylm/l"+std::to_string(l)+"/";
        write_corr_axial(h5,kp+"Apm",Cyl[b]); }
      std::cout<<"#   axial ylm l="<<l<<": conn(dt=4)="<<Cyl[0][4].real()<<"  ["<<timer.currentSeconds()<<" s]\n";
    }

    h5.createDataSet("complete",std::vector<int>{1});
    h5p.reset(); std::filesystem::rename(h5tmp,h5path);
    std::cout<<"# wrote "<<h5path<<"\n";
    cudaFree(d_L); cudaFree(d_R); cudaFree(d_C); cublasDestroy(cub);
    for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
    return 0;
  }
  // ======================================================================================

  // ---- single-insertion: K cache (RAW, kappa-in; shared with the vector exact1) ----
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
    MatPoly op_K; op_K.push_back(cplx(1.0), {&kop});   // D already updated above
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

  // K_ov_kappa = K / kappa (normalized current; same convention as the vector exact1)
  const double kappa_sp = kop.insertion_kappa(std::pair<int,BaseLink>{0, base.links[a.ins]});
  const double kappa_tp = kop.insertion_kappa(std::pair<int,Idx>{0, (Idx)a.ins});
  for(auto& z : K_sp) z /= kappa_sp;
  for(auto& z : K_tp) z /= kappa_tp;
  std::cout<<"# K_ov_kappa = K/kappa:  kappa_sp(link "<<a.ins<<")="<<kappa_sp
           <<"  kappa_t(site "<<a.ins<<")="<<kappa_tp<<"\n";

  std::vector<Complex> G, Gdag; build_or_load_G(op_oneMinusD, G); conj_transpose(G, Gdag);

  const std::string tag=std::string("free");
  const std::string esnid=tag+"_vmRe"+std::to_string(a.mass_re)+"vmIm"+std::to_string(a.mass_im);
  const std::string propdir="data_"+esnid+"/prop_deter_L"+std::to_string(Comp::N_REFINE)+"/";
  const std::string pfile=a.prop_file.empty()?propdir+"Dinv.0.h5":a.prop_file;
  if(!std::filesystem::exists(pfile)){ std::cout<<"# no propagator "<<pfile<<"\n"; return 1; }
  const std::string dname=a.out_tag.empty()?std::string("corr_deter_exact1_axial")
                                           :std::string("corr_deter_exact1_axial_")+a.out_tag;
  const std::string outdir="data_"+esnid+"/"+dname+"_L"+std::to_string(Comp::N_REFINE)+"/";
  std::filesystem::create_directories(outdir);
  const std::string h5path=outdir+"corr.0.h5";

  std::vector<Complex> P, Pdag; { HighFive::File f(pfile,HighFive::File::ReadOnly); load_mat(f,"Dm_inv",P); }
  conj_transpose(P, Pdag);
  std::cout<<"# loaded P "<<pfile<<"  (+ P^dag)\n";

  Timer timer;
  const std::string h5tmp=h5path+".tmp";
  auto h5p=std::make_unique<HighFive::File>(h5tmp,HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
  HighFive::File& h5=*h5p;
  h5.createDataSet("t0s",t0s); h5.createDataSet("n_t0",std::vector<int>{n_t0});
  h5.createDataSet("nhits",std::vector<int>{1}); h5.createDataSet("ins",std::vector<int>{a.ins});
  h5.createDataSet("kappa_sp",std::vector<double>{kappa_sp});
  h5.createDataSet("kappa_tp",std::vector<double>{kappa_tp});

  for(int which=0; which<2; which++){
    const std::string proj=(which==0)?"tp":"sp";
    const std::vector<Complex>& K=(which==0)?K_tp:K_sp;
    std::vector<Complex> cdt; Complex disc;
    axial_contract(cub, d_L, d_R, d_C, K, G, Gdag, P, Pdag, cdt, disc);
    std::vector<std::vector<Complex>> Cpp(n_t0,std::vector<Complex>(Nt,Complex(0,0)));
    std::vector<Complex> discvec(Nt, disc);
    for(int b=0;b<n_t0;b++) for(int t=0;t<Nt;t++){ const int dt=((t-t0s[b])%Nt+Nt)%Nt; Cpp[b][dt]=cdt[dt]; }
    for(int b=0;b<n_t0;b++){ const std::string kp="h0/t0_"+std::to_string(b)+"/";
      write_corr_axial(h5,kp+proj+"/Apm",Cpp[b]); }
    write_vec(h5,"h0/disc/"+proj+"/J",discvec);
    std::cout<<"#   "<<proj<<": disc(0)=("<<discvec[0].real()<<","<<discvec[0].imag()
             <<")  conn(dt=4)="<<Cpp[0][4].real()<<"  ["<<timer.currentSeconds()<<" s]\n";
  }
  h5.createDataSet("complete",std::vector<int>{1});
  h5p.reset(); std::filesystem::rename(h5tmp,h5path);
  std::cout<<"# wrote "<<h5path<<"\n";

  cudaFree(d_L); cudaFree(d_R); cudaFree(d_C); cublasDestroy(cub);
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
  return 0;
}
