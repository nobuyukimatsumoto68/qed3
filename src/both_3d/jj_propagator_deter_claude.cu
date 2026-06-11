// jj_propagator_deter_claude.cu
// -----------------------------------------------------------------------------
// EXACT all-to-all propagator builder (free field or per config), STAGE 1 of the exact
// current-correlator validation (jj_exact_freefield_impl_plan_claude.md).  Sibling of
// reweighting_R_claude.cu (left untouched): both build the dense massless overlap D_ov by applying
// it column-by-column to the N unit vectors (the dominant cost, N applies); here we then LU-invert
// (cuSOLVER getrf/getrs) the mass's operator(s) and SAVE the dense inverse(s) -- the all-to-all
// propagator -- for the contraction stage (jj_contract_deter_claude.cu).  No stochastic estimator.
//
// ONE mass per run (CLI --mass-re/--mass-im).  From the SINGLE dense D_ov:
//   massless (m=0):    save Dm_inv = D_ov^{-1}.
//   m_F (real m):      save Dm_inv = (D_ov + m)^{-1}.
//   m_P (imag m):      save Dm_inv = (D_ov + m)^{-1} AND Dtil_inv = (D_ov + m/(1-m))^{-1}.
// Since D_m = D_ov + c*I, every mass is a cheap shift + LU on the one dense D_ov.
// --with-R also runs cusolverDnXgeev on the same matrix and writes the Eq. 2.5 reweighting factor
// (the dense build is then shared with R -- the user's "combine with the R code", as a new sibling).
//
// L (= N_REFINE) is COMPILE-TIME (N is templated throughout): set -DN_REFINE_CLI=<L> (default 1);
// a run script builds one binary per L.  N = 2 * (10 L^2 + 2) * Nt.
//
// Config source (mirrors jj_corr): --ens-dir loops ckpoint_lat.k (k=0,ninter,...); OMIT => free
// field (single config k=1, U=1; the conformal test).  Output: data_<ESNID>/prop_deter_L<L>/
// Dinv.<k>.h5  (atomic write).   Build: nvcc (same flags as the other targets).
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
#include <ctime>
#include <getopt.h>
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

#ifndef N_REFINE_CLI
#define N_REFINE_CLI 1                 // L; override at compile time: nvcc -DN_REFINE_CLI=2 ...
#endif
  constexpr int N_REFINE=N_REFINE_CLI; // compile-time L (run script builds one binary per L)
  constexpr int NS=2;
#ifndef NT_CLI
#define NT_CLI 128                     // temporal extent; override: nvcc -DNT_CLI=16 ... (small-Nt checks)
#endif
  constexpr int Nt=NT_CLI;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;

  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;       // inner Zolotarev pole solves (for the D_ov apply)
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

#include "../../geometry/geodesic.h"

#include "sparse_matrix.h"
#include "dirac_base.h"
#include "dirac_simp.h"
#include "dirac_dual.h"
#include "dirac_ext.h"

// #include "sparse_dirac.h"
#include "sparse_dirac_claude.h"   // O(len) bucketing CSR build (was O(N*len)); -DCSR_VERIFY to check
#include "matpoly_claude.h"

#include "dirac_pf.h"
#include "overlap_wmass_claude.h"        // complex-mass overlap (massless at mass=0)

using BaseLink = std::array<Idx,2>;
using BaseFace = std::vector<Idx>;


// ---- write a length-m complex vector as key/{real,imag} ----
static void write_cvec(HighFive::File& h5, const std::string& key, const std::vector<Complex>& C){
  std::vector<double> re(C.size()), im(C.size());
  for(size_t k=0; k<C.size(); k++){ re[k] = C[k].real(); im[k] = C[k].imag(); }
  h5.createDataSet(key + "/real", re);
  h5.createDataSet(key + "/imag", im);
}

// ---- LU-solve (D_ov + c) X = I -> X = (D_ov+c)^{-1}; return ROW-MAJOR re/imag (length N^2). ----
// A0 = dense D_ov (column-major host); B_h = host inverse staging; I_h = identity RHS; d_* device.
static std::pair<std::vector<double>,std::vector<double> >
invert_shift(const Complex& c, int n, int lda, const CuC* A0, CuC* B_h, const std::vector<CuC>& I_h,
             CuC* d_A, CuC* d_B, int* d_ipiv, int* d_info, CuC* d_lu_work,
             cusolverDnHandle_t handle, cudaStream_t stream){
  std::vector<CuC> M(A0, A0 + (size_t)n*n);                          // copy dense D_ov
  for(int i=0; i<n; i++){
    CuC& d = M[(size_t)n*i + i];
    d = make_cuDoubleComplex(cuCreal(d) + c.real(), cuCimag(d) + c.imag());   // M = D_ov + c I
  }
  CUDA_CHECK(cudaMemcpy(d_A, M.data(),   CD*(size_t)n*n, H2D));
  CUDA_CHECK(cudaMemcpy(d_B, I_h.data(), CD*(size_t)n*n, H2D));
  int info = 0;
  CUSOLVER_CHECK( cusolverDnZgetrf(handle, n, n, d_A, lda, d_lu_work, d_ipiv, d_info) );
  CUDA_CHECK(cudaStreamSynchronize(stream));   // cusolver runs on `stream` (non-blocking); the
  CUDA_CHECK(cudaMemcpy(&info, d_info, sizeof(int), D2H)); assert(info==0);   // default-stream D2H
  CUSOLVER_CHECK( cusolverDnZgetrs(handle, CUBLAS_OP_N, n, n, d_A, lda, d_ipiv, d_B, n, d_info) );
  CUDA_CHECK(cudaStreamSynchronize(stream));   // does NOT wait for it -> MUST sync before reading.
  CUDA_CHECK(cudaMemcpy(&info, d_info, sizeof(int), D2H)); assert(info==0);
  CUDA_CHECK(cudaMemcpy(B_h, d_B, CD*(size_t)n*n, D2H));   // B_h = M^{-1}, column-major (col*n+row)
  std::vector<double> re((size_t)n*n), im((size_t)n*n);
  for(int col=0; col<n; col++)
    for(int row=0; row<n; row++){
      const CuC v = B_h[(size_t)n*col + row];
      re[(size_t)n*row + col] = cuCreal(v);   // -> row-major
      im[(size_t)n*row + col] = cuCimag(v);
    }
  return {std::move(re), std::move(im)};
}


// ---- CLI ----
void PrintHelp(){
  printf("Usage: ./a.out [options]   (exact all-to-all propagator D_m^{-1} via dense build + LU)\n");
  printf("  --mass-re <x>        parity mass Re (default: 0.0)\n");
  printf("  --mass-im <y>        parity mass Im (default: 0.0); pure parity = (0, y)\n");
  printf("  --ens-dir <path>     sea config directory; OMIT => free field (U=1) check\n");
  printf("  --ninter <N>         ensemble config stride: ckpoint_lat.k for k=0,N,2N,... (default: 10)\n");
  printf("  --nu0 <x>            (default 1.0)   --nu1 <y>  (default: nu0)\n");
  printf("  --gpu <id>           CUDA device (default 0)\n");
  printf("  --with-R             also diagonalize (geev) and write the Eq. 2.5 R (shares the build)\n");
  printf("  -h, --help\n");
}

void ParseArgs(int argc, char** argv, double& nu0, double& nu1,
               double& mass_re, double& mass_im, std::string& ens_dir, int& ninter, int& gpu,
               bool& with_R){
  static struct option long_opts[] = {
    {"nu0",     required_argument, nullptr, 'n'},
    {"nu1",     required_argument, nullptr, 'm'},
    {"mass-re", required_argument, nullptr, 'r'},
    {"mass-im", required_argument, nullptr, 'i'},
    {"ens-dir", required_argument, nullptr, 'e'},
    {"ninter",  required_argument, nullptr, 'I'},
    {"gpu",     required_argument, nullptr, 'G'},
    {"with-R",  no_argument,       nullptr, 'R'},
    {"help",    no_argument,       nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
  };
  int opt, idx;
  while((opt = getopt_long(argc, argv, "n:m:r:i:e:I:G:Rh", long_opts, &idx)) != -1){
    switch(opt){
      case 'n': nu0     = std::stod(optarg); break;
      case 'm': nu1     = std::stod(optarg); break;
      case 'r': mass_re = std::stod(optarg); break;
      case 'i': mass_im = std::stod(optarg); break;
      case 'e': ens_dir = optarg; break;
      case 'I': ninter  = std::stoi(optarg); break;
      case 'G': gpu     = std::stoi(optarg); break;
      case 'R': with_R  = true; break;
      case 'h': default: PrintHelp(); std::exit(0);
    }
  }
}

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  double nu0=1.0, nu1=-1.0;          // nu1<0 => default to nu0
  double mass_re=0.0, mass_im=0.0;
  std::string ens_dir="";            // empty => free-field mode
  int ninter=10, gpu=0;
  bool with_R=false;
  ParseArgs(argc, argv, nu0, nu1, mass_re, mass_im, ens_dir, ninter, gpu, with_R);
  if(nu1 < 0.0) nu1 = nu0;

  // GPU selection is via CUDA_VISIBLE_DEVICES (set by the run script), as in jj_corr/eig_wmass;
  // always use the first VISIBLE device.  (--gpu is the run-script env knob, NOT a device index here:
  // under CUDA_VISIBLE_DEVICES=<id> the chosen GPU is remapped to visible device 0.)
  (void)gpu;
  CUDA_CHECK(cudaSetDevice(0));
  cudaDeviceProp prop; cudaGetDeviceProperties(&prop, 0);
  std::cout << "# dev = " << prop.name << "  (visible device 0; CUDA_VISIBLE_DEVICES picks the GPU)" << std::endl;

  constexpr Idx N  = Comp::N;
  const Complex mP(mass_re, mass_im);
  const bool free_field = ens_dir.empty();
  std::cout << "# reweighting R: mP=("<<mass_re<<","<<mass_im<<")  N="<<N
            << "  ens_dir="<<(free_field?std::string("<free-field U=1>"):ens_dir) << std::endl;

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;
  using Fermion=OverlapWMass<WilsonDirac>;

  Base base(Comp::N_REFINE);
  std::cout << "# lattice set." << std::endl;
  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();

  const double r  = 1.0;
  const double M5 = -1.0;
  const double at = 0.2;
  WilsonDirac DW(base, 0.0, r, M5, at, nu1);

  Gauge U(base);                                 // free field (U=1) unless a config is read below

  // massless overlap D_ov (mass=0); its eigenvalues are the lam_i of Eq. 2.5.  npole=21.
  Fermion Dov(DW, Complex(0.0), 21);
  std::cout << "# D_ov set (massless; M5="<<M5<<")." << std::endl;

  // operator applied column-by-column to build the dense matrix: y = D_ov x
  // (use the multishift apply; the plain one is the reference: mult_deviceAsyncLaunch)
  // auto f_Op = std::bind(&Fermion::mult_deviceAsyncLaunch,    &Dov, std::placeholders::_1, std::placeholders::_2);
  auto f_Op = std::bind(&Fermion::mult_deviceAsyncLaunch_ms, &Dov, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_Op( f_Op );
  MatPoly Op;
  Op.push_back( cplx(1.0), {&M_Op} );

  // ---- output dir: data_<ESNID>/prop_deter_L<L>/  ----
  std::string ens_base = ens_dir;
  if(!ens_base.empty() && ens_base.back()=='/') ens_base.pop_back();
  { const auto slash = ens_base.find_last_of('/'); if(slash!=std::string::npos) ens_base = ens_base.substr(slash+1); }
  const std::string esnid = (free_field ? std::string("free") : ens_base)
                          + "_vmRe"+std::to_string(mass_re)+"vmIm"+std::to_string(mass_im);
  const std::string dir_out = "data_"+esnid+"/prop_deter_L"+std::to_string(Comp::N_REFINE)+"/";
  std::filesystem::create_directories(dir_out);
  std::cout << "# dir_out = " << dir_out << std::endl;

  const bool parity = (mass_im != 0.0);   // m_P => also need tilde D^{-1}
  const Complex mtil = parity ? (mP/(Complex(1.0,0.0) - mP)) : Complex(0.0,0.0);  // Eq. 3.64 shift

  // ---- cuSOLVER resources (n = N fixed; allocate once) ----
  const int n = (int)N, lda = (int)N;
  cusolverDnHandle_t handle = NULL; cudaStream_t stream = NULL; cusolverDnParams_t params = NULL;
  CUSOLVER_CHECK(cusolverDnCreate(&handle));
  CUDA_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
  CUSOLVER_CHECK(cusolverDnSetStream(handle, stream));
  CUSOLVER_CHECK(cusolverDnCreateParams(&params));

  // dense D_ov (column-major host A0), LU factor buffer d_A, RHS/solution d_B (becomes the inverse).
  CuC *A0    = (CuC*)malloc((size_t)n*n*CD);   // dense D_ov, column-major (j*n+i)=(row i,col j)
  CuC *B_h   = (CuC*)malloc((size_t)n*n*CD);   // host staging for the inverse (column-major)
  CuC *d_A=nullptr, *d_B=nullptr; int *d_ipiv=nullptr, *d_info=nullptr;
  CUDA_CHECK(cudaMalloc(&d_A,    CD*(size_t)n*n));
  CUDA_CHECK(cudaMalloc(&d_B,    CD*(size_t)n*n));
  CUDA_CHECK(cudaMalloc(&d_ipiv, sizeof(int)*(size_t)n));
  CUDA_CHECK(cudaMalloc(&d_info, sizeof(int)));
  int lwork_lu = 0;
  CUSOLVER_CHECK( cusolverDnZgetrf_bufferSize(handle, n, n, d_A, lda, &lwork_lu) );
  CuC *d_lu_work=nullptr; CUDA_CHECK(cudaMalloc(&d_lu_work, CD*(size_t)lwork_lu));

  // geev resources only when --with-R (eigenvalues for Eq. 2.5)
  CuC *W=nullptr, *d_W=nullptr, *d_VL=nullptr, *d_VR=nullptr; void *d_gwork=nullptr, *h_gwork=nullptr;
  size_t lwork_gd=0, lwork_gh=0;
  cusolverEigMode_t jobvl = CUSOLVER_EIG_MODE_NOVECTOR, jobvr = CUSOLVER_EIG_MODE_NOVECTOR;
  if(with_R){
    W = (CuC*)malloc((size_t)n*CD);
    CUDA_CHECK(cudaMalloc(&d_W,  CD*(size_t)n));
    CUDA_CHECK(cudaMalloc(&d_VL, CD*(size_t)n*n));
    CUDA_CHECK(cudaMalloc(&d_VR, CD*(size_t)n*n));
    CUSOLVER_CHECK( cusolverDnXgeev_bufferSize(handle, params, jobvl, jobvr, n,
                      CUDA_C_64F, d_A, lda, CUDA_C_64F, d_W, CUDA_C_64F, d_VL, n,
                      CUDA_C_64F, d_VR, n, CUDA_C_64F, &lwork_gd, &lwork_gh) );
    CUDA_CHECK(cudaMalloc(&d_gwork, lwork_gd));  h_gwork = malloc(lwork_gh);
  }

  // identity RHS (column-major == row-major identity), uploaded before each getrs.
  std::vector<CuC> I_h((size_t)n*n, make_cuDoubleComplex(0.0,0.0));
  for(int i=0;i<n;i++) I_h[(size_t)n*i+i] = make_cuDoubleComplex(1.0,0.0);

  // ---- config loop (mirrors jj_corr) ----
  const int k_ckpoint = free_field ? 1 : ninter;
  const int kmax      = free_field ? 0 : 1000000;
  Timer timer;
  for(int k=0; k<=kmax; k+=k_ckpoint){
    if(!free_field){
      const std::string str_lat = ens_dir + "ckpoint_lat." + std::to_string(k);
      if(!std::filesystem::exists(str_lat)){ if(k==0) continue; else break; }
      U.read(str_lat);
    }
    const std::string h5path = dir_out + "Dinv." + std::to_string(k) + ".h5";
    if(std::filesystem::exists(h5path)){
      bool complete=false;
      try { HighFive::File h5c(h5path, HighFive::File::ReadOnly); complete = h5c.exist("complete"); } catch(...) {}
      if(complete){ std::cout << "# skip k="<<k<<" (complete)" << std::endl; continue; }
    }
    Dov.update(U);
    std::cout << "# k="<<k<<(free_field?" (free field)":"")
              << "  lambda_min/max="<<Dov.lambda_min<<"/"<<Dov.lambda_max
              << "  building dense D_ov ("<<N<<" applies) ..." << std::endl;

    // build dense D_ov column-by-column: column i = D_ov e_i  (column-major A0[j*n+i] = (D_ov)_{j,i})
    {
      std::vector<Complex> e(N), De(N);
      for(Idx i=0; i<N; i++){
        std::fill(e.begin(), e.end(), Complex(0.0,0.0));
        e[i] = Complex(1.0,0.0);
        Op.from_cpu<N>( De.data(), e.data() );
        for(Idx j=0;j<N;j++) A0[(size_t)n*i+j] = cplx(De[j]);   // col i = De
        if((i % 256)==0) std::cout << "#   col "<<i<<"/"<<N<<" ["<<timer.currentSeconds()<<" s]" << std::endl;
      }
    }

    // optional R (Eq. 2.5) from the SAME dense matrix (geev destroys its input -> upload a copy)
    std::vector<Complex> lam; Complex Rval(1.0,0.0);
    if(with_R){
      CUDA_CHECK(cudaMemcpy(d_A, A0, CD*(size_t)n*n, H2D));
      CUDA_CHECK(cudaMemset(d_W, 0, CD*(size_t)n));
      int info=0;
      CUSOLVER_CHECK( cusolverDnXgeev(handle, params, jobvl, jobvr, n,
                        CUDA_C_64F, d_A, lda, CUDA_C_64F, d_W, CUDA_C_64F, d_VL, n,
                        CUDA_C_64F, d_VR, n, CUDA_C_64F, d_gwork, lwork_gd, h_gwork, lwork_gh, d_info) );
      CUDA_CHECK(cudaMemcpy(W, d_W, CD*(size_t)n, D2H));
      CUDA_CHECK(cudaMemcpy(&info, d_info, sizeof(int), D2H)); assert(info==0);
      lam.resize(N); Complex logR(0.0,0.0); const Complex one(1.0,0.0);
      for(Idx i=0;i<N;i++){ lam[i]=Complex(real(W[i]),imag(W[i]));
        logR += std::log((one-mP)*lam[i]+mP) - std::log(lam[i]+mP); }
      Rval = std::conj(std::exp(logR));
      std::cout << "#   R = ("<<Rval.real()<<","<<Rval.imag()<<")  |R|="<<std::abs(Rval)<<std::endl;
    }

    // LU-invert the mass's operator(s): Dm = D_ov + mP ; (parity) tilde = D_ov + mP/(1-mP)
    std::cout << "#   LU-inverting Dm" << (parity?" + Dtil":"") << " ..." << std::endl;
    auto Dm = invert_shift(mP, n, lda, A0, B_h, I_h, d_A, d_B, d_ipiv, d_info, d_lu_work, handle, stream);
    if(free_field){   // trivial self-check: || D_m . D_m^{-1} - I ||_F  (d_B holds D_m^{-1} from invert_shift)
      cublasHandle_t cub; cublasCreate(&cub);
      CuC* d_C=nullptr; CUDA_CHECK(cudaMalloc(&d_C, CD*(size_t)n*n));
      std::vector<CuC> Dm_mat(A0, A0+(size_t)n*n);                       // D_m = D_ov + mP (col-major)
      for(int i=0;i<n;i++){ CuC& d=Dm_mat[(size_t)n*i+i];
        d=make_cuDoubleComplex(cuCreal(d)+mP.real(), cuCimag(d)+mP.imag()); }
      CUDA_CHECK(cudaMemcpy(d_A, Dm_mat.data(), CD*(size_t)n*n, H2D));   // d_A = D_m ; d_B = D_m^{-1}
      const CuC c_one=make_cuDoubleComplex(1.0,0.0), c_zero=make_cuDoubleComplex(0.0,0.0);
      cublasZgemm(cub, CUBLAS_OP_N, CUBLAS_OP_N, n,n,n, &c_one, d_A,n, d_B,n, &c_zero, d_C,n);
      std::vector<CuC> C((size_t)n*n); CUDA_CHECK(cudaMemcpy(C.data(), d_C, CD*(size_t)n*n, D2H));
      double res=0.0;
      for(int col=0;col<n;col++) for(int row=0;row<n;row++){
        const CuC v=C[(size_t)n*col+row];
        const double dr=cuCreal(v)-(row==col?1.0:0.0), di=cuCimag(v); res+=dr*dr+di*di;
      }
      std::cout << "#   [check] || D_m . D_m^{-1} - I ||_F = " << std::sqrt(res) << std::endl;
      CUDA_CHECK(cudaFree(d_C)); cublasDestroy(cub);
    }
    std::pair<std::vector<double>,std::vector<double>> Dt;
    if(parity) Dt = invert_shift(mtil, n, lda, A0, B_h, I_h, d_A, d_B, d_ipiv, d_info, d_lu_work, handle, stream);
    std::cout << "#   inverted ["<<timer.currentSeconds()<<" s]" << std::endl;

    // ---- write data_<ESNID>/prop_deter_L<L>/Dinv.<k>.h5  (ATOMIC: .tmp + rename) ----
    // layout: Dm_inv/{real,imag} are ROW-MAJOR length-N^2: index = row*N + col, value=(D_m^{-1})_{row,col}.
    const std::string h5tmp = h5path + ".tmp";
    auto h5p = std::make_unique<HighFive::File>(h5tmp,
                 HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
    HighFive::File& h5 = *h5p;
    h5.createDataSet("N",        std::vector<int>{n});
    h5.createDataSet("N_REFINE", std::vector<int>{Comp::N_REFINE});
    h5.createDataSet("Nt",       std::vector<int>{Comp::Nt});
    h5.createDataSet("n_sites",  std::vector<int>{(int)Comp::N_SITES});
    h5.createDataSet("parity",   std::vector<int>{parity?1:0});
    write_cvec(h5, "mP", std::vector<Complex>{mP});
    h5.createDataSet("Dm_inv/real", Dm.first);   h5.createDataSet("Dm_inv/imag", Dm.second);
    if(parity){ h5.createDataSet("Dtil_inv/real", Dt.first); h5.createDataSet("Dtil_inv/imag", Dt.second); }
    {   // ALWAYS dump the dense D_ov (row-major Dov[r*N+c]=(D_ov)_{r,c}; A0 is col-major).  Needed by the
        // commutator-tp current jj_commut_tp (K^{tp}=[D_ov,Theta]).  NOTE disk: doubles the file (Dm_inv +
        // Dov, both N^2) -> L1 0.3 GB, L2 3.7 GB, L4 55 GB.
      std::vector<double> re((size_t)n*n), im((size_t)n*n);
      for(int r=0;r<n;r++) for(int c=0;c<n;c++){ const CuC v=A0[(size_t)n*c+r];
        re[(size_t)r*n+c]=cuCreal(v); im[(size_t)r*n+c]=cuCimag(v); }
      h5.createDataSet("Dov/real", re); h5.createDataSet("Dov/imag", im);
      std::cout << "#   dumped dense D_ov" << std::endl;
    }
    if(with_R){ write_cvec(h5,"lam",lam); write_cvec(h5,"R",std::vector<Complex>{Rval}); }
    h5.createDataSet("k", std::vector<int>{k});
    h5.createDataSet("complete", std::vector<int>{1});           // sentinel, written LAST
    h5p.reset();
    std::filesystem::rename(h5tmp, h5path);
    std::cout << "#   wrote " << h5path << std::endl;
  }

  // cleanup
  free(A0); free(B_h); if(W) free(W); if(h_gwork) free(h_gwork);
  CUDA_CHECK(cudaFree(d_A)); CUDA_CHECK(cudaFree(d_B)); CUDA_CHECK(cudaFree(d_ipiv));
  CUDA_CHECK(cudaFree(d_info)); CUDA_CHECK(cudaFree(d_lu_work));
  if(d_W){ CUDA_CHECK(cudaFree(d_W)); CUDA_CHECK(cudaFree(d_VL)); CUDA_CHECK(cudaFree(d_VR)); CUDA_CHECK(cudaFree(d_gwork)); }
  CUSOLVER_CHECK(cusolverDnDestroyParams(params));
  CUSOLVER_CHECK(cusolverDnDestroy(handle));
  CUDA_CHECK(cudaStreamDestroy(stream));
  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();
  return 0;
}
