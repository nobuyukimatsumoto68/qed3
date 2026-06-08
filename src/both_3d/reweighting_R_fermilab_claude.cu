// reweighting_R_fermilab_claude.cu
// -----------------------------------------------------------------------------
// FERMILAB build of reweighting_R_claude.cu: identical except the geometry paths are ABSOLUTE
// (the relative "../../geometry" did not resolve on the cluster).  Excluded from the local Makefile
// by the *_fermilab* filter; build/run on Fermilab.  Edit the GEOM root below if the repo path differs.
// -----------------------------------------------------------------------------
// Parity reweighting factor R (PDF qed3int_v2-10.pdf, Eq. 2.5), per gauge config:
//
//   R = conj( prod_i [ ((1-mP) lam_i + mP) / (lam_i + mP) ] ),
//
// where {lam_i} are the eigenvalues of the MASSLESS overlap D_ov (computed by dense
// diagonalization, as the PDF prescribes), and mP = Complex(mass_re, mass_im) is the
// (parity) mass from the CLI.  R -> 1 as mP -> 0.  Stage 3 of the analysis pipeline
// (jj_analysis_global_procedure_claude.md); multiplied into the correlators in the notebook.
//
// Method (mimics eig_wmass_claude.cu): build the dense N x N matrix of D_ov by applying it to
// each unit vector (one column per apply), then all eigenvalues via cuSOLVER cusolverDnXgeev.
// The eigenvalues are mP-INDEPENDENT, so the diagonalization runs ONCE per config; R for the CLI
// mP (and any other mP later) is a trivial product from the stored spectrum.
//
// Config source (mirrors jj_corr_claude.cu): --ens-dir loops ckpoint_lat.k (k=0,ninter,2*ninter,...);
// OMIT --ens-dir => free field (single config k=1, U=1; R->1, a check).  Output goes into the SAME
// obs dir as the correlators: data_<ESNID>/R/R.<k>.h5  (ESNID = ens + '_vmRe<re>vmIm<im>').
//
// Plan: reweighting_R_impl_plan_claude.md.   Build: nvcc (same flags as the other targets).
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

  constexpr int N_REFINE=1;            // match the correlator ensemble (L=1)
  constexpr int NS=2;
  constexpr int Nt=128;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;

  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;       // inner Zolotarev pole solves (for the D_ov apply)
  const double TOL_OUTER=1.0e-5;
}

// FERMILAB: absolute geometry path (edit this root if the repo lives elsewhere).
const std::string dir = "/project/qed3/qed3/geometry/data/";

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

#include "/project/qed3/qed3/geometry/geodesic.h"   // FERMILAB: absolute (matches `dir` above)

#include "sparse_matrix.h"
#include "dirac_base.h"
#include "dirac_simp.h"
#include "dirac_dual.h"
#include "dirac_ext.h"

#include "sparse_dirac.h"
#include "matpoly_claude.h"

#include "dirac_pf.h"
#include "overlap_wmass_claude.h"        // complex-mass overlap (massless at mass=0)

using BaseLink = std::array<Idx,2>;
using BaseFace = std::vector<Idx>;

// ---- CLI ----
void PrintHelp(){
  printf("Usage: ./a.out [options]   (reweighting factor R, Eq. 2.5, per config)\n");
  printf("  --mass-re <x>        parity mass Re (default: 0.0)\n");
  printf("  --mass-im <y>        parity mass Im (default: 0.0); pure parity = (0, y)\n");
  printf("  --ens-dir <path>     sea config directory; OMIT => free field (U=1) check\n");
  printf("  --ninter <N>         ensemble config stride: ckpoint_lat.k for k=0,N,2N,... (default: 10)\n");
  printf("  --nu0 <x>            (default 1.0)   --nu1 <y>  (default: nu0)\n");
  printf("  --gpu <id>           CUDA device (default 0)\n");
  printf("  -h, --help\n");
}

void ParseArgs(int argc, char** argv, double& nu0, double& nu1,
               double& mass_re, double& mass_im, std::string& ens_dir, int& ninter, int& gpu){
  static struct option long_opts[] = {
    {"nu0",     required_argument, nullptr, 'n'},
    {"nu1",     required_argument, nullptr, 'm'},
    {"mass-re", required_argument, nullptr, 'r'},
    {"mass-im", required_argument, nullptr, 'i'},
    {"ens-dir", required_argument, nullptr, 'e'},
    {"ninter",  required_argument, nullptr, 'I'},
    {"gpu",     required_argument, nullptr, 'G'},
    {"help",    no_argument,       nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
  };
  int opt, idx;
  while((opt = getopt_long(argc, argv, "n:m:r:i:e:I:G:h", long_opts, &idx)) != -1){
    switch(opt){
      case 'n': nu0     = std::stod(optarg); break;
      case 'm': nu1     = std::stod(optarg); break;
      case 'r': mass_re = std::stod(optarg); break;
      case 'i': mass_im = std::stod(optarg); break;
      case 'e': ens_dir = optarg; break;
      case 'I': ninter  = std::stoi(optarg); break;
      case 'G': gpu     = std::stoi(optarg); break;
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
  ParseArgs(argc, argv, nu0, nu1, mass_re, mass_im, ens_dir, ninter, gpu);
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

  // ---- output dir: data_<ESNID>/R/  (same obs dir as the correlators) ----
  std::string ens_base = ens_dir;
  if(!ens_base.empty() && ens_base.back()=='/') ens_base.pop_back();
  { const auto slash = ens_base.find_last_of('/'); if(slash!=std::string::npos) ens_base = ens_base.substr(slash+1); }
  const std::string esnid = (free_field ? std::string("free") : ens_base)
                          + "_vmRe"+std::to_string(mass_re)+"vmIm"+std::to_string(mass_im);
  const std::string dir_out = "data_"+esnid+"/R/";
  std::filesystem::create_directories(dir_out);
  std::cout << "# dir_out = " << dir_out << std::endl;

  // ---- cuSOLVER geev resources (n fixed = N; allocate once, reuse per config) ----
  const int n = (int)N, lda = (int)N, ldvl = (int)N, ldvr = (int)N;
  cusolverDnHandle_t handle = NULL; cudaStream_t stream = NULL; cusolverDnParams_t params = NULL;
  cusolverEigMode_t jobvl = CUSOLVER_EIG_MODE_NOVECTOR, jobvr = CUSOLVER_EIG_MODE_NOVECTOR;
  CUSOLVER_CHECK(cusolverDnCreate(&handle));
  CUDA_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
  CUSOLVER_CHECK(cusolverDnSetStream(handle, stream));
  CUSOLVER_CHECK(cusolverDnCreateParams(&params));

  CuC *A = (CuC*)malloc((size_t)n*n*CD);
  CuC *W = (CuC*)malloc((size_t)n*CD);
  CuC *d_A=nullptr, *d_W=nullptr, *d_VL=nullptr, *d_VR=nullptr; int *d_info=nullptr;
  CUDA_CHECK(cudaMalloc(&d_A,  CD*(size_t)n*n));
  CUDA_CHECK(cudaMalloc(&d_W,  CD*(size_t)n));
  CUDA_CHECK(cudaMalloc(&d_VL, CD*(size_t)n*n));
  CUDA_CHECK(cudaMalloc(&d_VR, CD*(size_t)n*n));
  CUDA_CHECK(cudaMalloc(&d_info, sizeof(int)));

  size_t lwork_d = 0, lwork_h = 0;
  CUSOLVER_CHECK( cusolverDnXgeev_bufferSize(handle, params, jobvl, jobvr, n,
                    CUDA_C_64F, d_A, lda, CUDA_C_64F, d_W,
                    CUDA_C_64F, d_VL, ldvl, CUDA_C_64F, d_VR, ldvr,
                    CUDA_C_64F, &lwork_d, &lwork_h) );
  void *d_work=nullptr, *h_work=nullptr;
  CUDA_CHECK(cudaMalloc(&d_work, lwork_d));
  h_work = malloc(lwork_h);

  // ---- helper: write a length-m complex vector as key/real + key/imag ----
  auto write_cvec = [](HighFive::File& h5, const std::string& key, const std::vector<Complex>& C){
    std::vector<double> re(C.size()), im(C.size());
    for(size_t k=0;k<C.size();k++){ re[k]=C[k].real(); im[k]=C[k].imag(); }
    h5.createDataSet(key+"/real", re);  h5.createDataSet(key+"/imag", im);
  };

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
    const std::string h5path = dir_out + "R." + std::to_string(k) + ".h5";
    if(std::filesystem::exists(h5path)){
      bool complete=false;
      try { HighFive::File h5c(h5path, HighFive::File::ReadOnly); complete = h5c.exist("complete"); } catch(...) {}
      if(complete){ std::cout << "# skip k="<<k<<" (complete)" << std::endl; continue; }
    }
    Dov.update(U);
    std::cout << "# k="<<k<<(free_field?" (free field)":"")
              << "  lambda_min/max="<<Dov.lambda_min<<"/"<<Dov.lambda_max
              << "  building dense D_ov ("<<N<<" applies) ..." << std::endl;

    // build dense matrix column-by-column: column i = D_ov e_i
    Eigen::MatrixXcd mat(N, N);
    std::vector<Complex> e(N), De(N);
    for(Idx i=0; i<N; i++){
      std::fill(e.begin(), e.end(), Complex(0.0,0.0));
      e[i] = Complex(1.0,0.0);
      Op.from_cpu<N>( De.data(), e.data() );
      for(Idx j=0;j<N;j++) mat(j,i) = De[j];
      if((i % 256)==0) std::cout << "#   col "<<i<<"/"<<N<<" ["<<timer.currentSeconds()<<" s]" << std::endl;
    }

    // diagonalize (all eigenvalues): A column-major = mat
    for(int j=0;j<n;j++) for(int i=0;i<n;i++) A[(size_t)n*j+i] = cplx(mat(i,j));
    CUDA_CHECK(cudaMemcpy(d_A, A, CD*(size_t)n*n, H2D));
    CUDA_CHECK(cudaMemset(d_W, 0, CD*(size_t)n));
    int info=0;
    CUSOLVER_CHECK( cusolverDnXgeev(handle, params, jobvl, jobvr, n,
                      CUDA_C_64F, d_A, lda, CUDA_C_64F, d_W,
                      CUDA_C_64F, d_VL, ldvl, CUDA_C_64F, d_VR, ldvr,
                      CUDA_C_64F, d_work, lwork_d, h_work, lwork_h, d_info) );
    CUDA_CHECK(cudaMemcpy(W, d_W, CD*(size_t)n, D2H));
    CUDA_CHECK(cudaMemcpy(&info, d_info, sizeof(int), D2H));
    assert(info==0);

    // R = conj( prod_i ((1-mP) lam_i + mP)/(lam_i + mP) ), accumulated in log space
    std::vector<Complex> lam(N);
    Complex logR(0.0,0.0);
    const Complex one(1.0,0.0);
    for(Idx i=0;i<N;i++){
      lam[i] = Complex(real(W[i]), imag(W[i]));
      const Complex num = (one - mP)*lam[i] + mP;
      const Complex den = lam[i] + mP;
      logR += std::log(num) - std::log(den);
    }
    const Complex Rval = std::conj( std::exp(logR) );
    std::cout << "# k="<<k<<"  R = ("<<Rval.real()<<", "<<Rval.imag()<<")  |R|="<<std::abs(Rval)
              << "  ["<<timer.currentSeconds()<<" s]" << std::endl;

    // write data_<ESNID>/R/R.<k>.h5
    HighFive::File h5(h5path, HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
    write_cvec(h5, "lam", lam);                                   // eigenvalues of D_ov (mP-independent)
    write_cvec(h5, "R",   std::vector<Complex>{Rval});            // the reweighting factor
    write_cvec(h5, "mP",  std::vector<Complex>{mP});              // mass label
    h5.createDataSet("k", std::vector<int>{k});
    h5.createDataSet("complete", std::vector<int>{1});           // sentinel, written LAST
    std::cout << "#   wrote " << h5path << std::endl;
  }

  // cleanup
  free(A); free(W); free(h_work);
  CUDA_CHECK(cudaFree(d_A)); CUDA_CHECK(cudaFree(d_W)); CUDA_CHECK(cudaFree(d_VL));
  CUDA_CHECK(cudaFree(d_VR)); CUDA_CHECK(cudaFree(d_info)); CUDA_CHECK(cudaFree(d_work));
  CUSOLVER_CHECK(cusolverDnDestroyParams(params));
  CUSOLVER_CHECK(cusolverDnDestroy(handle));
  CUDA_CHECK(cudaStreamDestroy(stream));
  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();
  return 0;
}
