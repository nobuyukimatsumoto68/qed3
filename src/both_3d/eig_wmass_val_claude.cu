// eig_wmass_val_claude.cu
// _claude: DENSE reference for validating the IRL (eig_lanczos_claude.cu).  Builds A = (D_ov+m)^\dagger(D_ov+m)
// column-by-column (bind DHD_deviceAsyncLaunch, NOT mult -- so we compare the SAME Hermitian A the IRL sees,
// with NO normality assumption) and diagonalizes it with cusolverDnXgeev.  A is Hermitian PSD so its geev
// eigenvalues are real; we sort ascending and write the smallest -- to be diffed 1:1 against the IRL evals.
// Cloned from saved_scripts_claude/eig_wmass_claude.cu (which geev'd D_m itself via mult).  Free-field
// (cold gauge) at matched L/Nt via -DLREF/-DNT_CLI; identical M5/r/T/at/mass to eig_lanczos_claude.cu.
// Usage: ./bin [mass_re] [mass_im]      Output: eig_wmass_L<L>_nt<Nt>_claude.dat  (i  eval, ascending).

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <filesystem>
#include <cstdint>
#include <complex>
#include <array>
#include <vector>
#include <map>
#include <string>
#include <cmath>
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
  constexpr int NPARALLEL=NPARALLEL_DUPDATE;      // match the known-good bench Comp (uses _claude includes)
  constexpr int NSTREAMS=4;
  constexpr int NPARALLEL_GAUGE=NPARALLEL_DUPDATE;
  constexpr int NPARALLEL_SORT=NPARALLEL_DUPDATE;

#ifndef LREF
#define LREF 1
#endif
#ifndef NT_CLI
#define NT_CLI 32
#endif
  constexpr int N_REFINE=LREF;
  constexpr int NS=2;
  constexpr int Nt=NT_CLI;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;
  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-8;
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

#include "overlap_wmass_claude.h"
#include "frozen_window_claude.h"       // FIXED Zolotarev window (lmin,lmax) per L -- same source of truth as HMC

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(12);
  std::clog << std::scientific << std::setprecision(12);

  int device;
  CUDA_CHECK(cudaGetDeviceCount(&device));
  cudaDeviceProp dp;
  cudaGetDeviceProperties(&dp, 0);
  std::cout << "# dev = " << dp.name << std::endl;
  CUDA_CHECK(cudaSetDevice(0));

  constexpr Idx N = Comp::N;

  double mass_re = 0.0;
  double mass_im = 0.0;
  if(argc>1) mass_re = atof(argv[1]);
  if(argc>2) mass_im = atof(argv[2]);
  const Complex mass = Complex(mass_re, mass_im);
  std::cout << "# L(N_REFINE)="<<Comp::N_REFINE<<" Nt="<<Comp::Nt<<" N="<<N<<" mass="<<mass<<std::endl;

  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;
  using Fermion=OverlapWMass<WilsonDirac>;

  Base base(Comp::N_REFINE);
  const double r = 1.0;
  const double M5 = -1.0;
  const double at = 0.2;         // MUST match the ensemble (dir "at0.200000") and production (jj/hmc at=0.2)
  const double nu0 = 1.0;
  WilsonDirac DW(base, 0.0, r, M5, at, nu0);

  Gauge U(base);                 // cold (free field) by default; CONFIG_LAT env overrides with a real config
  const char* cfg = std::getenv("CONFIG_LAT");
  if(cfg != nullptr && std::string(cfg).size() > 0){
    assert(std::filesystem::exists(cfg) && "CONFIG_LAT path not found");
    U.read(cfg);
    std::cout << "# read gauge config: " << cfg << std::endl;
  }
  else{
    std::cout << "# free field (cold gauge)" << std::endl;
  }
  Fermion Dov(DW, mass, 21);
  double lmin_frozen, lmax_frozen;
  frozen_window(Comp::N_REFINE, lmin_frozen, lmax_frozen);   // FIXED Zolotarev window per L (matches HMC)
  Dov.set_lambda(lmin_frozen, lmax_frozen);
  Dov.update(U);
  std::cout << "# Dov (FROZEN window): lambda_min/max = " << Dov.lambda_min << " / " << Dov.lambda_max
            << "  Delta = " << Dov.Delta() << std::endl;

  // Bind the HERMITIAN normal operator A = (D+m)^dag (D+m) (NOT the bare D_m) so the dense spectrum is
  // directly the IRL target -- no normality assumption in the comparison.
  MatPoly Op;
#ifdef SPECTRUM_DW
  // Bare Wilson D_W (the overlap's public CSR).  D_ov is a UNITARY PROJECTOR of D_W (D_ov = 1 + V, V = sign/
  // unitary part), so it keeps only the PHASE of each D_W eigenvalue -> plot D_W's complex spectrum to see it.
  Op.push_back( cplx(1.0), {&Dov.M_DW} );
#elif defined(SPECTRUM_DOV)
  // Build the DENSE overlap operator D_ov ITSELF -> complex spectrum on the GW circle |z-1|=1.
  auto f_Op = std::bind(&Fermion::mult_deviceAsyncLaunch, &Dov, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_Op( f_Op );
  Op.push_back( cplx(1.0), {&M_Op} );
#else
  auto f_Op = std::bind(&Fermion::DHD_deviceAsyncLaunch, &Dov, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_Op( f_Op );
  Op.push_back( cplx(1.0), {&M_Op} );
#endif

  Eigen::MatrixXcd mat(N, N);
  {
    Timer tb;
    for(Idx i=0; i<N; i++){
      Eigen::VectorXcd e = Eigen::VectorXcd::Zero(N);
      e(i) = 1.0;
      std::vector<Complex> xi(e.data(), e.data()+N);
      std::vector<Complex> Dxi(N);
      Op.from_cpu<N>( Dxi.data(), xi.data() );
      mat.block(0,i,N,1) = Eigen::Map<Eigen::MatrixXcd>(Dxi.data(), N, 1);
      if(i%64==0) std::cout << "# col " << i << "/" << N << "  (" << tb.currentSeconds() << " s)" << std::endl;
    }
  }

  // =========================================
  // cusolver geev (A is Hermitian -> eigenvalues real up to roundoff).
  cusolverDnHandle_t handle = NULL;
  cudaStream_t stream = NULL;
  cusolverDnParams_t params = NULL;

  const int n = mat.cols();
  const int lda = n;

  CuC *A, *W;
  A = (CuC*)malloc((size_t)n*n*CD);
  W = (CuC*)malloc((size_t)n*CD);
  for(int j=0; j<n; j++) for(int i=0; i<n; i++) A[(size_t)n*j+i] = cplx(mat(i,j));
  for(int i=0; i<n; i++) W[i] = cplx(0.);

  CuC *d_A, *d_W, *d_VL, *d_VR;
  cusolverEigMode_t jobvl = CUSOLVER_EIG_MODE_NOVECTOR;
  cusolverEigMode_t jobvr = CUSOLVER_EIG_MODE_NOVECTOR;
  int ldvl = n;
  int ldvr = n;
  int info = 0;
  int *d_info = nullptr;

  size_t workspaceInBytesOnDevice = 0;
  void *d_work = nullptr;
  size_t workspaceInBytesOnHost = 0;
  void *h_work = nullptr;

  CUSOLVER_CHECK(cusolverDnCreate(&handle));
  CUDA_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
  CUSOLVER_CHECK(cusolverDnSetStream(handle, stream));
  CUSOLVER_CHECK(cusolverDnCreateParams(&params));

  CUDA_CHECK(cudaMalloc( &d_A, CD * (size_t)n*n ));
  CUDA_CHECK(cudaMalloc( &d_W, CD * (size_t)n ));
  CUDA_CHECK(cudaMalloc( &d_VL, CD * (size_t)n*n ));
  CUDA_CHECK(cudaMalloc( &d_VR, CD * (size_t)n*n ));
  CUDA_CHECK(cudaMalloc( &d_info, sizeof(int)));

  CUDA_CHECK( cudaMemcpy(d_A, A, CD*(size_t)n*n, H2D) );
  CUDA_CHECK( cudaMemset(d_W, 0, CD * (size_t)n) );
  CUDA_CHECK( cudaMemset(d_VL, 0, CD * (size_t)n*n) );
  CUDA_CHECK( cudaMemset(d_VR, 0, CD * (size_t)n*n) );

  CUSOLVER_CHECK( cusolverDnXgeev_bufferSize( handle, params, jobvl, jobvr, n,
                                              CUDA_C_64F, d_A, lda,
                                              CUDA_C_64F, d_W,
                                              CUDA_C_64F, d_VL, ldvl,
                                              CUDA_C_64F, d_VR, ldvr,
                                              CUDA_C_64F,
                                              &workspaceInBytesOnDevice, &workspaceInBytesOnHost) );

  CUDA_CHECK(cudaMalloc( &d_work, workspaceInBytesOnDevice ));
  h_work = malloc(workspaceInBytesOnHost);

  CUSOLVER_CHECK( cusolverDnXgeev( handle, params, jobvl, jobvr, n,
                                   CUDA_C_64F, d_A, lda,
                                   CUDA_C_64F, d_W,
                                   CUDA_C_64F, d_VL, ldvl,
                                   CUDA_C_64F, d_VR, ldvr,
                                   CUDA_C_64F,
                                   d_work, workspaceInBytesOnDevice,
                                   h_work, workspaceInBytesOnHost, d_info) );

  CUDA_CHECK(cudaMemcpy( W, d_W, CD*(size_t)n, D2H) );
  CUDA_CHECK(cudaMemcpy( &info, d_info, sizeof(int), D2H ));
  std::cout << "# geev info (0=success) = " << info << std::endl;
  assert( info==0 );

#if defined(SPECTRUM_DOV) || defined(SPECTRUM_DW)
  // COMPLEX spectrum: write Re, Im, |z| for every eigenvalue (unsorted geev order).  Plot on the complex
  // plane.  D_ov -> GW circle |z-1|=1; D_W -> the Wilson eigenvalue cloud that D_ov projects onto that circle.
  const char* out_tag_env = std::getenv("OUT_TAG");                 // e.g. "_gsq2" to keep runs separate
  const std::string out_tag = (out_tag_env != nullptr) ? std::string(out_tag_env) : std::string("");
#ifdef SPECTRUM_DW
  const std::string opname = "dw";
#else
  const std::string opname = "dov";
#endif
  const std::string outname = "eig_"+opname+"_spectrum_L"+std::to_string(Comp::N_REFINE)+out_tag
                            + "_nt"+std::to_string(Comp::Nt)+"_claude.dat";
  std::ofstream out(outname);
  out << std::scientific << std::setprecision(14);
  out << "# DENSE geev COMPLEX spectrum of "<<opname<<"   L="<<Comp::N_REFINE<<" Nt="<<Comp::Nt
      << " mass="<<mass<<"  (frozen window)  tag="<<out_tag << std::endl;
  out << "# i   Re(z)   Im(z)   |z|" << std::endl;
  for(int i=0;i<n;i++){
    const double re = cuCreal(W[i]);
    const double im = cuCimag(W[i]);
    out << i << "   " << re << "   " << im << "   " << std::sqrt(re*re + im*im) << std::endl;
  }
  out.close();
  std::cout << "# wrote " << outname << " (" << n << " complex eigenvalues of " << opname << ")" << std::endl;
#else
  // Real eigenvalues of A, sorted ascending; report max |imag| as a Hermiticity check.
  std::vector<double> ev(n);
  double max_imag = 0.0;
  for(int i=0;i<n;i++){
    ev[i] = cuCreal(W[i]);
    max_imag = std::max(max_imag, std::fabs(cuCimag(W[i])));
  }
  std::sort(ev.begin(), ev.end());
  std::cout << "# max|imag(eval)| = " << max_imag << "  (should be ~roundoff; A Hermitian)" << std::endl;

  const std::string outname = "eig_wmass_L"+std::to_string(Comp::N_REFINE)
                            + "_nt"+std::to_string(Comp::Nt)+"_claude.dat";
  std::ofstream out(outname);
  out << std::scientific << std::setprecision(14);
  out << "# DENSE geev spectrum of A=(D_ov+m)^dag(D_ov+m)   L="<<Comp::N_REFINE<<" Nt="<<Comp::Nt
      << " mass="<<mass<<"  max|imag|="<<max_imag<<std::endl;
  out << "# i   eval(A, ascending)   sqrt(eval)=singular_value" << std::endl;
  for(int i=0;i<n;i++){
    out << i << "   " << ev[i] << "   " << std::sqrt(std::max(0.0, ev[i])) << std::endl;
  }
  out.close();
  std::cout << "# wrote " << outname << " (" << n << " evals)" << std::endl;
#endif

  free(A);
  free(W);
  free(h_work);
  CUDA_CHECK(cudaFree(d_A));
  CUDA_CHECK(cudaFree(d_W));
  CUDA_CHECK(cudaFree(d_VL));
  CUDA_CHECK(cudaFree(d_VR));
  CUDA_CHECK(cudaFree(d_info));
  CUDA_CHECK(cudaFree(d_work));

  CUSOLVER_CHECK(cusolverDnDestroyParams(params));
  CUSOLVER_CHECK(cusolverDnDestroy(handle));
  CUDA_CHECK(cudaStreamDestroy(stream));

  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
  return 0;
}
