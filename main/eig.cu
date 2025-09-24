#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <cstdint>
#include <complex>
#include <array>
#include <vector>
#include <map>
#include <Eigen/Dense>

using Double = double;
using Idx = std::int32_t;
using Complex = std::complex<double>;

using Link = std::array<Idx,2>; // <int,int>;
using Face = std::vector<Idx>;


using MS=Eigen::Matrix2cd;
using VD=Eigen::Vector2d;
using VE=Eigen::Vector3d;
using VC=Eigen::VectorXcd;

static constexpr int NS = 2;
static constexpr int DIM = 2;
static constexpr Complex I = Complex(0.0, 1.0);

// #define GAUGE_TRSF

// #define IS_OVERLAP

// #define IsVerbose
// #define IsVerbose2
// #define InfoDelta

namespace Comp{
  constexpr bool is_compact=false;

  // d_DW.update() is always done independently
#ifdef IS_OVERLAP
  constexpr int NPARALLEL_DUPDATE=1;
  constexpr int NPARALLEL=12; // 12
  constexpr int NSTREAMS=4; // 4
#else
  constexpr int NPARALLEL_DUPDATE=12;
  constexpr int NPARALLEL=1; // 12
  constexpr int NSTREAMS=12; // for grad loop
#endif
  constexpr int NPARALLEL_GAUGE=12; // 12
  constexpr int NPARALLEL_SORT=12; // 12

  constexpr int N_REFINE=4;
  constexpr int NS=2;

  // constexpr int Nt=16;
  constexpr int Nt=24;
  // constexpr int Nt=48;
  // constexpr int Nt=96;
  // constexpr int Nt=2;
  // constexpr int Nt=1;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE; // 30, 120, 480

  constexpr Idx Nx=NS*N_SITES; // matrix size of DW
  constexpr Idx N=Nx*Nt; // matrix size of DW

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-8;
}

const std::string dir = "../dats/";


#include "timer.h"

#include "s2n_simp.h"
#include "rng.h"
#include "valence.h"
#include "gauge_ext.h"
#include "action_ext.h"

#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cublas_api.h>
#include <cusolverDn.h>
using CuC = cuDoubleComplex;
#include "gpu_header.h"

// ======================================

#include "../connection/geodesic.h"

#include "sparse_matrix.h"
#include "dirac_base.h"
#include "dirac_simp.h"
#include "dirac_ext.h"

#include "sparse_dirac.h"
#include "matpoly.h"

#include "dirac_pf.h"

#include "overlap.h"


using BaseLink = std::array<Idx,2>; // <int,int>;
using BaseFace = std::vector<Idx>;



int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  int device;
  CUDA_CHECK(cudaGetDeviceCount(&device));
  cudaDeviceProp device_prop[device];
  cudaGetDeviceProperties(&device_prop[0], 0);
  std::cout << "# dev = " << device_prop[0].name << std::endl;
  CUDA_CHECK(cudaSetDevice(0));// "TITAN V"
  std::cout << "# (GPU device is set.)" << std::endl;

  // ---------------------------------------
  constexpr Idx N = Comp::N;
  constexpr int Nt = Comp::Nt;

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;

  using Force=GaugeExt<Base,Nt,Comp::is_compact>;
  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  using Action=U1WilsonExt<Base>;

  using Rng=ParallelRngExt<Base,Nt>;
  using Overlap=Overlap<WilsonDirac>;
  using Fermion=DiracPf<WilsonDirac>;

  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;

  // ----------------------


  // U.random_gauge_trsf(rng);



#ifdef IS_OVERLAP
  const double r = 1.0;
  const double M5 = -1.0;
#else // if not overlap
  const double r = 1.0;
  const double M5 = 0.0;
#endif
  // const double at = base.mean_ell * 1.0;
  const double T = 4.0;
  // const double T = 16.0;
  const double at = T/Comp::Nt;
  // const double at = base.mean_ell * 3.0/4.0;

  assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);

  WilsonDirac DW(base, 0.0, r, M5, at);

  std::cout << "# alat = " << base.mean_ell << std::endl;

  Gauge U(base);
  // srand( time(NULL) );
  // Rng rng(base, rand());
  // U.gaussian( rng, 0.2 );
  std::cout << "# U set" << std::endl;

  Fermion D(DW);
  D.update( U );

  std::cout << "# D updated" << std::endl;

  COO<N> gmfourth;
  DW.volume_matrix( gmfourth.en, -0.5 );
  gmfourth.do_it();


  MatPoly Op;
#ifdef IS_OVERLAP
  Overlap Dov(DW, 51);
  Dov.update(U);
  std::cout << "# Dov set; M5 = " << M5 << std::endl;
  std::cout << "# min max ratio: "
            << Dov.lambda_min << " "
            << Dov.lambda_max << " "
            << Dov.lambda_min/Dov.lambda_max << std::endl;
  std::cout << "# delta = " << Dov.Delta() << std::endl;

  auto f_Op = std::bind(&Overlap::mult_deviceAsyncLaunch, &Dov, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_Op( f_Op );
  Op.push_back ( cplx(1.0), {&gmfourth, &M_Op, &gmfourth} );
#else
  auto f_Op = std::bind(&Fermion::mult_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_Op( f_Op );
  Op.push_back ( cplx(1.0), {&gmfourth, &M_Op, &gmfourth} );
#endif

  Eigen::MatrixXcd mat(N, N);
  {
    for(Idx i=0; i<N; i++){
      Eigen::VectorXcd e = Eigen::VectorXcd::Zero(N);
      e(i) = 1.0;
      std::vector<Complex> xi(e.data(), e.data()+N);
      std::vector<Complex> Dxi(N);

      Op.from_cpu<N>( Dxi, xi );

      mat.block(0,i,N,1) = Eigen::Map<Eigen::MatrixXcd>(Dxi.data(), N, 1);
      std::cout << "# i = " << i << " finished." << std::endl;
    }
  }


  // {
  //   // Eigen::IOFormat fmt(Eigen::FullPrecision, 0, ", ", ",\n", "{", "}", "{", "}");
  //   // Eigen::IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " << ", ";");
  //   std::clog << mat.real() << std::endl;
  //   std::clog << mat.imag() << std::endl;
  //   // std::clog << mat.real().format(fmt) << std::endl;
  //   // std::clog << mat.imag().format(fmt) << std::endl;
  //   return 0;
  // }



  // =========================================
  // cusolver
  cusolverDnHandle_t handle = NULL;
  cudaStream_t stream = NULL;
  cusolverDnParams_t params = NULL;

  const int n = mat.cols(); // Number of rows (or columns) of matrix A.
  const int lda = n;

  CuC *A, *W;
  A = (CuC*)malloc(n*n*CD);
  W = (CuC*)malloc(n*CD);
  for(int j=0; j<n; j++) for(int i=0; i<n; i++) A[n*j+i] = cplx(mat(i,j));
  for(int i=0; i<n; i++) W[i] = cplx(0.);

  CuC *d_A, *d_W, *d_VL, *d_VR;

  cusolverEigMode_t jobvl = CUSOLVER_EIG_MODE_NOVECTOR;
  cusolverEigMode_t jobvr = CUSOLVER_EIG_MODE_NOVECTOR;
  int ldvl = n;
  int ldvr = n;
  //
  int info = 0;
  int *d_info = nullptr;

  size_t workspaceInBytesOnDevice = 0; /* size of workspace */
  void *d_work = nullptr;              /* device workspace */
  size_t workspaceInBytesOnHost = 0;   /* size of workspace */
  void *h_work = nullptr;              /* host workspace for */

  /* step 1: create cusolver handle, bind a stream */
  CUSOLVER_CHECK(cusolverDnCreate(&handle));
  CUDA_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
  CUSOLVER_CHECK(cusolverDnSetStream(handle, stream));
  CUSOLVER_CHECK(cusolverDnCreateParams(&params));

  CUDA_CHECK(cudaMalloc( &d_A, CD * n*n ));
  CUDA_CHECK(cudaMalloc( &d_W, CD * n ));
  CUDA_CHECK(cudaMalloc( &d_VL, CD * n*n ));
  CUDA_CHECK(cudaMalloc( &d_VR, CD * n*n ));
  // CUDA_CHECK(cudaMalloc( &d_VL, CD * 0 ));
  // CUDA_CHECK(cudaMalloc( &d_VR, CD * 0 ));
  CUDA_CHECK(cudaMalloc( &d_info, sizeof(int)));

  CUDA_CHECK( cudaMemcpy(d_A, A, CD*n*n, H2D) );
  CUDA_CHECK( cudaMemset(d_W, 0, CD * n) );
  CUDA_CHECK( cudaMemset(d_VL, 0, CD * n*n) );
  CUDA_CHECK( cudaMemset(d_VR, 0, CD * n*n) );

  // step 3: query working space of syevd
  // cusolverEigMode_t jobvl = CUSOLVER_EIG_MODE_NOVECTOR;
  // cusolverEigMode_t jobvr = CUSOLVER_EIG_MODE_VECTOR;
  cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;

  CUSOLVER_CHECK( cusolverDnXgeev_bufferSize( handle,
        				      params,
        				      jobvl,
        				      jobvr,
        				      n,
        				      CUDA_C_64F,
        				      d_A, // device
        				      lda,
        				      CUDA_C_64F,
        				      d_W, // Array holding the computed eigenvalues of A
        				      CUDA_C_64F,
        				      d_VL,
        				      ldvl,
        				      CUDA_C_64F,
        				      d_VR,
        				      ldvr,
        				      CUDA_C_64F,
        				      &workspaceInBytesOnDevice,
        				      &workspaceInBytesOnHost)
        	  );

  CUDA_CHECK(cudaMalloc( &d_work, workspaceInBytesOnDevice ) );
  h_work = malloc(workspaceInBytesOnHost);

  // step 4: compute spectrum
  CUSOLVER_CHECK( cusolverDnXgeev( handle,
        			   params,
        			   jobvl,
        			   jobvr,
        			   n,
        			   CUDA_C_64F,
        			   d_A,
        			   lda,
        			   CUDA_C_64F,
        			   d_W,
        			   CUDA_C_64F,
        			   d_VL,
        			   ldvl,
        			   CUDA_C_64F,
        			   d_VR,
        			   ldvr,
        			   CUDA_C_64F,
        			   d_work, // void *bufferOnDevice,
        			   workspaceInBytesOnDevice,
        			   h_work, // void *bufferOnHost,
        			   workspaceInBytesOnHost,
        			   d_info)
        	  );

  // ---------------------------------------------

  CUDA_CHECK(cudaMemcpy( W, d_W, CD*n, D2H) );
  CUDA_CHECK(cudaMemcpy( &info, d_info, sizeof(int), D2H ));

  std::vector<std::complex<double>> vr(n*n);
  for(Idx i=0; i<N; i++) gmfourth(d_VL+i*N, d_VR+i*N);
  CUDA_CHECK(cudaMemcpy( reinterpret_cast<CuC*>(vr.data()), d_VL, CD * n*n, D2H ));

  std::cout << "# info (0=success) = " << info << std::endl;
  assert( info==0 );

  for(int i=0; i<n; i++) std::clog << i << " " << real(W[i]) << " " << imag(W[i]) << " " << abs(W[i]) << std::endl;

  /* free resources */
  free(A);
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


  return 0; // EXIT_SUCCESS;

}

