#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <algorithm>

#include <cstdint>
#include <complex>
using Idx = std::int32_t;
using Complex = std::complex<double>;

namespace Comp{
  constexpr int NPARALLEL=12;
  constexpr int NSTREAMS=4;

  constexpr int N_REFINE=16;
  constexpr int NS=2;
  constexpr Idx N_SITES=20*N_REFINE*N_REFINE;
  constexpr Idx N=NS*N_SITES; // matrix size of DW

  // const double TOL=1.0e-9;
  const double TOL_INNER=1.0e-10;
  const double TOL_OUTER=1.0e-9;
}

// #define IsVerbose
#define IsVerbose2
// #define InfoForce
#define InfoDelta

#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cublas_api.h>
#include <cusolverDn.h>
using CuC = cuDoubleComplex;
#include "gpu_header.h"

// ======================================

// #include "s2n.h"
#include "s2.h"
#include "rng.h"
// #include "gauge.h"
// #include "action.h"

#include "sparse_matrix.h"
// #include "pseudofermion.h"
// #include "dirac.h"
#include "dirac_simp.h"

#include "sparse_dirac.h"
#include "matpoly.h"

#include "overlap.h"

// #include "hmc.h"
// #include "dirac_s2_dual.h"
// #include "header_cusolver.hpp"


// TODO: Cusparse for SparseMatrix::act_gpu, probably defining handle in matpoly.h
// all the operation on GPU in Overlap::operator()
// gradient of Dov (Overlap class, in parallel to Dirac)
// pseudofermion
// 3d
// __m256 to vectorize with AVX2


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

  using Gauge=U1onS2<false>;
  // using Force=U1onS2<false>;
  // using Action=U1Wilson;
  // using Fermion=Dirac1fonS2;
  // using HMC=HMC<Force,Gauge,Action,Fermion>;
  // using Rng=ParallelRng;
  using Lattice=S2Trivalent;
  using Rng=ParallelRng<Lattice>;

  Lattice lattice(Comp::N_REFINE);
  // Dirac1fonS2 D(lattice, 0.0, 1.0);

  using WilsonDirac=Dirac1fonS2;
  // using Overlap=OverlapPseudoFermion;

  Gauge U(lattice);
  Rng rng(lattice);
  // U.gaussian( rng, 0.2 );

  const double M5 = -1.8;
  // const double M5 = 0.0;
  // const double M5 = -2.5;
  WilsonDirac DW(lattice, M5, 1.0/3.0);
  // Overlap Dov(DW);
  // Overlap Dov(DW, 1.0e-4, 21);
  Overlap Dov(DW, 31);
  // Dov.compute(U);
  Dov.update(U);
  std::cout << "# min max ratio: "
            << Dov.lambda_min << " "
            << Dov.lambda_max << " "
            << Dov.lambda_min/Dov.lambda_max << std::endl;
  std::cout << "# delta = " << Dov.Delta() << std::endl;


  // MatPoly Op;
  // Op.push_back ( cplx(1.0), {&(Dov.M_DW), &(Dov.M_DWH)} );
  // auto f_Op = std::bind(&Overlap::sq_device, &Dov, std::placeholders::_1, std::placeholders::_2);
  // auto f_Op = std::bind(&Overlap::sq_device, &Dov, std::placeholders::_1, std::placeholders::_2);
  // auto f_Op = std::bind(&Overlap::mult_device, &Dov, std::placeholders::_1, std::placeholders::_2);
  auto f_Op = std::bind(&Overlap::mult_deviceAsyncLaunch, &Dov, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_Op( f_Op );

  MatPoly Op;
  Op.push_back ( cplx(1.0), {&M_Op} );
  // Op.push_back ( cplx(1.0), {&Dov.M_DW} );

  constexpr Idx N = Comp::N;
  Eigen::MatrixXcd mat(N, N);
  {
    for(Idx i=0; i<N; i++){
      Eigen::VectorXcd e = Eigen::VectorXcd::Zero(N);
      e(i) = 1.0;
      std::vector<Complex> xi(e.data(), e.data()+N);
      std::vector<Complex> Dxi(N);

      //   // Op.solve<N>( d_eta, d_xi );
      Op.from_cpu<N>( Dxi, xi );

      // for(Idx j=0; j<N; j++) Dxi[j] -= xi[j];
      // for(Idx j=0; j<N; j++) Dxi[j] -= M5*xi[j];
      // std::cout << "debug. i=" << i << std::endl;
      // Op.from_cpu<N>( Dxi, xi );
      mat.block(0,i,N,1) = Eigen::Map<Eigen::MatrixXcd>(Dxi.data(), N, 1);
      std::clog << "i = " << i << " finished." << std::endl;
    }
  }

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
  // for(int j=0; j<n; j++) for(int i=0; i<n; i++) A[n*j+i] = reinterpret_cast<CuC*>(&mat(i,j));
  // for(int j=0; j<n; j++) for(int i=0; i<n; i++) A[n*j+i] = cplxmat(i,j));
  for(int i=0; i<n; i++) W[i] = cplx(0.);

  CuC *d_A, *d_W, *d_VL, *d_VR;
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
  CUDA_CHECK(cudaMalloc( &d_info, sizeof(int)));

  CUDA_CHECK( cudaMemcpy(d_A, A, CD*n*n, H2D) );

  // step 3: query working space of syevd
  cusolverEigMode_t jobvl = CUSOLVER_EIG_MODE_NOVECTOR;
  cusolverEigMode_t jobvr = CUSOLVER_EIG_MODE_NOVECTOR;
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

  std::cout << "# info (0=success) = " << info << std::endl;
  assert( info==0 );

  // std::vector<double> res(n);
  // for(int i=0; i<n; i++) res[i] = real(W[i]);
  // std::sort(res.begin(), res.end());
  // for(int i=0; i<n; i++) std::cout << i << " "
  // 				   << res[i] << " "
  // 				   << Dov.sgn(res[i]) << std::endl;

  for(int i=0; i<n; i++) std::cout << real(W[i]) << " " << imag(W[i]) << " " << abs(W[i]) << std::endl;

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

  // CUDA_CHECK(cudaDeviceReset());



  // // 2.4.5.7. cusolverDnXgeev()
  // cusolverStatus_t
  //   cusolverDnXgeev_bufferSize(
  // 			       cusolverDnHandle_t handle,
  // 			       cusolverDnParams_t params,
  // 			       cusolverEigMode_t jobvl,
  // 			       cusolverEigMode_t jobvr,
  // 			       int64_t n,
  // 			       cudaDataType dataTypeA,
  // 			       const void *A,
  // 			       int64_t lda,
  // 			       cudaDataType dataTypeW,
  // 			       const void *W,
  // 			       cudaDataType dataTypeVL,
  // 			       const void *VL,
  // 			       int64_t ldvl,
  // 			       cudaDataType dataTypeVR,
  // 			       const void *VR,
  // 			       int64_t ldvr,
  // 			       cudaDataType computeType,
  // 			       size_t *workspaceInBytesOnDevice,
  // 			       size_t *workspaceInBytesOnHost);

  // // ss. 2.5.2.5. cusolverSp<t>csreigvsi()
  // cusolverSpZcsreigvsi(cusolverSpHandle_t handle,
  // 		       int m,
  // 		       int nnz,
  // 		       const cusparseMatDescr_t descrA,
  // 		       const cuDoubleCuC *csrValA,
  // 		       const int *csrRowPtrA,
  // 		       const int *csrColIndA,
  // 		       cuDoubleCuC mu0,
  // 		       const cuDoubleCuC *x0,
  // 		       int maxite,
  // 		       double tol,
  // 		       cuDoubleCuC *mu,
  // 		       cuDoubleCuC *x);

  // Eigen::CuCEigenSolver<Eigen::MatrixXcd> solver( mat );
  // const Eigen::MatrixXcd evec = solver.eigenvectors();
  // Eigen::VectorXcd ev = solver.eigenvalues();
  // for(int i=0; i<evec.rows(); i++){
  //   const Eigen::VectorXcd check1 = sq * evec.col(i);
  //   const Eigen::VectorXcd check2 = eval[i] * evec.col(i);
  //   assert( (check1-check2).norm() < 1.0e-8 );

  //   const Eigen::VectorXcd MV = mat * evec.col(i);
  //   std::cout << ( MV.array() / evec.col(i).array() - 1.0).abs().maxCoeff() << std::endl;
  // }

  // auto ev = mat.eigenvalues();
  // for(int i=0; i<ev.size(); i++){
  //   std::cout << ev[i].real() << " " << ev[i].imag() << std::endl;
  // }

  // ----------------------------------

    // return 0;
}



  // for(int ix=0; ix<lattice.n_sites; ix++){
  //   for(int jj=0; jj<lattice.sites[ix].nn; jj++){
  //     const int iy = lattice.sites[ix].neighbors[jj];
  //     auto mat1 = ( D.sigma[0] - D.gamma(ix, iy) ) * D.Omega(ix, iy);
  //     auto mat2 = D.Omega(ix, iy) * ( D.sigma[0] - D.gamma(iy, ix, M_PI) );
  //     std::cout << mat1-mat2 << std::endl;
  //   }}

  // ----------------------------------
