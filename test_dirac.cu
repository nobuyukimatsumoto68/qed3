#include <iostream>
#include <fstream>
#include <cstdlib>

#include "dirac_s2.h"

#include "header_cusolver.hpp"


int main(int argc, char* argv[]){

  int device;
  cudacheck(cudaGetDeviceCount(&device));
  cudaDeviceProp device_prop[device];
  cudaGetDeviceProperties(&device_prop[0], 0);
  std::cout << "dev = " << device_prop[0].name << std::endl;
  cudacheck(cudaSetDevice(0));// "TITAN V"
  std::cout << "(GPU device is set.)" << std::endl;

  // ---------------------------------------

  const int q=5; // icosahedron
  const int n_refine=3; // no refinement


  QfeLatticeS2 lattice(q, n_refine);
  Dirac1fonS2 D(lattice, n_refine, 0.0, 1.0);
  auto mat = D.matrix_form();

  // =========================================
  // cusolver
  cusolverDnHandle_t handle = NULL;
  cudaStream_t stream = NULL;
  cusolverDnParams_t params = NULL;

  const int n = mat.cols(); // Number of rows (or columns) of matrix A.
  const int lda = n;

  complex *A, *W;
  A = (complex*)malloc(n*n*CD);
  W = (complex*)malloc(n*CD);
  for(int j=0; j<n; j++) for(int i=0; i<n; i++) A[n*j+i] = cplx(mat(i,j));
  for(int i=0; i<n; i++) W[i] = cplx(0.);

  complex *d_A, *d_W, *d_VL, *d_VR;
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
  cudacheck(cusolverDnCreate(&handle));
  cudacheck(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
  cudacheck(cusolverDnSetStream(handle, stream));
  cudacheck(cusolverDnCreateParams(&params));

  cudacheck(cudaMalloc( &d_A, CD * n*n ));
  cudacheck(cudaMalloc( &d_W, CD * n ));
  cudacheck(cudaMalloc( &d_VL, CD * n*n ));
  cudacheck(cudaMalloc( &d_VR, CD * n*n ));
  cudacheck(cudaMalloc( &d_info, sizeof(int)));

  cudacheck( cudaMemcpy(d_A, A, CD*n*n, H2D) );

  // step 3: query working space of syevd
  cusolverEigMode_t jobvl = CUSOLVER_EIG_MODE_NOVECTOR;
  cusolverEigMode_t jobvr = CUSOLVER_EIG_MODE_NOVECTOR;
  cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;

  cudacheck( cusolverDnXgeev_bufferSize( handle,
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

  cudacheck(cudaMalloc( &d_work, workspaceInBytesOnDevice ) );
  h_work = malloc(workspaceInBytesOnHost);

  // step 4: compute spectrum
  cudacheck( cusolverDnXgeev( handle,
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

  cudacheck(cudaMemcpy( W, d_W, CD*n, D2H) );
  cudacheck(cudaMemcpy( &info, d_info, sizeof(int), D2H ));

  std::cout << "info (0=success) = " << info << std::endl;
  for(int i=0; i<n; i++){
    std::clog << real(W[i]) << " " << imag(W[i]) << std::endl;
  }

  /* free resources */
  free(A);
  free(h_work);

  cudacheck(cudaFree(d_A));
  cudacheck(cudaFree(d_W));
  cudacheck(cudaFree(d_VL));
  cudacheck(cudaFree(d_VR));
  cudacheck(cudaFree(d_info));
  cudacheck(cudaFree(d_work));

  cudacheck(cusolverDnDestroyParams(params));
  cudacheck(cusolverDnDestroy(handle));
  cudacheck(cudaStreamDestroy(stream));

  cudacheck(cudaDeviceReset());

  return 0; // EXIT_SUCCESS;














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
  // 		       const cuDoubleComplex *csrValA,
  // 		       const int *csrRowPtrA,
  // 		       const int *csrColIndA,
  // 		       cuDoubleComplex mu0,
  // 		       const cuDoubleComplex *x0,
  // 		       int maxite,
  // 		       double tol,
  // 		       cuDoubleComplex *mu,
  // 		       cuDoubleComplex *x);

  // Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver( mat );
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
