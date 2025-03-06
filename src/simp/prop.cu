#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>

#include <cstdint>
#include <complex>
using Idx = std::int32_t;

using Double = double;
using Complex = std::complex<Double>;

#include "s2n.h"
#include "dirac_s2.h"



#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cublas_api.h>
#include <cusolverDn.h>
using CuC = cuDoubleComplex;
#include "header_cusolver.hpp"


#include "../../integrator/geodesic.h"

int main(int argc, char* argv[]){

  int device;
  cudacheck(cudaGetDeviceCount(&device));
  cudaDeviceProp device_prop[device];
  cudaGetDeviceProperties(&device_prop[0], 0);
  std::cout << "1. dev = " << device_prop[0].name << std::endl;
  cudacheck(cudaSetDevice(0));// "TITAN V"
  std::cout << "(GPU device is set.)" << std::endl;

  // ---------------------------------------

  const int q=5; // icosahedron
  int n_refine=1; // no refinement

  if (argc>1) n_refine = atoi(argv[1]);

  using Lattice=S2Simp;
  using WilsonDirac=DiracS2Simp<Lattice>;

  // std::cout << "debug.1" << std::endl;
  Lattice lattice(n_refine);
  // std::cout << "debug.2" << std::endl;


  WilsonDirac D(lattice, n_refine, 0.0, 1.0);
  // std::cout << "debug.3" << std::endl;
  auto A = D.matrix_form();
  Eigen::VectorXcd B = Eigen::VectorXcd::Zero(A.cols());
  B[0] = 1.0;
  // std::cout << "debug.4" << std::endl;

  {
    cusolverDnHandle_t cusolverH = NULL;
    cudaStream_t stream = NULL;

    const int m = A.cols();
    const int lda = m;
    const int ldb = m;

    // const std::vector<double> A = ; // {1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 10.0};
    // const std::vector<double> B = {1.0, 2.0, 3.0};
    std::vector<Complex> X(m, 0);
    std::vector<Complex> LU(lda * m, 0);
    std::vector<int> Ipiv(m, 0);
    int info = 0;

    CuC *d_A = nullptr; /* device copy of A */
    CuC *d_B = nullptr; /* device copy of B */
    int *d_Ipiv = nullptr; /* pivoting sequence */
    int *d_info = nullptr; /* error info */

    int lwork = 0;            /* size of workspace */
    CuC *d_work = nullptr; /* device workspace for getrf */

    const int pivot_on = 1;

    if (pivot_on) {
        printf("pivot is on : compute P*A = L*U \n");
    } else {
        printf("pivot is off: compute A = L*U (not numerically stable)\n");
    }

    printf("A = (matlab base-1)\n");
    //    print_matrix(m, m, A.data(), lda);
    printf("=====\n");

    printf("B = (matlab base-1)\n");
    // print_matrix(m, 1, B.data(), ldb);
    printf("=====\n");

    /* step 1: create cusolver handle, bind a stream */
    CUSOLVER_CHECK(cusolverDnCreate(&cusolverH));

    CUDA_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    CUSOLVER_CHECK(cusolverDnSetStream(cusolverH, stream));

    /* step 2: copy A to device */
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_A), CD * A.size()));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_B), CD * B.size()));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_Ipiv), sizeof(int) * Ipiv.size()));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_info), sizeof(int)));

    CUDA_CHECK(
               cudaMemcpyAsync(d_A, reinterpret_cast<CuC*>(A.data()), CD * A.size(), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(
               cudaMemcpyAsync(d_B, reinterpret_cast<CuC*>(B.data()), CD * B.size(), cudaMemcpyHostToDevice, stream));

    /* step 3: query working space of getrf */
    CUSOLVER_CHECK(cusolverDnZgetrf_bufferSize(cusolverH, m, m, d_A, lda, &lwork));

    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_work), CD * lwork));

    /* step 4: LU factorization */
    if (pivot_on) {
        CUSOLVER_CHECK(cusolverDnZgetrf(cusolverH, m, m, d_A, lda, d_work, d_Ipiv, d_info));
    } else {
        CUSOLVER_CHECK(cusolverDnZgetrf(cusolverH, m, m, d_A, lda, d_work, NULL, d_info));
    }

    if (pivot_on) {
        CUDA_CHECK(cudaMemcpyAsync(Ipiv.data(), d_Ipiv, sizeof(int) * Ipiv.size(),
                                   cudaMemcpyDeviceToHost, stream));
    }
    CUDA_CHECK(
        cudaMemcpyAsync(LU.data(), d_A, CD * A.size(), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost, stream));

    CUDA_CHECK(cudaStreamSynchronize(stream));

    if (0 > info) {
        printf("%d-th parameter is wrong \n", -info);
        exit(1);
    }
    if (pivot_on) {
        printf("pivoting sequence, matlab base-1\n");
        for (int j = 0; j < m; j++) {
            printf("Ipiv(%d) = %d\n", j + 1, Ipiv[j]);
        }
    }
    printf("L and U = (matlab base-1)\n");
    // print_matrix(m, m, LU.data(), lda);
    printf("=====\n");

    /*
     * step 5: solve A*X = B
     *       | 1 |       | -0.3333 |
     *   B = | 2 |,  X = |  0.6667 |
     *       | 3 |       |  0      |
     *
     */
    if (pivot_on) {
        CUSOLVER_CHECK(cusolverDnZgetrs(cusolverH, CUBLAS_OP_N, m, 1, /* nrhs */
                                        d_A, lda, d_Ipiv, d_B, ldb, d_info));
    } else {
        CUSOLVER_CHECK(cusolverDnZgetrs(cusolverH, CUBLAS_OP_N, m, 1, /* nrhs */
                                        d_A, lda, NULL, d_B, ldb, d_info));
    }

    CUDA_CHECK(
        cudaMemcpyAsync(X.data(), d_B, CD * X.size(), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));

    std::vector<double> lengths;
    {
      const auto x0 = lattice.r[0];
      for(int ix=0; ix<lattice.n_sites; ix++){
        const auto x1 = lattice.r[ix];
        double len = Geodesic::geodesicLength(Geodesic::Pt(x0), Geodesic::Pt(x1));
        // std::cout << "len = " << len << std::endl;
        lengths.push_back(len);
      }
    }


    printf("X = (matlab base-1)\n");
    // print_matrix(m, 1, X.data(), ldb);
    int counter=0;
    for(auto elem : X){
      std::clog << std::setprecision(16)
                << std::setw(25) << lengths[int(counter/2)] << " "
                << std::setw(25) << 1.0/D.a * real(elem) << " "
                << std::setw(25) << 1.0/D.a * imag(elem) << std::endl;
      counter++;
    }
    printf("=====\n");

    /* free resources */
    CUDA_CHECK(cudaFree(d_A));
    CUDA_CHECK(cudaFree(d_B));
    CUDA_CHECK(cudaFree(d_Ipiv));
    CUDA_CHECK(cudaFree(d_info));
    CUDA_CHECK(cudaFree(d_work));

    CUSOLVER_CHECK(cusolverDnDestroy(cusolverH));

    CUDA_CHECK(cudaStreamDestroy(stream));

    CUDA_CHECK(cudaDeviceReset());
  }

}



  // for(int ix=0; ix<lattice.n_sites; ix++){
  //   for(int jj=0; jj<lattice.sites[ix].nn; jj++){
  //     const int iy = lattice.sites[ix].neighbors[jj];
  //     auto mat1 = ( D.sigma[0] - D.gamma(ix, iy) ) * D.Omega(ix, iy);
  //     auto mat2 = D.Omega(ix, iy) * ( D.sigma[0] - D.gamma(iy, ix, M_PI) );
  //     std::cout << mat1-mat2 << std::endl;
  //   }}

  // ----------------------------------
