#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdint>
#include <complex>

using Double = double;
using Idx = std::int32_t;
using Complex = std::complex<double>;

// #define IS_DUAL

namespace Comp{
  constexpr int NPARALLEL=12;
  constexpr int NSTREAMS=4;

  constexpr int N_REFINE=4;
  constexpr int NS=2;

#ifdef IS_DUAL
  constexpr Idx N_SITES=20*N_REFINE*N_REFINE;
#else
  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
#endif

  constexpr Idx N=NS*N_SITES; // matrix size of DW

  // const double TOL=1.0e-9;
  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-8;
}

const std::string dir = "/mnt/hdd_barracuda/qed3/dats/";

#include "s2n_dual.h"
#include "s2n_simp.h"
#include "rng.h"
#include "gauge.h"
#include "action.h"





#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cublas_api.h>
#include <cusolverDn.h>
using CuC = cuDoubleComplex;
#include "gpu_header.h"
// #include "header_cusolver.hpp"

#include "sparse_matrix.h"
// #include "pseudofermion.h"
#include "dirac_dual.h"
#include "dirac_simp.h"

#include "sparse_dirac.h"
#include "matpoly.h"

#include "overlap.h"


#include "../../integrator/geodesic.h"

int main(int argc, char* argv[]){

  int device;
  CUDA_CHECK(cudaGetDeviceCount(&device));
  cudaDeviceProp device_prop[device];
  cudaGetDeviceProperties(&device_prop[0], 0);
  std::cout << "# dev = " << device_prop[0].name << std::endl;
  CUDA_CHECK(cudaSetDevice(0));// "TITAN V"
  std::cout << "# (GPU device is set.)" << std::endl;

  // ---------------------------------------

#ifdef IS_DUAL
  using Lattice=S2Trivalent;
  using Gauge=U1onS2<Lattice,false>;
  using WilsonDirac=DiracS2Dual<Gauge>;
#else
  using Lattice=S2Simp;
  using Gauge=U1onS2<Lattice,false>;
  using WilsonDirac=DiracS2Simp<Gauge>;
#endif
  using Rng=ParallelRng<Lattice>;
  using Overlap=Overlap<Gauge,WilsonDirac,Lattice>;

  Lattice lattice(Comp::N_REFINE);

  // const double r = 1.0;
  // const double M5 = 0.0;

#ifdef IS_DUAL
  const double r = 1.0/3.0;
  const double M5 = -1.8;
#else
  const double r = 1.0;
  const double M5 = -1.8;
#endif
  WilsonDirac DW(lattice, 0.0, r, M5);

  Gauge U(lattice);
  Rng rng(lattice);
  // U.gaussian( rng, 0.2 );

  DWDevice<WilsonDirac,Lattice> d_DW(DW); // actual data used in M_DW, M_DWH
  CSR M_DW;
  CSR M_DWH;
  d_DW.associateCSR( M_DW, false );
  d_DW.associateCSR( M_DWH, true );
  d_DW.update( U );

  // Overlap Dov(DW, 31);
  // Dov.update(U);
  // std::cout << "# min max ratio: "
  //           << Dov.lambda_min << " "
  //           << Dov.lambda_max << " "
  //           << Dov.lambda_min/Dov.lambda_max << std::endl;
  // std::cout << "# delta = " << Dov.Delta() << std::endl;


  // // auto A = DW.matrix_form();

  // auto f_Op = std::bind(&Overlap::mult_deviceAsyncLaunch, &Dov, std::placeholders::_1, std::placeholders::_2);
  // LinOpWrapper M_Op( f_Op );

  MatPoly Op;
  // Op.push_back ( cplx(1.0), {&M_Op} );
  // Op.push_back ( cplx(1.0), {&Dov.M_DW} );
  Op.push_back ( cplx(1.0), {&M_DW} );

  constexpr Idx N = Comp::N;
  Eigen::MatrixXcd A(N, N);
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
      A.block(0,i,N,1) = Eigen::Map<Eigen::MatrixXcd>(Dxi.data(), N, 1);
      std::cout << "i = " << i << " finished." << std::endl;
    }
  }





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



#ifdef IS_DUAL
    std::vector<double> lengths;
  {
    std::string dir = "/mnt/hdd_barracuda/qed3/dats/";
    std::vector<Geodesic::V3> sites;
    {
      std::ifstream file(dir+"pts_dual_n"+std::to_string(Comp::N_REFINE)+"_singlepatch.dat");

      std::string str;
      while (std::getline(file, str)){
        std::istringstream iss(str);
        double v1, v2, v3;
        iss >> v1;
        iss >> v2;
        iss >> v3;
        sites.push_back( Geodesic::V3(v1, v2, v3) );
      }
    }
    const auto x0 = sites[0];
    for(const auto& elem : sites){
      double len = Geodesic::geodesicLength(Geodesic::Pt(x0), Geodesic::Pt(elem));
      // std::cout << "len = " << len << std::endl;
      lengths.push_back(len);
    }
  }
  double alat;
  {
    std::string dir = "/mnt/hdd_barracuda/qed3/dats/";
    std::ifstream file(dir+"alat_n"+std::to_string(Comp::N_REFINE)+"_singlepatch.dat");

    std::string str;
    std::getline(file, str);
    std::istringstream iss(str);
    iss >> alat;
  }


        printf("X = (matlab base-1)\n");
    // print_matrix(m, 1, X.data(), ldb);
    int counter=0;
    for(auto elem : X){
      std::clog << std::setprecision(16)
                << std::setw(25) << lengths[int(counter/2)] << " "
                << std::setw(25) << 1.0/alat * real(elem) << " "
                << std::setw(25) << 1.0/alat * imag(elem) << std::endl;
      counter++;
    }
#else
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
              << std::setw(25) << 1.0/DW.a * real(elem) << " "
              << std::setw(25) << 1.0/DW.a * imag(elem) << std::endl;
    counter++;
  }
#endif

    printf("=====\n");

    /* free resources */
    CUDA_CHECK(cudaFree(d_A));
    CUDA_CHECK(cudaFree(d_B));
    CUDA_CHECK(cudaFree(d_Ipiv));
    CUDA_CHECK(cudaFree(d_info));
    CUDA_CHECK(cudaFree(d_work));

    CUSOLVER_CHECK(cusolverDnDestroy(cusolverH));

    CUDA_CHECK(cudaStreamDestroy(stream));
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
