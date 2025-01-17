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

namespace CompilationConst{
  constexpr int NPARALLEL=10;

  constexpr int N_REFINE=2;
  constexpr int NS=2;
  constexpr Idx N_SITES=20*N_REFINE*N_REFINE;
  constexpr Idx N=NS*N_SITES; // matrix size of DW
}

#define IsVerbose

#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cublas_api.h>
#include <cusolverDn.h>
using CuC = cuDoubleComplex;
#include "gpu_header.h"

// ======================================

#include "s2n.h"
#include "rng.h"
#include "gauge.h"
#include "action.h"

#include "sparse_matrix.h"
// #include "pseudofermion.h"

#include "dirac.h"

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

  using Gauge=U1onS2;
  // using Force=U1onS2;
  // using Action=U1Wilson;
  // using Fermion=Dirac1fonS2;
  // using HMC=HMC<Force,Gauge,Action,Fermion>;
  using Rng=ParallelRng;

  Lattice lattice(CompilationConst::N_REFINE);
  // Dirac1fonS2 D(lattice, 0.0, 1.0);

  using WilsonDirac=Dirac1fonS2;
  // using Overlap=OverlapPseudoFermion;

  Gauge U(lattice);
  Rng rng(lattice);
  // U.gaussian( rng );

  const double M5 = -2.0;
  WilsonDirac DW(lattice, M5);

  constexpr Idx N = CompilationConst::N;

  using Link = std::array<Idx,2>; // <int,int>;

  // const Idx ix=2;
  // const Idx iy=lattice.nns[ix][0];
  // const Link ell{ix,iy};

  // COO coo;
  // DW.d_coo_format(coo.en, U, ell);
  // coo.do_it();








  // {
  //   // PseudoFermion phi( D, U, rng );
  //   // Force exact = phi.dS( U );

  //   Overlap Dov(DW);
  //   Dov.compute(U);

  //   std::vector<Complex> xi(N);
  //   for(int i=0; i<N; i++) xi[i] = rng.gaussian();

  //   // -------------

  //   std::vector<Complex> eta(N);
  //   {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0), {&(Dov.M_DW), &(Dov.M_DWH)} );
  //     Op.solve<N>( eta, xi );
  //   }

  //   double Sf = 0.0;
  //   for(Idx i=0; i<N; i++) Sf += std::real( std::conj(xi[i]) * eta[i] );

  //   // -------------

  //   Idx il=0;

  //   {
  //     double dSf = Dov.grad( lattice.links[il], U, eta );
  //     std::cout << "grad = " << dSf << std::endl;
  //   }

  //   // -------------

  //   {
  //     double Sfp = 0.0, Sfm = 0.0;
  //     const double eps = 1.0e-5;

  //     {
  //       Gauge UP(U);
  //       UP[il] += eps;
  //       Dov.compute(UP);

  //       MatPoly Op;
  //       Op.push_back ( cplx(1.0), {&(Dov.M_DW), &(Dov.M_DWH)} );
  //       Op.solve<N>( eta, xi );
  //       for(Idx i=0; i<N; i++) Sfp += std::real( std::conj(xi[i]) * eta[i] );
  //     }

  //     {
  //       Gauge UM(U);
  //       UM[il] -= eps;
  //       Dov.compute(UM);

  //       MatPoly Op;
  //       Op.push_back ( cplx(1.0), {&(Dov.M_DW), &(Dov.M_DWH)} );
  //       Op.solve<N>( eta, xi );
  //       for(Idx i=0; i<N; i++) Sfm += std::real( std::conj(xi[i]) * eta[i] );
  //     }

  //     std::cout << (Sfp-Sfm)/(2.0*eps) << std::endl;
  //   }

  // }



















    // -------------

    // std::vector<Complex> deriv(N);
    // // for(int il=0; il<U.lattice.n_links; il++){
    // {
    //   Idx il=0;

    //   const double eps = 1.0e-5;
    //   Gauge UP(U);
    //   Gauge UM(U);

    //   UP[il] += eps;
    //   UM[il] -= eps;

    //   std::vector<Complex> DxiP(N), DxiM(N);
    //   Dov.compute(UP);
    //   Dov( DxiP, xi );
    //   Dov.compute(UM);
    //   Dov( DxiM, xi );

    //   for(int i=0; i<N; i++) deriv[i] = (DxiP[i]-DxiM[i])/(2.0*eps);

    //   // double numeric = ( phi.S(UP) - phi.S(UM) ) / (2.0*eps);
    //   for(int i=0; i<N; i++) {
    //     std::cout << deriv[i] << std::endl;
    //   }
    // }






  // Eigen::MatrixXcd mat(N, N);
  // {
  //   for(Idx i=0; i<N; i++){
  //     Eigen::VectorXcd e = Eigen::VectorXcd::Zero(N);
  //     e(i) = 1.0;
  //     std::vector<Complex> xi(e.data(), e.data()+N);
  //     std::vector<Complex> Dxi(N);
  //     Dov( Dxi, xi );
  //     mat.block(0,i,N,1) = Eigen::Map<Eigen::MatrixXcd>(Dxi.data(), N, 1);
  //   }
  // }
  // std::cout << Dov.lambda_max << std::endl;


  // return 0; // EXIT_SUCCESS;

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
