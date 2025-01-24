#include <typeinfo>
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
  constexpr int NPARALLEL=10;

  constexpr int N_REFINE=2;
  constexpr int NS=2;
  constexpr Idx N_SITES=20*N_REFINE*N_REFINE;
  constexpr Idx N=NS*N_SITES; // matrix size of DW
}

// #define IsVerbose
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

#include "s2n.h"
#include "rng.h"
#include "gauge.h"
#include "force.h"
#include "action.h"
#include "sparse_matrix.h"
#include "dirac.h"
#include "sparse_dirac.h"
#include "matpoly.h"
#include "overlap.h"
#include "pseudofermion.h"

# include "integrator.h"
#include "hmc.h"
// #include "dirac_s2_dual.h"
// #include "header_cusolver.hpp"


// TODO: Cusparse for SparseMatrix::act_gpu, probably defining handle in matpoly.h
// make 2 streams in V Vdag in square in Overlap
// all the operation on GPU in Overlap::operator()
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

  using Lattice=S2Trivalent;
  using Gauge=U1onS2<false>;
  using Force=U1onS2<false>;
  using Action=U1Wilson;
  using Fermion=Overlap;
  using Rng=ParallelRng<Lattice>;
  using WilsonDirac=Dirac1fonS2;

  using Link = std::array<Idx,2>; // <int,int>;
  constexpr Idx N = Comp::N;

  // ----------------------

  Lattice lattice(Comp::N_REFINE);
  Gauge U(lattice);
  Rng rng(lattice);
  U.gaussian( rng, 0.2 );

  // ------------------

  const double gR = 0.4;
  const double beta = 1.0/(gR*gR);
  Action SW(beta);

  // -----------------

  const double M5 = -1.8;
  WilsonDirac DW(lattice, M5);

  Fermion Dov(DW, 15);

  const auto f_DHDov = std::bind(&Overlap::sq_device, &Dov,
                                 std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_DHDov( f_DHDov );
  MatPoly Op_DHDov; Op_DHDov.push_back ( cplx(1.0), {&M_DHDov} );
  auto f_DHov = std::bind(&Overlap::adj_device, &Dov,
                          std::placeholders::_1, std::placeholders::_2);
  auto f_mgrad_DHDov = std::bind(&Overlap::grad_device, &Dov,
                                 std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);

  // // ------------------

  // Idx il=2;
  // Link ell = lattice.links[il];

  // const double eps = 1.0e-5;
  // Gauge UP(U);
  // UP[il] += eps;
  // Gauge UM(U);
  // UM[il] -= eps;

  // Dov.update(U);
  // pf.gen( rng );
  // Force grad(lattice);
  // pf.get_force( grad, U );

  // Dov.update(UP);
  // pf.update_eta();
  // double sfp = pf.S();

  // Dov.update(UM);
  // pf.update_eta();
  // double sfm = pf.S();

  // double chck = (sfp-sfm)/(2.0*eps);
  // std::cout << "grad = " << grad[il] << std::endl;
  // std::cout << "check = " << chck << std::endl;

  // // -----------------

  Force pi( lattice );
  pi.gaussian( rng );
  Force pi0=pi;

  Gauge U0=U;
  Dov.update(U);

  PseudoFermion pf( Op_DHDov, f_DHov, f_mgrad_DHDov );
  pf.gen( rng );

  double tmax = 1.0; // 0.1
  // for(int nsteps=1; nsteps<=1; nsteps+=1){
  const int nsteps=4;
  // ExplicitLeapfrog integrator( tmax, nsteps );
  ExplicitLeapfrogML integrator( tmax, nsteps, 100 );
  HMC hmc(rng, &SW, &Dov, U, pi, &pf, &integrator);
  // pi = pi0;
  // U = U0;
  // Dov.update( U ); pf.update_eta();
  //     const double h0 = hmc.H();
  //   hmc.integrate();
  //   const double h1 = hmc.H();
  //   double dH = h1-h0;
  //   std::cout << tmax/nsteps << " " << dH << std::endl;
  // }

  double r, dH;
  bool is_accept;
  for(int k=0; k<20; k++){
    hmc.run( r, dH, is_accept, true);
    std::cout << "# dH : " << dH
              << " is_accept : " << is_accept << std::endl;
  }
  for(int k=0; k<20; k++){
    hmc.run( r, dH, is_accept);
    std::cout << "# dH : " << dH
              << " is_accept : " << is_accept << std::endl;
  }

  // CUDA_CHECK(cudaDeviceReset());
  return 0;

}

