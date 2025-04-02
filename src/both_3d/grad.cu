#include <typeinfo>
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

// BaseLink
// using Link = std::array<Idx,2>; // <int,int>;
// using Face = std::vector<Idx>;

using MS=Eigen::Matrix2cd;
using VD=Eigen::Vector2d;
using VE=Eigen::Vector3d;
using VC=Eigen::VectorXcd;

static constexpr int NS = 2;
static constexpr int DIM = 2;
static constexpr Complex I = Complex(0.0, 1.0);


#define IS_DUAL
// #define IS_OVERLAP

// #define IsVerbose
// #define InfoForce
// #define InfoDelta

namespace Comp{
  constexpr bool is_compact=false;

#ifdef IS_OVERLAP
  constexpr int NPARALLEL=12; // 12
  constexpr int NPARALLEL2=1; // 12
  constexpr int NSTREAMS=4; // 4
#else
  constexpr int NPARALLEL=1; // 12
  constexpr int NPARALLEL2=12; // 12
  constexpr int NSTREAMS=12; // 4
#endif
  constexpr int NPARALLEL3=10; // 12

  constexpr int N_REFINE=2;
  constexpr int NS=2;

  constexpr int Nt=10;

#ifdef IS_DUAL
  constexpr Idx N_SITES=20*N_REFINE*N_REFINE;
#else
  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
#endif

  constexpr Idx Nx=NS*N_SITES; // matrix size of DW
  constexpr Idx N=Nx*Nt; // matrix size of DW

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-8;
}

const std::string dir = "/mnt/hdd_barracuda/qed3/dats/";

#include "timer.h"

#include "s2n_simp.h"
#include "s2n_dual.h"
#include "rng.h"
// #include "gauge.h"
#include "gauge_ext.h"
// #include "action.h"
#include "action_ext.h"

#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cublas_api.h>
#include <cusolverDn.h>
using CuC = cuDoubleComplex;
#include "gpu_header.h"

// ======================================

#include "sparse_matrix.h"
#include "dirac_simp.h"
#include "sparse_dirac.h"
#include "matpoly.h"
#include "dirac_pf.h"
#include "overlap.h"
#include "pseudofermion.h"

# include "integrator.h"
#include "hmc.h"

#include "obs.h"



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
  // using Link = std::array<Idx,2>; // <int,int>;
  constexpr Idx N = Comp::N;
  constexpr int Nt = Comp::Nt;

#ifdef IS_DUAL
  using Base=S2Trivalent;
  using WilsonDirac=DiracExt<Base, DiracS2Dual>;
#else
  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
#endif

  using Force=GaugeExt<Base,Nt,Comp::is_compact>;
  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  // using Action=U1WilsonExt;

  using Rng=ParallelRngExt<Base,Nt>;
  using Overlap=Overlap<WilsonDirac>;

  using Fermion=DiracPf<WilsonDirac>;

  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;

  // ----------------------

#ifdef IS_OVERLAP
  const double r = 1.0;
  const double M5 = -1.6/2.0 * 0.5*(1.0 + std::sqrt( 5.0 + 2.0*std::sqrt(2.0) ));
#else
  const double r = 1.0;
  const double M5 = 0.0;
#endif
  // WilsonDirac DW(lattice, 0.0, r, M5);
  const double c = 1.0;
  WilsonDirac DW(base, 0.0, 1.0, M5, c);


  std::cout << "# DW set" << std::endl;

  Gauge U(lattice);
  Rng rng(lattice);
  // U.gaussian( rng, 0.2 );

  // ---------------------

#ifdef IS_OVERLAP
  // Fermion D(DW, 21);
  // std::cout << "# Dov set; M5 = " << M5 << std::endl;
  // D.update(U);
  // std::cout << "# min max ratio: "
  //           << D.lambda_min << " "
  //           << D.lambda_max << " "
  //           << D.lambda_min/D.lambda_max << std::endl;
  // std::cout << "# delta = " << D.Delta() << std::endl;

  // auto f_DHD = std::bind(&Fermion::sq_deviceAsyncLaunch, &D,
  //                        std::placeholders::_1, std::placeholders::_2);
  // auto f_DH = std::bind(&Fermion::adj_deviceAsyncLaunch, &D,
  //                       std::placeholders::_1, std::placeholders::_2);
  // LinOpWrapper M_DHD( f_DHD );
  // // LinOpWrapper M_DH( f_DH );

  // // MatPoly DHD;
  // // DHD.push_back ( cplx(1.0), {&M_DHD} );
  // //
  // // MatPoly DH;
  // // DH.push_back ( cplx(1.0), {&M_DH} );
  // MatPoly Op_DHD; Op_DHD.push_back ( cplx(1.0), {&M_DHD} );
  // auto f_mgrad_DHD = std::bind(&Fermion::grad_deviceAsyncLaunch, &D,
  //                              std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
#else
  Fermion D(DW);
  // DWDevice<WilsonDirac,Lattice> d_DW(DW); // actual data used in M_DW, M_DWH
  // CSR M_DW;
  // CSR M_DWH;
  // d_DW.associateCSR( M_DW, false );
  // d_DW.associateCSR( M_DWH, true );
  D.update( U );

  auto f_DHD = std::bind(&Fermion::sq_deviceAsyncLaunch, &D,
                         std::placeholders::_1, std::placeholders::_2);
  auto f_DH = std::bind(&Fermion::adj_deviceAsyncLaunch, &D,
                        std::placeholders::_1, std::placeholders::_2);

  LinOpWrapper M_DHD( f_DHD );
  MatPoly Op_DHD; Op_DHD.push_back ( cplx(1.0), {&M_DHD} );
  // MatPoly Op_DHD;
  // Op_DHD.push_back ( cplx(1.0), {&D.M_DW, &D.M_DWH} );
  //
  // MatPoly DH;
  // DH.push_back ( cplx(1.0), {&D.M_DWH} );
  // auto f_DHD = std::bind(&Fermion::sq_deviceAsyncLaunch, &Dov,
  //                        std::placeholders::_1, std::placeholders::_2);
  auto f_mgrad_DHD = std::bind(&Fermion::grad_deviceAsyncLaunch, &D,
                               std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);

#endif


  // -----------------------------------------------------------

  const double gR = 0.4;
  const double beta = 1.0/(gR*gR);
  Action SW(beta);


  PseudoFermion pf( Op_DHD, f_DH, f_mgrad_DHD, lattice );

  Timer timer;

  // ------------------

  Idx il=1;
  Link ell = lattice.links[il];
  std::cout << "debug. ell = " << ell[0] << " " << ell[1] << std::endl;

  const double eps = 1.0e-5;
  Gauge UP(U);
  UP[il] += eps;
  Gauge UM(U);
  UM[il] -= eps;

  std::cout << " --- Dov.update : " << timer.currentSeconds() << std::endl;
  D.update(U);
  std::cout << " --- pf.gen : " << timer.currentSeconds() << std::endl;
  pf.gen( rng );

  std::cout << " --- grad constructor : " << timer.currentSeconds() << std::endl;
  Force grad(lattice);

  std::cout << " --- pre calc : " << timer.currentSeconds() << std::endl;
  D.precalc_grad_deviceAsyncLaunch( U, pf.d_eta );
  std::cout << " --- get force : " << timer.currentSeconds() << std::endl;
  pf.get_force( grad, U );

  std::cout << " --- fin : " << timer.currentSeconds() << std::endl;

  std::cout << "grad = " << grad[il] << std::endl;
  D.update(UP);
  pf.update_eta();
  double sfp = pf.S();

  D.update(UM);
  pf.update_eta();
  double sfm = pf.S();

  double chck = (sfp-sfm)/(2.0*eps);
  std::cout << "check = " << chck << std::endl;

  // -----------------


  // const double eps = 1.0e-5;

  // std::cout << " --- Dov.update : " << timer.currentSeconds() << std::endl;
  // Dov.update(U);
  // std::cout << " --- pf.gen : " << timer.currentSeconds() << std::endl;
  // pf.gen( rng );

  // std::cout << " --- grad constructor : " << timer.currentSeconds() << std::endl;
  // Force dSf(lattice);
  // std::cout << " --- pre calc : " << timer.currentSeconds() << std::endl;
  // Dov.precalc_grad_deviceAsyncLaunch( U, pf.d_eta );
  // std::cout << " --- get force : " << timer.currentSeconds() << std::endl;
  // pf.get_force( dSf, U );
  // std::cout << " --- fin : " << timer.currentSeconds() << std::endl;


  // for(Idx il=0; il<lattice.n_links; il++) std::cout << "grad = " << il << " " << dSf[il] << std::endl;

  // const double tmax = 0.5; // 1.0; // 0.1
  // const int nsteps=5;
  // ExplicitLeapfrogML integrator( tmax, nsteps, 10 );


  // Force pi( lattice );
  // pi.gaussian( rng );
  // // Force pi0=pi;

  // for(Idx il=0; il<lattice.n_links; il++){
  //   //   Idx il=3;
  //   // Link ell = lattice.links[il];

  //   Gauge UP(U);
  //   UP[il] += eps;
  //   Gauge UM(U);
  //   UM[il] -= eps;


  //   double Hp, Hm;
  //   {
  //     HMC hmc(rng, &SW, &Dov, UP, pi, &pf, &integrator);
  //     Dov.update(UP);
  //     pf.update_eta();
  //     Hp = hmc.H();
  //   }

  //   {
  //     HMC hmc(rng, &SW, &Dov, UM, pi, &pf, &integrator);
  //     Dov.update(UM);
  //     pf.update_eta();
  //     Hm = hmc.H();
  //   }

  //   double chck = (Hp-Hm)/(2.0*eps);
  //   std::cout << "check = " << il << " " << chck << std::endl;
  // }
  // // -----------------


  // Force pi( lattice );
  // pi.gaussian( rng );
  // Force pi0=pi;

  // Gauge U0=U;
  // D.update(U);
  // pf.gen( rng );
  // D.precalc_grad_deviceAsyncLaunch( U, pf.d_eta );

  // const double tmax = 0.2; // 1.0; // 0.1
  // for(int nsteps=1; nsteps<=5; nsteps+=1){
  //   // const int nsteps=5;
  //   ExplicitLeapfrogML integrator( tmax, nsteps, 10 );
  //   // ExplicitLeapfrogML integrator( tmax, nsteps, 100 );
  //   pi = pi0;
  //   U = U0;
  //   HMC hmc(rng, &SW, &D, U, pi, &pf, &integrator);
  //   D.update( U ); pf.update_eta();
  //   D.precalc_grad_deviceAsyncLaunch( U, pf.d_eta );
  //   const double h0 = hmc.H();
  //   hmc.integrate();
  //   const double h1 = hmc.H();
  //   double dH = h1-h0;
  //   std::cout << tmax/nsteps << " " << dH << std::endl;
  // }



  // Force pi( lattice );
  // pi.gaussian( rng );
  // Force pi0=pi;

  // Gauge U0=U;
  // D.update(U);
  // pf.gen( rng );
  // D.precalc_grad_deviceAsyncLaunch( U, pf.d_eta );

  // const double tmax = 1.0; // 0.1
  // // for(int nsteps=1; nsteps<=5; nsteps+=1){
  // const int nsteps=8;
  // ExplicitLeapfrogML integrator( tmax, nsteps, 20 );
  // // ExplicitLeapfrogML integrator( tmax, nsteps, 100 );
  // pi = pi0;
  // U = U0;
  // HMC hmc(rng, &SW, &D, U, pi, &pf, &integrator);
  // D.update( U ); pf.update_eta();
  // D.precalc_grad_deviceAsyncLaunch( U, pf.d_eta );

  // double rate, dH;
  // bool is_accept;
  // for(int k=0; k<10; k++){
  //   Timer timer;
  //   hmc.run( rate, dH, is_accept, true);
  //   std::cout << "# dH : " << dH
  //             << " is_accept : " << is_accept << std::endl;
  //   std::cout << "# HMC : " << timer.currentSeconds() << " sec" << std::endl;
  // }

  // double r_mean;
  // const int kmax=50;
  // for(int k=0; k<kmax; k++){
  //   Timer timer;
  //   hmc.run( rate, dH, is_accept);
  //   std::cout << "# dH : " << dH
  //             << " is_accept : " << is_accept << std::endl;
  //   r_mean += rate;
  //   std::cout << "# HMC : " << timer.currentSeconds() << " sec" << std::endl;
  // }
  // r_mean /= kmax;
  // std::cout << "# r_mean = " << r_mean << std::endl;



  // CUDA_CHECK(cudaDeviceReset());
  return 0;

}

