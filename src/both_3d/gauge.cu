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

using Link = std::array<Idx,2>; // <int,int>;
using Face = std::vector<Idx>;

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
  constexpr int NPARALLEL3=12; // 12

  constexpr int N_REFINE=2;
  constexpr int NS=2;

  constexpr int Nt=12;

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

// #include "sparse_matrix.h"
// #include "dirac_simp.h"
// #include "sparse_dirac.h"
// #include "matpoly.h"
// #include "dirac_pf.h"
// #include "overlap.h"
// #include "pseudofermion.h"

# include "integrator.h"
#include "hmc.h"

#include "obs.h"


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

  std::cout << "# Nx = " << Comp::Nx << std::endl;

  // --------------------
  using Link = std::array<Idx,2>; // <int,int>;
  constexpr Idx Nx = Comp::Nx;
  constexpr int Nt = Comp::Nt;

#ifdef IS_DUAL
  using Base=S2Trivalent;
#else
  using Base=S2Simp;
#endif

  using Force=GaugeExt<Base,Nt,Comp::is_compact>;
  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  using Action=U1WilsonExt;

  using Rng=ParallelRngExt<Base,Nt>;

  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;

  // ----------------------

  // const double gR = 10.0;
  double beta = 24.0; // 1.0/(gR*gR);
  Action SW(beta, beta);

  Gauge U(base);
  srand( time(NULL) );
  Rng rng(base, rand());
  U.gaussian( rng, 0.2 );


  // //--------------------------------

  // int s = 0;
  // const Idx ix = 0;
  // const int yy=0;
  // Link link{ix, base.nns[ix][yy]};
  // std::cout << U.sp( s, link ) << std::endl
  //           << U.tp( s, link[1] ) << std::endl
  //           << U.sp( s+1, link[1] ) << std::endl
  //           << U.tp( s, link[0] ) << std::endl;
  // // std::cout << "angle = " << U.plaquette_angle( s, link ) << std::endl;

  // // std::cout << "debug. tp = " << std::endl;
  // // for(const auto& v : U.temporal) {
  // //   for(const double& elem : v) {
  // //     std::cout << elem << " ";
  // //   }
  // //   std::cout << std::endl;
  // // }

  // std::cout << "S(U) = " << SW(U) << std::endl;

  // Force grad(base);
  // SW.get_force( grad, U );


  // {  // sp
  //   for(int sl=0; sl<Comp::Nt; sl++){ // const int sl=0;
  //     for(Idx il=0; il<base.n_links; il++){ //Idx il=0;
  //       Link ell = base.links[il];
  //       // std::cout << "dS = " << grad.sp(sl, il) << std::endl;

  //       const double eps = 1.0e-5;

  //       Gauge UP(U);
  //       UP.sp(sl,il) += eps;

  //       Gauge UM(U);
  //       UM.sp(sl,il) -= eps;

  //       double sfp = SW(UP);
  //       double sfm = SW(UM);

  //       double chck = (sfp-sfm)/(2.0*eps);
  //       // std::cout << "check = " << chck << std::endl;
  //       std::cout << "diff = " << grad.sp(sl, il)-chck << std::endl;
  //       assert( std::abs(grad.sp(sl, il)-chck) < 1.0e-5 );
  //     }}
  // }


  // {  // tp
  //   // const int sl=0;
  //   // Idx ix=0;
  //   for(int sl=0; sl<Comp::Nt; sl++){ // const int sl=0;
  //     for(Idx ix=0; ix<base.n_sites; ix++){ //Idx ix=0;
  //       // std::cout << "dS = " << grad.tp(sl, ix) << std::endl;

  //       const double eps = 1.0e-5;

  //       Gauge UP(U);
  //       UP.tp(sl,ix) += eps;
  //       Gauge UM(U);
  //       UM.tp(sl,ix) -= eps;

  //       double sfp = SW(UP);
  //       double sfm = SW(UM);

  //       double chck = (sfp-sfm)/(2.0*eps);
  //       // std::cout << "check = " << chck << std::endl;
  //       std::cout << "diff = " << grad.tp(sl, ix)-chck << std::endl;
  //       assert( std::abs(grad.tp(sl, ix)-chck) < 1.0e-5 );
  //     }}
  // }


  // --------------------------------

  // Force pi(base);

  // HMCPureGauge hmc(rng, &SW, U, pi, 1.0, 10);

  // double rate, dH;
  // bool is_accept;
  // for(int k=0; k<10; k++){
  //   Timer timer;
  //   hmc.run( rate, dH, is_accept, true);
  //   std::cout << "# dH : " << dH
  //             << " is_accept : " << is_accept << std::endl;
  //   std::cout << "# HMC : " << timer.currentSeconds() << " sec" << std::endl;
  // }

  // pi.gaussian( rng );
  // Force pi0=pi;
  // Gauge U0=U;

  // const double tmax = 0.2; // 1.0; // 0.1
  // for(int nsteps=4; nsteps<=10; nsteps+=1){
  //   pi = pi0;
  //   U = U0;
  //   HMCPureGauge hmc(rng, &SW, U, pi, tmax, nsteps);
  //   const double h0 = hmc.H();
  //   hmc.integrate();
  //   const double h1 = hmc.H();
  //   double dH = h1-h0;
  //   std::cout << tmax/nsteps << " " << dH << std::endl;
  // }




  // --------------------------------


  //   using WilsonDirac=DiracS2Simp<Gauge>;

// #ifdef IS_OVERLAP
//   using Fermion=Overlap<Gauge,WilsonDirac,Lattice>;
// #else
//   using Fermion=DiracPf<Gauge,WilsonDirac,Lattice>;
// #endif

// #ifdef IS_OVERLAP
//   const double r = 1.0;
//   const double M5 = -1.6/2.0 * 0.5*(1.0 + std::sqrt( 5.0 + 2.0*std::sqrt(2.0) ));
// #else
//   const double r = 1.0;
//   const double M5 = 0.0;
// #endif
//   WilsonDirac DW(lattice, 0.0, r, M5);

//   std::cout << "# DW set" << std::endl;

//   Gauge U(lattice);
//   Rng rng(lattice);
//   // U.gaussian( rng, 0.2 );

//   // ---------------------

// #ifdef IS_OVERLAP
//   Fermion D(DW, 21);
//   std::cout << "# Dov set; M5 = " << M5 << std::endl;
//   D.update(U);
//   std::cout << "# min max ratio: "
//             << D.lambda_min << " "
//             << D.lambda_max << " "
//             << D.lambda_min/D.lambda_max << std::endl;
//   std::cout << "# delta = " << D.Delta() << std::endl;

//   auto f_DHD = std::bind(&Fermion::sq_deviceAsyncLaunch, &D,
//                          std::placeholders::_1, std::placeholders::_2);
//   auto f_DH = std::bind(&Fermion::adj_deviceAsyncLaunch, &D,
//                         std::placeholders::_1, std::placeholders::_2);
//   LinOpWrapper M_DHD( f_DHD );
//   // LinOpWrapper M_DH( f_DH );

//   // MatPoly DHD;
//   // DHD.push_back ( cplx(1.0), {&M_DHD} );
//   //
//   // MatPoly DH;
//   // DH.push_back ( cplx(1.0), {&M_DH} );
//   MatPoly Op_DHD; Op_DHD.push_back ( cplx(1.0), {&M_DHD} );
//   auto f_mgrad_DHD = std::bind(&Fermion::grad_deviceAsyncLaunch, &D,
//                                std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);

// #else
//   Fermion D(DW);
//   // DWDevice<WilsonDirac,Lattice> d_DW(DW); // actual data used in M_DW, M_DWH
//   // CSR M_DW;
//   // CSR M_DWH;
//   // d_DW.associateCSR( M_DW, false );
//   // d_DW.associateCSR( M_DWH, true );
//   D.update( U );

//   auto f_DHD = std::bind(&Fermion::sq_deviceAsyncLaunch, &D,
//                          std::placeholders::_1, std::placeholders::_2);
//   auto f_DH = std::bind(&Fermion::adj_deviceAsyncLaunch, &D,
//                         std::placeholders::_1, std::placeholders::_2);

//   LinOpWrapper M_DHD( f_DHD );
//   MatPoly Op_DHD; Op_DHD.push_back ( cplx(1.0), {&M_DHD} );
//   // MatPoly Op_DHD;
//   // Op_DHD.push_back ( cplx(1.0), {&D.M_DW, &D.M_DWH} );
//   //
//   // MatPoly DH;
//   // DH.push_back ( cplx(1.0), {&D.M_DWH} );
//   // auto f_DHD = std::bind(&Fermion::sq_deviceAsyncLaunch, &Dov,
//   //                        std::placeholders::_1, std::placeholders::_2);
//   auto f_mgrad_DHD = std::bind(&Fermion::grad_deviceAsyncLaunch, &D,
//                                std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);

// #endif


  // -----------------------------------------------------------

  // const double gR = 0.4;
  // const double beta = 1.0/(gR*gR);
  // Action SW(beta);


  // PseudoFermion pf( Op_DHD, f_DH, f_mgrad_DHD, lattice );

  // Timer timer;

  // ------------------

  // Idx il=1;
  // Link ell = lattice.links[il];
  // std::cout << "debug. ell = " << ell[0] << " " << ell[1] << std::endl;

  // const double eps = 1.0e-5;
  // Gauge UP(U);
  // UP[il] += eps;
  // Gauge UM(U);
  // UM[il] -= eps;

  // std::cout << " --- Dov.update : " << timer.currentSeconds() << std::endl;
  // D.update(U);
  // std::cout << " --- pf.gen : " << timer.currentSeconds() << std::endl;
  // pf.gen( rng );

  // std::cout << " --- grad constructor : " << timer.currentSeconds() << std::endl;
  // Force grad(lattice);

  // std::cout << " --- pre calc : " << timer.currentSeconds() << std::endl;
  // D.precalc_grad_deviceAsyncLaunch( U, pf.d_eta );
  // std::cout << " --- get force : " << timer.currentSeconds() << std::endl;
  // pf.get_force( grad, U );

  // std::cout << " --- fin : " << timer.currentSeconds() << std::endl;

  // std::cout << "grad = " << grad[il] << std::endl;
  // D.update(UP);
  // pf.update_eta();
  // double sfp = pf.S();

  // D.update(UM);
  // pf.update_eta();
  // double sfm = pf.S();

  // double chck = (sfp-sfm)/(2.0*eps);
  // std::cout << "check = " << chck << std::endl;

  // -----------------


//   // const double eps = 1.0e-5;

//   // std::cout << " --- Dov.update : " << timer.currentSeconds() << std::endl;
//   // Dov.update(U);
//   // std::cout << " --- pf.gen : " << timer.currentSeconds() << std::endl;
//   // pf.gen( rng );

//   // std::cout << " --- grad constructor : " << timer.currentSeconds() << std::endl;
//   // Force dSf(lattice);
//   // std::cout << " --- pre calc : " << timer.currentSeconds() << std::endl;
//   // Dov.precalc_grad_deviceAsyncLaunch( U, pf.d_eta );
//   // std::cout << " --- get force : " << timer.currentSeconds() << std::endl;
//   // pf.get_force( dSf, U );
//   // std::cout << " --- fin : " << timer.currentSeconds() << std::endl;


//   // for(Idx il=0; il<lattice.n_links; il++) std::cout << "grad = " << il << " " << dSf[il] << std::endl;

//   // const double tmax = 0.5; // 1.0; // 0.1
//   // const int nsteps=5;
//   // ExplicitLeapfrogML integrator( tmax, nsteps, 10 );


//   // Force pi( lattice );
//   // pi.gaussian( rng );
//   // // Force pi0=pi;

//   // for(Idx il=0; il<lattice.n_links; il++){
//   //   //   Idx il=3;
//   //   // Link ell = lattice.links[il];

//   //   Gauge UP(U);
//   //   UP[il] += eps;
//   //   Gauge UM(U);
//   //   UM[il] -= eps;


//   //   double Hp, Hm;
//   //   {
//   //     HMC hmc(rng, &SW, &Dov, UP, pi, &pf, &integrator);
//   //     Dov.update(UP);
//   //     pf.update_eta();
//   //     Hp = hmc.H();
//   //   }

//   //   {
//   //     HMC hmc(rng, &SW, &Dov, UM, pi, &pf, &integrator);
//   //     Dov.update(UM);
//   //     pf.update_eta();
//   //     Hm = hmc.H();
//   //   }

//   //   double chck = (Hp-Hm)/(2.0*eps);
//   //   std::cout << "check = " << il << " " << chck << std::endl;
//   // }
//   // // -----------------


//   // Force pi( lattice );
//   // pi.gaussian( rng );
//   // Force pi0=pi;

//   // Gauge U0=U;
//   // D.update(U);
//   // pf.gen( rng );
//   // D.precalc_grad_deviceAsyncLaunch( U, pf.d_eta );

//   // const double tmax = 0.2; // 1.0; // 0.1
//   // for(int nsteps=1; nsteps<=5; nsteps+=1){
//   //   // const int nsteps=5;
//   //   ExplicitLeapfrogML integrator( tmax, nsteps, 10 );
//   //   // ExplicitLeapfrogML integrator( tmax, nsteps, 100 );
//   //   pi = pi0;
//   //   U = U0;
//   //   HMC hmc(rng, &SW, &D, U, pi, &pf, &integrator);
//   //   D.update( U ); pf.update_eta();
//   //   D.precalc_grad_deviceAsyncLaunch( U, pf.d_eta );
//   //   const double h0 = hmc.H();
//   //   hmc.integrate();
//   //   const double h1 = hmc.H();
//   //   double dH = h1-h0;
//   //   std::cout << tmax/nsteps << " " << dH << std::endl;
//   // }



  Force pi(base);
  pi.gaussian( rng );
  Force pi0=pi;
  Gauge U0=U;

  const double tmax = 1.0; // 0.1
  const int nsteps=50;
  pi = pi0;
  U = U0;
  HMCPureGauge hmc(rng, &SW, U, pi, tmax, nsteps);

  double rate, dH;
  bool is_accept;
  for(int k=0; k<10; k++){
    Timer timer;
    hmc.run( rate, dH, is_accept, true );
    if constexpr(Comp::is_compact) U.project();
    std::cout << "# dH : " << dH
              << " is_accept : " << is_accept << std::endl;
    // std::cout << "# HMC : " << timer.currentSeconds() << " sec" << std::endl;
  }
  const int kmax=4e6;
  for(int k=0; k<1e4; k++){
  // const int kmax=4000;
  // for(int k=0; k<100; k++){
    Timer timer;
    hmc.run( rate, dH, is_accept );
    if constexpr(Comp::is_compact) U.project();
    std::cout << "# dH : " << dH
              << " is_accept : " << is_accept << std::endl;
    // std::cout << "# HMC : " << timer.currentSeconds() << " sec" << std::endl;
  }



  std::vector<std::vector<double>> plaq_s0(Comp::Nt);
  std::vector<std::vector<double>> plaq_s1(Comp::Nt);
  std::vector<double> polyakov_s0(Comp::Nt);
  std::vector<double> polyakov_s1(Comp::Nt);
  // std::vector<double> plaq_t0;
  // const int ix0 = 0;
  const int il0 = 1;

#ifdef IS_DUAL
  int iface0 = 2; // 2,3
  int iface1 = 3; // 2,3
#else
  int iface0 = 2; // 2,3
  int iface1 = 3; // 2,3
#endif
  // if(argc==2) iface0 = atoi( argv[1] );

  // std::vector<double> 1;
  // const int s=0;

  double r_mean;
  const int interval=20;

  for(int k=0; k<kmax; k++){
    Timer timer;
    hmc.run( rate, dH, is_accept);
    if constexpr(Comp::is_compact) U.project();
    std::cout << "# dH : " << dH
              << " is_accept : " << is_accept << std::endl;
    r_mean += rate;
    // std::cout << "# HMC : " << timer.currentSeconds() << " sec" << std::endl;

    if(k%interval==0){
      // std::vector<double> tmp1(Comp::Nt, 0.0);
#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL3)
#endif
      for(int t=0; t<Comp::Nt; t++){
        int counter = 0;
        double tmp1 = 0.0;
        for(int s=0; s<Comp::Nt; s++){
          for(int i_face=0; i_face<base.n_faces; i_face++){
            if( std::abs( base.vols[i_face]-base.vols[iface0] )>1.0e-12 ) continue;
            tmp1 += std::pow( U.plaquette_angle(s, U.lattice.faces[i_face]), 1) * std::pow( U.plaquette_angle(s+t, U.lattice.faces[i_face]), 1);
            counter++;
            tmp1 += std::pow( U.plaquette_angle(s, U.lattice.faces[i_face]), 1) * std::pow( U.plaquette_angle(s-t, U.lattice.faces[i_face]), 1);
            counter++;
          }
        }
        tmp1 /= counter;
        plaq_s0[t].push_back( tmp1 );
      }
      polyakov_s0.push_back( std::real( get_polyakov(U, 0) ) );

#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL3)
#endif
      for(int t=0; t<Comp::Nt; t++){
        int counter = 0;
        double tmp1 = 0.0;
        for(int s=0; s<Comp::Nt; s++){
          for(int i_face=0; i_face<base.n_faces; i_face++){
            if( std::abs( base.vols[i_face]-base.vols[iface1] )>1.0e-12 ) continue;
            tmp1 += std::pow( U.plaquette_angle(s, U.lattice.faces[i_face]), 1) * std::pow( U.plaquette_angle(s+t, U.lattice.faces[i_face]), 1);
            counter++;
            tmp1 += std::pow( U.plaquette_angle(s, U.lattice.faces[i_face]), 1) * std::pow( U.plaquette_angle(s-t, U.lattice.faces[i_face]), 1);
            counter++;
          }
        }
        tmp1 /= counter;
        plaq_s1[t].push_back( tmp1 );
      }
      polyakov_s1.push_back( std::real( get_polyakov(U, 1) ) );

    }
    if(k%100==0){
      std::cout << "# k = " << k << std::endl;
    }
  }
  r_mean /= kmax;
  std::cout << "# r_mean = " << r_mean << std::endl;

  // int counter1 = 0;
  // for(int i_face=0; i_face<base.n_faces; i_face++){
  //   if( std::abs( base.vols[i_face]-base.vols[iface0] )>1.0e-12 ) continue;
  //   counter1++;
  // }
  {
    std::vector<double> mean(Comp::Nt, 0.0), var(Comp::Nt, 0.0);
#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL3)
#endif
    for(int t=0; t<Comp::Nt; t++){
      for(const double elem : plaq_s0[t]) mean[t] += elem;
      mean[t] /= plaq_s0[t].size();
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL3)
#endif
    for(int t=0; t<Comp::Nt; t++){
      for(const double elem : plaq_s1[t]) var[t] += std::pow( elem-mean[t], 2);
      var[t] /= std::pow( plaq_s0[t].size(), 2 );
    }

    std::string path = "plaq_s0.dat";
#ifdef IS_DUAL
    path = "dual_"+path;
#endif
    std::ofstream ofs(path);
    ofs << "# kmax = " << kmax << std::endl;
    for(int t=0; t<Comp::Nt; t++){
      ofs << mean[t] << " " << std::sqrt(var[t]) << " " << base.vols[iface0] << std::endl;
    }
  }
  {
    std::vector<double> mean(Comp::Nt, 0.0), var(Comp::Nt, 0.0);
#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL3)
#endif
    for(int t=0; t<Comp::Nt; t++){
      for(const double elem : plaq_s1[t]) mean[t] += elem;
      mean[t] /= plaq_s1[t].size();
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL3)
#endif
    for(int t=0; t<Comp::Nt; t++){
      for(const double elem : plaq_s1[t]) var[t] += std::pow( elem-mean[t], 2);
      var[t] /= std::pow( plaq_s1[t].size(), 2 );
    }

    std::string path = "plaq_s1.dat";
#ifdef IS_DUAL
    path = "dual_"+path;
#endif
    std::ofstream ofs(path);
    ofs << "# kmax = " << kmax << std::endl;
    for(int t=0; t<Comp::Nt; t++){
      ofs << mean[t] << " " << std::sqrt(var[t]) << " " << base.vols[iface1] << std::endl;
    }
  }

  {
    double mn=0.0, vr=0.0;
    for(const double elem : polyakov_s0) mn += elem;
    mn /= polyakov_s0.size();
    for(const double elem : polyakov_s0) vr += std::pow(elem-mn, 2);
    vr /= std::pow( polyakov_s0.size(), 2 );
    std::cout << mn << " " << std::sqrt(vr) << " " << vr << std::endl;
  }

  {
    double mn=0.0, vr=0.0;
    for(const double elem : polyakov_s1) mn += elem;
    mn /= polyakov_s1.size();
    for(const double elem : polyakov_s1) vr += std::pow(elem-mn, 2);
    vr /= std::pow( polyakov_s1.size(), 2 );
    std::cout << mn << " " << std::sqrt(vr) << " " << vr << std::endl;
  }

  // CUDA_CHECK(cudaDeviceReset());
  return 0;

}

