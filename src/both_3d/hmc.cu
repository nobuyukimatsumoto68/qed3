#include <typeinfo>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>

#include <algorithm>
#include <filesystem>
#include <chrono>


#include <cstdint>
#include <complex>

#include <array>
#include <vector>
#include <map>
#include <Eigen/Dense>

using Double = double;
using Idx = std::int32_t;
using Complex = std::complex<double>;

using MS=Eigen::Matrix2cd;
using VD=Eigen::Vector2d;
using VE=Eigen::Vector3d;
using VC=Eigen::VectorXcd;

static constexpr int NS = 2;
static constexpr int DIM = 2;
static constexpr Complex I = Complex(0.0, 1.0);

// #define Nf2
#define Nf4

// #define IS_DUAL
#define IS_OVERLAP

// #define IsVerbose
// #define InfoForce
#define InfoDelta


namespace Comp{
  constexpr bool is_compact=false;

  // d_DW.update() is always done independently
#ifdef IS_OVERLAP
  constexpr int NPARALLEL_DUPDATE=1;
  constexpr int NPARALLEL=16; // 12
  constexpr int NSTREAMS=2; // 4
#else
  constexpr int NPARALLEL_DUPDATE=12;
  constexpr int NPARALLEL=1; // 12
  constexpr int NSTREAMS=12; // for grad loop
#endif
  constexpr int NPARALLEL_GAUGE=16; // 12
  constexpr int NPARALLEL_SORT=16; // 12

  constexpr int N_REFINE=1;
  constexpr int NS=2;

  constexpr int Nt=24;

  // constexpr int Nf=4; // even

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

const std::string dir = "../../dats/";

#include "timer.h"

#include "s2n_simp.h"
#include "s2n_dual.h"
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

#include "sparse_matrix.h"

#include "dirac_simp.h"
#include "dirac_dual.h"
#include "dirac_ext.h"

#include "sparse_dirac.h"
#include "matpoly.h"
#include "dirac_pf.h"
#include "overlap.h"
#include "pseudofermion.h"

# include "integrator.h"
#include "hmc.h"

// #include "obs.h" // to be developed



int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  int device;
  CUDA_CHECK(cudaGetDeviceCount(&device));
  cudaDeviceProp device_prop[device];
  cudaGetDeviceProperties(&device_prop[0], 0);
  std::cout << "# dev = " << device_prop[0].name << std::endl;
  // CUDA_CHECK(cudaSetDevice(0));// "TITAN V"
  std::cout << "# (GPU device is set.)" << std::endl;

  // ---------------------------------------
  using BaseLink = std::array<Idx,2>; // <int,int>;
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
  using Action=U1WilsonExt<Base>;
  using Rng=ParallelRngExt<Base,Nt>;


  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;

  // ----------------------

#ifdef IS_OVERLAP
  const double r = 1.0;
  const double M5 = -1.8; // -1.6/2.0 * 0.5*(1.0 + std::sqrt( 5.0 + 2.0*std::sqrt(2.0) ));
  using Fermion=Overlap<WilsonDirac>;
#else
  const double r = 1.0;
  const double M5 = 0.0;
  using Fermion=DiracPf<WilsonDirac>;
#endif
  // const double c = 1.0;
  // double at = 0.05; // base.mean_ell * 0.125 * ratio;
  const double T = 16;
  const double at = T/Comp::Nt;
  assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);
  WilsonDirac DW(base, 0.0, 1.0, M5, at);


  std::cout << "# DW set" << std::endl;

  Gauge U(base);
  Rng rng(base);
  // U.gaussian( rng, 0.01 );

  // ---------------------

#ifdef IS_OVERLAP
  Fermion D(DW, 31);
  std::cout << "# Dov set; M5 = " << M5 << std::endl;
  D.update(U);
  std::cout << "# min max ratio: "
            << D.lambda_min << " "
            << D.lambda_max << " "
            << D.lambda_min/D.lambda_max << std::endl;
  std::cout << "# delta = " << D.Delta() << std::endl;

#else
  Fermion D(DW);
  D.update( U );

#endif


  // -----------------------------------------------------------

  const double gsq = 0.1;
  // const double beta = 1.0/(gR*gR);
  Action SW( gsq, at, base );
  std::cout << "# alat = " << base.mean_ell << std::endl;

  std::vector<PseudoFermion<Fermion>*> pfs;
  //if(Comp::Nf>=2){
#ifdef Nf2
  PseudoFermion pf1(D);
  pfs.push_back(&pf1);
#else
#ifdef Nf4
  PseudoFermion pf1(D);
  pfs.push_back(&pf1);
  PseudoFermion pf2(D);
  pfs.push_back(&pf2);
#else
  static_assert(false, "inappropriate Nf");
#endif
#endif
  // }
  // if(Comp::Nf>=4){
  // }

  Timer timer;

  // ------------------

  // {
  //   int s=Nt-1;
  //   Idx il=4;
  //   BaseLink ell = base.links[il];
  //   std::cout << "debug. ell = " << ell[0] << " " << ell[1] << std::endl;

  //   const double eps = 1.0e-5;
  //   Gauge UP(U);
  //   UP.sp(s,il) += eps;
  //   Gauge UM(U);
  //   UM.sp(s,il) -= eps;

  //   std::cout << " --- Dov.update : " << timer.currentSeconds() << std::endl;
  //   D.update(U);
  //   std::cout << " --- pf.gen : " << timer.currentSeconds() << std::endl;
  //   pf.gen( rng );

  //   std::cout << " --- grad constructor : " << timer.currentSeconds() << std::endl;
  //   Force grad(base);

  //   std::cout << " --- pre calc : " << timer.currentSeconds() << std::endl;
  //   D.precalc_grad_deviceAsyncLaunch( U, pf.d_eta );
  //   std::cout << " --- get force : " << timer.currentSeconds() << std::endl;
  //   pf.get_force( grad, U );

  //   std::cout << " --- fin : " << timer.currentSeconds() << std::endl;

  //   std::cout << "grad = " << grad.sp(s,il) << std::endl;
  //   D.update(UP);
  //   pf.update_eta();
  //   double sfp = pf.S();

  //   D.update(UM);
  //   pf.update_eta();
  //   double sfm = pf.S();

  //   double chck = (sfp-sfm)/(2.0*eps);
  //   std::cout << "check = " << chck << std::endl;
  // }

  // // // -----------------

  // {
  //   int s=Nt-1;
  //   Idx ix=4;

  //   const double eps = 1.0e-5;
  //   Gauge UP(U);
  //   UP.tp(s,ix) += eps;
  //   Gauge UM(U);
  //   UM.tp(s,ix) -= eps;

  //   std::cout << " --- Dov.update : " << timer.currentSeconds() << std::endl;
  //   D.update(U);
  //   std::cout << " --- pf.gen : " << timer.currentSeconds() << std::endl;
  //   pf.gen( rng );

  //   std::cout << " --- grad constructor : " << timer.currentSeconds() << std::endl;
  //   Force grad(base);

  //   std::cout << " --- pre calc : " << timer.currentSeconds() << std::endl;
  //   D.precalc_grad_deviceAsyncLaunch( U, pf.d_eta );
  //   std::cout << " --- get force : " << timer.currentSeconds() << std::endl;
  //   pf.get_force( grad, U );

  //   std::cout << " --- fin : " << timer.currentSeconds() << std::endl;

  //   std::cout << "grad = " << grad.tp(s, ix) << std::endl;
  //   D.update(UP);
  //   pf.update_eta();
  //   double sfp = pf.S();

  //   D.update(UM);
  //   pf.update_eta();
  //   double sfm = pf.S();

  //   double chck = (sfp-sfm)/(2.0*eps);
  //   std::cout << "check = " << chck << std::endl;
  // }

  // -----------------


  // const double eps = 1.0e-5;

  // std::cout << " --- D.update : " << timer.currentSeconds() << std::endl;
  // D.update(U);
  // std::cout << " --- pf.gen : " << timer.currentSeconds() << std::endl;
  // pf.gen( rng );

  // std::cout << " --- grad constructor : " << timer.currentSeconds() << std::endl;
  // Force dSf(base);
  // std::cout << " --- pre calc : " << timer.currentSeconds() << std::endl;
  // D.precalc_grad_deviceAsyncLaunch( U, pf.d_eta );
  // std::cout << " --- get force : " << timer.currentSeconds() << std::endl;
  // pf.get_force( dSf, U );
  // std::cout << " --- fin : " << timer.currentSeconds() << std::endl;


  // const double tmax = 0.5; // 1.0; // 0.1
  // const int nsteps=5;
  // ExplicitLeapfrogML integrator( tmax, nsteps, 10 );

  // // for(Idx il=0; il<base.n_links; il++) std::cout << "grad = " << il << " " << dSf[il] << std::endl;

  // Force pi( base );
  // pi.gaussian( rng );
  // // Force pi0=pi;

  // int s=0;
  // for(Idx il=0; il<base.n_links; il++){
  //   //   Idx il=3;
  //   // Link ell = base.links[il];

  //   Gauge UP(U);
  //   UP.sp(s,il) += eps;
  //   Gauge UM(U);
  //   UM.sp(s,il) -= eps;

  //   double Hp, Hm;
  //   {
  //     HMC hmc(rng, &SW, &D, UP, pi, &pf, &integrator);
  //     D.update(UP);
  //     pf.update_eta();
  //     Hp = hmc.H();
  //   }

  //   {
  //     HMC hmc(rng, &SW, &D, UM, pi, &pf, &integrator);
  //     D.update(UM);
  //     pf.update_eta();
  //     Hm = hmc.H();
  //   }

  //   double chck = (Hp-Hm)/(2.0*eps);
  //   std::cout << "check = " << il << " " << chck << " " << dSf.sp(s,il) << std::endl;
  // }
  // return 1;

  // -----------------


  // Force pi( base );
  // pi.gaussian( rng );
  // Force pi0=pi;

  // Gauge U0=U;
  // D.update(U);

  // for(PseudoFermion<Fermion>* pf : pfs){
  //   // pf.gen( rng );
  //   // D.precalc_grad_deviceAsyncLaunch( U, pf.d_eta );
  //   pf->gen( rng );
  //   D.precalc_grad_deviceAsyncLaunch( U, pf->d_eta );
  // }

  // const double tmax = 1.0; // 1.0; // 0.1
  // for(int nsteps=12; nsteps<=24; nsteps+=4){
  //   // const int nsteps=5;
  //   // ExplicitLeapfrogML integrator( tmax, nsteps, 20 );
  //   ExplicitLeapfrogML integrator( tmax, nsteps, 100 );
  //   pi = pi0;
  //   U = U0;
  //   HMC2 hmc(rng, &SW, &D, U, pi, pfs, &integrator);

  //   D.update( U );
  //   for(PseudoFermion<Fermion>* pf : pfs){
  //     // pf.update_eta();
  //     // D.precalc_grad_deviceAsyncLaunch( U, pf.d_eta );
  //     pf->update_eta();
  //     D.precalc_grad_deviceAsyncLaunch( U, pf->d_eta );
  //   }

  //   const double h0 = hmc.H();
  //   hmc.integrate();
  //   const double h1 = hmc.H();
  //   double dH = h1-h0;
  //   std::cout << tmax/nsteps << " " << dH << std::endl;
  //   std::cout << " --- hmc : " << timer.currentSeconds() << std::endl;
  // }




  // -----------------// -----------------// -----------------// -----------------// -----------------

  std::string dir3;
#ifdef Nf2
  dir3="Nf2_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
#else
#ifdef Nf4
  dir3="Nf4_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
#endif
#endif
  std::filesystem::create_directory(dir3);
  const int k_ckpoint=10;
  const int kmax=1e5;


  {
    int k_tmp=0;
    for(k_tmp=k_ckpoint; k_tmp<=kmax; k_tmp+=k_ckpoint ){
      const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k_tmp);
      const std::string str_rng=dir3+"ckpoint_rng."+std::to_string(k_tmp);

      const bool bool_lat = std::filesystem::exists(str_lat);
      const bool bool_rng = std::filesystem::exists(str_rng);

      if(!(bool_lat&&bool_rng)) break;
    }
    k_tmp -= k_ckpoint;

    if(k_tmp>0){ // from existing
      std::cout << "read from k_tmp = " << k_tmp << std::endl;
      const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k_tmp);
      const std::string str_rng=dir3+"ckpoint_rng."+std::to_string(k_tmp);
      U.read( str_lat );
      rng.read( str_rng );
    }
  }


  Force pi( base );
  pi.gaussian( rng );
  Force pi0=pi;

  Gauge U0=U;
  D.update(U);

  for(PseudoFermion<Fermion>* pf : pfs){
    pf->gen( rng );
    D.precalc_grad_deviceAsyncLaunch( U, pf->d_eta );
  }

  const double tmax = 1.0; // 0.1
  const int nsteps=20;
  ExplicitLeapfrogML integrator( tmax, nsteps, 100 );
  pi = pi0;
  U = U0;
  HMC2 hmc(rng, &SW, &D, U, pi, pfs, &integrator);
  D.update( U );

  for(PseudoFermion<Fermion>* pf : pfs){
    pf->gen( rng );
    D.precalc_grad_deviceAsyncLaunch( U, pf->d_eta );
  }

  double rate, dH;
  bool is_accept;
  for(int k=0; k<10; k++){
    Timer timer;
    hmc.run( rate, dH, is_accept, true);
    std::cout << "# dH : " << dH
              << " is_accept : " << is_accept
              << " rate : " << rate << std::endl;
    std::cout << "# HMC : " << timer.currentSeconds() << " sec" << std::endl;
  }

  double r_mean;
  for(int k=0; k<kmax; k++){
    Timer timer;
    hmc.run( rate, dH, is_accept);
    std::cout << "# dH : " << dH
              << " is_accept : " << is_accept
              << " rate : " << rate << std::endl;
    r_mean += rate;
    std::cout << "# HMC : " << timer.currentSeconds() << " sec" << std::endl;

    if(k%100==0){
      std::cout << "# k = " << k << std::endl;
    }

    if(k%k_ckpoint==0){
      const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k);
      const std::string str_rng=dir3+"ckpoint_rng."+std::to_string(k);
      U.ckpoint( str_lat );
      rng.ckpoint( str_rng );
    }
  }
  r_mean /= kmax;
  std::cout << "# r_mean = " << r_mean << std::endl;



  // CUDA_CHECK(cudaDeviceReset());
  return 0;

}

