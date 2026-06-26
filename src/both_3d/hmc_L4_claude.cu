#include <typeinfo>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>

#include <algorithm>
#include <filesystem>
#include <thread>
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

using Face = std::vector<Idx>;


using MS=Eigen::Matrix2cd;
using VD=Eigen::Vector2d;
using VE=Eigen::Vector3d;
using VC=Eigen::VectorXcd;

static constexpr int NS = 2;
static constexpr int DIM = 2;
static constexpr Complex I = Complex(0.0, 1.0);


// #define IS_DUAL
#define IS_OVERLAP

// #define IsVerbose
// #define IsVerbose2
// #define InfoForce
#define InfoDelta


namespace Comp{
  constexpr bool is_compact=false;
  // constexpr bool is_compact=true;

  // d_DW.update() is always done independently
#ifdef IS_OVERLAP
  constexpr int NPARALLEL_DUPDATE=1;
  constexpr int NPARALLEL=NPARALLEL_DUPDATE; // 12
  constexpr int NSTREAMS=NPARALLEL_DUPDATE; // 4
#else
  constexpr int NPARALLEL_DUPDATE=12;
  constexpr int NPARALLEL=1; // 12
  constexpr int NSTREAMS=12; // for grad loop
#endif
  constexpr int NPARALLEL_GAUGE=1; // 12
  constexpr int NPARALLEL_SORT=1; // 12

  constexpr int N_REFINE=4; // L=4 massless SEA run (copy of hmc_claude.cu)
  constexpr int NS=2;

  // constexpr int Nt=96; // @@@
  constexpr int Nt=128; // @@@
  // constexpr int Nt=16;

  // constexpr int Nf=4; // even

#ifdef IS_DUAL
  constexpr Idx N_SITES=20*N_REFINE*N_REFINE;
#else
  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE; // 30, 120, 480
#endif

  constexpr Idx Nx=NS*N_SITES; // matrix size of DW
  constexpr Idx N=Nx*Nt; // matrix size of DW

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-8;
}

// const std::string dir = "../../dats/";
const std::string dir = "../../geometry/data/";
// #include "../../integrator/geodesic.h"
#include "../../geometry/geodesic.h"

#include "timer.h"

#include "s2n_simp.h"
// #include "s2n_dual.h"
#include "rng.h"

#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cublas_api.h>
#include <cusolverDn.h>
using CuC = cuDoubleComplex;
#include "gpu_header.h"


#include "valence.h"
#include "gauge_ext.h"
#include "action_ext.h"


// ======================================

// #include "sparse_matrix.h"             // _ms: multishift copy below
#include "sparse_matrix_claude.h"

#include "dirac_simp.h"
// #include "dirac_dual.h"
#include "dirac_ext.h"

// #include "sparse_dirac.h"
#include "sparse_dirac_claude.h"   // O(len) bucketing CSR build (was O(N*len)); -DCSR_VERIFY to check
// #include "matpoly.h"
#include "matpoly_claude.h"
// #include "dirac_pf.h"
// #include "overlap.h"                    // _ms: OverlapWMass copy below
#define GRAD_L4   // HMC force opt L1+L2+L4 (hoist X Z_m/X Y_m + block poles + skip do_it); force==ref ~1e-16 (test PASS, ~3.4x grad)
#include "overlap_wmass_claude.h"
// #include "pseudofermion.h"              // _ms: multishift copy below
#include "pseudofermion_claude.h"

# include "integrator.h"
#include "hmc.h"

// #include "obs.h" // to be developed


int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) == "-h") {
      printf("Usage: ./a.out [gsq] [Nf] [nu0]\n");
      printf("  gsq  Wilson coupling squared (default: 8.0)\n");
      printf("  Nf   number of fermion flavors (default: 2)\n");
      printf("  nu0  mass parameter (default: 1.0)\n");
      return 0;
    }
  }

  // double gsq = 8.0;
  double gsq = 8.0;
  if(argc>1) gsq = atof(argv[1]);
  int Nf = 2;
  if(argc>2) Nf = atoi(argv[2]);
  // double nu0=2.0;
  double nu0=1.0;
  if(argc>3) nu0 = atof(argv[3]);
  std::cout << "# gsq = " << gsq << " Nf = " << Nf << " nu0 = " << nu0 << std::endl;


  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();

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
  const double M5 = -1.0; // -1.6/2.0 * 0.5*(1.0 + std::sqrt( 5.0 + 2.0*std::sqrt(2.0) ));
  // using Fermion=Overlap<WilsonDirac>;             // _ms: OverlapWMass below
  using Fermion=OverlapWMass<WilsonDirac>;           // massless (mass=0), enables _ms solver
#else
  const double r = 1.0;
  const double M5 = 0.0;
  using Fermion=DiracPf<WilsonDirac>;
#endif
  // const double c = 1.0;
  // double at = 0.05; // base.mean_ell * 0.125 * ratio;
  // const double T = 24;
  const double at = 0.2; // T/Comp::Nt;
  assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);


  std::cout << "# DW set" << std::endl;

  Gauge U(base);
  Rng rng(base);
  // U.gaussian( rng, 0.01 );

  // ---------------------

  // HERE
#ifdef IS_OVERLAP
  const Complex mass = Complex(0.0, 0.0); // massless
  // Fermion D(DW, 11); // 21                         // _ms: OverlapWMass ctor below
  Fermion D(DW, mass, 11);      // n=11 (5 poles) -- reverted per request (hmc_claude.cu ONLY; jj keeps n=21)
  // Fermion D(DW, mass, 21);   // n=21 (10 poles): unified pole count (multishift-validated, freeze-safe)
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


  // const double beta = 1.0/(gR*gR);
  Action SW( gsq, at, base );
  std::cout << "# alat = " << base.mean_ell << std::endl;

  std::vector<std::shared_ptr<PseudoFermion<Fermion>>> pfs;
  assert(Nf%2==0);
  for(int f=0; f<Nf/2; f++) pfs.push_back( std::shared_ptr<PseudoFermion<Fermion>>( new PseudoFermion<Fermion>(D) ) );

  Timer timer;





  // -----------------// -----------------// -----------------// -----------------// -----------------

  std::string dir3;
  dir3="Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
  // dir3="Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"C"+"/";
  std::filesystem::create_directory(dir3);
  // const int k_ckpoint=1;
  const int k_ckpoint=1;
  const int k_ckpoint_rng=1; // 1000; // L=4: keep rng every conf (match hmc_fermilab_wmass_L4; 1000=rolling-latest blocked clean rollback)
  const int kmax=4e3; // @@@@
  // const int kmax=2;

  int k_tmp=0;
  {
    for(int k_scan=k_ckpoint; k_scan<=kmax; k_scan+=k_ckpoint ){
      const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k_scan);
      const std::string str_rng=dir3+"ckpoint_rng."+std::to_string(k_scan);

      const bool bool_lat = std::filesystem::exists(str_lat);
      const bool bool_rng = std::filesystem::exists(str_rng);

      if(!bool_lat) break; // no lat: stop scanning
      if(bool_lat && bool_rng) k_tmp = k_scan; // both present: candidate
    }

    if(k_tmp>0){ // from existing
      std::cout << "read from k_tmp = " << k_tmp << std::endl;
      const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k_tmp);
      const std::string str_rng=dir3+"ckpoint_rng."+std::to_string(k_tmp);
      U.read( str_lat );
      rng.read( str_rng );
    }
  }
  std::cout << "#starting from k_tmp = " << k_tmp << std::endl;


  Force pi( base );
  const double tmax = 1.0; // 1.9; // 0.1 -- reduced for L>1 integrator stability (Nf2 step was 0.2375 -> stuck)
  int nsteps;
  // 2026-06-02 15:03: bumped +3 (Nf=2: 4->7, Nf=4,6: 5->8) to reduce discretization error after Nf=4,6 runs stuck at 100% rejection
  // 2026-06-04 10:42: bumped to 2x the original (Nf=2: 4->8, Nf=4,6: 5->10) to reduce discretization error after Nf=4,6 runs stuck at 100% rejection
  // 2026-06-20: L=4 fixed nsteps=6 (tmax=1.0 -> step 0.1667) per request; old Nf-conditional below
  // if(Nf==2) nsteps = 8;
  // else if(Nf==4) nsteps = 10;
  // else if(Nf==6) nsteps = 10;
  // else nsteps = 10;
  // nsteps = 6;
  nsteps = 8; // 2026-06-23: match hmc_fermilab_wmass_L4_claude.cu (tmax=1.0 -> step 0.125)
  std::cout << "# tmax = " << tmax << std::endl
            << "# nsteps = " << nsteps << std::endl;

  MinimumNorm2 integrator( tmax, nsteps, 100 );
  HMC2 hmc(rng, &SW, &D, U, pi, pfs, &integrator);
  D.update( U );

  double rate, dH;
  bool is_accept;

  double r_mean;
  for(int k=k_tmp+1; k<kmax; k++){
    Timer timer;
    hmc.run( rate, dH, is_accept);
    std::cout << "# dH : " << dH
              << " is_accept : " << is_accept
              << " rate : " << rate << std::endl;
    r_mean += rate;
    std::cout << "# HMC : " << timer.currentSeconds() << " sec" << std::endl;

    if(k%20==0) D.is_update = false;
    if(k%100==0){
      std::cout << "# k = " << k << std::endl;
    }

    if(k%k_ckpoint==0){
      const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k);
      const std::string str_rng=dir3+"ckpoint_rng."+std::to_string(k);
      U.ckpoint( str_lat );
      rng.ckpoint( str_rng );
      int k_prev = k - k_ckpoint;
      if(k_prev > 0 && k%k_ckpoint_rng != 0){
        std::error_code ec;
        std::filesystem::remove(dir3+"ckpoint_rng."+std::to_string(k_prev), ec);
        if(ec){
          std::cout << "# error removing ckpoint_rng." << k_prev << ": " << ec.message() << std::endl;
          assert(false);
        }
      }
    }
  }
  r_mean /= kmax;
  std::cout << "# r_mean = " << r_mean << std::endl;

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();

  return 0;

}





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
  // Force pi1=pi;

  // Gauge U0=U;
  // D.update(U);

  // for(auto pf : pfs){
  //   pf->gen( rng );
  //   D.precalc_grad_deviceAsyncLaunch( U, pf->d_eta );
  // }

  // const double tmax = 1.0; // 1.0; // 0.1
  // for(int nsteps=4; nsteps<=24; nsteps+=4){
  //   // for(int nsteps=8; nsteps<=40; nsteps+=8){
  //   MinimumNorm2 integrator( tmax, nsteps, 100 );
  //   // MinimumNorm2 integrator( tmax, nsteps, 1 );
  //   // ExplicitLeapfrogML integrator( tmax, nsteps, 1 );
  //   HMC2 hmc(rng, &SW, &D, U, pi, pfs, &integrator);

  //   pi = pi0;
  //   U = U0;
  //   D.update( U );
  //   for(auto pf : pfs){
  //     pf->update_eta(); // fixed phi; update eta according to U
  //     D.precalc_grad_deviceAsyncLaunch( U, pf->d_eta );
  //   }

  //   const double h0 = hmc.H();
  //   hmc.integrate();
  //   const double h1 = hmc.H();

  //   double dH = h1-h0;
  //   std::cout << tmax/nsteps << " " << dH << std::endl;
  //   std::cout << "# --- hmc : " << timer.currentSeconds() << std::endl;

  //   // reversibility
  //   // pi *= -1.0;
  //   // hmc.integrate();
  //   // std::cout << "# rev. check = " << (pi+pi0).norm() << std::endl;
  //   // U *= -1.0;
  //   // std::cout << "# rev. check = " << (U+U0).norm() << std::endl;
  // }
