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


// #define IsVerbose
// #define IsVerbose2
// #define InfoForce
#define InfoDelta


namespace Comp{
  constexpr bool is_compact=false;

  // overlap only (no IS_DUAL)
  constexpr int NPARALLEL_DUPDATE=4;
  constexpr int NPARALLEL=NPARALLEL_DUPDATE;
  constexpr int NSTREAMS=NPARALLEL_DUPDATE;

  constexpr int NPARALLEL_GAUGE=16;
  constexpr int NPARALLEL_SORT=16;

  constexpr int N_REFINE=1;
  constexpr int NS=2;

  // constexpr int Nt=96; // @@@
  constexpr int Nt=128; // @@@
  // constexpr int Nt=16;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;

  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-8;
}


const std::string dir = "/project/affine/nmatsum/qed3/geometry/data/";
#include "/project/affine/nmatsum/qed3/geometry/geodesic.h"

#include "timer.h"

#include "s2n_simp.h"
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

#include "sparse_matrix.h"

#include "dirac_simp.h"
#include "dirac_ext.h"

#include "sparse_dirac.h"
// #include "matpoly.h"
#include "matpoly_claude.h"
#include "includes/overlap_wmass_claude.h"
#include "pseudofermion.h"

#include "integrator.h"
#include "hmc.h"


int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) == "-h") {
      printf("Usage: ./a.out [gsq] [Nf] [nu0] [mass_re] [mass_im]\n");
      printf("  gsq      Wilson coupling squared (default: 8.0)\n");
      printf("  Nf       number of fermion flavors (default: 2)\n");
      printf("  nu0      mass parameter (default: 1.0)\n");
      printf("  mass_re  real part of additive mass (default: 0.0)\n");
      printf("  mass_im  imaginary part of additive mass (default: 0.0)\n");
      return 0;
    }
  }

  double gsq = 8.0;
  if(argc>1) gsq = atof(argv[1]);
  int Nf = 2;
  if(argc>2) Nf = atoi(argv[2]);
  double nu0 = 1.0;
  if(argc>3) nu0 = atof(argv[3]);
  double mass_re = 0.0;
  if(argc>4) mass_re = atof(argv[4]);
  double mass_im = 0.0;
  if(argc>5) mass_im = atof(argv[5]);
  Complex mass = Complex(mass_re, mass_im);
  std::cout << "# gsq = " << gsq << " Nf = " << Nf << " nu0 = " << nu0 << " mass = " << mass << std::endl;


  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();

  // ---------------------------------------
  using BaseLink = std::array<Idx,2>;
  constexpr Idx N = Comp::N;
  constexpr int Nt = Comp::Nt;

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Force=GaugeExt<Base,Nt,Comp::is_compact>;
  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  using Action=U1WilsonExt<Base>;
  using Rng=ParallelRngExt<Base,Nt>;
  using Fermion=OverlapWMass<WilsonDirac>;

  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;

  // ----------------------

  const double r = 1.0;
  const double M5 = -1.0;
  const double at = 0.2;
  assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);

  std::cout << "# DW set" << std::endl;

  Gauge U(base);
  Rng rng(base);

  // ---------------------

  Fermion D(DW, mass, 21);
  std::cout << "# Dov set; M5 = " << M5 << std::endl;
  D.update(U);
  std::cout << "# min max ratio: "
            << D.lambda_min << " "
            << D.lambda_max << " "
            << D.lambda_min/D.lambda_max << std::endl;
  std::cout << "# delta = " << D.Delta() << std::endl;

  // -----------------------------------------------------------

  Action SW( gsq, at, base );
  std::cout << "# alat = " << base.mean_ell << std::endl;

  std::vector<std::shared_ptr<PseudoFermion<Fermion>>> pfs;
  assert(Nf%2==0);
  for(int f=0; f<Nf/2; f++) pfs.push_back( std::shared_ptr<PseudoFermion<Fermion>>( new PseudoFermion<Fermion>(D) ) );

  Timer timer;


  // -----------------

  std::string dir3;
  dir3="Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"mRe"+std::to_string(mass.real())+"mIm"+std::to_string(mass.imag())+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
  std::filesystem::create_directory(dir3);
  const int k_ckpoint=1;

  const int k_ckpoint_rng=100;
  const int kmax=1000;

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

    if(k_tmp>0){
      std::cout << "read from k_tmp = " << k_tmp << std::endl;
      const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k_tmp);
      const std::string str_rng=dir3+"ckpoint_rng."+std::to_string(k_tmp);
      U.read( str_lat );
      rng.read( str_rng );
    }
  }
  std::cout << "#starting from k_tmp = " << k_tmp << std::endl;


  Force pi( base );
  const double tmax = 1.9;
  int nsteps;
  // 2026-06-02 15:03: bumped +3 (Nf=2: 4->7, Nf=4,6: 5->8) to reduce discretization error after Nf=4,6 runs stuck at 100% rejection
  // 2026-06-04 10:42: bumped again +3 (Nf=2: 7->8, Nf=4,6: 5->10) to reduce discretization error after Nf=4,6 runs stuck at 100% rejection
  if(Nf==2) nsteps = 8;
  else if(Nf==4) nsteps = 10;
  else if(Nf==6) nsteps = 10;
  else nsteps = 10;
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
