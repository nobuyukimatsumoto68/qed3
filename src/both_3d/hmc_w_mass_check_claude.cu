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

  constexpr int NPARALLEL_DUPDATE=1;
  constexpr int NPARALLEL=4;
  constexpr int NSTREAMS=4;

  constexpr int NPARALLEL_GAUGE=4;
  constexpr int NPARALLEL_SORT=4;

  constexpr int N_REFINE=1;
  constexpr int NS=2;

  constexpr int Nt=16;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;

  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-8;
}

const std::string dir = "../../geometry/data/";
#include "../../geometry/geodesic.h"

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
#include "matpoly.h"
#include "includes/overlap_wmass_claude.h"
#include "pseudofermion.h"


int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

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
  const Complex mass = Complex(mass_re, mass_im);
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
  std::cout << "# lattice set." << std::endl;

  // ----------------------

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

  Action SW(gsq, at, base);
  std::cout << "# alat = " << base.mean_ell << std::endl;

  PseudoFermion<Fermion> pf(D);

  Timer timer;

  // --- spatial link force check: grad.sp vs finite difference of S_PF ---
  {
    int s = Nt-1;
    Idx il = 4;
    BaseLink ell = base.links[il];
    std::cout << "debug. ell = " << ell[0] << " " << ell[1] << std::endl;

    const double eps = 1.0e-5;
    Gauge UP(U);
    UP.sp(s,il) += eps;
    Gauge UM(U);
    UM.sp(s,il) -= eps;

    std::cout << " --- Dov.update : " << timer.currentSeconds() << std::endl;
    D.update(U);
    std::cout << " --- pf.gen : " << timer.currentSeconds() << std::endl;
    pf.gen(rng);

    std::cout << " --- grad constructor : " << timer.currentSeconds() << std::endl;
    Force grad(base);

    std::cout << " --- pre calc : " << timer.currentSeconds() << std::endl;
    D.precalc_grad_deviceAsyncLaunch(U, pf.d_eta);
    std::cout << " --- get force : " << timer.currentSeconds() << std::endl;
    pf.get_force(grad, U);
    std::cout << " --- fin : " << timer.currentSeconds() << std::endl;

    std::cout << "grad = " << grad.sp(s,il) << std::endl;

    D.update(UP);
    pf.update_eta();
    double sfp = pf.S();

    D.update(UM);
    pf.update_eta();
    double sfm = pf.S();

    double chck = (sfp-sfm)/(2.0*eps);
    std::cout << "check = " << chck << std::endl;
  }

  // --- temporal link force check: grad.tp vs finite difference of S_PF ---
  {
    int s = Nt-1;
    Idx ix = 4;

    const double eps = 1.0e-5;
    Gauge UP(U);
    UP.tp(s,ix) += eps;
    Gauge UM(U);
    UM.tp(s,ix) -= eps;

    std::cout << " --- Dov.update : " << timer.currentSeconds() << std::endl;
    D.update(U);
    std::cout << " --- pf.gen : " << timer.currentSeconds() << std::endl;
    pf.gen(rng);

    std::cout << " --- grad constructor : " << timer.currentSeconds() << std::endl;
    Force grad(base);

    std::cout << " --- pre calc : " << timer.currentSeconds() << std::endl;
    D.precalc_grad_deviceAsyncLaunch(U, pf.d_eta);
    std::cout << " --- get force : " << timer.currentSeconds() << std::endl;
    pf.get_force(grad, U);
    std::cout << " --- fin : " << timer.currentSeconds() << std::endl;

    std::cout << "grad = " << grad.tp(s,ix) << std::endl;

    D.update(UP);
    pf.update_eta();
    double sfp = pf.S();

    D.update(UM);
    pf.update_eta();
    double sfm = pf.S();

    double chck = (sfp-sfm)/(2.0*eps);
    std::cout << "check = " << chck << std::endl;
  }

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();

  return 0;
}
