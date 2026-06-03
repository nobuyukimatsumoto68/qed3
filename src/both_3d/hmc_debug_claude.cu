// Debug binary: loads the checkpoint pair (ckpoint_lat.1269 + ckpoint_rng.1269)
// that immediately precedes the k=1271 acceptance crash and re-runs a small
// number of trajectories to confirm reproducibility.
// Writes NO files (no checkpoints, no rng saves, no directory creation).
// D.is_update is forced false after setup to match the original run state
// (was permanently latched false from k=1220 in the crashed run).

#include <typeinfo>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>

#include <algorithm>
#include <filesystem>

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

#define IS_OVERLAP
#define IsVerbose
#define IsVerbose2
#define InfoForce
#define InfoDelta
#define PrintH
#define PrintHinner
#define PrintEtaPhi

namespace Comp{
  constexpr bool is_compact=false;

#ifdef IS_OVERLAP
  constexpr int NPARALLEL_DUPDATE=1;
  constexpr int NPARALLEL=NPARALLEL_DUPDATE;
  constexpr int NSTREAMS=NPARALLEL_DUPDATE;
#else
  constexpr int NPARALLEL_DUPDATE=12;
  constexpr int NPARALLEL=1;
  constexpr int NSTREAMS=12;
#endif
  constexpr int NPARALLEL_GAUGE=1;
  constexpr int NPARALLEL_SORT=1;

  constexpr int N_REFINE=1;
  constexpr int NS=2;

  constexpr int Nt=128;

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

#include "sparse_matrix.h"

#include "dirac_simp.h"
#include "dirac_ext.h"

#include "sparse_dirac.h"
// #include "matpoly.h"
#include "matpoly_claude.h"
#include "overlap.h"
#include "pseudofermion.h"

// #include "integrator.h"
#include "integrator_claude.h"
#include "hmc.h"


int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  // --- fixed parameters matching the crashed run ---
  const double gsq = 8.0;
  const int    Nf  = 4;
  const double nu0 = 1.0;
  std::cout << "# gsq = " << gsq << " Nf = " << Nf << " nu0 = " << nu0 << std::endl;

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();

  constexpr Idx N  = Comp::N;
  constexpr int Nt = Comp::Nt;

  using Base       = S2Simp;
  using WilsonDirac= DiracExt<Base, DiracS2Simp>;
  using Force      = GaugeExt<Base,Nt,Comp::is_compact>;
  using Gauge      = GaugeExt<Base,Nt,Comp::is_compact>;
  using Action     = U1WilsonExt<Base>;
  using Rng        = ParallelRngExt<Base,Nt>;
  using Fermion    = Overlap<WilsonDirac>;

  Base base(Comp::N_REFINE);
  std::cout << "# lattice set." << std::endl;

  const double r  = 1.0;
  const double M5 = -1.0;
  const double at = 0.2;
  assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);
  std::cout << "# DW set" << std::endl;

  Gauge U(base);
  Rng   rng(base);

  Fermion D(DW, 11);
  std::cout << "# Dov set; M5 = " << M5 << std::endl;

  Action SW( gsq, at, base );

  std::vector<std::shared_ptr<PseudoFermion<Fermion>>> pfs;
  assert(Nf%2==0);
  for(int f=0; f<Nf/2; f++)
    pfs.push_back( std::shared_ptr<PseudoFermion<Fermion>>(
                     new PseudoFermion<Fermion>(D) ) );

  // --- hardcoded checkpoint: the last complete pair before the k=1271 crash ---
  // ckpoint_lat.1270 + ckpoint_rng.1269 are both present.
  // (ckpoint_rng.1270 was deleted during the original run at k=1271.)
  // Loading lat.1270 and rng.1269 reproduces the exact (U, RNG) state from
  // which trajectory k=1270 was originally run, so k=1270 and k=1271 should
  // repeat identically.
  const int k_debug_start = 1269; // RNG checkpoint index
  const int k_lat_start   = 1270; // gauge field checkpoint index
  const int k_run_upto    = k_lat_start + 10; // run 10 trajectories

  const std::string dir3 =
    "Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)
    +"at"+std::to_string(at)+"nu0"+std::to_string(nu0)
    +"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";

  const std::string str_lat = dir3+"ckpoint_lat."+std::to_string(k_lat_start);
  const std::string str_rng = dir3+"ckpoint_rng."+std::to_string(k_debug_start);

  if( !std::filesystem::exists(str_lat) ){
    std::cerr << "# ERROR: " << str_lat << " not found" << std::endl;
    return 1;
  }
  if( !std::filesystem::exists(str_rng) ){
    std::cerr << "# ERROR: " << str_rng << " not found" << std::endl;
    return 1;
  }

  U.read( str_lat );
  rng.read( str_rng );
  std::cout << "# loaded lat from k=" << k_lat_start
            << "  rng from k=" << k_debug_start << std::endl;

  D.update( U );
  std::cout << "# lambda_min = " << D.lambda_min
            << "  lambda_max = " << D.lambda_max
            << "  ratio = " << D.lambda_min/D.lambda_max
            << "  Delta = " << D.Delta() << std::endl;

  // Force is_update false to match the original run state.
  // In the crashed run D.is_update was permanently latched false from k=1220.
  D.is_update = false;
  std::cout << "# D.is_update forced false (mimicking original run state post k=1220)"
            << std::endl;

  Force pi( base );
  const double tmax  = 1.9;
  const int    nsteps = 5; // Nf=4 setting from hmc_claude.cu
  MinimumNorm2 integrator( tmax, nsteps, 100 );
  HMC2 hmc(rng, &SW, &D, U, pi, pfs, &integrator);
  D.update( U );

  double rate, dH;
  bool   is_accept;

  // Start the loop from k_lat_start (one step past the loaded gauge field).
  // The RNG loaded is from k_debug_start = k_lat_start - 1, which is
  // exactly the state the original run used to begin trajectory k_lat_start.
  const int k_crash = 1272;
  const int nsteps_normal = nsteps;
  const int nsteps_crash  = 7;

  for(int k = k_lat_start+1; k <= k_run_upto; k++){
    integrator.nsteps = (k == k_crash) ? nsteps_crash : nsteps_normal;
    std::clog << "# nsteps = " << integrator.nsteps
              << "  tau = " << integrator.stot / integrator.nsteps << std::endl;
    Timer timer;
    hmc.run( rate, dH, is_accept );
    std::cout << "# k = " << k
              << "  dH = " << dH
              << "  is_accept = " << is_accept
              << "  rate = " << rate << std::endl;
    std::cout << "# lambda_min = " << D.lambda_min
              << "  lambda_max = " << D.lambda_max
              << "  ratio = " << D.lambda_min/D.lambda_max
              << "  Delta = " << D.Delta() << std::endl;
    std::cout << "# HMC : " << timer.currentSeconds() << " sec" << std::endl;
    // no checkpoint writing
  }

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();
  return 0;
}
