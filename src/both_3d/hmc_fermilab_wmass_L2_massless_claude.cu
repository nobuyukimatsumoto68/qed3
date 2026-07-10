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
  constexpr int NPARALLEL_DUPDATE=1;   // was 4 (set 2026-06-16); -> NPARALLEL=NSTREAMS=1 via deps below: single CUDA stream for MPS packing (2 clients/GPU)
  constexpr int NPARALLEL=NPARALLEL_DUPDATE;
  constexpr int NSTREAMS=NPARALLEL_DUPDATE;

  constexpr int NPARALLEL_GAUGE=16;
  constexpr int NPARALLEL_SORT=16;

  constexpr int N_REFINE=2;
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

// #include "sparse_matrix.h"            // C4b -> multishift copy below
#include "sparse_matrix_claude.h"

#include "dirac_simp.h"
#include "dirac_ext.h"

// #include "sparse_dirac.h"
#include "sparse_dirac_claude.h"   // O(len) bucketing CSR build (was O(N*len)); -DCSR_VERIFY to check
// #include "matpoly.h"
#include "matpoly_claude.h"
#define GRAD_L4   // HMC force opt L1+L2+L4 (hoist + block poles + skip do_it); force==ref ~1e-16 (~3.4x grad)
#include "includes/overlap_wmass_claude.h"
// #include "pseudofermion.h"            // C4b -> multishift copy below
#include "pseudofermion_claude.h"

#include "integrator.h"
#include "hmc.h"


int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) == "-h") {
      printf("Usage: ./a.out [gsq] [Nf] [nu0] [mass_re] [mass_im] [max_sec]\n");
      printf("  gsq      Wilson coupling squared (default: 8.0)\n");
      printf("  Nf       number of fermion flavors (default: 2)\n");
      printf("  nu0      mass parameter (default: 1.0)\n");
      printf("  mass_re  real part of PHYSICAL mass m, R=1 units (diagonal m_L = m*A_y/abar_s built internally) (default: 0.0)\n");
      printf("  mass_im  imaginary part of physical mass m (default: 0.0)\n");
      printf("  max_sec  wall-time budget in seconds, 0 = unlimited (default: 0.0)\n");
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
  double max_sec = 0.0;   // wall-time budget [s] (0 = unlimited); stop before a traj that would overrun (set 2026-06-16)
  if(argc>6) max_sec = atof(argv[6]);
  Complex mass = Complex(mass_re, mass_im);
  std::cout << "# gsq = " << gsq << " Nf = " << Nf << " nu0 = " << nu0 << " physical_m = " << mass << " (R=1 units; diagonal m_L = m*A_y/abar_s)" << std::endl;
  std::cout << "# max_sec = " << max_sec << " (wall-time budget; 0 = unlimited)" << std::endl;
  Timer wall_timer;   // elapsed since program start (includes structure build); drives the graceful wall-time stop


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
  std::cout << "# mass_coeff = physical_m * mean_dual_area/mean_ell = " << mass*base.mean_dual_area/base.mean_ell
            << "  (uniform-measure equivalent; at L=1 equals the old bare mass)" << std::endl;

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

  // gsq=12 study (created 2026-07-02): override to the safe WIDE Zolotarev window n=21/k=0.001.
  // Stronger coupling -> rougher gauge fields -> lower Wilson kernel low-modes; the wide fixed
  // window guards sign(H_W) against the zero-crossing spikes seen at L=4 (A.D.Kennedy hep-lat/0402038).
  // Original campaign L2 value (n=17, k_=0.01 default) preserved below, commented, per convention.
  // Fermion D(DW, mass, 17);   // npole=17 (set 2026-06-16, was 21): kernel ratio lambda_min/lambda_max=0.149 sits well inside the k_=0.01 Zolotarev window, so 17 poles suffice. projected delta ~3e-6 (log-linear fit -0.357 dec/pole off n=21 @ 1.2e-7); confirm printed "# delta" at startup
  Fermion D(DW, mass, 21, 0.001);   // n=21, k=0.001 (gsq=12 wide-window override)
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

  // const int k_ckpoint_rng=100;
  const int k_ckpoint_rng=5;   // L=2: keep rng every 5 conf (set 2026-06-17 after dup-chain incident; was 100 = rolling-latest, blocked clean rollback)
  // const int kmax=1200;
  // const int kmax=600;   // L=2 max conf (set 2026-06-15)
  // const int kmax=120;   // L=2 max conf (set 2026-06-21)
  // const int kmax=420;   // L=2 max conf bumped 120 -> 420 (set 2026-06-22; L2 reached 119 cap)
  // const int kmax=320;   // L=2 max conf (recapped 420 -> 320 on 2026-06-23; realistic target vs alloc)
  // const int kmax=620;   // 2026-07-03 MASSLESS strong-coupling study (gsq=12,16): L=2 cap 620 (NM)
  const int kmax=1200;  // 2026-07-10 MASSLESS study L=2 cap bumped 620 -> 1200 (NM)

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
  // const double tmax = 1.9;
  const double tmax = 1.0;   // 2026-06-20: shortened trajectory tmax 1.9 -> 1.0
  int nsteps;
  // 2026-06-02 15:03: bumped +3 (Nf=2: 4->7, Nf=4,6: 5->8) to reduce discretization error after Nf=4,6 runs stuck at 100% rejection
  // 2026-06-04 10:42: bumped to 2x the original (Nf=2: 4->8, Nf=4,6: 5->10) to reduce discretization error after Nf=4,6 runs stuck at 100% rejection
  // if(Nf==2) nsteps = 10;     // L2 originals (more steps for finer lattice)
  // else if(Nf==4) nsteps = 12;
  // else if(Nf==6) nsteps = 12;
  // else nsteps = 12;
  // if(Nf==2) nsteps = 9;       // benchmark: |dH|~0.2 at nsteps=8 (100% accept); 9 for margin (set 2026-06-16)
  // else if(Nf==4) nsteps = 10;
  // else if(Nf==6) nsteps = 10;
  // else nsteps = 10;
  // 2026-06-20: UNIFY nsteps across Nf to the Nf=2 value (L2: 9). Same integrator resolution
  // for all flavors; Nf=4/6 were at 10. Acceptance already ~100% at L2 so 9 is safe.
  // if(Nf==2) nsteps = 9;
  // else if(Nf==4) nsteps = 9;
  // else if(Nf==6) nsteps = 9;
  // else nsteps = 9;
  // 2026-06-20: tmax 1.9 -> 1.0, nsteps L2 -> 5 (all Nf)
  // if(Nf==2) nsteps = 5;
  // else if(Nf==4) nsteps = 5;
  // else if(Nf==6) nsteps = 5;
  // else nsteps = 5;
  // 2026-06-23: nsteps L2 5 -> 6 (all Nf), tmax=1.0 -> dt=1/6; finer integration
  // if(Nf==2) nsteps = 6;
  // else if(Nf==4) nsteps = 6;
  // else if(Nf==6) nsteps = 6;
  // else nsteps = 6;
  // 2026-07-03: gsq12 nsteps 6 -> 8 (all Nf). At g^2=12 (stronger coupling -> rougher fields ->
  // larger fermion force) nsteps=6 gave dH~36 all-rejected on a hard Nf6 heavy mRe0.4229 config
  // (NOT a Zolotarev-window issue: n=21/k=0.001 already, delta=1.5e-5, no "eval below" warning).
  if(Nf==2) nsteps = 8;
  else if(Nf==4) nsteps = 8;
  else if(Nf==6) nsteps = 8;
  else nsteps = 8;
  std::cout << "# tmax = " << tmax << std::endl
            << "# nsteps = " << nsteps << std::endl;

  MinimumNorm2 integrator( tmax, nsteps, 100 );
  HMC2 hmc(rng, &SW, &D, U, pi, pfs, &integrator);
  D.update( U );

  double rate, dH;
  bool is_accept;

  double r_mean;
  double last_traj_sec = 0.0;   // wall time of the previous trajectory (drives the budget estimate)
  for(int k=k_tmp+1; k<kmax; k++){
    // graceful wall-time stop: never START a trajectory we cannot finish (+checkpoint) within max_sec.
    // The 1.3x margin covers per-traj variance; the first traj always runs (last_traj_sec==0).
    if(max_sec > 0.0 && last_traj_sec > 0.0){
      const double elapsed = wall_timer.currentSeconds();
      const double est = 1.3*last_traj_sec;
      if(elapsed + est > max_sec){
        std::cout << "# wall budget reached: stopping before traj " << k
                  << " (elapsed " << elapsed << "s + est " << est << "s > budget " << max_sec << "s)" << std::endl;
        break;
      }
    }
    Timer timer;
    hmc.run( rate, dH, is_accept);
    std::cout << "# dH : " << dH
              << " is_accept : " << is_accept
              << " rate : " << rate << std::endl;
    r_mean += rate;
    last_traj_sec = timer.currentSeconds();
    std::cout << "# HMC : " << last_traj_sec << " sec" << std::endl;

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
