// hmc_Nfblocked_claude.cu -- UNIFIED Nf-block HMC driver (action-inversion mrhs blocking).
// One source for all (L, Nf) targets: the Nf/2 pseudo-fermion ACTION inversions (heat-bath, update_eta,
// S) are packed into ONE N*(Nf/2) block CG via BlockedMat (blocked_mat_claude.h) instead of the serial
// per-flavor loop; the FORCE stays serial per flavor. Replaces the seven serial hmc_(fermilab_wmass_)*
// drivers for Nf=4/6 (Nf=2 keeps the serial originals + MPS). Build-time knobs (all have defaults):
//   -DLREFINE=<1|2|4>   lattice refinement (sizes Comp::N)
//   -DNFPF=<4|6>        Nf baked in; NSTACK = NFPF/2 (block width). argv Nf must match (asserted).
//   -DGEODESIC_H="..."  geodesic.h include path (local vs fermilab); -DGEOM_DATA="..." geometry data dir
//   -DDIR_NO_MASS       output dir omits mRe/mIm (local L1 convention; assumes massless)  [default: include]
//   -DKMAX=<n> -DKRNG=<n>  trajectory cap / rng-checkpoint keep-stride (per-campaign policy)
// Zolotarev fixed at n=21, window k=0.001 for all L (current standard). tmax/nsteps are L-keyed below.
//
// Reference (algorithm): B. Jegerlehner, hep-lat/9612014 (inner multishift) + block CG (outer).

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
#include <cmath>

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


// ---- Nf-block build-time knobs (defaults; the build script overrides per (L,Nf,environment)) ----
#ifndef LREFINE
#define LREFINE 2
#endif
#ifndef NFPF
#define NFPF 6
#endif
#ifndef GEODESIC_H
#define GEODESIC_H "../../geometry/geodesic.h"
#endif
#ifndef GEOM_DATA
#define GEOM_DATA "../../geometry/data/"
#endif
#ifndef KMAX
#define KMAX 320
#endif
#ifndef KRNG
#define KRNG 5
#endif

// #define IsVerbose
// #define IsVerbose2
// #define InfoForce
#define InfoDelta


namespace Comp{
  constexpr bool is_compact=false;

  // overlap only: single CUDA stream (MPS packing done at the process level, not here)
  constexpr int NPARALLEL_DUPDATE=1;
  constexpr int NPARALLEL=NPARALLEL_DUPDATE;
  constexpr int NSTREAMS=NPARALLEL_DUPDATE;

  constexpr int NPARALLEL_GAUGE=16;
  constexpr int NPARALLEL_SORT=16;

  constexpr int N_REFINE=LREFINE;
  constexpr int NS=2;

  constexpr int Nt=128;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;

  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-8;
}


const std::string dir = GEOM_DATA;
#include GEODESIC_H

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

#include "sparse_matrix_claude.h"

#include "dirac_simp.h"
#include "dirac_ext.h"

#include "sparse_dirac_claude.h"   // O(len) bucketing CSR build
#include "matpoly_claude.h"        // pulls multishift_block_kernels_claude.h (block kernels)
#define GRAD_L4                    // HMC force opt L1+L2+L4 (serial force; unchanged here)
#include "overlap_wmass_claude.h"

// ---- action-inversion path ----
// Default: Nf-block (mrhs) path. -DSERIAL_REF selects the ORIGINAL serial per-flavor path
// (pfs + HMC2 + MinimumNorm2) for the validation harness (byte-same params, so dH/timing compare
// cleanly). Production block builds never define SERIAL_REF.
#ifdef SERIAL_REF
#include "pseudofermion_claude.h"
#include "integrator.h"
#include "hmc.h"
#else
#include "blocked_mat_claude.h"
#ifdef BLOCK_FORCE
#include "overlap_force_Nfblock_claude.h"   // Phase 2: BlockedForce (must precede PseudoFermionBlock)
#endif
#include "pseudofermion_Nfblocked_claude.h"
#include "integrator_Nfblocked_claude.h"
#include "hmc_Nfblocked_claude.h"
#endif


int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) == "-h") {
      printf("Usage: ./a.out [gsq] [Nf] [nu0] [mass_re] [mass_im] [max_sec]\n");
      printf("  gsq      Wilson coupling squared (default: 8.0)\n");
      printf("  Nf       number of fermion flavors -- MUST equal the compiled NFPF=%d\n", NFPF);
      printf("  nu0      mass parameter (default: 1.0)\n");
      printf("  mass_re  real part of PHYSICAL mass m, R=1 units (diagonal m_L = m*A_y/abar_s) (default: 0.0)\n");
      printf("  mass_im  imaginary part of physical mass m (default: 0.0)\n");
      printf("  max_sec  wall-time budget in seconds, 0 = unlimited (default: 0.0)\n");
      return 0;
    }
  }

  double gsq = 8.0;
  if(argc>1) gsq = atof(argv[1]);
  int Nf = NFPF;
  if(argc>2) Nf = atoi(argv[2]);
  double nu0 = 1.0;
  if(argc>3) nu0 = atof(argv[3]);
  double mass_re = 0.0;
  if(argc>4) mass_re = atof(argv[4]);
  double mass_im = 0.0;
  if(argc>5) mass_im = atof(argv[5]);
  double max_sec = 0.0;   // wall-time budget [s] (0 = unlimited)
  if(argc>6) max_sec = atof(argv[6]);
  Complex mass = Complex(mass_re, mass_im);

  constexpr int NSTACK_PF = NFPF/2;                 // block width = Nf/2 (compile-time)
  assert(Nf==NFPF && "this binary is compiled for a fixed Nf (=NFPF); pass a matching Nf");
#ifdef DIR_NO_MASS
  assert(mass==Complex(0.0,0.0) && "DIR_NO_MASS build assumes massless (dir omits mRe/mIm)");
#endif

  std::cout << "# gsq = " << gsq << " Nf = " << Nf << " (NFPF=" << NFPF << ", NSTACK=" << NSTACK_PF << ")"
            << " nu0 = " << nu0 << " physical_m = " << mass << " (R=1 units; diagonal m_L = m*A_y/abar_s)" << std::endl;
  std::cout << "# max_sec = " << max_sec << " (wall-time budget; 0 = unlimited)" << std::endl;
  Timer wall_timer;


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

  // Zolotarev: n=21, fixed window k=0.001 for ALL L (current standard; guards sign(H_W) at strong coupling).
  Fermion D(DW, mass, 21, 0.001);
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

#ifdef SERIAL_REF
  // ---- serial per-flavor pseudo-fermions (validation reference) ----
  std::vector<std::shared_ptr<PseudoFermion<Fermion>>> pfs;
  for(int f=0; f<NSTACK_PF; f++) pfs.push_back( std::shared_ptr<PseudoFermion<Fermion>>( new PseudoFermion<Fermion>(D) ) );
  std::cout << "# SERIAL_REF path: " << NSTACK_PF << " serial PseudoFermions" << std::endl;
#else
  // ---- Nf-block pseudo-fermion manager (+ device-memory check) ----
  size_t mfree0=0, mtot=0;
  CUDA_CHECK(cudaMemGetInfo(&mfree0, &mtot));
  PseudoFermionBlock<Fermion,NSTACK_PF> bpf(D);
  size_t mfree1=0;
  CUDA_CHECK(cudaMemGetInfo(&mfree1, &mtot));
  size_t est = PseudoFermionBlock<Fermion,NSTACK_PF>::block_bytes(D.size-1);
#ifdef BLOCK_FORCE
  est += bpf.bforce.bytes();   // Phase-2 force output blocks
#endif
  std::cout << "# Nf-block: NSTACK = " << NSTACK_PF << " npole = " << (D.size-1)
            << " ; scratch est " << (est>>20) << " MiB, used " << ((mfree0-mfree1)>>20)
            << " MiB, free " << (mfree1>>20) << " / total " << (mtot>>20) << " MiB" << std::endl;
  assert(est < mtot && "Nf-block scratch exceeds device memory");
#endif

  Timer timer;

  // -----------------

  std::string dir3;
#ifdef DIR_NO_MASS
  dir3="Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
#else
  dir3="Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"mRe"+std::to_string(mass.real())+"mIm"+std::to_string(mass.imag())+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
#endif
  std::filesystem::create_directory(dir3);
  const int k_ckpoint=1;
  const int k_ckpoint_rng=KRNG;   // keep rng checkpoint every this many trajectories
  const int kmax=KMAX;

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

  // L-keyed trajectory tuning (per user 2026-07-09): L1 tmax=1.9/nsteps=10 ; L2,L4 tmax=1.0/nsteps=8.
#if LREFINE==1
  const double tmax = 1.9;
  const int nsteps = 10;
#else
  const double tmax = 1.0;
  const int nsteps = 8;
#endif
  std::cout << "# tmax = " << tmax << std::endl
            << "# nsteps = " << nsteps << std::endl;

#ifdef SERIAL_REF
  MinimumNorm2 integrator( tmax, nsteps, 100 );
  HMC2 hmc(rng, &SW, &D, U, pi, pfs, &integrator);
#elif defined(BLOCK_FORCE)
  MinimumNorm2BlockF integrator( tmax, nsteps, 100 );   // action + FORCE blocked
  HMC2Block hmc(rng, &SW, &D, U, pi, bpf, &integrator);
#else
  MinimumNorm2Block integrator( tmax, nsteps, 100 );    // action-only block; serial force
  HMC2Block hmc(rng, &SW, &D, U, pi, bpf, &integrator);
#endif
  D.update( U );

  double rate, dH;
  bool is_accept;

  double r_mean;
  double last_traj_sec = 0.0;   // wall time of the previous trajectory (drives the budget estimate)
  for(int k=k_tmp+1; k<kmax; k++){
    // graceful wall-time stop: never START a trajectory we cannot finish (+checkpoint) within max_sec.
    if(max_sec > 0.0 && last_traj_sec > 0.0){
      const double elapsed = wall_timer.currentSeconds();
      const double est_sec = 1.3*last_traj_sec;
      if(elapsed + est_sec > max_sec){
        std::cout << "# wall budget reached: stopping before traj " << k
                  << " (elapsed " << elapsed << "s + est " << est_sec << "s > budget " << max_sec << "s)" << std::endl;
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
