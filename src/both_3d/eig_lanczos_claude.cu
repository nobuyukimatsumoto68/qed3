// eig_lanczos_claude.cu
// Entry point for the IRL + Chebyshev eigensolver.
// Computes low-lying eigenvalues of A = (D_ov + m)^\dagger (D_ov + m)
// using the Implicit Restarted Lanczos algorithm with Chebyshev
// polynomial acceleration.
//
// Usage:
//   ./a.out mass_re mass_im alpha beta [cheb_degree]
//
// argv[1] = mass_re   (real part of mass, default 0)
// argv[2] = mass_im   (imaginary part of mass, default 0)
// argv[3] = alpha     (upper edge of target window, i.e. smaller value)
// argv[4] = beta      (spectral upper bound, i.e. larger value)
// argv[5] = cheb_degree (Chebyshev degree, default 12)
//
// Output: eigenvalues written to std::clog in the format
//   i  eval  0.0  eval
// one per line, for i=0..Nk-1.

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

using Link = std::array<Idx,2>;
using Face = std::vector<Idx>;

using MS=Eigen::Matrix2cd;
using VD=Eigen::Vector2d;
using VE=Eigen::Vector3d;
using VC=Eigen::VectorXcd;

static constexpr int NS = 2;
static constexpr int DIM = 2;
static constexpr Complex I = Complex(0.0, 1.0);


namespace Comp{
  constexpr bool is_compact=false;

  // IS_OVERLAP always active; no IS_DUAL
  constexpr int NPARALLEL_DUPDATE=1;
  constexpr int NPARALLEL=12;
  constexpr int NSTREAMS=4;

  constexpr int NPARALLEL_GAUGE=12;
  constexpr int NPARALLEL_SORT=12;

  constexpr int N_REFINE=1;
  constexpr int NS=2;

  constexpr int Nt=128;

  // non-dual (S2Simp) site count
  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;

  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-8;
}

// const std::string dir = "/mnt/hdd_barracuda/qed3/dats/";
const std::string dir = "../../geometry/data/";


#include "timer.h"

#include "s2n_simp.h"
#include "s2n_dual.h"
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

#include "../../geometry/geodesic.h"

#include "sparse_matrix.h"
#include "dirac_base.h"
#include "dirac_simp.h"
#include "dirac_dual.h"
#include "dirac_ext.h"

#include "sparse_dirac.h"
#include "matpoly.h"

#include "dirac_pf.h"

// Always use OverlapWMass (supports complex mass).
#include "includes/overlap_wmass_claude.h"

#include "includes/lanczos_claude.h"


using BaseLink = std::array<Idx,2>;
using BaseFace = std::vector<Idx>;



int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  int device;
  CUDA_CHECK(cudaGetDeviceCount(&device));
  cudaDeviceProp device_prop[device];
  cudaGetDeviceProperties(&device_prop[0], 0);
  std::cout << "# dev = " << device_prop[0].name << std::endl;
  CUDA_CHECK(cudaSetDevice(0));
  std::cout << "# (GPU device is set.)" << std::endl;

  // ---------------------------------------
  constexpr Idx N = Comp::N;
  constexpr int Nt = Comp::Nt;

  double mass_re = 0.0;
  if(argc>1) mass_re = atof(argv[1]);
  double mass_im = 0.0;
  if(argc>2) mass_im = atof(argv[2]);
  Complex mass = Complex(mass_re, mass_im);
  std::cout << "# mass = " << mass << std::endl;

  // Chebyshev filter parameters.
  double alpha_cheb = 0.0;  // upper edge of target window
  double beta_cheb  = 1.0;  // spectral upper bound
  int cheb_degree   = 12;   // polynomial degree

  if(argc>3) alpha_cheb = atof(argv[3]);
  if(argc>4) beta_cheb  = atof(argv[4]);
  if(argc>5) cheb_degree = atoi(argv[5]);

  std::cout << "# alpha (target window upper edge) = " << alpha_cheb << std::endl;
  std::cout << "# beta  (spectral upper bound)     = " << beta_cheb  << std::endl;
  std::cout << "# Chebyshev degree                 = " << cheb_degree << std::endl;

  // always S2Simp (non-dual)
  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;

  using Force=GaugeExt<Base,Nt,Comp::is_compact>;
  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  using Action=U1WilsonExt<Base>;

  using Rng=ParallelRngExt<Base,Nt>;
  using Fermion=OverlapWMass<WilsonDirac>;

  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();

  // ----------------------

  const double r = 1.0;
  const double M5 = -1.0;

  const double T = 4.0;
  const double at = T/Comp::Nt;

  WilsonDirac DW(base, 0.0, r, M5, at, 1.0);

  const double gsq = 0.05;
  Action SW( gsq, at, base );
  std::cout << "# alat = " << base.mean_ell << std::endl;

  Gauge U(base);
  srand( time(NULL) );
  Rng rng(base, rand());

  std::cout << "SW = " << SW(U) << std::endl;
  std::cout << "ch.= " << SW(U) << std::endl;

  Fermion Dov(DW, mass, 21);
  Dov.update(U);
  std::cout << "# Dov set; M5 = " << M5 << std::endl;
  std::cout << "# min max ratio: "
            << Dov.lambda_min << " "
            << Dov.lambda_max << " "
            << Dov.lambda_min/Dov.lambda_max << std::endl;
  std::cout << "# delta = " << Dov.Delta() << std::endl;

  // Build HermOp = LinOpDHDWrapper (applies A = (D+m)^\dagger(D+m)).
  LinOpDHDWrapper<Fermion> HermOp(Dov);

  // Build Chebyshev-filtered operator.
  ChebyshevFilter ChebOp(HermOp, cheb_degree, alpha_cheb, beta_cheb);

  // Build IRL state: Nk=50 converged eigenpairs, Krylov dim Nm=150.
  constexpr int Nk_irl = 50;
  constexpr int Nm_irl = 150;
  IRLState s(Nk_irl, Nm_irl);

  // Build random initial vector on CPU, copy to device.
  CuC* d_src;
  CUDA_CHECK(cudaMalloc(&d_src, N * CD));
  {
    std::vector<Complex> h_src(N);
    srand( time(NULL) + 1 );
    for(Idx i = 0; i < N; ++i){
      h_src[i] = Complex( (double)rand()/RAND_MAX - 0.5,
                          (double)rand()/RAND_MAX - 0.5 );
    }
    CUDA_CHECK(cudaMemcpy(d_src,
                          reinterpret_cast<const CuC*>(h_src.data()),
                          N * CD, H2D));
  }

  // Create cuBLAS handle.
  cublasHandle_t cublas_handle;
  CUBLAS_CHECK(cublasCreate(&cublas_handle));

  // Run IRL eigensolver.
  calc(ChebOp, HermOp, s, d_src, 1.0e-8, 200, cublas_handle);

  // Clean up.
  CUBLAS_CHECK(cublasDestroy(cublas_handle));
  CUDA_CHECK(cudaFree(d_src));

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();

  return 0;
}
