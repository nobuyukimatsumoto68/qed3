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


#define IS_OVERLAP
// #define IS_DAGGER
// #undef _OPENMP

// #define GAUGE_TRSF


namespace Comp{
  constexpr bool is_compact=false;

  // d_DW.update() is always done independently
#ifdef IS_OVERLAP
  constexpr int NPARALLEL_DUPDATE=1;
  constexpr int NPARALLEL=12; // 12
  constexpr int NSTREAMS=4; // 4
#else
  constexpr int NPARALLEL_DUPDATE=12;
  constexpr int NPARALLEL=1; // 12
  constexpr int NSTREAMS=12; // for grad loop
#endif
  constexpr int NPARALLEL_GAUGE=12; // 12
  constexpr int NPARALLEL_SORT=12; // 12

  constexpr int N_REFINE=8;
  constexpr int NS=2;

  // constexpr int Nt=24;
  // constexpr int Nt=48; // add 2
  //
  // constexpr int Nt=64;
  // constexpr int Nt=96; // add 4
  // constexpr int Nt=120;
  // constexpr int Nt=144; // add 8
  constexpr int Nt=168;

  constexpr double T = 12.;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE; // 30, 120, 480

  constexpr Idx Nx=NS*N_SITES; // matrix size of DW
  constexpr Idx N=Nx*Nt; // matrix size of DW

  const double TOL_INNER=1.0e-15;
  const double TOL_OUTER=1.0e-14;
}

const std::string dir = "../geometry/data/";

// // #define IsVerbose
// #define IsVerbose2
// // #define InfoForce
#define InfoDelta

#include "timer.h"

#include "../geometry/geodesic.h"

#include "s2n_simp.h"
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
#include "dirac_base.h"
#include "dirac_simp.h"
#include "dirac_ext.h"

#include "sparse_dirac.h"
#include "matpoly.h"

#include "dirac_pf.h"
#include "overlap.h"



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

  constexpr Idx N = Comp::N;
  constexpr int Nt = Comp::Nt;

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;

  using Force=GaugeExt<Base,Nt,Comp::is_compact>;
  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  using Action=U1WilsonExt<Base>;


  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;

  // ----------------------
  const double at = Comp::T/Comp::Nt;
  assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);


  Gauge U(base);
  // srand( time(NULL) );
  // Rng rng(base);


#ifdef IS_OVERLAP
  const double M5 = -1.0;

  WilsonDirac DW(base, 0.0, 1.0, M5, at);
  std::cout << "# DW set. " << std::endl;

  using Fermion=Overlap<WilsonDirac>;
  Fermion D(DW, 51);
  std::cout << "# D set. " << std::endl;
#else
  const double M5 = 0.0;
  WilsonDirac DW(base, 0.0, 1.0, M5, at);
  std::cout << "# DW set. " << std::endl;

  using Fermion=DiracPf<WilsonDirac>;
  Fermion D(DW);
  std::cout << "# D set. " << std::endl;
#endif

  D.update( U );
  std::cout << "# D updated. " << std::endl;

#ifdef IS_DAGGER
  auto f_pre = std::bind(&Fermion::mult_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);
  auto f_sq = std::bind(&Fermion::DDH_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);
#else
  auto f_pre = std::bind(&Fermion::adj_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);
  auto f_sq = std::bind(&Fermion::DHD_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);
#endif
  LinOpWrapper M_pre( f_pre );
  MatPoly pre; pre.push_back ( cplx(1.0), {&M_pre} );
  LinOpWrapper M_sq( f_sq );
  MatPoly sq; sq.push_back ( cplx(1.0), {&M_sq} );

  // ---------------------

  std::cout << "# calculating src " << std::endl;

#ifdef GAUGE_TRSF
  FermionVector gauge; // (base, Nt, rng);
  gauge.set_random_gauge(rng);
#endif

  FermionVector src1; // (base, Nt, rng);
  FermionVector src; // (base, Nt, rng);
  src1.set_pt_source(0, 0, 0);
  // src1.set_pt_source(Comp::Nt/4, 0, 1);

#ifdef GAUGE_TRSF
  DW.bd.random_gauge_trsf(base,rng);

  src1.gauge_trsf( gauge );
  U.gauge_trsf( gauge );
  D.update( U );
#endif

  pre.from_cpu<N>( src.field, src1.field );

  FermionVector sink; // (base, Nt, rng);

  std::cout << "# calculating sink" << std::endl;

  sq.solve<N>( sink.field, src.field );

  std::cout << "# done" << std::endl;

#ifdef GAUGE_TRSF
  sink.gauge_trsf( gauge, -1 );
#endif

  const double width = 0.05;

  double factor = at*base.mean_ell;
  if(Comp::Nt==1) factor = base.mean_ell;

  {
    std::string path = "prop_temporal_L"+std::to_string(Comp::N_REFINE)+"_Nt"+std::to_string(Nt)+"_T"+std::to_string(Comp::T)+".dat";
#ifdef IS_DUAL
    path = "dual_"+path;
#endif
#ifdef IS_OVERLAP
    path = "ov_"+path;
#endif
#ifdef IS_DAGGER
    path = "dagger_"+path;
#endif
    std::ofstream ofs("./data/"+path);
    ofs << std::scientific << std::setprecision(15);

    // Idx counter=0;
    for(Idx s=0; s<Comp::Nt; s++) {
      {
        const auto elem = sink(s,0,0);
        ofs << std::setw(25) << at*s << " "
          // ofs << std::setw(25) << s << " "
            // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.real() << " "
            // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.imag() << std::endl;
            << std::setw(25) << 1.0 * elem.real() / factor << " "
            << std::setw(25) << 1.0 * elem.imag() / factor << std::endl;
      }
      {
        const auto elem = sink(s,0,1);
        ofs << std::setw(25) << at*s << " "
          // ofs << std::setw(25) << s << " "
          // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.real() << " "
            // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.imag() << std::endl;
            << std::setw(25) << 1.0 * elem.real() / factor << " "
            << std::setw(25) << 1.0 * elem.imag() / factor << std::endl;
      }
      // counter++;
    }
  }





  // ------------------


  return 0;

}

