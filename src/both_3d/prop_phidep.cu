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


// #define IS_DUAL
#define IS_OVERLAP
// #define IS_DAGGER
// #undef _OPENMP

// #define IS_FLAT

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

  constexpr int N_REFINE=4;
  constexpr int NS=2;

  // constexpr int Nt=24;
  // constexpr int Nt=48; // add 2
  // constexpr int Nt=64;
  constexpr int Nt=96; // add 4
  // constexpr int Nt=120;
  // constexpr int Nt=128; // @@@
  // constexpr int Nt=144; // add 8
  // constexpr int Nt=168;

  // constexpr int Nt=24;
  // constexpr int Nt=192;
  // constexpr int Nt=1;
  // constexpr int Nt=16;

#ifdef IS_DUAL
  constexpr Idx N_SITES=20*N_REFINE*N_REFINE;
#else
  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
#endif
  constexpr int N_LINKS=30*N_REFINE*N_REFINE; // 30, 120, 480

  constexpr Idx Nx=NS*N_SITES; // matrix size of DW
  constexpr Idx N=Nx*Nt; // matrix size of DW

  const double TOL_INNER=1.0e-15;
  const double TOL_OUTER=1.0e-14;
}

const std::string dir = "../../dats/";

// // #define IsVerbose
// #define IsVerbose2
// // #define InfoForce
// #define InfoDelta

#include "timer.h"

#include "../../geometry/geodesic.h"

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

#include "sparse_matrix.h"
#include "dirac_base.h"
#include "dirac_simp.h"
#include "dirac_dual.h"
#include "dirac_ext.h"
// // #include "pseudofermion.h"
// #include "dirac.h"

#include "sparse_dirac.h"
#include "matpoly.h"

#include "dirac_pf.h"
#include "overlap.h"


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
  // using FermionVector = FermionVector<Base>;


  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;

  // ----------------------
  // const double at = 0.5;
  // const double T = 0.2;
  // const double T = 16;
  // const double at = T/Comp::Nt;
  const double at = M_PI/20; // 0.1;
  if(Nt!=1) assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);

  double nu0 = 1.0;
  if(argc>1) nu0 = atof(argv[1]);
  std::cout << "# debug. nu0 = " << nu0 << std::endl;

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();


  const double gsq = 0.05;
  // double at = 0.05; // base.mean_ell * 0.125 * ratio;
  // if(Comp::Nt==1) at=0.;
  Action SW( gsq, at, base );
  std::cout << "# alat = " << base.mean_ell << std::endl;


  Gauge U(base);
  srand( time(NULL) );
  Rng rng(base, rand());


  const double M5 = -1.0;

  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);
  std::cout << "# DW set. " << std::endl;

  using Fermion=Overlap<WilsonDirac>;
  Fermion D(DW, 31);
  std::cout << "# D set. " << std::endl;

  D.update( U );
  std::cout << "# D updated. " << std::endl;

  auto f_pre = std::bind(&Fermion::adj_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);
  auto f_sq = std::bind(&Fermion::DHD_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_pre( f_pre );
  MatPoly pre; pre.push_back ( cplx(1.0), {&M_pre} );
  LinOpWrapper M_sq( f_sq );
  MatPoly sq; sq.push_back ( cplx(1.0), {&M_sq} );

  // ---------------------

  std::cout << "# calculating src " << std::endl;


  FermionVector src1; // (base, Nt, rng);
  FermionVector src; // (base, Nt, rng);
  const Idx ix0=15;
  src1.set_pt_source(0, ix0, 0);
  // src1.set_pt_source(Comp::Nt/4, 0, 1);


  pre.from_cpu<N>( src.field, src1.field );
  FermionVector sink; // (base, Nt, rng);

  std::cout << "# calculating sink" << std::endl;
  sq.solve<N>( sink.field, src.field );
  std::cout << "# done" << std::endl;


  std::vector<double> thetas;
  std::vector<double> phis;
  // std::vector<double> lengths;
  {
    // const auto x0 = base.sites[0];
    for(int ix=0; ix<base.n_sites; ix++){
      const auto x1 = base.sites[ix];
      // double len = Geodesic::geodesicLength(Geodesic::Pt(x0), Geodesic::Pt(x1));
      // std::cout << "len = " << len << std::endl;
      // lengths.push_back(len);
      thetas.push_back( Geodesic::projectionS2(x1)[0] );
      phis.push_back( Geodesic::projectionS2(x1)[1] );
    }
  }
  std::cout << "theta0 = " << thetas[ix0] << std::endl;
  std::cout << "phi0 = " << phis[ix0] << std::endl;

  const double width = 0.05;
  double factor = at*base.mean_ell;
  // if(Comp::Nt==1) factor = base.mean_ell;
  std::string path = "prop_phi_L"+std::to_string(Comp::N_REFINE)+"_Nt"+std::to_string(Nt)+"_at"+std::to_string(at)+"_nu0"+std::to_string(nu0)+"_ix0"+std::to_string(ix0)+".dat";
  path = "ov_"+path;
  std::ofstream ofs(path);
  ofs << std::scientific << std::setprecision(15);
  ofs << std::scientific << std::setprecision(15);

  // Idx counter=0;
  const Idx s=4;
  // for(Idx s=0; s<Comp::Nt; s++) {
  for(int ix=0; ix<base.n_sites; ix++){
    if(std::abs(thetas[ix]-thetas[ix0]-0.5*M_PI)<width){
      auto elem = sink(s,ix,0);
      //ofs << std::setw(25) << Geodesic::Mod(phis[ix]-phis[ix0]) << " "
      ofs << std::setw(25) << Geodesic::Mod(phis[ix]) << " "
        // ofs << std::setw(25) << s << " "
        // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.real() << " "
        // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.imag() << std::endl;
          << std::setw(25) << 1.0 * elem.real() / factor << " "
          << std::setw(25) << 1.0 * elem.imag() / factor << std::endl;

      elem = sink(s,ix,1);
      ofs << std::setw(25) << Geodesic::Mod(phis[ix]) << " "
        // ofs << std::setw(25) << s << " "
        // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.real() << " "
        // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.imag() << std::endl;
          << std::setw(25) << 1.0 * elem.real() / factor << " "
          << std::setw(25) << 1.0 * elem.imag() / factor << std::endl;
    }
  }





  // ------------------

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();

  return 0;

}

