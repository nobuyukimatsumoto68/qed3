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

  constexpr int N_REFINE=1;
  constexpr int NS=2;

  // constexpr int Nt=24;
  // constexpr int Nt=48; // add 2
  // constexpr int Nt=64;
  constexpr int Nt=96; // add 4
  // constexpr int Nt=120;
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

#include "../../integrator/geodesic.h"

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
  const double at = 0.2;
  if(Nt!=1) assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();


  double gsq = 4.0;
  if(argc>1) gsq = atof(argv[1]);
  int Nf = 6;
  if(argc>2) Nf = atoi(argv[2]);
  double nu0=1.0;
  if(argc>3) nu0 = atof(argv[3]);
  std::cout << "# gsq = " << gsq << " Nf = " << Nf << " nu0 = " << nu0 << std::endl;



  
  // const double gsq = 4.0;
  // double at = 0.05; // base.mean_ell * 0.125 * ratio;
  // if(Comp::Nt==1) at=0.;
  Action SW( gsq, at, base );
  std::cout << "# alat = " << base.mean_ell << std::endl;


  Gauge U(base);
  srand( time(NULL) );
  Rng rng(base, rand());

#ifdef GAUGE_TRSF
  // std::cout << "SW = " << SW(U) << std::endl;
  // U.random_gauge_trsf(rng);
  // std::cout << "ch.= " << SW(U) << std::endl;
#endif



#ifdef IS_OVERLAP
  // Overlap Dov(DW, 31);
  // Dov.update(U);
  // std::cout << "# Dov set; M5 = " << M5 << std::endl;
  // std::cout << "# min max ratio: "
  //           << Dov.lambda_min << " "
  //           << Dov.lambda_max << " "
  //           << Dov.lambda_min/Dov.lambda_max << std::endl;
  // std::cout << "# delta = " << Dov.Delta() << std::endl;

  // auto f_Op = std::bind(&Overlap::mult_deviceAsyncLaunch, &Dov, std::placeholders::_1, std::placeholders::_2);
  // LinOpWrapper M_Op( f_Op );
  // Op.push_back ( cplx(1.0), {&M_Op} );

#ifdef IS_DUAL
  const double M5 = -1.0;
  // const double M5 = -1.0;
#else
  const double M5 = -1.0;
#endif

  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);
  std::cout << "# DW set. " << std::endl;

  using Fermion=Overlap<WilsonDirac>;
  Fermion D(DW, 51);
  std::cout << "# D set. " << std::endl;
#else
  const double M5 = 0.0;
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);
  std::cout << "# DW set. " << std::endl;

  using Fermion=DiracPf<WilsonDirac>;
  Fermion D(DW);
  std::cout << "# D set. " << std::endl;
#endif




  std::string dir3;
  dir3="Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
  std::filesystem::create_directory(dir3);
  // const int k_ckpoint=1;


  int k_tmp = 100;
  const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k_tmp);
  const std::string str_rng=dir3+"ckpoint_rng."+std::to_string(k_tmp);
  U.read( str_lat );
  rng.read( str_rng );
  std::cout << "# U read. " << std::endl;















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


//   std::string path = "prop_temporal_L"+std::to_string(Comp::N_REFINE)+"_Nt"+std::to_string(Nt)+"_at"+std::to_string(at)+"_nu0"+std::to_string(nu0)+".dat";
// #ifdef IS_DUAL
//     path = "dual_"+path;
// #endif
// #ifdef IS_OVERLAP
//     path = "ov_"+path;
// #endif
// #ifdef IS_DAGGER
//     path = "dagger_"+path;
// #endif
// #ifdef IS_FLAT
//     path = "flat_"+path;
// #endif
//     std::ofstream ofs(path);
  // ofs << std::scientific << std::setprecision(15);
  // ofs << std::scientific << std::setprecision(15);

  // Idx counter=0;
  for(Idx s=0; s<Comp::Nt; s++) {
    {
      const auto elem = sink(s,0,0);
      std::cout << std::setw(25) << at*s << " "
        // ofs << std::setw(25) << s << " "
        // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.real() << " "
        // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.imag() << std::endl;
          << std::setw(25) << 1.0 * elem.real()  << " "
          << std::setw(25) << 1.0 * elem.imag()  << std::endl;
    }
    {
      const auto elem = sink(s,0,1);
      std::cout << std::setw(25) << at*s << " "
        // ofs << std::setw(25) << s << " "
        // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.real() << " "
        // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.imag() << std::endl;
          << std::setw(25) << 1.0 * elem.real()  << " "
          << std::setw(25) << 1.0 * elem.imag()  << std::endl;
    }
    // counter++;
  }
  // }





  // ------------------

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();

  return 0;

}

