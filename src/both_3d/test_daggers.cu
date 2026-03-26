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



#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Easy.hpp>

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
  constexpr int Nt=96; // add 4

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

// #include "../../integrator/geodesic.h"
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

//------------------------------------------
// https://gist.github.com/ashwin/d88184923c7161d368a9
#include <getopt.h>

void PrintHelp()
{
  std::cout <<
    "--ell1 <ell1>"
    "--em1 <em1>"
    "--ell2 <ell2>"
    "--em2 <em2>"
    "--a1 <a>"
    "--a2 <b>";
  exit(1);
}

void parse_args(int argc, char** argv,
                int& ell1,
                int& em1,
                int& ell2,
                int& em2,
                int& a,
                int& b
                ){
  const char* const short_opts = "lmpqabh";
  const option long_opts[] = {
    {"ell1", required_argument, nullptr, 'l'},
    {"em1", required_argument, nullptr, 'm'},
    {"ell2", required_argument, nullptr, 'p'},
    {"em2", required_argument, nullptr, 'q'},
    {"a1", required_argument, nullptr, 'a'},
    {"a2", required_argument, nullptr, 'b'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, no_argument, nullptr, 0}
  };

  while (true){
    const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

    if (-1 == opt) break;

    switch (opt) {
    case 'l':
      ell1 = std::stoi(optarg);
      // std::cout << "Num set to: " << num << std::endl;
      break;

    case 'm':
      em1 = std::stoi(optarg);
      // std::cout << "Beep is set to true\n";
      break;

    case 'p':
      ell2 = std::stoi(optarg);
      // std::cout << "Num set to: " << num << std::endl;
      break;

    case 'q':
      em2 = std::stoi(optarg);
      // std::cout << "Beep is set to true\n";
      break;

    case 'a':
      a = std::stoi(optarg);
      // std::cout << "Num set to: " << num << std::endl;
      break;

    case 'b':
      b = std::stoi(optarg);
      // std::cout << "Beep is set to true\n";
      break;

    case 'h': // -h or --help
    case '?': // Unrecognized option
    default: PrintHelp();
      break;
    }
  }
}
//------------------------------------------


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




  constexpr Idx N = Comp::N;
  constexpr int Nt = Comp::Nt;

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;

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

  double nu0 = 1.0;
  if(argc>1) nu0 = atof(argv[1]);
  std::cout << "# debug. nu0 = " << nu0 << std::endl;

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();

  Gauge U(base);
  Rng rng(base, rand());
  U.gaussian( rng, 0.1 );

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
  //
  auto f_pre_adj = std::bind(&Fermion::mult_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);
  auto f_sq_adj = std::bind(&Fermion::DDH_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);

  LinOpWrapper M_pre( f_pre );
  MatPoly pre; pre.push_back ( cplx(1.0), {&M_pre} );
  LinOpWrapper M_sq( f_sq );
  MatPoly sq; sq.push_back ( cplx(1.0), {&M_sq} );

  LinOpWrapper M_pre_adj( f_pre_adj );
  MatPoly pre_adj; pre_adj.push_back ( cplx(1.0), {&M_pre_adj} );
  LinOpWrapper M_sq_adj( f_sq_adj );
  MatPoly sq_adj; sq_adj.push_back ( cplx(1.0), {&M_sq_adj} );


  // ----------------------

  std::cout << "# calculating src " << std::endl;

  FermionMatrix eta;
  FermionMatrix eta_rev;
  FermionMatrix eta_adj;

  const int s1 = 0;
  const int s2 = 10;

  const int ix1 = 4;
  const int ix2 = 31;

  std::cout << "# ix1 = " << ix1 << std::endl;
  std::cout << "# ix2 = " << ix2 << std::endl;

  for(int spin=0; spin<2; spin++){
    {
      FermionVector src1; // (base, Nt, rng);
      FermionVector src; // (base, Nt, rng);
      src1.set_pt_source(s2, ix2, spin);

      pre.from_cpu<N>( src.field, src1.field );
      FermionVector& sink = eta[spin]; // (base, Nt, rng);
      std::cout << "# calculating sink" << std::endl;
      sq.solve<N>( sink.field, src.field );
      std::cout << "# done" << std::endl;
    }
    {
      FermionVector src1; // (base, Nt, rng);
      FermionVector src; // (base, Nt, rng);
      src1.set_pt_source(s2, ix2, spin);

      pre_adj.from_cpu<N>( src.field, src1.field );
      FermionVector& sink = eta_adj[spin]; // (base, Nt, rng);
      std::cout << "# calculating sink" << std::endl;
      sq_adj.solve<N>( sink.field, src.field );
      std::cout << "# done" << std::endl;
    }
    {
      FermionVector src1; // (base, Nt, rng);
      FermionVector src; // (base, Nt, rng);
      src1.set_pt_source(s1, ix1, spin);

      pre.from_cpu<N>( src.field, src1.field );
      FermionVector& sink = eta_rev[spin]; // (base, Nt, rng);
      std::cout << "# calculating sink" << std::endl;
      sq.solve<N>( sink.field, src.field );
      std::cout << "# done" << std::endl;
    }
  }


  const MS orig = eta.get_spinmatrix(s1, ix1);
  const MS madj = -eta_adj.get_spinmatrix(s1, ix1);
  const MS mrevadj = -(eta_rev.get_spinmatrix(s2, ix2)).adjoint();
  std::cout << "orig = " << std::endl
            << orig << std::endl;
  std::cout << "madj = " << std::endl
            << madj << std::endl;
  std::cout << "mrevadj = " << std::endl
            << mrevadj << std::endl;
  // std::cout << "diff = " << std::endl
  //           << ( eta.get_spinmatrix(s2, ix2)+(eta_adj.get_spinmatrix(s1, ix1)).adjoint() ).norm() << std::endl;


  // ------------------
  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();



  return 0;
}
