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
// #define IS_OVERLAP
// #define IS_DAGGER
// #undef _OPENMP


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
  constexpr int NPARALLEL_SORT=16; // 12

  constexpr int N_REFINE=4;
  constexpr int NS=2;

  constexpr int Nt=24;
  // constexpr int Nt=1;
  // constexpr int Nt=16;

#ifdef IS_DUAL
  constexpr Idx N_SITES=20*N_REFINE*N_REFINE;
#else
  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
#endif

  constexpr Idx Nx=NS*N_SITES; // matrix size of DW
  constexpr Idx N=Nx*Nt; // matrix size of DW

  const double TOL_INNER=1.0e-15;
  const double TOL_OUTER=1.0e-14;
}

const std::string dir = "/mnt/hdd_barracuda/qed3/dats/";

// // #define IsVerbose
// #define IsVerbose2
// // #define InfoForce
// #define InfoDelta

#include "timer.h"

#include "s2n_simp.h"
#include "s2n_dual.h"
#include "rng.h"
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
#include "dirac_dual.h"
#include "dirac_ext.h"
// // #include "pseudofermion.h"
// #include "dirac.h"

#include "sparse_dirac.h"
#include "matpoly.h"

#include "dirac_pf.h"
#include "overlap.h"

#include "valence.h"

#include "../../integrator/geodesic.h"


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

  using Rng=ParallelRngExt<Base,Nt>;


  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;

  // ----------------------


  Gauge U(base);
  srand( time(NULL) );
  Rng rng(base, rand());

  // const double at = 0.5;
  const double T = 0.2;
  const double at = T/Comp::Nt;
  assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);




  using Action=U1WilsonExt;
  const double gR = 10.0;
  double beta = 4.0; // 1.0/(gR*gR);
  Action SW(beta, at);

  U.gaussian( rng, 0.001 );

  std::cout << SW(U) << std::endl;

  U.random_gauge_trsf(rng, 2.0);
  U.random_gauge_trsf(rng, 2.0);
  U.random_gauge_trsf(rng, 2.0);

  std::cout << SW(U) << std::endl;







  return 1;





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
  const double M5 = -1.5;
#else
  const double M5 = -1.0;
#endif

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

  FermionVector src1; // (base, Nt, rng);
  FermionVector src; // (base, Nt, rng);
  src1.set_pt_source(0, 0, 0);
  // src1.set_pt_source(Comp::Nt/4, 0, 1);
  pre.from_cpu<N>( src.field, src1.field );

  FermionVector sink; // (base, Nt, rng);

  std::cout << "# calculating sink" << std::endl;

  sq.solve<N>( sink.field, src.field );

  std::cout << "# done" << std::endl;



  std::vector<double> thetas;
  std::vector<double> phis;
  std::vector<double> lengths;
#ifdef IS_DUAL
  {
    std::string dir = "/mnt/hdd_barracuda/qed3/dats/";
    std::vector<Geodesic::V3> sites;
    {
      std::ifstream file(dir+"pts_dual_n"+std::to_string(Comp::N_REFINE)+"_singlepatch.dat");

      std::string str;
      while (std::getline(file, str)){
        std::istringstream iss(str);
        double v1, v2, v3;
        iss >> v1;
        iss >> v2;
        iss >> v3;
        sites.push_back( Geodesic::V3(v1, v2, v3) );
      }
    }
    const auto x0 = sites[0];
    for(const auto& elem : sites){
      double len = Geodesic::geodesicLength(Geodesic::Pt(x0), Geodesic::Pt(elem));
      // std::cout << "len = " << len << std::endl;
      lengths.push_back(len);
      thetas.push_back( Geodesic::projectionS2(elem)[0] );
      phis.push_back( Geodesic::projectionS2(elem)[1] );
    }
  }
  // double alat;
  // {
  //   std::string dir = "/mnt/hdd_barracuda/qed3/dats/";
  //   std::ifstream file(dir+"alat_n"+std::to_string(Comp::N_REFINE)+"_singlepatch.dat");

  //   std::string str;
  //   std::getline(file, str);
  //   std::istringstream iss(str);
  //   iss >> alat;
  // }
#else
  {
    const auto x0 = base.sites[0];
    for(int ix=0; ix<base.n_sites; ix++){
      const auto x1 = base.sites[ix];
      double len = Geodesic::geodesicLength(Geodesic::Pt(x0), Geodesic::Pt(x1));
      // std::cout << "len = " << len << std::endl;
      lengths.push_back(len);
      thetas.push_back( Geodesic::projectionS2(x1)[0] );
      phis.push_back( Geodesic::projectionS2(x1)[1] );
    }
  }
#endif

  const double width = 0.05;

  double factor = at*base.mean_ell;
  if(Comp::Nt==1) factor = base.mean_ell;

  {
    std::string path = "prop_spacial_L"+std::to_string(Comp::N_REFINE)+"_Nt"+std::to_string(Nt)+".dat1";
#ifdef IS_DUAL
    path = "dual_"+path;
#endif
#ifdef IS_OVERLAP
    path = "ov_"+path;
#endif
#ifdef IS_DAGGER
    path = "dagger_"+path;
#endif
    std::ofstream ofs(path);

    // Idx counter=0;
    for(Idx ix=0; ix<base.n_sites; ix++) {
      if( phis[ix]>width || phis[ix]<0. ) continue;
      {
        const auto elem = sink(0,ix,0);
        ofs << std::setw(25) << thetas[ix] << " "
          // ofs << std::setw(25) << lengths[ix] << " "
            // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.real() << " "
            // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.imag() << std::endl;
        << std::setw(25) << 1.0 * elem.real() / factor << " "
        << std::setw(25) << 1.0 * elem.imag() / factor << std::endl;
      }
      {
        const auto elem = sink(0,ix,1);
        ofs << std::setw(25) << thetas[ix] << " "
          // ofs << std::setw(25) << lengths[ix] << " "
            // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.real() << " "
            // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.imag() << std::endl;
            << std::setw(25) << 1.0 * elem.real() / factor << " "
            << std::setw(25) << 1.0 * elem.imag() / factor << std::endl;
      }
      // counter++;
    }
  }
  {
    std::string path = "prop_temporal_L"+std::to_string(Comp::N_REFINE)+"_Nt"+std::to_string(Nt)+".dat1";
#ifdef IS_DUAL
    path = "dual_"+path;
#endif
#ifdef IS_OVERLAP
    path = "ov_"+path;
#endif
#ifdef IS_DAGGER
    path = "dagger_"+path;
#endif
    std::ofstream ofs(path);

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

