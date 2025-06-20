#include <typeinfo>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <cstdint>
#include <complex>

using Double = double;
using Idx = std::int32_t;
using Complex = std::complex<double>;

// #define IS_DUAL
// #define IS_OVERLAP

namespace Comp{
  constexpr int NPARALLEL=12; // 12
  constexpr int NSTREAMS=4; // 4

  constexpr int N_REFINE=8;
  constexpr int NS=2;

#ifdef IS_DUAL
  constexpr Idx N_SITES=20*N_REFINE*N_REFINE;
#else
  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
#endif

  constexpr Idx N=NS*N_SITES; // matrix size of DW

  // const double TOL=1.0e-9;
  const double TOL_INNER=1.0e-15;
  const double TOL_OUTER=1.0e-14;
}

// #define IsVerbose
#define IsVerbose2
// #define InfoForce
#define InfoDelta

const std::string dir = "/mnt/hdd_barracuda/qed3/dats/";

#include "s2n_dual.h"
#include "s2n_simp.h"
#include "rng.h"
#include "gauge.h"
#include "action.h"


#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cublas_api.h>
#include <cusolverDn.h>
using CuC = cuDoubleComplex;
#include "gpu_header.h"

// ======================================

// #include "timer.h"

// #include "s2n_dual.h"
// #include "s2n_simp.h"
// #include "rng.h"
// #include "gauge.h"
// #include "force.h"
// #include "action.h"
// #include "sparse_matrix.h"
// #include "dirac_dual.h"
// #include "dirac_simp.h"
// #include "sparse_dirac.h"
// #include "matpoly.h"
// #include "overlap.h"
// #include "pseudofermion.h"

// #include "integrator.h"
// #include "hmc.h"

#include "sparse_matrix.h"
// #include "pseudofermion.h"
#include "dirac_dual.h"
#include "dirac_simp.h"

#include "sparse_dirac.h"
#include "matpoly.h"

#include "overlap.h"


#include "valence.h"

using Double = double;
#include "../../integrator/geodesic.h"
// #include "dirac_s2_dual.h"
// #include "header_cusolver.hpp"


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

#ifdef IS_DUAL
  using Lattice=S2Trivalent;
  using Gauge=U1onS2<Lattice,false>;
  using WilsonDirac=DiracS2Dual<Gauge>;
#else
  using Lattice=S2Simp;
  using Gauge=U1onS2<Lattice,false>;
  using WilsonDirac=DiracS2Simp<Gauge>;
#endif
  using Rng=ParallelRng<Lattice>;
  using Overlap=Overlap<Gauge,WilsonDirac,Lattice>;

  Lattice lattice(Comp::N_REFINE);

  using Link = std::array<Idx,2>; // <int,int>;
  constexpr Idx N = Comp::N;

  // ----------------------

  std::cout << "# lattice set. " << std::endl;

#ifdef IS_OVERLAP
  const double r = 1.0;
#ifdef IS_DUAL
  const double M5 = -1.6/2.0 * 0.5*3.0/2.0;
#else
  const double M5 = -1.6/2.0 * 0.5*(1.0 + std::sqrt( 5.0 + 2.0*std::sqrt(2.0) ));
#endif
#else
  const double r = 1.0;
  const double M5 = 0.0;
#endif
  WilsonDirac DW(lattice, 0.0, r, M5);

  std::cout << "# DW set" << std::endl;

  Gauge U(lattice);
  Rng rng(lattice);
  // U.gaussian( rng, 0.2 );

  // ---------------------

  Overlap Dov(DW, 31);
  std::cout << "# Dov set; M5 = " << M5 << std::endl;
  Dov.update(U);
  std::cout << "# min max ratio: "
            << Dov.lambda_min << " "
            << Dov.lambda_max << " "
            << Dov.lambda_min/Dov.lambda_max << std::endl;
  std::cout << "# delta = " << Dov.Delta() << std::endl;

#ifdef IS_OVERLAP
  auto f_DHD = std::bind(&Overlap::sq_deviceAsyncLaunch, &Dov,
                         std::placeholders::_1, std::placeholders::_2);
  auto f_DH = std::bind(&Overlap::adj_deviceAsyncLaunch, &Dov,
                        std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_DHD( f_DHD );
  LinOpWrapper M_DH( f_DH );

  MatPoly DHD;
  DHD.push_back ( cplx(1.0), {&M_DHD} );
  //
  MatPoly DH;
  DH.push_back ( cplx(1.0), {&M_DH} );
#else
  MatPoly DHD;
  DHD.push_back ( cplx(1.0), {&Dov.M_DW, &Dov.M_DWH} );
  //
  MatPoly DH;
  DH.push_back ( cplx(1.0), {&Dov.M_DWH} );
#endif

  // ---------------------


  FermionVector src1(lattice, rng);
  FermionVector src(lattice, rng);
  src1.set_pt_source(0, 0);
  DH.from_cpu<N>( src.field, src1.field );

  


  FermionVector sink(lattice, rng);
  // Op.from_cpu<N>( sink.field, src.field );
  // FermionVector rc(lattice, rng);
  // rc.set_random();

  // DHD.bicgstab<N>( sink.field, src.field, rc.field, 1.0e-3, 1e8, 1.0e-8 );
  DHD.solve<N>( sink.field, src.field );



#ifdef IS_DUAL
  std::vector<double> lengths;
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
  std::vector<double> lengths;
  {
    const auto x0 = lattice.sites[0];
    for(int ix=0; ix<lattice.n_sites; ix++){
      const auto x1 = lattice.sites[ix];
      double len = Geodesic::geodesicLength(Geodesic::Pt(x0), Geodesic::Pt(x1));
      // std::cout << "len = " << len << std::endl;
      lengths.push_back(len);
    }
  }
#endif

  Idx counter=0;
  for(auto& elem : sink) {
    std::cout << std::setw(25) << lengths[int(counter/2)] << " "
              << std::setw(25) << 1.0/lattice.alat * elem.real() << " "
              << std::setw(25) << 1.0/lattice.alat * elem.imag() << std::endl;
    counter++;
  }



  // ------------------


  return 0;

}

