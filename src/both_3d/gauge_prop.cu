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


  constexpr int Nt=24;
  // constexpr int Nt=192;
  // constexpr int Nt=1;
  // constexpr int Nt=16;

#ifdef IS_DUAL
  constexpr Idx N_SITES=20*N_REFINE*N_REFINE;
#else
  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE; // 30, 120, 480
#endif


  constexpr Idx Nx=NS*N_SITES; // matrix size of DW
  constexpr Idx N=Nx*Nt; // matrix size of DW

  const double TOL_INNER=1.0e-15;
  const double TOL_OUTER=1.0e-14;
}

const std::string dir = "../../dats/";

#define IsVerbose
#define IsVerbose2
// // #define InfoForce
// #define InfoDelta

#include "timer.h"

#include "../../integrator/geodesic.h"

#include "s2n_simp.h"
#include "s2n_dual.h"
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
  const double T = 12;
  const double at = T/Comp::Nt;
  // assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);

  const double gsq = 0.05;
  Action SW( gsq, at, base );
  std::cout << "# alat = " << base.mean_ell << std::endl;

  Gauge U(base);
  srand( time(NULL) );
  Rng rng(base, rand());
  U.gaussian( rng, 0.3 );

  // std::cout << "debug. pt2" << std::endl;
  std::vector<Idx> is;
  std::vector<Idx> js;
  std::vector<Complex> vs;
  SW.coo_format(is, js, vs, U);

  const Idx len = vs.size();

  std::vector<COOEntry> coo;
  for(Idx i=0; i<len; i++) coo.push_back( COOEntry( vs[i], is[i], js[i] ) );

  COO op; op.en = coo; op.do_it();
  MatPoly mat; mat.push_back ( cplx(1.0), {&op} );

  constexpr Idx Ng = Comp::Nt * (Comp::N_LINKS + Comp::N_SITES);
  std::vector<Complex> res(Ng), v(Ng);
  U.vectorize( v );
  matmulcoo( reinterpret_cast<CuC*>(res.data()), reinterpret_cast<CuC*>(v.data()), coo, v.size() );

  Complex action = 0.0;
  for(Idx i=0; i<res.size(); i++) action += res[i] * v[i];
  action *= 0.5;
  std::cout << "action = " << action.real() << std::endl;
  std::cout << "SW = " << SW(U) << std::endl;

  // LinOpWrapper M_pre( f_pre );
  // MatPoly mat;
  // mat.push_back ( cplx(1.0), {&op} );
  // @@@@ DEBUG HERE !!!!!
  {
    // std::vector<Complex> res(Ng), v(Ng);
    U.gaussian( rng, 1.2 );
    U.vectorize( v );

    mat.from_cpu<Ng>( res, v );
    Complex action2 = 0.0;
    for(Idx i=0; i<res.size(); i++) action2 += res[i] * v[i];
    action2 *= 0.5;
    std::cout << "action = " << action2 << std::endl;
    std::cout << "SW = " << SW(U) << std::endl;
  }

  // const Face& face = base.faces[0];
  // std::vector<double> res(Comp::Nt, 0.0);

  // for(Idx i=0; i<face.size(); i++) {
  //   const Idx ix = face[i];
  //   const Idx iy = face[(i+1)%face.size()];
  //   const Link li{ix,iy};

  //   std::cout << "# prepare source" << std::endl;
  //   GaugeVector sink(U), source(U);
  //   source.set_zero();
  //   sink.set_zero();
  //   source.sp(0, li) = 1.0 * base.map2sign.at(li);

  //   // project
  //   mat.from_cpu<Ng>( sink.field, source.field );
  //   source.set_zero();
  //   mat.solve<Ng>( source.field, sink.field );

  //   std::cout << "# calculate sink" << std::endl;
  //   // invert
  //   sink.set_zero();
  //   mat.solve<Ng>( sink.field, source.field );

  //   // std::cout << "sink = " << std::endl;
  //   // for(const auto elem : sink) std::cout << elem.real() << std::endl;

  //   for(Idx j=0; j<face.size(); j++) {
  //     const Idx jx = face[j];
  //     const Idx jy = face[(j+1)%face.size()];
  //     const Link lj{jx,jy};
  //     for(int t=0; t<Comp::Nt; t++) res[t] += base.map2sign.at(lj) * sink.sp(t, lj).real();
  //   }
  // }

  // std::cout << "res = " << std::endl;
  // for(const auto elem : res) std::cout << elem << std::endl;

  // std::cout << "# done" << std::endl;




  // ------------------


  return 0;

}

