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



int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);
  constexpr Idx N = Comp::N;
  constexpr int Nt = Comp::Nt;
  using Base=S2Simp;
  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;
  const double at = 0.2;
  if(Nt!=1) assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);
  double nu0 = 1.0;
  // if(argc>1) nu0 = atof(argv[1]);

  // ----------------------

  std::string path = "A2A_L"+std::to_string(Comp::N_REFINE)+"_Nt"+std::to_string(Nt)+"_at"+std::to_string(at)+"_nu0"+std::to_string(nu0)+".dat";
  path = "ov_"+path;
  const HighFive::File f(path.c_str(), HighFive::File::ReadOnly);

  DiracBase sigma;

  const int ell1=0, em1=0, a=1;
  const int ell2=0, em2=0, b=1;
  std::vector<Complex> Cab(Comp::Nt, 0.0);

  for(Idx x2=0; x2<base.n_sites; x2++){
    FermionMatrix eta;
    const auto r2 = base.sites[x2];
    // std::cout << "debug. Ylm2 = " << Ylm_real( ell2, em2, r2 ) << std::endl;

    for(int spin=0; spin<2; spin++){
      std::vector<Complex> vector = f.getDataSet("ix"+std::to_string(x2)+"/s"+std::to_string(spin)).read<std::vector<Complex>>();

      // for(auto elem : vector){
      //   std::cout << "# debug. " << std::abs(elem) << std::endl;
      // }
      // std::cout << "debug. size = " << vector.size() << std::endl;
      memcpy( eta[spin].data(), vector.data(), vector.size()*CD );
    }

    for(Idx x1=0; x1<Comp::N_SITES; x1++){
      const auto r1 = base.sites[x1];
      for(int s=0; s<Nt; s++){
        const MS Delta = eta.get_spinmatrix(s, x1);
        const Complex Csp = ( sigma[a] * Delta * sigma[b] * Delta.adjoint() ).trace();
        // std::cout << "# debug. " << s << " " << std::abs(Csp) << std::endl;
        Cab[s] += Ylm_real( ell1, em1, r1 )*Ylm_real( ell2, em2, r2 ) * Csp;
      }
    }
  }

  for(int s=0; s<Nt; s++){
    std::cout << std::setw(5) << s << " "
              << std::setw(15) << Cab[s].real() << " "
              << std::setw(15) << Cab[s].imag() << " "
              << std::endl;
  }


  return 0;
}
