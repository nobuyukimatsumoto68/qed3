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


#define IS_DUAL
// #define IS_OVERLAP
#undef _OPENMP


namespace Comp{
  constexpr bool is_compact=false;

#ifdef IS_OVERLAP
  constexpr int NPARALLEL=4; // 12
  // constexpr int NPARALLEL2=1; // 12
  constexpr int NSTREAMS=4; // 4
#else
  constexpr int NPARALLEL=1; // 12
  // constexpr int NPARALLEL2=4; // 12
  constexpr int NSTREAMS=1; // 4
#endif
  constexpr int NPARALLEL3=1; // 12

  constexpr int N_REFINE=1;
  constexpr int NS=2;

  constexpr int Nt=128;

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
// #include "gauge.h"
#include "gauge_ext.h"
// #include "action.h"
// #include "action_ext.h"

#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cublas_api.h>
#include <cusolverDn.h>
using CuC = cuDoubleComplex;
#include "gpu_header.h"

// #include <cuComplex.h>
// #include <cuda_runtime.h>
// #include <cublas_v2.h>
// #include <cublas_api.h>
// #include <cusolverDn.h>
// using CuC = cuDoubleComplex;
// #include "gpu_header.h"

// ======================================

// #include "s2n.h"
// #include "rng.h"
// #include "gauge.h"
// #include "action.h"

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
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <cstdlib>
// #include <cassert>
// #include <algorithm>
// #include <cstdint>
// #include <complex>

// using Double = double;
// using Idx = std::int32_t;
// using Complex = std::complex<double>;

// // #define IS_DUAL
// #define IS_OVERLAP

// namespace Comp{
//   constexpr int NPARALLEL=12; // 12
//   constexpr int NSTREAMS=4; // 4

//   constexpr int N_REFINE=1;
//   constexpr int NS=2;

// #ifdef IS_DUAL
//   constexpr Idx N_SITES=20*N_REFINE*N_REFINE;
// #else
//   constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
// #endif

//   constexpr Idx N=NS*N_SITES; // matrix size of DW

//   // const double TOL=1.0e-9;
//   const double TOL_INNER=1.0e-15;
//   const double TOL_OUTER=1.0e-14;
// }

// // #define IsVerbose
// #define IsVerbose2
// // #define InfoForce
// #define InfoDelta

// const std::string dir = "/mnt/hdd_barracuda/qed3/dats/";

// #include "s2n_dual.h"
// #include "s2n_simp.h"
// #include "rng.h"
// #include "gauge.h"
// #include "action.h"


// #include <cuComplex.h>
// #include <cuda_runtime.h>
// #include <cublas_v2.h>
// #include <cublas_api.h>
// #include <cusolverDn.h>
// using CuC = cuDoubleComplex;
// #include "gpu_header.h"

// // ======================================

// // #include "timer.h"

// // #include "s2n_dual.h"
// // #include "s2n_simp.h"
// // #include "rng.h"
// // #include "gauge.h"
// // #include "force.h"
// // #include "action.h"
// // #include "sparse_matrix.h"
// // #include "dirac_dual.h"
// // #include "dirac_simp.h"
// // #include "sparse_dirac.h"
// // #include "matpoly.h"
// // #include "overlap.h"
// // #include "pseudofermion.h"

// // #include "integrator.h"
// // #include "hmc.h"

// #include "sparse_matrix.h"
// // #include "pseudofermion.h"
// #include "dirac_dual.h"
// #include "dirac_simp.h"

// #include "sparse_dirac.h"
// #include "matpoly.h"

// #include "overlap.h"


#include "valence.h"

// using Double = double;
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

// #ifdef IS_DUAL
//   using Lattice=S2Trivalent;
//   using Gauge=U1onS2<Lattice,false>;
//   using WilsonDirac=DiracS2Dual<Gauge>;
// #else
//   using Lattice=S2Simp;
//   using Gauge=U1onS2<Lattice,false>;
//   using WilsonDirac=DiracS2Simp<Gauge>;
// #endif
//   using Rng=ParallelRng<Lattice>;
//   using Overlap=Overlap<Gauge,WilsonDirac,Lattice>;

//   Lattice lattice(Comp::N_REFINE);

//   using Link = std::array<Idx,2>; // <int,int>;
//   constexpr Idx N = Comp::N;

//   // ----------------------

//   std::cout << "# lattice set. " << std::endl;

// #ifdef IS_OVERLAP
//   const double r = 1.0;
// #ifdef IS_DUAL
//   const double M5 = -1.6/2.0 * 0.5*3.0/2.0;
// #else
//   const double M5 = -1.6/2.0 * 0.5*(1.0 + std::sqrt( 5.0 + 2.0*std::sqrt(2.0) ));
// #endif
// #else
//   const double r = 1.0;
//   const double M5 = 0.0;
// #endif
//   WilsonDirac DW(lattice, 0.0, r, M5);

//   std::cout << "# DW set" << std::endl;

//   Gauge U(lattice);
//   Rng rng(lattice);
//   // U.gaussian( rng, 0.2 );

//   // ---------------------

//   Overlap Dov(DW, 31);
//   std::cout << "# Dov set; M5 = " << M5 << std::endl;
//   Dov.update(U);
//   std::cout << "# min max ratio: "
//             << Dov.lambda_min << " "
//             << Dov.lambda_max << " "
//             << Dov.lambda_min/Dov.lambda_max << std::endl;
//   std::cout << "# delta = " << Dov.Delta() << std::endl;

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
  // using Action=U1WilsonExt;

  using Rng=ParallelRngExt<Base,Nt>;

  // using WilsonDirac=DiracExt<Base>;
  using Fermion=DiracPf<WilsonDirac>;

  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;

  // using Gauge=U1onS2<false>;
  // // using Force=U1onS2<false>;
  // // using Action=U1Wilson;
  // // using Fermion=Dirac1fonS2;
  // // using HMC=HMC<Force,Gauge,Action,Fermion>;
  // // using Rng=ParallelRng;
  // using Lattice=S2Trivalent;
  // using Rng=ParallelRng<Lattice>;

  // ----------------------

  // const double gR = 10.0;
  // double beta = 4.0; // 1.0/(gR*gR);
  // Action SW(beta, beta);

  Gauge U(base);
  srand( time(NULL) );
  Rng rng(base, rand());
  // U.gaussian( rng, 0.2 );
  // const double M5 = -1.8;
  const double M5 = 0.0;
  const double c = 1.0;

  WilsonDirac DW(base, 0.0, 1.0, M5, c);
  std::cout << "# DW set. " << std::endl;

  Fermion D(DW);
  std::cout << "# D set. " << std::endl;

  D.update( U );
  std::cout << "# D updated. " << std::endl;

  // std::cout << "kappa" << std::endl;
  // for(double elem : DW.bd.kappa){
  //   std::cout << elem << std::endl;
  // }
  // std::cout << "ell" << std::endl;
  // for(double elem : base.ell){
  //   std::cout << elem << std::endl;
  // }
  // std::cout << "link_volume" << std::endl;
  // for(double elem : base.link_volume){
  //   std::cout << elem << std::endl;
  // }
  // return 0;


// #ifdef IS_OVERLAP
//   auto f_DHD = std::bind(&Overlap::sq_deviceAsyncLaunch, &Dov,
//                          std::placeholders::_1, std::placeholders::_2);
//   auto f_DH = std::bind(&Overlap::adj_deviceAsyncLaunch, &Dov,
//                         std::placeholders::_1, std::placeholders::_2);
//   LinOpWrapper M_DHD( f_DHD );
//   LinOpWrapper M_DH( f_DH );

//   MatPoly DHD;
//   DHD.push_back ( cplx(1.0), {&M_DHD} );
//   //
//   MatPoly DH;
//   DH.push_back ( cplx(1.0), {&M_DH} );
// #else
  // MatPoly DHD;
  // DHD.push_back ( cplx(1.0), {&Dov.M_DW, &Dov.M_DWH} );
  // //
  // MatPoly DH;
  // DH.push_back ( cplx(1.0), {&Dov.M_DWH} );
// #endif
  // auto f_DH = std::bind(&Fermion::adj_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);
  // @@@debug
  auto f_DH = std::bind(&Fermion::adj_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);
  //
  auto f_DHD = std::bind(&Fermion::sq_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_DH( f_DH );
  MatPoly DH; DH.push_back ( cplx(1.0), {&M_DH} );
  LinOpWrapper M_DHD( f_DHD );
  MatPoly DHD; DHD.push_back ( cplx(1.0), {&M_DHD} );

  // ---------------------

  std::cout << "# calculating src " << std::endl;

  // std::cout << "debug. n = " << base.n_sites*NS << std::endl;
  // FermionVector src1(base, rng);
  // FermionVector src(base, rng);
  FermionVector src1; // (base, Nt, rng);
  FermionVector src; // (base, Nt, rng);
  src1.set_pt_source(0, 0, 0);
  // DH.from_cpu<N>( src.field, src1.field );
  DH.from_cpu<N>( src.field, src1.field );

//   const Idx len = D.d_DW.is.size();

//   std::vector<Complex> v_coo;
//   v_coo.resize(len);
//   DW.coo_format( v_coo, U );
//   for(Idx i=0; i<len; i++){
//     // if(std::abs(v_coo[i])<1.0e-14) continue;
//     std::cout << "COO "
//               << std::setw(5) << D.d_DW.is[i] << " "
//               << std::setw(5) << D.d_DW.js[i] << " "
//               << std::setw(35) << v_coo[i] << " "
//               << std::endl;
//   }


//   std::cout << "csr." << std::endl;
//   D.d_DW.coo2csr_csrH( D.d_DW.v_csr, D.d_DW.v_csrH, v_coo );
//   {
//     Idx count = 0;
//     for(Idx i=0; i<len; i++){
//       // if(std::abs(D.d_DW.v_csr[i])<1.0e-14) continue;
//       std::cout << "CSR "
//                 << std::setw(5) << count << " "
//                 << std::setw(5) << D.d_DW.ell2em[i] << " "
//                 << std::setw(35) << D.d_DW.v_csr[i] << " "
//                 << std::endl;
//       count++;
//     }
//   }
//   {
//     Idx count = 0;
//     std::cout << "csrH." << std::endl;
//     for(Idx i=0; i<len; i++){
//       // if(std::abs(D.d_DW.v_csrH[i])<1.0e-14) continue;
//       std::cout << "CSRH " << std::setw(5) << count << " "
//                 << std::setw(5) << D.d_DW.ell2emT[i] << " "
//                 << std::setw(35) << D.d_DW.v_csrH[i] << " "
//                 << std::endl;
//       count++;
//     }
//   }

//   for(const auto elem : src){
//     std::cout << elem << std::endl;
//   }
// //   void coo2csr_csrH( std::vector<Complex>& v_csr,
// // 		     std::vector<Complex>& v_csrH,
// // 		     const std::vector<Complex>& v_coo) const {
// // #ifdef _OPENMP
// // #pragma omp parallel for num_threads(Comp::NPARALLEL)
// // #endif
// //     for(Idx ell=0; ell<len; ell++) {
// //       v_csr[ ell2em[ell] ] = v_coo[ell];
// //       v_csrH[ ell2emT[ell] ] = std::conj( v_coo[ell] );
// //     }
// //   }

//   return 1;





  // FermionVector sink(base, Nt, rng);
  FermionVector sink; // (base, Nt, rng);
  // Op.from_cpu<N>( sink.field, src.field );
  // FermionVector rc(base, rng);
  // rc.set_random();

  std::cout << "# calculating sink" << std::endl;

  // DHD.bicgstab<N>( sink.field, src.field, rc.field, 1.0e-3, 1e8, 1.0e-8 );
  DHD.solve<N>( sink.field, src.field );

  std::cout << "# done" << std::endl;



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
    const auto x0 = base.sites[0];
    for(int ix=0; ix<base.n_sites; ix++){
      const auto x1 = base.sites[ix];
      double len = Geodesic::geodesicLength(Geodesic::Pt(x0), Geodesic::Pt(x1));
      // std::cout << "len = " << len << std::endl;
      lengths.push_back(len);
    }
  }
#endif

  {
    std::string path = "prop_spacial_L"+std::to_string(Comp::N_REFINE)+"_Nt"+std::to_string(Nt)+".dat2";
#ifdef IS_DUAL
    path = "dual_"+path;
#endif
    std::ofstream ofs(path);

    // Idx counter=0;
    for(Idx ix=0; ix<base.n_sites; ix++) {
      {
        const auto elem = sink(0,ix,0);
        ofs << std::setw(25) << lengths[ix] << " "
            << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.real() << " "
            << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.imag() << std::endl;
        // << std::setw(25) << 1.0 * elem.real() << " "
        // << std::setw(25) << 1.0 * elem.imag() << std::endl;
      }
      {
        const auto elem = sink(0,ix,1);
        ofs << std::setw(25) << lengths[ix] << " "
            << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.real() << " "
            << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.imag() << std::endl;
        // << std::setw(25) << 1.0 * elem.real() << " "
        // << std::setw(25) << 1.0 * elem.imag() << std::endl;
      }
      // counter++;
    }
  }
  {
    std::string path = "prop_temporal_L"+std::to_string(Comp::N_REFINE)+"_Nt"+std::to_string(Nt)+".dat2";
#ifdef IS_DUAL
    path = "dual_"+path;
#endif
    std::ofstream ofs(path);

    // Idx counter=0;
    for(Idx s=0; s<Comp::Nt; s++) {
      {
        const auto elem = sink(s,0,0);
        ofs << std::setw(25) << s << " "
            << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.real() << " "
            << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.imag() << std::endl;
        // << std::setw(25) << 1.0 * elem.real() << " "
        // << std::setw(25) << 1.0 * elem.imag() << std::endl;
      }
      {
        const auto elem = sink(s,0,1);
        ofs << std::setw(25) << s << " "
            << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.real() << " "
            << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.imag() << std::endl;
        // << std::setw(25) << 1.0 * elem.real() << " "
        // << std::setw(25) << 1.0 * elem.imag() << std::endl;
      }
      // counter++;
    }
  }


  // ------------------


  return 0;

}

