#include <typeinfo>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <filesystem>
#include <chrono>
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

// #define GAUGE_TRSF


namespace Comp{
  constexpr bool is_compact=false;

  // d_DW.update() is always done independently
#ifdef IS_OVERLAP
  constexpr int NPARALLEL_DUPDATE=1;
  // constexpr int NPARALLEL=12; // 12
  // constexpr int NSTREAMS=4; // 4
  constexpr int NPARALLEL=1; // 12
  constexpr int NSTREAMS=1; // 4
#else
  constexpr int NPARALLEL_DUPDATE=12;
  constexpr int NPARALLEL=1; // 12
  constexpr int NSTREAMS=12; // for grad loop
#endif
  constexpr int NPARALLEL_GAUGE=1; // 12
  constexpr int NPARALLEL_SORT=1; // 12

  constexpr int N_REFINE=1;
  constexpr int NS=2;

  // constexpr int Nt=24;
  // constexpr int Nt=48; // add 2
  // constexpr int Nt=64;
  // constexpr int Nt=96; // add 4
  // constexpr int Nt=120;
  // constexpr int Nt=144; // add 8
  // constexpr int Nt=168;

  // constexpr int Nt=24;
  // constexpr int Nt=192;
  constexpr int Nt=96;
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

// #include "dirac_pf.h"
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

  double gsq = 1.0;
  if(argc>1) gsq = atof(argv[1]);
  int Nf = 2;
  if(argc>2) Nf = atoi(argv[2]);
  // std::cout << "# gsq = " << gsq << " Nf = " << Nf << std::endl;
  double nu0 = 1.0;
  if(argc>3) nu0 = atof(argv[3]);
  std::cout << "# gsq = " << gsq << " Nf = " << Nf << " nu0 = " << nu0 << std::endl;

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();


  // int igam=0;
  // if(argc>3) igam = atoi(argv[3]);
  // std::cout << "# igam = " << igam << " where 0=id., i=sigma_3" << std::endl;



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
  // const double T = 24;
  const double at = 0.2; // T/Comp::Nt;
  if(Nt!=1) assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);


  // const double gsq = 0.05;
  // double at = 0.05; // base.mean_ell * 0.125 * ratio;
  // if(Comp::Nt==1) at=0.;
  Action SW( gsq, at, base );
  std::cout << "# alat = " << base.mean_ell << std::endl;

  // std::string dir3, dir4;
  // // #ifdef Nf2
  // dir3="Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
  // dir4="data_Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
  std::string dir3, dir4;
  // #ifdef Nf2
  dir3="Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
  // dir3="Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
  dir4="data_Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
  // dir4="data_Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";


  std::filesystem::create_directory(dir4);


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
  const double M5 = -1.5;
#else
  const double M5 = -1.0;
#endif

  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);
  std::cout << "# DW set. " << std::endl;

  using Fermion=Overlap<WilsonDirac>;
  Fermion D(DW, 21);
  std::cout << "# D set. " << std::endl;
#else
  const double M5 = 0.0;
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);
  std::cout << "# DW set. " << std::endl;

  using Fermion=DiracPf<WilsonDirac>;
  Fermion D(DW);
  std::cout << "# D set. " << std::endl;
#endif

  D.update( U );
  std::cout << "# D updated. " << std::endl;

  auto f_D = std::bind(&Fermion::mult_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);
  auto f_DDH = std::bind(&Fermion::DDH_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);

  auto f_DH = std::bind(&Fermion::adj_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);
  auto f_DHD = std::bind(&Fermion::DHD_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);

  LinOpWrapper M_D( f_D );
  MatPoly op_D; op_D.push_back ( cplx(1.0), {&M_D} );
  LinOpWrapper M_DDH( f_DDH );
  MatPoly op_DDH; op_DDH.push_back ( cplx(1.0), {&M_DDH} );

  LinOpWrapper M_DH( f_DH );
  MatPoly op_DH; op_DH.push_back ( cplx(1.0), {&M_DH} );
  LinOpWrapper M_DHD( f_DHD );
  MatPoly op_DHD; op_DHD.push_back ( cplx(1.0), {&M_DHD} );

  // ---------------------

  std::cout << "# calculating src " << std::endl;

#ifdef GAUGE_TRSF
  FermionVector gauge; // (base, Nt, rng);
  gauge.set_random_gauge(rng);
#endif

  FermionVector src0; // (base, Nt, rng);
  FermionVector DHsrc0; // (base, Nt, rng);
  // FermionVector Dsrc0; // (base, Nt, rng);

  FermionVector src1; // (base, Nt, rng);
  FermionVector DHsrc1; // (base, Nt, rng);
  // FermionVector Dsrc1; // (base, Nt, rng);

  // pt source
  const Idx iy = 0;
  src0.set_pt_source(0, iy, 0);
  src1.set_pt_source(0, iy, 1);

#ifdef GAUGE_TRSF
  DW.bd.random_gauge_trsf(base,rng);

  src1.gauge_trsf( gauge );
  U.gauge_trsf( gauge );
  D.update( U );
#endif


  // op_D.from_cpu<N>( Dsrc0.field, src0.field );
  // op_D.from_cpu<N>( Dsrc1.field, src1.field );

  FermionVector Dinv0; // (base, Nt, rng);
  FermionVector Dinv1; // (base, Nt, rng);
  // FermionVector DHinv0; // (base, Nt, rng);
  // FermionVector DHinv1; // (base, Nt, rng);

  // FermionVector DHDinv0; // (base, Nt, rng);
  // FermionVector DHDinv1; // (base, Nt, rng);

  std::cout << "# calculating sink" << std::endl;


  // free
  op_DH.from_cpu<N>( DHsrc0.field, src0.field );
  op_DH.from_cpu<N>( DHsrc1.field, src1.field );
  op_DHD.solve<N>( Dinv0.field, DHsrc0.field );
  op_DHD.solve<N>( Dinv1.field, DHsrc1.field );

#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL_SORT) schedule(static)
#endif
  for(int igam=0; igam<=3; igam++){
    const std::string path=dir4+"meson_"+std::to_string(igam)+"_corr.free";
    // const std::string path=dir4+"meson_free_corr."+std::to_string(k);
    std::ofstream ofs(path);
    ofs << std::scientific << std::setprecision(15);

    for(Idx s=0; s<Comp::Nt; s++) {
      const auto Dinv_s0_O0 = Dinv0(s,iy,0);
      const auto Dinv_s1_O0 = Dinv0(s,iy,1);

      const auto Dinv_s0_O1 = Dinv1(s,iy,0);
      const auto Dinv_s1_O1 = Dinv1(s,iy,1);

      MS Dinv_sO; // << inits in row major way
      Dinv_sO << Dinv_s0_O0, Dinv_s0_O1, Dinv_s1_O0, Dinv_s1_O1;

      MS DinvH_sO = Dinv_sO.adjoint(); // << inits in row major way
      Complex tr = (DinvH_sO * DW.sigma[igam] * Dinv_sO * DW.sigma[igam] ).trace();

      // std::cout << s << " " <<  tr.real() << " " << tr.imag() << std::endl;
      ofs << s << " " <<  tr.real() << " " << tr.imag() << std::endl;
    }
  }




  const int k_ckpoint=10;
  const int kmax=1e5;

  int k_tmp=0;
  {
    for(k_tmp=k_ckpoint; k_tmp<=kmax; k_tmp+=k_ckpoint ){
      const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k_tmp);
      const bool bool_lat = std::filesystem::exists(str_lat);
      if(!bool_lat) break;
    }
    k_tmp -= k_ckpoint;
  }


  for(int k=k_ckpoint; k<=k_tmp; k+=k_ckpoint ){
    std::cout << "# read from k = " << k << std::endl;
    const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k);
    U.read( str_lat );

    D.update( U );

    op_DH.from_cpu<N>( DHsrc0.field, src0.field );
    op_DH.from_cpu<N>( DHsrc1.field, src1.field );
    op_DHD.solve<N>( Dinv0.field, DHsrc0.field );
    op_DHD.solve<N>( Dinv1.field, DHsrc1.field );

    std::cout << "# done" << std::endl;

#ifdef GAUGE_TRSF
    sink.gauge_trsf( gauge, -1 );
#endif


#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL_SORT) schedule(static)
#endif
    for(int igam=0; igam<=3; igam++){

      const std::string path=dir4+"meson_"+std::to_string(igam)+"_corr."+std::to_string(k);
      std::ofstream ofs(path);
      ofs << std::scientific << std::setprecision(15);

      for(Idx s=0; s<Comp::Nt; s++) {
        const auto Dinv_s0_O0 = Dinv0(s,iy,0);
        const auto Dinv_s1_O0 = Dinv0(s,iy,1);

        const auto Dinv_s0_O1 = Dinv1(s,iy,0);
        const auto Dinv_s1_O1 = Dinv1(s,iy,1);

        MS Dinv_sO; // << inits in row major way
        Dinv_sO << Dinv_s0_O0, Dinv_s0_O1, Dinv_s1_O0, Dinv_s1_O1;

        MS DinvH_sO = Dinv_sO.adjoint(); // << inits in row major way
        Complex tr = (DinvH_sO * DW.sigma[igam] * Dinv_sO * DW.sigma[igam] ).trace();

        // std::cout << s << " " <<  tr.real() << " " << tr.imag() << std::endl;
        ofs << s << " " <<  tr.real() << " " << tr.imag() << std::endl;
      }
    }

  } // end for k
















  //   std::vector<double> thetas;
//   std::vector<double> phis;
//   std::vector<double> lengths;
// #ifdef IS_DUAL
//   {
//     std::string dir = "/mnt/hdd_barracuda/qed3/dats/";
//     std::vector<Geodesic::V3> sites;
//     {
//       std::ifstream file(dir+"pts_dual_n"+std::to_string(Comp::N_REFINE)+"_singlepatch.dat");

//       std::string str;
//       while (std::getline(file, str)){
//         std::istringstream iss(str);
//         double v1, v2, v3;
//         iss >> v1;
//         iss >> v2;
//         iss >> v3;
//         sites.push_back( Geodesic::V3(v1, v2, v3) );
//       }
//     }
//     const auto x0 = sites[0];
//     for(const auto& elem : sites){
//       double len = Geodesic::geodesicLength(Geodesic::Pt(x0), Geodesic::Pt(elem));
//       // std::cout << "len = " << len << std::endl;
//       lengths.push_back(len);
//       thetas.push_back( Geodesic::projectionS2(elem)[0] );
//       phis.push_back( Geodesic::projectionS2(elem)[1] );
//     }
//   }
//   // double alat;
//   // {
//   //   std::string dir = "/mnt/hdd_barracuda/qed3/dats/";
//   //   std::ifstream file(dir+"alat_n"+std::to_string(Comp::N_REFINE)+"_singlepatch.dat");

//   //   std::string str;
//   //   std::getline(file, str);
//   //   std::istringstream iss(str);
//   //   iss >> alat;
//   // }
// #else
//   {
//     const auto x0 = base.sites[0];
//     for(int ix=0; ix<base.n_sites; ix++){
//       const auto x1 = base.sites[ix];
//       double len = Geodesic::geodesicLength(Geodesic::Pt(x0), Geodesic::Pt(x1));
//       // std::cout << "len = " << len << std::endl;
//       lengths.push_back(len);
//       thetas.push_back( Geodesic::projectionS2(x1)[0] );
//       phis.push_back( Geodesic::projectionS2(x1)[1] );
//     }
//   }
// #endif

//   const double width = 0.05;

//   double factor = at*base.mean_ell;
// #ifdef IS_DUAL
//   if(Comp::Nt==1) factor = 1.0; // /base.mean_ell;
// #else
//   if(Comp::Nt==1) factor = base.mean_ell;
// #endif

//   {
//     std::string path = "prop_spacial_L"+std::to_string(Comp::N_REFINE)+"_Nt"+std::to_string(Nt)+"_T"+std::to_string(T)+".dat";
// #ifdef IS_DUAL
//     path = "dual_"+path;
// #endif
// #ifdef IS_OVERLAP
//     path = "ov_"+path;
// #endif
// #ifdef IS_DAGGER
//     path = "dagger_"+path;
// #endif
//     std::ofstream ofs(path);

//     // Idx counter=0;
//     for(Idx ix=0; ix<base.n_sites; ix++) {
//       if( phis[ix]>width || phis[ix]<0. ) continue;
//       {
//         const auto elem = sink(0,ix,0);
//         ofs << std::setw(25) << thetas[ix] << " "
//           // ofs << std::setw(25) << lengths[ix] << " "
//             // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.real() << " "
//             // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.imag() << std::endl;
//         << std::setw(25) << 1.0 * elem.real() / factor << " "
//         << std::setw(25) << 1.0 * elem.imag() / factor << std::endl;
//       }
//       {
//         const auto elem = sink(0,ix,1);
//         ofs << std::setw(25) << thetas[ix] << " "
//           // ofs << std::setw(25) << lengths[ix] << " "
//             // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.real() << " "
//             // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.imag() << std::endl;
//             << std::setw(25) << 1.0 * elem.real() / factor << " "
//             << std::setw(25) << 1.0 * elem.imag() / factor << std::endl;
//       }
//       // counter++;
//     }
//   }
// //   {
// //     std::string path = "prop_temporal_L"+std::to_string(Comp::N_REFINE)+"_Nt"+std::to_string(Nt)+"_T"+std::to_string(T)+".dat";
// // #ifdef IS_DUAL
// //     path = "dual_"+path;
// // #endif
// // #ifdef IS_OVERLAP
// //     path = "ov_"+path;
// // #endif
// // #ifdef IS_DAGGER
// //     path = "dagger_"+path;
// // #endif
// //     std::ofstream ofs(path);

// //     // Idx counter=0;
// //     for(Idx s=0; s<Comp::Nt; s++) {
// //       {
// //         const auto elem = sink(s,0,0);
// //         ofs << std::setw(25) << at*s << " "
// //           // ofs << std::setw(25) << s << " "
// //             // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.real() << " "
// //             // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.imag() << std::endl;
// //         << std::setw(25) << 1.0 * elem.real() / factor << " "
// //         << std::setw(25) << 1.0 * elem.imag() / factor << std::endl;
// //       }
// //       {
// //         const auto elem = sink(s,0,1);
// //         ofs << std::setw(25) << at*s << " "
// //           // ofs << std::setw(25) << s << " "
// //           // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.real() << " "
// //             // << std::setw(25) << 1.0/std::pow(base.mean_ell,2) * elem.imag() << std::endl;
// //         << std::setw(25) << 1.0 * elem.real() / factor << " "
// //         << std::setw(25) << 1.0 * elem.imag() / factor << std::endl;
// //       }
// //       // counter++;
// //     }
// //   }





  // ------------------

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();


  return 0;

}

