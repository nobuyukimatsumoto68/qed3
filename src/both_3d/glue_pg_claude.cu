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


// #define IS_FLOW

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

#include "dirac_pf.h"
#include "overlap.h"

#include "flow.h"


// TODO: Cusparse for SparseMatrix::act_gpu, probably defining handle in matpoly.h
// make 2 streams in V Vdag in square in Overlap
// all the operation on GPU in Overlap::operator()
// pseudofermion
// 3d
// __m256 to vectorize with AVX2


int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) == "-h") {
      printf("Usage: ./a.out [gsq]\n");
      printf("  gsq  Wilson coupling squared (default: 12.0)\n");
      return 0;
    }
  }

  double gsq = 12.0;
  if(argc>1) gsq = atof(argv[1]);
  // int Nf = 2;
  // if(argc>2) Nf = atoi(argv[2]);
  // double nu0 = 1.0;
  // if(argc>3) nu0 = atof(argv[3]);
  std::cout << "# gsq = " << gsq << std::endl;


  // int igam=0;
  // if(argc>3) igam = atoi(argv[3]);
  // std::cout << "# igam = " << igam << " where 0=id., i=sigma_3" << std::endl;



  // int device;
  // CUDA_CHECK(cudaGetDeviceCount(&device));
  // cudaDeviceProp device_prop[device];
  // cudaGetDeviceProperties(&device_prop[0], 0);
  // std::cout << "# dev = " << device_prop[0].name << std::endl;
  // CUDA_CHECK(cudaSetDevice(0));// "TITAN V"
  // std::cout << "# (GPU device is set.)" << std::endl;

  // ---------------------------------------

  constexpr Idx N = Comp::N;
  constexpr int Nt = Comp::Nt;

#ifdef IS_DUAL
  using Base=S2Trivalent;
  // using WilsonDirac=DiracExt<Base, DiracS2Dual>;
#else
  using Base=S2Simp;
  // using WilsonDirac=DiracExt<Base, DiracS2Simp>;
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

  std::string dir3, dir4;
  // #ifdef Nf2
  // std::string dir2="gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
  // std::filesystem::create_directory(dir2);
  dir3="gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
  // dir3="Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
  dir4="data_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
  // dir4="data_Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";

  std::filesystem::create_directory(dir4);


  Gauge U(base);

  const int k_ckpoint=1;
  const int kmax=1e5;
  std::cout << "# dir3 = " << dir3 << std::endl;

  int k_tmp=0;
  {
    for(k_tmp=k_ckpoint; k_tmp<=kmax; k_tmp+=k_ckpoint ){
      const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k_tmp);
      const bool bool_lat = std::filesystem::exists(str_lat);
      if(!bool_lat) break;
    }
    k_tmp -= k_ckpoint;
  }

  // #ifdef IS_FLOW
  Flow flow(&SW);
  // #endif

#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL_GAUGE)
#endif
  for(int k=k_ckpoint; k<=k_tmp; k+=k_ckpoint ){
    std::cout << "# read from k = " << k << std::endl;
    const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k);
    U.read( str_lat );

// #ifdef IS_FLOW
    Gauge Uflow = U;
    flow(Uflow);
// #endif

    std::vector<double> plaq_avg(Comp::Nt);
    for(int t=0; t<Comp::Nt; t++) plaq_avg[t] = U.plaquette_angle_avg(t);

    std::vector<double> flow_plaq_avg(Comp::Nt);
    for(int t=0; t<Comp::Nt; t++) flow_plaq_avg[t] = Uflow.plaquette_angle_avg(t);

    std::vector<double> chair_avg(Comp::Nt);
    for(int t=0; t<Comp::Nt; t++) chair_avg[t] = U.chair_angle_avg(t);

    std::vector<std::vector<double>*> obs_ptrs;
    obs_ptrs.push_back( &plaq_avg );
    obs_ptrs.push_back( &flow_plaq_avg );
    obs_ptrs.push_back( &chair_avg );
    const int nops = obs_ptrs.size();

    {
      const std::string path=dir4+"F_corr."+std::to_string(k);
      std::ofstream ofs(path);
      ofs << std::scientific << std::setprecision(15);

      for(int dt=0; dt<Comp::Nt; dt++){
        // double cdt_avg1 = 0.0;
        // double cdt_avg2 = 0.0;
        // double cdt_avg3 = 0.0;
        // double cdt_avg4 = 0.0;
        Eigen::MatrixXd cdt_avg = Eigen::MatrixXd::Zero( nops, nops );

        for(int t=0; t<Comp::Nt; t++) {
          for(int i=0; i<nops; i++){
            for(int j=0; j<nops; j++){
              cdt_avg(i,j) += (*obs_ptrs[i])[t] * (*obs_ptrs[j])[(t+dt)%Comp::Nt];
            }
          }
          // cdt_avg1 += plaq_avg[t]*plaq_avg[(t+dt)%Comp::Nt];
          // cdt_avg2 += plaq_avg[t]*flow_plaq_avg[(t+dt)%Comp::Nt];
          // cdt_avg3 += flow_plaq_avg[t]*plaq_avg[(t+dt)%Comp::Nt];
          // cdt_avg4 += flow_plaq_avg[t]*flow_plaq_avg[(t+dt)%Comp::Nt];
        }
        for(int i=0; i<nops; i++){
          for(int j=0; j<nops; j++){
            cdt_avg(i,j) / Comp::Nt;
          }}
        // cdt_avg1 /= Comp::Nt;
        // cdt_avg2 /= Comp::Nt;
        // cdt_avg3 /= Comp::Nt;
        // cdt_avg4 /= Comp::Nt;
        ofs << dt << " ";
        for(int i=0; i<nops; i++){
          for(int j=0; j<nops; j++){
            ofs << cdt_avg(i,j) << " ";
          }}
        ofs << std::endl;
        // << cdt_avg2 << " " << cdt_avg3 << " " << cdt_avg4 << std::endl;
        // std::cout << dt << " " << cdt_avg << std::endl;
      }
    }

    {
      const std::string path=dir4+"F."+std::to_string(k);
      std::ofstream ofs(path);
      ofs << std::scientific << std::setprecision(15);

      // double avg1 = 0.0;
      // double avg2 = 0.0;
      Eigen::VectorXd avg = Eigen::VectorXd::Zero( nops );
      for(int t=0; t<Comp::Nt; t++) {
        for(int i=0; i<nops; i++) avg(i) += (*obs_ptrs[i])[t];
        // avg1 += plaq_avg[t];
        // avg2 += flow_plaq_avg[t];
      }
      for(int i=0; i<nops; i++) avg(i) /= Comp::Nt;
      // avg1 /= Comp::Nt;
      // avg2 /= Comp::Nt;
      // ofs << avg1 << " " << avg2 << std::endl;
      for(int i=0; i<nops; i++) ofs << avg(i) << " ";
      ofs << std::endl;
      // std::cout << avg << std::endl;
    }
  } // end for k (omp)

  // ------------------


  return 0;

}

