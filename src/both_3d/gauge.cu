#include <typeinfo>

#include <iomanip>
#include <fstream>
#include <iostream>
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
// using Idx = long int;
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

// #define IsVerbose
// #define InfoForce
// #define InfoDelta

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
  constexpr int N_REFINE=1;
  constexpr int NPARALLEL_GAUGE=12; // 12
  constexpr int NPARALLEL_SORT=NPARALLEL_GAUGE; // 12

  constexpr int Nt=96; // 10

  constexpr int NS=2;

#ifdef IS_DUAL
  constexpr Idx N_SITES=20*N_REFINE*N_REFINE;
#else
  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
#endif

  constexpr Idx Nx=NS*N_SITES; // matrix size of DW
  constexpr Idx N=Nx*Nt; // matrix size of DW

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-8;
}

const std::string dir = "../../dats/";

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

// # include "integrator.h"
#include "hmc.h"

#include "obs.h"


// #include "../../integrator/geodesic.h"
// #include "s2.h"

// struct MultTableG {
//   const QfeLatticeS2& lattice;
//   std::vector<std::vector<int>> table;
//   const double TOL;

//   MultTableG( const QfeLatticeS2& lattice, const double TOL=1.0e-10 ) // TODO: TOL should be set from given geometry
//     : lattice(lattice)
//     , table( lattice.n_sites, std::vector<int>(lattice.G.size(), -1) )
//     , TOL(TOL)
//   {
//     for(int ir=0; ir<lattice.n_sites; ir++){ const Vec3 r = lattice.r[ir];
//       for(int g=0; g<lattice.G.size(); g++){ const Vec3 gr = lattice.G[g] * r;
//         for(int irp=0; irp<lattice.n_sites; irp++){ const Vec3 rp = lattice.r[irp];
//           if( (gr-rp).norm()<TOL ) { table[ir][g] = irp; continue; }}}}

//     // check
//     for(int ir=0; ir<lattice.n_sites; ir++){
//       for(int g=0; g<lattice.G.size(); g++){
//         if( table[ir][g] < 0 ) assert( false );
//       }}}

//   int operator()(const int g, const int ir) const { return table[ir][g]; }
// };



// TODO: Cusparse for SparseMatrix::act_gpu, probably defining handle in matpoly.h
// make 2 streams in V Vdag in square in Overlap
// all the operation on GPU in Overlap::operator()
// pseudofermion
// 3d
// __m256 to vectorize with AVX2


int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  // int device;
  // CUDA_CHECK(cudaGetDeviceCount(&device));
  // cudaDeviceProp device_prop[device];
  // cudaGetDeviceProperties(&device_prop[0], 0);
  // std::cout << "# dev = " << device_prop[0].name << std::endl;
  // CUDA_CHECK(cudaSetDevice(0));// "TITAN V"
  // std::cout << "# (GPU device is set.)" << std::endl;

  // ---------------------------------------

  std::cout << "# Nx = " << Comp::Nx << std::endl;

  // --------------------
  using Link = std::array<Idx,2>; // <int,int>;
  constexpr Idx Nx = Comp::Nx;
  constexpr int Nt = Comp::Nt;

#ifdef IS_DUAL
  using Base=S2Trivalent;
#else
  using Base=S2Simp;
#endif

  using Force=GaugeExt<Base,Nt,Comp::is_compact>;
  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  using Action=U1WilsonExt<Base>;

  using Rng=ParallelRngExt<Base,Nt>;

  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;

  // ----------------------

  // const double gsqR = 0.02;
  const double gsq = 0.05;
  // const double gsqR = 0.2;
  // double beta = 28.0; // 1.0/(gR*gR);
  // double beta = 1.0/gsqR; // 1.0/(gR*gR);
  // double at = base.mean_ell * 0.125;
  // double ratio = 1.0/2.0;
  double at = 0.05; // base.mean_ell * 0.125 * ratio;
  // double beta_t = beta_s; // 1.0/(gR*gR);
  if(Comp::Nt==1) at=0.;
  Action SW( gsq, at, base );
  std::cout << "# alat = " << base.mean_ell << std::endl;

  Gauge U(base);
  srand( time(NULL) );
  Rng rng(base, rand());
  // Rng rng(base);
  U.gaussian( rng, 0.2 );

  // std::string dir2="beta"+std::to_string(beta)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"ratio"+std::to_string(ratio)+"/";
  std::string dir2="gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"_v3/";
  std::filesystem::create_directory(dir2);


  // //--------------------------------


  const double tmax = 1.0; // 0.1
  // const int nsteps=250;
  const int nsteps=40;


  double rate, dH;
  bool is_accept;

  // const int kmax=2e7;
  const int kmax=1e2;
  const int interval=10;
  const int k_ckpoint=10;
  const int k_therm=1e1;

  Force pi(base);
  pi.gaussian( rng );


  int k_tmp=0;
  for(k_tmp=k_ckpoint; k_tmp<kmax; k_tmp+=k_ckpoint ){
    const std::string str_lat=dir2+"ckpoint_lat."+std::to_string(k_tmp);
    const std::string str_rng=dir2+"ckpoint_rng."+std::to_string(k_tmp);

    const bool bool_lat = std::filesystem::exists(str_lat);
    const bool bool_rng = std::filesystem::exists(str_rng);

    if(!(bool_lat&&bool_rng)) break;
  }
  k_tmp -= k_ckpoint;

  if(k_tmp>0){ // from existing
    std::cout << "read from k_tmp = " << k_tmp << std::endl;
    const std::string str_lat=dir2+"ckpoint_lat."+std::to_string(k_tmp);
    const std::string str_rng=dir2+"ckpoint_rng."+std::to_string(k_tmp);
    U.read( str_lat );
    rng.read( str_rng );
  }


  Force pi0=pi;
  Gauge U0=U;

  pi = pi0;
  U = U0;
  HMCPureGauge hmc(rng, &SW, U, pi, tmax, nsteps);


  if(k_tmp<=0){ // thermalize
    std::cout << "thermalize" << std::endl;

    for(int k=0; k<10; k++){
      Timer timer;
      hmc.run( rate, dH, is_accept, true );
      // if constexpr(Comp::is_compact) U.project();
      std::cout << "# dH : " << dH
                << " is_accept : " << is_accept << std::endl;
      // std::cout << "# HMC : " << timer.currentSeconds() << " sec" << std::endl;
    }

    for(int k=0; k<k_therm; k++){
      Timer timer;
      hmc.run( rate, dH, is_accept );
      std::cout << "# dH : " << dH
                << " is_accept : " << is_accept << std::endl;
    }
    k_tmp = 0;
  }


  std::vector<double> plaq_s0(Comp::Nt);
  double r_mean;

  for(int k=k_tmp+1; k<kmax; k++){
    Timer timer;
    hmc.run( rate, dH, is_accept);
    // if constexpr(Comp::is_compact) U.project();
    std::cout << "# dH : " << dH
              << " is_accept : " << is_accept << std::endl;
    r_mean += rate;
    // std::cout << "# HMC : " << timer.currentSeconds() << " sec" << std::endl;

    if(k%interval==0){

      std::string path = "plaq_ss_t_"+std::to_string(k)+".dat";
      std::ofstream ofs(dir2+path);
      // ofs << "# kmax = " << kmax << std::endl;

      // -----------------------------------------------------------------

#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL_GAUGE)
#endif
      for(int t=0; t<Comp::Nt; t++){
        int counter = 0;
        double tmp1 = 0.0;
        for(int s=0; s<Comp::Nt; s++){
          for(int i_face=0; i_face<base.n_faces; i_face++){
            tmp1 += U.plaquette_angle(s, U.lattice.faces[i_face]) * U.plaquette_angle(s+t, U.lattice.faces[i_face]) / std::pow(base.vols[i_face],2);
            counter++;
            tmp1 += U.plaquette_angle(s, U.lattice.faces[i_face]) * U.plaquette_angle(s-t, U.lattice.faces[i_face]) / std::pow(base.vols[i_face],2);
            counter++;
          }
        }
        tmp1 /= counter;
        plaq_s0[t] = tmp1;
      }

      for(int t=0; t<Comp::Nt; t++){
        ofs << plaq_s0[t] << std::endl;
      }

    }
    if(k%100==0){
      std::cout << "# k = " << k << std::endl;
    }

    if(k%k_ckpoint==0){
      const std::string str_lat=dir2+"ckpoint_lat."+std::to_string(k);
      const std::string str_rng=dir2+"ckpoint_rng."+std::to_string(k);
      U.ckpoint( str_lat );
      rng.ckpoint( str_rng );
    }

  }
  r_mean /= kmax;
  std::cout << "# r_mean = " << r_mean << std::endl;


  return 0;
}

