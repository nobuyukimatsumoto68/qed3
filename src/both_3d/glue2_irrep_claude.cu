// ============================================================================
// glue2_irrep_claude.cu -- item 6 staging copy of glue2_msm_claude.cu.
//
// PLAN (see glue_msm_irrep_impl_plan_claude.md): replace the continuum Ylm
// weight in the operator with icosahedral-group irrep-adapted weights. The
// operator is already O = (1/N_face) sum_face theta(face) w(face); an irrep
// operator is just a different weight vector w^{Gamma,k}(face). Concretely:
//   (1) add Gauge::plaquette_angle_avg_weight(s, w) that dots the per-face
//       plaquette angles with an arbitrary precomputed weight table w[face];
//   (2) generate irrep weight tables offline (60 rotations of I as a face
//       permutation rep -> projectors P^Gamma -> orthonormal weight basis).
// For the current l<=3 data the only genuine within-l mixing is l=3 -> T_2 + G,
// which can also be done as a notebook recombination of the 7 l=3 operators.
// NOT a massive change to this file; the work is the offline weight tables.
// NOT YET IMPLEMENTED -- this copy is identical to glue2_msm for now.
// ============================================================================
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
  constexpr int Nt=128;
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

// const std::string dir = "../../dats/";
const std::string dir = "../../geometry/data/";


// // #define IsVerbose
// #define IsVerbose2
// // #define InfoForce
// #define InfoDelta

#include "timer.h"

// #include "../../integrator/geodesic.h"

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

// #include "sparse_dirac.h"
#include "sparse_dirac_claude.h"   // O(len) bucketing CSR build (was O(N*len)); -DCSR_VERIFY to check
// #include "matpoly.h"
#include "matpoly_claude.h"

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
      printf("Usage: ./a.out [gsq] [Nf] [nu0]\n");
      printf("  gsq  Wilson coupling squared (default: 8.0)\n");
      printf("  Nf   number of fermion flavors (default: 2)\n");
      printf("  nu0  mass parameter (default: 1.0)\n");
      return 0;
    }
  }

  double gsq = 8.0;
  if(argc>1) gsq = atof(argv[1]);
  int Nf = 2;
  if(argc>2) Nf = atoi(argv[2]);
  double nu0 = 1.0;
  if(argc>3) nu0 = atof(argv[3]);
  std::cout << "# gsq = " << gsq << " Nf = " << Nf << " nu0 = " << nu0 << std::endl;


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
  if(Nf==0){
    dir3="gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
    // dir4 gets a "_msm" suffix so the N_FLOW*n_lm operator output does not clobber
    // the 16-op data of glue2_claude.cu (dir3 input ckpoints are shared, unchanged).
    dir4="data_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"_msm/";
    std::cout << "dir3 = " << dir3 << std::endl;
  }
  else{
    dir3="Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
    dir4="data_Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"_msm/";
  }
  std::filesystem::create_directory(dir4);


  Gauge U(base);

  const int k_ckpoint=1;
  // const int kinit=10;
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

  // std::cout << "# kinit = " << kinit << " k_tmp = " << k_tmp << std::endl;

  // ---- item 7: multi-smearing variational basis via cumulative Wilson flow ----
  // Same Ylm channels measured at N_FLOW cumulative flow times
  //   {1*FLOW_INCR, 2*FLOW_INCR, ..., N_FLOW*FLOW_INCR}.
  // Cumulative: the running Uflow is advanced by FLOW_INCR each checkpoint
  // (FLOW_NSTEP integrator steps), never restarted from U. Gradient-flow analog
  // of multi-level APE/HYP smearing. Ref: Morningstar & Peardon, PRD 60 (1999) 034509.
  constexpr int N_FLOW = 3;
  constexpr double FLOW_INCR = 1.0;
  constexpr int FLOW_NSTEP = 100;

  // (ell, em) channel list; must match the analysis opset in glue_ylms3_claude.ipynb.
  const std::vector<std::array<int,2>> lm_set = {
    {0,0},
    {1,-1},{1,0},{1,1},
    {2,-2},{2,-1},{2,0},{2,1},{2,2},
    {3,-3},{3,-2},{3,-1},{3,0},{3,1},{3,2},{3,3},
  };
  const int n_lm = lm_set.size();
  const int nops = N_FLOW * n_lm; // operator index op = iflow*n_lm + ilm

  // #ifdef IS_FLOW
  // Flow flow(&SW, 1.5, 100);
  Flow flow(&SW, FLOW_INCR, FLOW_NSTEP);
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
// #endif

    // item 7: measure the n_lm Ylm channels at N_FLOW cumulative flow times.
    // obs[op][t], op = iflow*n_lm + ilm; flow(Uflow) advances Uflow in place.
    std::vector<std::vector<double>> obs( nops, std::vector<double>(Comp::Nt, 0.0) );
    for(int iflow=0; iflow<N_FLOW; iflow++){
      flow(Uflow);

      for(int ilm=0; ilm<n_lm; ilm++){
        const int ell = lm_set[ilm][0];
        const int em  = lm_set[ilm][1];
        const int op  = iflow*n_lm + ilm;

        for(int t=0; t<Comp::Nt; t++){
          obs[op][t] = Uflow.plaquette_angle_avg_Ylm_real(t, ell, em);
        }
      }
    }

    // std::vector<double> plaq_avg(Comp::Nt);
    // for(int t=0; t<Comp::Nt; t++) plaq_avg[t] = U.plaquette_angle_avg(t);

    // for(int t=0; t<Comp::Nt; t++) flow_plaq_avg[t] = Uflow.plaquette_angle_avg(t);
    // std::vector<double> plaq_avg_00(Comp::Nt);
    // std::vector<double> plaq_avg_1m1(Comp::Nt);
    // std::vector<double> plaq_avg_10(Comp::Nt);
    // std::vector<double> plaq_avg_11(Comp::Nt);
    // std::vector<double> plaq_avg_2m2(Comp::Nt);
    // std::vector<double> plaq_avg_2m1(Comp::Nt);
    // std::vector<double> plaq_avg_20(Comp::Nt);
    // std::vector<double> plaq_avg_21(Comp::Nt);
    // std::vector<double> plaq_avg_22(Comp::Nt);
    // std::vector<double> plaq_avg_3m3(Comp::Nt);
    // std::vector<double> plaq_avg_3m2(Comp::Nt);
    // std::vector<double> plaq_avg_3m1(Comp::Nt);
    // std::vector<double> plaq_avg_30(Comp::Nt);
    // std::vector<double> plaq_avg_31(Comp::Nt);
    // std::vector<double> plaq_avg_32(Comp::Nt);
    // std::vector<double> plaq_avg_33(Comp::Nt);
    //
    // std::vector<double> plaq_avg_temporal_00(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_1m1(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_10(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_11(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_2m2(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_2m1(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_20(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_21(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_22(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_3m3(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_3m2(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_3m1(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_30(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_31(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_32(Comp::Nt);
    // std::vector<double> plaq_avg_temporal_33(Comp::Nt);
    //
    // (named flow_plaq_avg_* vectors + obs_ptrs of the original now replaced by
    //  the obs[][] flow-checkpoint loop above; pristine version in glue2_claude.cu)

    {
      const std::string path=dir4+"F_corr."+std::to_string(k);
      std::ofstream ofs(path);
      ofs << std::scientific << std::setprecision(15);

      for(int dt=0; dt<Comp::Nt; dt++){
        Eigen::MatrixXd cdt_avg = Eigen::MatrixXd::Zero( nops, nops );

        for(int t=0; t<Comp::Nt; t++) {
          for(int i=0; i<nops; i++){
            for(int j=0; j<nops; j++){
              cdt_avg(i,j) += obs[i][t] * obs[j][(t+dt)%Comp::Nt];
            }
          }
        }
        for(int i=0; i<nops; i++){
          for(int j=0; j<nops; j++){
            cdt_avg(i,j) /= Comp::Nt; // item 2 fix: was a no-op `/`; needed so F_corr and F share the <..> scale for item 1 vacuum subtraction
          }}
        // ofs << dt << " ";
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
        for(int i=0; i<nops; i++) avg(i) += obs[i][t];
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

