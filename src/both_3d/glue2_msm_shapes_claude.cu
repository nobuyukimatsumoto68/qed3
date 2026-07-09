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
#include <highfive/H5File.hpp>   // per-config HDF5 output (fast binary read for analysis)

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

  // Refinement level L is a compile-time flag: -DN_REFINE_CLI=2 (or 4) for L=2,4.
#ifdef N_REFINE_CLI
  constexpr int N_REFINE=N_REFINE_CLI;
#else
  constexpr int N_REFINE=1;
#endif
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

#include "icos_orbits_claude.h"    // native Ih orbit table on the lattice
#include "wilson_shapes_claude.h"  // generic shape (triangle/rectangle) orbit operators


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
  // config selection (added _claude): argv[4]=kmax_run cap, argv[5]=kmin, argv[6]=stride.
  // Disjoint (kmin,stride) across processes -> collision-free packing (per-k files).
  int kmax_run = 100000;
  if(argc>4) kmax_run = atoi(argv[4]);
  int kmin = 1;
  if(argc>5) kmin = atoi(argv[5]);
  int stride = 1;
  if(argc>6) stride = atoi(argv[6]);
  std::cout << "# gsq = " << gsq << " Nf = " << Nf << " nu0 = " << nu0
            << " kmax_run = " << kmax_run << " kmin = " << kmin << " stride = " << stride << std::endl;


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
    dir4="data_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";  // shared data_ dir (same as mesonic correlators)
    std::cout << "dir3 = " << dir3 << std::endl;
  }
  else{
    dir3="Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
    dir4="data_Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";  // shared data_ dir (same as mesonic correlators)
  }
  std::filesystem::create_directory(dir4);


  Gauge U(base);

  const int k_ckpoint=1;
  // const int kinit=10;
  const int kmax=1e5;

  int k_tmp=0;
  {
    // find the highest available config, scanning from kmin in steps of `stride` so ensembles that
    // are NOT numbered contiguously from 1 are detected (e.g. the free gauge set numbered 10,20,30,...).
    // (Previously scanned 1,2,3,... in steps of k_ckpoint=1 and stopped at the first gap -> k_tmp=0
    //  whenever ckpoint_lat.1 was absent, which silently measured nothing.)
    for(k_tmp=kmin; k_tmp<=kmax; k_tmp+=stride ){
      const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k_tmp);
      const bool bool_lat = std::filesystem::exists(str_lat);
      if(!bool_lat) break;
    }
    k_tmp -= stride;
  }

  // std::cout << "# kinit = " << kinit << " k_tmp = " << k_tmp << std::endl;

  // ---- single Wilson-flow smearing (matches non-_claude glue2.cu: tmax=1.0, 100 steps) ----
  constexpr double FLOW_TMAX = 1.0;
  constexpr int FLOW_NSTEP = 100;
  Flow flow(&SW, FLOW_TMAX, FLOW_NSTEP);

  // ---- shape operator basis: icosahedral orbits of spatial Wilson-loop shapes ----
  // LINEAR (F_12) operators; Y_lm tower ell=0..3. l=0 = orbit-summed total flux = monopole/
  // topological mode (~constant within a topological sector); it is NOT a propagating glueball,
  // but with vacuum subtraction C_ij -= <O_i><O_j> in the GEVP analysis its constant part is
  // removed, so it is re-included as an extra variational operator (near-singular directions are
  // rtol-pruned in inv_sqrt_sym). Analysis must be run with drop_l0=0 to keep these columns.
  // {0,0} MUST come first so the analysis (l=0..3, op%n_lm=ilm) indexing matches.
  IcosOrbits<Base> orb( base );
  WilsonShapes<Base> shp( base, orb );
  // -DNO_FACE_SIGN: measure the RAW-orientation operator (matches glue_ylms3's plaquette_angle_avg_Ylm_real),
  // written to a DISTINCT prefix so it does not clobber the production (face-sign) h5.
#ifdef NO_FACE_SIGN
  shp.use_face_sign = false;
  const std::string H5PREFIX = "glue_msm_shapes_nofs";
#else
  const std::string H5PREFIX = "glue_msm_shapes";
#endif
  using Inst = typename WilsonShapes<Base>::Instance;
  std::vector<std::vector<Inst>> orbits;
  {
    // five shape types: triangle + {rect, twisted rect, figure-8, twisted figure-8}
    // TWISTED shapes REMOVED: the twist Phi0-Phi1 is a flux DIFFERENCE that cancels the leading
    // (smooth) F_12 mode, so its GEVP ground is a spurious sub-sqrt(2) short-distance lattice
    // artifact (free L=1: eigenvector ~90% twisted, ground 1.05 vs true sqrt2). Keep triangle +
    // untwisted rect + figure-8 (Phi0+Phi1), which reproduce the free-limit sqrt(2) cleanly.
    std::vector<std::vector<Inst>> all[4] = {
      shp.orbits_from( shp.triangles() ),
      shp.orbits_from( shp.rectangles() ),
      // shp.orbits_from( shp.twisted_rectangles() ),   // removed (twisted artifact)
      shp.orbits_from( shp.figure8s() ),
      // shp.orbits_from( shp.twisted_figure8s() ),      // removed (twisted artifact)
      shp.orbits_from( shp.three_triangles() ),         // NEW: central triangle + 2 edge-neighbors
    };
    for(int is=0; is<4; is++) for(auto& o : all[is]) orbits.push_back(std::move(o));
  }
  const int n_orbits = (int)orbits.size();
  const bool SQUARED = false;

  const std::vector<std::array<int,2>> lm_set = {
    {0,0},
    {1,-1},{1,0},{1,1},
    {2,-2},{2,-1},{2,0},{2,1},{2,2},
    // {3,-3},{3,-2},{3,-1},{3,0},{3,1},{3,2},{3,3},   // l=3 DROPPED for production (disk): F ground is l=1, first excited l=2
  };
  const int n_lm = lm_set.size();
  const int nops = n_orbits * n_lm; // op = iorbit*n_lm + ilm
  std::cout << "# n_orbits = " << n_orbits << " n_lm = " << n_lm << " nops = " << nops << std::endl;

  if(k_tmp > kmax_run) k_tmp = kmax_run;
  // serial over configs k; parallelism is ensemble-level (one process per Nf).
  for(int k=kmin; k<=k_tmp; k+=stride ){
    // resume-safe: skip configs already fully measured (h5 with "complete" flag)
    const std::string h5path = dir4+H5PREFIX+"."+std::to_string(k)+".h5"; // distinct prefix in shared dir
    {
      bool done=false;
      if(std::filesystem::exists(h5path)){
        try { HighFive::File f(h5path, HighFive::File::ReadOnly); done = f.exist("complete"); } catch(...) {}
      }
      if(done) continue;
    }
    std::cout << "# read from k = " << k << std::endl;
    const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k);
    U.read( str_lat );

    Gauge Uflow = U;
    flow(Uflow);   // single flow (measure on the flowed field)

    // obs[op][t], op = iorbit*n_lm + ilm ; linear shape operators (F_12)
    std::vector<std::vector<double>> obs( nops, std::vector<double>(Comp::Nt, 0.0) );
    for(int iorbit=0; iorbit<n_orbits; iorbit++){
      for(int ilm=0; ilm<n_lm; ilm++){
        const int ell = lm_set[ilm][0];
        const int em  = lm_set[ilm][1];
        const int op  = iorbit*n_lm + ilm;
        for(int t=0; t<Comp::Nt; t++){
          obs[op][t] = shp.op( Uflow, t, orbits[iorbit], ell, em, SQUARED );
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

    // correlator matrix C(dt)[i][j] = (1/Nt) sum_t obs[i][t] obs[j][t+dt], and one-point <O_i>.
    // Only the small-separation window dt = 0..TMAX_CORR is stored: the GEVP uses dt up to tcut(~5),
    // and the backward fold slice is recovered LOSSLESSLY from the transpose of a stored forward
    // slice via the periodicity identity  C_ij(Nt-d) = C_ji(d)  (see glue_gevp_analysis_claude.cu
    // fold). This cuts both the correlator cost (Nt^2 -> Nt*TMAX) and the h5 size (~Nt/TMAX).
    constexpr int TMAX_CORR = 16;
    const int nsep = std::min(TMAX_CORR + 1, Comp::Nt); // stored separations dt = 0..nsep-1
    std::vector<std::vector<double>> Fcorr( nsep, std::vector<double>(nops*nops, 0.0) );
    for(int dt=0; dt<nsep; dt++){
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
          Fcorr[dt][i*nops+j] = cdt_avg(i,j) / Comp::Nt; // F_corr and F share the <..> scale for vacuum subtraction
        }}
    }

    std::vector<double> F1( nops, 0.0 );
    for(int t=0; t<Comp::Nt; t++) {
      for(int i=0; i<nops; i++) F1[i] += obs[i][t];
    }
    for(int i=0; i<nops; i++) F1[i] /= Comp::Nt;

    // write per-config HDF5 (F_corr: Nt x nops^2, F: nops); "complete" flag written LAST
    HighFive::File h5( h5path, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate );
    h5.createDataSet( "F_corr", Fcorr );
    h5.createDataSet( "F", F1 );
    h5.createDataSet( "n_lm", std::vector<int>{n_lm} );  // (l,m) count so the analysis auto-adapts to the l-tower
    h5.createDataSet( "complete", std::vector<int>{1} );
  } // end for k

  // ------------------


  return 0;

}

