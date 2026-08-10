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

// glue_f2_claude.cu
//
// F_{\mu\nu}^2 action-density (0^{++} scalar-glueball) measurement.
// Copy of glue2_msm_claude.cu; the per-timeslice interpolator is swapped from the
// LINEAR spatial plaquette angle (F_{12}, plaquette_angle_avg_Ylm_real) to the
// QUADRATIC local action density, split into two channels:
//   ch 0 = magnetic B^2 = F_{12}^2   (SW.density_Ylm_spatial, over faces)
//   ch 1 = electric E^2 = F_{0i}^2   (SW.density_Ylm_temporal, over links)
// each measured at N_FLOW cumulative Wilson-flow times and Y_{lm}-projected.
// Operator index op = iflow*(NCH*n_lm) + ich*n_lm + ilm; must match glue_f2_claude.ipynb.
// The density is the local Wilson-action integrand (action_ext_claude.h), so B^2_{00}+E^2_{00}
// is proportional to the total action density; VEV subtraction is done in analysis.
// Ref: N. Matsumoto et al., qed3_v2-6.pdf, Sec. IV, Eqs. (IV.24)-(IV.35).

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
  // Same #ifndef-guard convention as the L2/L4 stochastic Y_lm drivers.
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
#include "action_ext_claude.h"   // adds density_Ylm_spatial (B^2) / density_Ylm_temporal (E^2)

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

// #include "flow.h"
#include "flow_claude.h"   // adds -DFLOW_FULL switch (full 3D flow via get_force vs spatial-only)

#include "icos_orbits_claude.h"    // native Ih orbit table on the lattice
#include "wilson_shapes_claude.h"  // generic shape (triangle/rectangle) orbit operators


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
  // Optional cap on the number of configs processed (default: all). Handy for a
  // quick smoke test / flow-time tuning on a handful of configs.
  int kmax_run = 100000;
  if(argc>4) kmax_run = atoi(argv[4]);
  // Config subsampling / process-level packing: process k = kmin, kmin+stride, ...
  // Disjoint (kmin,stride) across processes -> collision-free parallelism (each writes
  // per-k files). Default kmin=1 stride=1 = all configs.
  int kmin = 1;
  if(argc>5) kmin = atoi(argv[5]);
  int stride = 1;
  if(argc>6) stride = atoi(argv[6]);
  // arg7 = explicit ensemble dir (redo production `..._hb*` names, which the legacy dir3
  // construction below cannot form; plan redo_obs_measurement_impl_plan_claude.md Sec 5)
  std::string ens_dir = "";
  if(argc>7) ens_dir = argv[7];
  std::cout << "# gsq = " << gsq << " Nf = " << Nf << " nu0 = " << nu0
            << " kmax_run = " << kmax_run << " kmin = " << kmin << " stride = " << stride
            << " ens_dir = " << (ens_dir.empty() ? std::string("<legacy dir3>") : ens_dir) << std::endl;


  // ---------------------------------------

  constexpr Idx N = Comp::N;
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
  const double at = 0.2; // production a_t
  if(Nt!=1) assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);

  Action SW( gsq, at, base );
  std::cout << "# alat = " << base.mean_ell << std::endl;

  std::string dir3, dir4;
  if(!ens_dir.empty()){
    // explicit ensemble dir; output h5 shares data_<basename>/ with the fermionic obs
    dir3 = ens_dir;
    if(dir3.back()!='/') dir3 += "/";
    std::string bn = dir3.substr(0, dir3.size()-1);
    const std::size_t pos = bn.find_last_of('/');
    if(pos!=std::string::npos) bn = bn.substr(pos+1);
    dir4 = "data_"+bn+"/";
    std::cout << "dir3 = " << dir3 << std::endl;
  }
  else if(Nf==0){
    dir3="gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
    // "_f2" suffix so the action-density output does not clobber the _msm (linear F_{12}) data;
    // dir3 input ckpoints are shared, unchanged.
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
  const int kmax=1e5;

  int k_tmp=0;
  {
    // scan from kmin in steps of `stride` so non-contiguously-numbered ensembles are detected
    // (e.g. the free gauge set numbered 10,20,30,...). Previously scanned 1,2,3,... and stopped at
    // the first gap -> k_tmp=0 whenever ckpoint_lat.1 was absent (silently measured nothing).
    for(k_tmp=kmin; k_tmp<=kmax; k_tmp+=stride ){
      const std::string str_lat=dir3+"ckpoint_lat."+std::to_string(k_tmp);
      const bool bool_lat = std::filesystem::exists(str_lat);
      if(!bool_lat) break;
    }
    k_tmp -= stride;
  }
  if(k_tmp > kmax_run) k_tmp = kmax_run; // cap for smoke test / flow tuning

  // ---- MULTI-FLOW smearing levels: the SAME shapes measured at N_FLOW cumulative flow times ----
  // Multi-smearing variational basis: C. Morningstar & M. Peardon, arXiv:hep-lat/9901004
  // (Phys. Rev. D 60, 034509 (1999)).  Wilson flow as the smearing: M. Luescher, arXiv:1006.4518.
  // The flow acts IN PLACE, so ONE trajectory 0 -> 0.5 -> 1.0 -> 2.0 -> 3.0 yields every level:
  // only the (cheap) shape evaluations multiply; the integration cost is just the longest time.
  // dt is UNCHANGED at 0.01, so the t=2.0 level reproduces the old single-flow data exactly
  // (same total step count, same integrator) -- used as the correctness check.
  // OLD single level (kept for reference):
  // constexpr double FLOW_TMAX = 2.0;   // bumped 1.0 -> 2.0
  // constexpr int FLOW_NSTEP = 200;     // keep dt=0.01 (100 -> 200)
  // Flow flow(&SW, FLOW_TMAX, FLOW_NSTEP);
  // MULTI-FLOW SET BACK TO ONE LEVEL (t=2.0 = the production smearing): the 4-level basis was tested
  // on L3 Nf2 g1.5 (202 cfg) and gave NO improvement -- at 14 ops the error was unchanged
  // (F l=1 0.0453 -> 0.0449) and at 21/28 ops the GEVP collapsed to nan, because the levels are the
  // same 7 shapes merely smoothed differently (redundant -> near-singular metric).  Kept as N_FLOW=1
  // so the machinery stays available; FLOW_T={2.0} reproduces the old single flow EXACTLY.
  constexpr double FLOW_DT = 0.01;                          // UNCHANGED dt
  constexpr int N_FLOW = 1;
  const double FLOW_T[N_FLOW] = { 2.0 };                    // CUMULATIVE flow times

  // ---- shape operator basis: icosahedral orbits of spatial Wilson-loop shapes ----
  // SQUARED (F^2 / 0++) operators; full Y_lm tower ell=0..3 (l=0 KEPT = scalar signal).
  IcosOrbits<Base> orb( base );
  WilsonShapes<Base> shp( base, orb );
  using Inst = typename WilsonShapes<Base>::Instance;
  std::vector<std::vector<Inst>> orbits;
  std::array<int,9> shape_sizes{};   // # icosahedral orbits per shape type (for consolidation)
  {
    // 9 shape types (TWISTED shapes removed). Order fixed: triangle, rect, fig8, three-tri, star,
    // trio (star-center), five-six (site contour), four-five, three-four.
    // The last two are the site contour with 1 / 2 non-adjacent faces removed (NM 2026-08-08).
    // They are SQUARED-ONLY: for the linear operator they collapse to multiples of the five-six
    // (see wilson_shapes_claude.h), but squared they add sum_k theta_k^2 and sum_{non-adj} th_k th_m.
    std::vector<std::vector<Inst>> all[9] = {
      shp.orbits_from( shp.triangles() ),
      shp.orbits_from( shp.rectangles() ),
      shp.orbits_from( shp.figure8s() ),
      shp.orbits_from( shp.three_triangles() ),         // central triangle + 2 edge-neighbors
      shp.orbits_from( shp.four_triangles() ),          // 4-triangle STAR (central + all 3 edge-nbrs)
      shp.orbits_from( shp.trios() ),                   // star minus center (3 edge-neighbors)
      shp.orbits_from( shp.site_contours() ),           // five-six contour around each site
      shp.orbits_from( shp.site_contours_minus1() ),    // NEW: four-five (contour minus 1 face)
      shp.orbits_from( shp.site_contours_minus2() ),    // NEW: three-four (minus 2 non-adjacent)
    };
    for(int is=0; is<9; is++) shape_sizes[is] = (int)all[is].size();
    for(int is=0; is<9; is++) for(auto& o : all[is]) orbits.push_back(std::move(o));
  }
  // CONSOLIDATE raw icosahedral orbits into ONE operator per shape (equal weight per orbit) -> 7-shape
  // basis at every L (F_corr_blk = n_lm*n_shapes^2).  9 shapes now (four-five, three-four added).
  const int n_orbits_raw = (int)orbits.size();   // 7 (L1), grows with L
  // const int n_shapes = 7;                     // OLD: single flow level
  const int n_shapes_geom = 9;                   // geometric shape types (7 + four-five + three-four)
  const int n_shapes = N_FLOW * n_shapes_geom;   // 28 = flow level FOLDED INTO the shape axis, so the
                                                 // analysis (op = ishape*n_lm + ilm, n_shapes read
                                                 // from h5) needs NO change; select flow subsets at
                                                 // analysis time via its orbits_arg operator subset.
  const bool SQUARED = true;

  // F^2: physical 0++ is l=0; keep l=1 too (vector F^2 channel). l=2 DROPPED for disk.
  // The exact (l,m) list is written to h5 (ell/em datasets).
  const std::vector<std::array<int,2>> lm_set = {
    {0,0},
    {1,-1},{1,0},{1,1},
  };
  const int n_lm = lm_set.size();
  const int nops = n_shapes * n_lm;        // consolidated 5-shape ops: op = ishape*n_lm + ilm
  const int nops_raw = n_orbits_raw * n_lm; // raw per-orbit ops (computed, then summed by shape)
  std::cout << "# n_shapes = " << n_shapes << " (n_flow = " << N_FLOW << " x n_shapes_geom = "
            << n_shapes_geom << ") n_lm = " << n_lm << " nops = " << nops << std::endl;

  // serial over configs k; parallelism is ensemble-level (one process per Nf).
  for(int k=kmin; k<=k_tmp; k+=stride ){
    // resume-safe: skip configs already fully measured (h5 with "complete" flag)
    // const std::string h5path = dir4+"glue_f2_shapes."+std::to_string(k)+".h5"; // distinct prefix in shared dir
#ifdef FLOW_FULL
    const std::string h5path = dir4+"glue_f2_shapes_fullflow."+std::to_string(k)+".h5"; // full 3D flow variant
#else
    // "_s9" = 9-SHAPE basis (7 + four-five + three-four): DISTINCT prefix so it never clobbers (nor is
    // skipped by the per-config "complete" gate of) the existing 7-shape glue_f2_shapes production data.
    const std::string h5path = dir4+"glue_f2_shapes_s9."+std::to_string(k)+".h5";
#endif
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

    // obs[(iflow*n_shapes_geom + is)*n_lm + ilm][t] -- ONE in-place flow trajectory, measured at each
    // cumulative level FLOW_T[iflow].  Per level: obs_raw[op][t] per RAW orbit (op = iorbit*n_lm+ilm),
    // then CONSOLIDATE (equal weight per orbit) into the shape slot.  SQUARED operators (F^2 / 0++).
    std::vector<std::vector<double>> obs( nops, std::vector<double>(Comp::Nt, 0.0) );
    std::vector<std::vector<double>> obs_raw( nops_raw, std::vector<double>(Comp::Nt, 0.0) );
    for(int iflow=0; iflow<N_FLOW; iflow++){
      // advance the SAME field from the previous level to FLOW_T[iflow] (cumulative, in place)
      const double t_prev = (iflow==0) ? 0.0 : FLOW_T[iflow-1];
      const double seg = FLOW_T[iflow] - t_prev;
      const int nstep_seg = (int)(seg/FLOW_DT + 0.5);
      Flow flow_seg(&SW, seg, nstep_seg);
      flow_seg(Uflow);

      for(int iorbit=0; iorbit<n_orbits_raw; iorbit++){
        for(int ilm=0; ilm<n_lm; ilm++){
          const int ell = lm_set[ilm][0];
          const int em  = lm_set[ilm][1];
          const int op  = iorbit*n_lm + ilm;
          for(int t=0; t<Comp::Nt; t++){
            obs_raw[op][t] = shp.op( Uflow, t, orbits[iorbit], ell, em, SQUARED );
          }
        }
      }
      for(int ilm=0; ilm<n_lm; ilm++){
        int o0=0;
        for(int is=0; is<n_shapes_geom; is++){
          const int ishape = iflow*n_shapes_geom + is;
          for(int oo=0; oo<shape_sizes[is]; oo++){
            const int orbit = o0+oo;
            for(int t=0; t<Comp::Nt; t++) obs[ishape*n_lm+ilm][t] += obs_raw[orbit*n_lm+ilm][t];
          }
          o0 += shape_sizes[is];
        }
      }
    }

    // correlator matrix C(dt)[i][j] = (1/Nt) sum_t obs[i][t] obs[j][t+dt], and one-point <O_i>.
    // Only the small-separation window dt = 0..TMAX_CORR is stored: the GEVP uses dt up to tcut(~5),
    // and the backward fold slice is recovered LOSSLESSLY from the transpose of a stored forward
    // slice via the periodicity identity  C_ij(Nt-d) = C_ji(d)  (see glue_gevp_analysis_claude.cu
    // fold). This cuts both the correlator cost (Nt^2 -> Nt*TMAX) and the h5 size (~Nt/TMAX).
    constexpr int TMAX_CORR = 16;
    const int nsep = std::min(TMAX_CORR + 1, Comp::Nt); // stored separations dt = 0..nsep-1
    // Store ONLY the symmetry-allowed per-(l,m) blocks (see glue2_msm_shapes_claude.cu for the layout):
    // blk[dt][ilm*nob*nob + a*nob + b] = C(a*n_lm+ilm, b*n_lm+ilm). Cross-(l,m) vanish by symmetry.
    const int nob = n_shapes;
    std::vector<std::vector<double>> Fcorr( nsep, std::vector<double>(n_lm*nob*nob, 0.0) );
    for(int dt=0; dt<nsep; dt++){
      for(int t=0; t<Comp::Nt; t++) {
        const int tj = (t+dt)%Comp::Nt;
        for(int ilm=0; ilm<n_lm; ilm++){
          for(int a=0; a<nob; a++){
            const double oa = obs[a*n_lm+ilm][t];
            for(int b=0; b<nob; b++){
              Fcorr[dt][ilm*nob*nob + a*nob + b] += oa * obs[b*n_lm+ilm][tj];
            }
          }
        }
      }
      for(double& x : Fcorr[dt]) x /= Comp::Nt; // F_corr and F share the <..> scale for vacuum subtraction
    }

    std::vector<double> F1( nops, 0.0 );
    for(int t=0; t<Comp::Nt; t++) {
      for(int i=0; i<nops; i++) F1[i] += obs[i][t];
    }
    for(int i=0; i<nops; i++) F1[i] /= Comp::Nt;

    // write per-config HDF5 (F_corr_blk: nsep x n_lm*n_shapes^2 block-diagonal, F: nops); "complete" LAST
    HighFive::File h5( h5path, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate );
    h5.createDataSet( "F_corr_blk", Fcorr );  // per-(l,m) shape blocks (symmetry-allowed only)
    h5.createDataSet( "F", F1 );
    h5.createDataSet( "n_lm", std::vector<int>{n_lm} );  // (l,m) count so the analysis auto-adapts to the l-tower
    h5.createDataSet( "n_shapes", std::vector<int>{n_shapes} );  // shapes/orbit count -> nops = n_shapes*n_lm
    // multi-flow provenance: n_shapes = n_flow * n_shapes_geom, ishape = iflow*n_shapes_geom + is
    h5.createDataSet( "n_flow", std::vector<int>{N_FLOW} );
    h5.createDataSet( "n_shapes_geom", std::vector<int>{n_shapes_geom} );
    h5.createDataSet( "flow_times", std::vector<double>( FLOW_T, FLOW_T + N_FLOW ) );
    // explicit channel labels: ell[ilm], em[ilm] (op = iorbit*n_lm + ilm) so the analysis reads the
    // EXACT (l,m) list -- required now that the driver saves a reduced l-tower (F^2 saves l=0,1).
    std::vector<int> ell_v( n_lm ), em_v( n_lm );
    for(int i=0; i<n_lm; i++){ ell_v[i] = lm_set[i][0]; em_v[i] = lm_set[i][1]; }
    h5.createDataSet( "ell", ell_v );
    h5.createDataSet( "em", em_v );
    h5.createDataSet( "complete", std::vector<int>{1} );
  } // end for k

  return 0;

}
