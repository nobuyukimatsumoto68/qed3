// gR_flux_measure_claude.cu
//
// Bare (UNFLOWED) magnetic flux harmonic on every config of one ensemble:
//   b_{lm}(t) = \sum_f \epsilon_f \theta_f(t) Y_{lm}(n_f),   \epsilon_f = face_signs[f],
// l=1,2, all m, all t. Numerator of the renormalized coupling g^2(R)R
// (final/analysis_coupling/gR_flux_impl_plan_claude.md, chunk 2). Denominator = gR_flux_gauss_claude.cu.
// Method: finite-volume coupling, cf. Luscher-Weisz-Wolff Nucl. Phys. B359 (1991) 221.
// Gauge-only: no flow, no solver, CPU only.
//
// Usage: ./a.out <ens_dir> <out_h5>      (L via -DN_REFINE_CLI, Nt = 128)
// Reads ckpoint_lat.1, .2, ... until the first missing index (configs are contiguous from 1).

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
#include <filesystem>
#include <sstream>
#include <Eigen/Dense>
#include <highfive/H5File.hpp>

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


namespace Comp{
  constexpr bool is_compact=false;

  constexpr int NPARALLEL_DUPDATE=12;
  constexpr int NPARALLEL=1;
  constexpr int NSTREAMS=12;
  constexpr int NPARALLEL_GAUGE=4; // <= 4 cores
  constexpr int NPARALLEL_SORT=4;

#ifdef N_REFINE_CLI
  constexpr int N_REFINE=N_REFINE_CLI;
#else
  constexpr int N_REFINE=1;
#endif

  constexpr int NS=2;
  constexpr int Nt=128;   // production ensembles

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE; // 30, 120, 480

  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-15;
  const double TOL_OUTER=1.0e-14;
}


const std::string dir = "../../geometry/data/";   // production tree

#include "timer.h"

#include "../../geometry/geodesic.h"

#include "s2n_simp.h"
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

#include "sparse_matrix.h"
#include "dirac_base.h"
#include "dirac_simp.h"
#include "dirac_ext.h"

#include "sparse_dirac.h"
#include "matpoly.h"

#include "dirac_pf.h"
#include "overlap.h"




int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  if(argc<3){
    std::cout << "Usage: ./a.out <ens_dir> <out_h5>" << std::endl;
    return 1;
  }
  std::string ens_dir = argv[1];
  if(ens_dir.back()!='/') ens_dir += "/";
  const std::string out_h5 = argv[2];

  // L of the binary must match the ensemble
  const std::string tagL = "L"+std::to_string(Comp::N_REFINE)+"_";
  if(ens_dir.find(tagL)==std::string::npos){
    std::cout << "# ERROR: ensemble " << ens_dir << " does not contain " << tagL << std::endl;
    return 1;
  }

  // resume at ensemble level
  if(std::filesystem::exists(out_h5)){
    bool done = false;
    try {
      HighFive::File f(out_h5, HighFive::File::ReadOnly);
      done = f.exist("complete");
    } catch(...) {}
    if(done){
      std::cout << "# " << out_h5 << " complete, skip" << std::endl;
      return 0;
    }
  }

  constexpr int Nt = Comp::Nt;

  using Base=S2Simp;
  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;

  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. L = " << Comp::N_REFINE << " Nt = " << Nt << std::endl;
  std::cout << "# ens_dir = " << ens_dir << std::endl;

  const std::vector<std::array<int,2>> lm_set = {
    {1,-1},{1,0},{1,1},
    {2,-2},{2,-1},{2,0},{2,1},{2,2},
  };
  const int n_lm = lm_set.size();

  // face weights c_f = \epsilon_f Y_{lm}(n_f), identical to gR_flux_gauss_claude.cu
  std::vector<std::vector<double>> cf( n_lm, std::vector<double>(base.n_faces) );
  for(int ilm=0; ilm<n_lm; ilm++){
    for(Idx f=0; f<base.n_faces; f++){
      cf[ilm][f] = base.face_signs[f] * Ylm_real( lm_set[ilm][0], lm_set[ilm][1], base.dual_sites[f] );
    }
  }

  Gauge U(base);

  std::vector<int> k_list;
  std::vector<std::vector<std::vector<double>>> b; // [ik][ilm][t]

  for(int k=1; ; k++){
    const std::string str_lat = ens_dir+"ckpoint_lat."+std::to_string(k);
    if(!std::filesystem::exists(str_lat)) break;
    U.read( str_lat );

    std::vector<std::vector<double>> bk( n_lm, std::vector<double>(Nt, 0.0) );
    for(int t=0; t<Nt; t++){
      for(Idx f=0; f<base.n_faces; f++){
        const double theta = U.plaquette_angle( t, base.faces[f] );
        for(int ilm=0; ilm<n_lm; ilm++) bk[ilm][t] += cf[ilm][f] * theta;
      }
    }
    b.push_back( bk );
    k_list.push_back( k );
    if(k%200==0) std::cout << "# k = " << k << std::endl;
  }

  const int nk = k_list.size();
  std::cout << "# n_cfg = " << nk << std::endl;
  if(nk==0){
    std::cout << "# ERROR: no configs in " << ens_dir << std::endl;
    return 1;
  }

  // quick look: <b^2> averaged over cfg, m, t, per l
  for(int ell=1; ell<=2; ell++){
    double sum = 0.0;
    int cnt = 0;
    for(int ik=0; ik<nk; ik++){
      for(int ilm=0; ilm<n_lm; ilm++){
        if(lm_set[ilm][0]!=ell) continue;
        for(int t=0; t<Nt; t++){
          sum += b[ik][ilm][t] * b[ik][ilm][t];
          cnt++;
        }
      }
    }
    std::cout << "# <b^2>(l=" << ell << ") all cfg = " << sum/cnt << std::endl;
  }

  std::vector<int> ell_v;
  std::vector<int> em_v;
  for(int ilm=0; ilm<n_lm; ilm++){
    ell_v.push_back( lm_set[ilm][0] );
    em_v.push_back( lm_set[ilm][1] );
  }

  {
    std::filesystem::path pp(out_h5);
    if(pp.has_parent_path()) std::filesystem::create_directories( pp.parent_path() );
    HighFive::File h5( out_h5, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate );
    h5.createDataSet( "b", b );          // [cfg][lm][t]
    h5.createDataSet( "k", k_list );
    h5.createDataSet( "ell", ell_v );
    h5.createDataSet( "em", em_v );
    h5.createDataSet( "Nt", std::vector<int>{Nt} );
    h5.createDataSet( "L", std::vector<int>{Comp::N_REFINE} );
    h5.createDataSet( "complete", std::vector<int>{1} );
  }
  std::cout << "# wrote " << out_h5 << std::endl;
  std::cout << "# done" << std::endl;

  return 0;
}
