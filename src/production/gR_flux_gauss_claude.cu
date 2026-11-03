// gR_flux_gauss_claude.cu
//
// Exact free (Gaussian, non-compact) propagator of the magnetic flux harmonic
//   b_{lm}(t) = \sum_f \epsilon_f \theta_f(t) Y_{lm}(n_f),   \epsilon_f = face_signs[f],
// i.e. <b_{lm}(t) b_{lm}(0)>_Gauss for all t, l=1,2, all m, at gsq=1 (result is exactly \propto g^2).
// Denominator of the renormalized coupling g^2(R)R (final/analysis_coupling/gR_flux_impl_plan_claude.md).
//
// Adapted from N. Matsumoto, gauge_prop.cu (github nobuyukimatsumoto68/qed3, release/main;
// local copy src/both_3d/saved_scripts/gauge_prop.cu): Hessian M of the Gaussian action via
// U1WilsonExt::coo_format, source projected onto range(M) (apply M, then CG), then CG for M^+ s.
// Original single-face (face 0, t-correlator) block kept commented below for A/B.
//
// Usage: ./a.out [at]      (Nt = 128, L via -DN_REFINE_CLI)

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
  constexpr int NPARALLEL=1; // 12
  constexpr int NSTREAMS=12; // for grad loop
  constexpr int NPARALLEL_GAUGE=4; // <= 4 cores
  constexpr int NPARALLEL_SORT=4;

  // constexpr int N_REFINE=8;
#ifdef N_REFINE_CLI
  constexpr int N_REFINE=N_REFINE_CLI;
#else
  constexpr int N_REFINE=1;
#endif

  constexpr int NS=2;

  // constexpr int Nt=120;
  constexpr int Nt=128;   // production ensembles

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE; // 30, 120, 480

  constexpr Idx Nx=NS*N_SITES; // matrix size of DW
  constexpr Idx N=Nx*Nt; // matrix size of DW

  const double TOL_INNER=1.0e-15;
  const double TOL_OUTER=1.0e-14;
}


// const std::string dir = "../geometry/data/";
const std::string dir = "../../geometry/data/";   // production tree

#include "timer.h"

#include "../../geometry/geodesic.h"

#include "s2n_simp.h"
#include "rng.h"
// #include "valence.h"     // needs CD / CUDA_CHECK from gpu_header.h -> moved below
// #include "gauge_ext.h"
// #include "action_ext.h"

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
#include "dirac_ext.h"

#include "sparse_dirac.h"
#include "matpoly.h"

#include "dirac_pf.h"
#include "overlap.h"




int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  double at = 0.2;
  if(argc>1) at = atof(argv[1]);

  int device;
  CUDA_CHECK(cudaGetDeviceCount(&device));
  cudaDeviceProp device_prop[device];
  cudaGetDeviceProperties(&device_prop[0], 0);
  std::cout << "# dev = " << device_prop[0].name << std::endl;
  CUDA_CHECK(cudaSetDevice(0));// "TITAN V"
  std::cout << "# (GPU device is set.)" << std::endl;

  // ---------------------------------------

  constexpr int Nt = Comp::Nt;

  using Base=S2Simp;

  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  using Action=U1WilsonExt<Base>;

  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;

  // ----------------------
  // const double T = Comp::T;
  // const double at = T/Comp::Nt;
  std::cout << "# L = " << Comp::N_REFINE << " Nt = " << Nt << " at = " << at << std::endl;

  // const double gsq = 0.05;
  const double gsq = 1.0;   // result is exactly \propto gsq
  Action SW( gsq, at, base );
  std::cout << "# alat = " << base.mean_ell << std::endl;

  // R = 1 check: total face area vs 4\pi
  double area_tot = 0.0;
  for(Idx i=0; i<base.n_faces; i++) area_tot += base.vols[i];
  std::cout << "# sum_f vols = " << area_tot << " 4pi = " << 4.0*M_PI << std::endl;

  Gauge U(base);

  std::vector<COOEntry> coo;
  {
    std::vector<Idx> is;
    std::vector<Idx> js;
    std::vector<Complex> vs;
    SW.coo_format(is, js, vs, U);
    const Idx len = vs.size();
    for(Idx i=0; i<len; i++) coo.push_back( COOEntry( vs[i], is[i], js[i] ) );
  }

  constexpr Idx Ng = Comp::Nt * (Comp::N_LINKS + Comp::N_SITES);
  COO<Ng> op;
  op.en = coo;
  op.do_it();
  MatPoly mat;
  mat.push_back ( cplx(1.0), {&op} );

  // ---- ORIGINAL (gauge_prop.cu): plaquette-plaquette t-correlator of face 0 ----
  // const Idx iface=0;
  // const Face& face = base.faces[iface];
  // std::vector<double> res(Comp::Nt, 0.0);
  // for(Idx i=0; i<face.size(); i++) {
  //   const Idx ix = face[i];
  //   const Idx iy = face[(i+1)%face.size()];
  //   const Link li{ix,iy};
  //   GaugeVector sink(U), source(U);
  //   source.set_zero();
  //   sink.set_zero();
  //   source.sp(0, li) = 1.0 * base.map2sign.at(li);
  //   mat.from_cpu<Ng>( sink.field, source.field );
  //   source.set_zero();
  //   mat.solve<Ng>( source.field, sink.field, 1.0e-11 );
  //   sink.set_zero();
  //   mat.solve<Ng>( sink.field, source.field, 1.0e-11 );
  //   for(Idx j=0; j<face.size(); j++) {
  //     const Idx jx = face[j];
  //     const Idx jy = face[(j+1)%face.size()];
  //     const Link lj{jx,jy};
  //     for(int t=0; t<Comp::Nt; t++) res[t] += base.map2sign.at(lj) * sink.sp(t, lj).real();
  //   }
  // }
  // const double factor = base.vols[iface];
  // for(int t=0; t<Comp::Nt; t++) res[t] /= factor*factor;

  // ---- NEW: flux-harmonic source b_{lm} at t=0, sink b_{lm}(t) ----
  const std::vector<std::array<int,2>> lm_set = {
    {1,-1},{1,0},{1,1},
    {2,-2},{2,-1},{2,0},{2,1},{2,2},
  };
  const int n_lm = lm_set.size();

  std::vector<std::vector<double>> corr( n_lm, std::vector<double>(Nt, 0.0) );

  for(int ilm=0; ilm<n_lm; ilm++){
    const int ell = lm_set[ilm][0];
    const int em = lm_set[ilm][1];
    std::cout << "# (l,m) = (" << ell << "," << em << ")" << std::endl;

    // face weights c_f = \epsilon_f Y_{lm}(n_f)
    std::vector<double> cf(base.n_faces);
    for(Idx f=0; f<base.n_faces; f++) cf[f] = base.face_signs[f] * Ylm_real( ell, em, base.dual_sites[f] );

    // source s = \sum_f c_f \partial\theta_f(0)/\partial A
    GaugeVector sink(U), source(U);
    source.set_zero();
    sink.set_zero();
    for(Idx f=0; f<base.n_faces; f++){
      const Face& face = base.faces[f];
      for(Idx i=0; i<face.size(); i++) {
        const Idx ix = face[i];
        const Idx iy = face[(i+1)%face.size()];
        const Link li{ix,iy};
        source.sp(0, li) += cf[f] * base.map2sign.at(li);
      }
    }

    // project onto range(M), as in gauge_prop.cu
    mat.from_cpu<Ng>( sink.field.data(), source.field.data() );
    source.set_zero();
    mat.solve<Ng>( source.field.data(), sink.field.data(), 1.0e-11 );

    // invert: sink = M^+ s
    sink.set_zero();
    mat.solve<Ng>( sink.field.data(), source.field.data(), 1.0e-11 );

    // b_{lm}(t) of the response = <b_{lm}(t) b_{lm}(0)>_Gauss
    for(Idx f=0; f<base.n_faces; f++){
      const Face& face = base.faces[f];
      for(Idx j=0; j<face.size(); j++) {
        const Idx jx = face[j];
        const Idx jy = face[(j+1)%face.size()];
        const Link lj{jx,jy};
        for(int t=0; t<Nt; t++) corr[ilm][t] += cf[f] * base.map2sign.at(lj) * sink.sp(t, lj).real();
      }
    }

    const double cont = 0.5*std::sqrt( ell*(ell+1.0) );
    std::cout << "# equal-time <b^2> = " << corr[ilm][0] << "  continuum sqrt(l(l+1))/2 = " << cont
              << "  ratio = " << corr[ilm][0]/cont << std::endl;
  }

  std::cout << "# writing" << std::endl;
  {
    const std::string outdir = "final/analysis_coupling/gauss/";
    std::filesystem::create_directories(outdir);
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6) << at;
    const std::string path = outdir+"gR_flux_gauss_L"+std::to_string(Comp::N_REFINE)+"_nt"+std::to_string(Nt)+"_at"+oss.str()+"_claude.dat";
    std::ofstream ofs(path);
    ofs << std::scientific << std::setprecision(15);
    ofs << "# <b_lm(t) b_lm(0)>_Gauss at gsq=1; b_lm = sum_f face_sign Y_lm(n_f) theta_f; cols: ell em t corr" << std::endl;
    ofs << "# L = " << Comp::N_REFINE << " Nt = " << Nt << " at = " << at << " sum_f vols = " << area_tot << std::endl;
    for(int ilm=0; ilm<n_lm; ilm++){
      for(int t=0; t<Nt; t++){
        ofs << lm_set[ilm][0] << " " << lm_set[ilm][1] << " " << t << " " << corr[ilm][t] << std::endl;
      }
    }
    std::cout << "# wrote " << path << std::endl;
  }
  std::cout << "# done" << std::endl;

  return 0;
}
