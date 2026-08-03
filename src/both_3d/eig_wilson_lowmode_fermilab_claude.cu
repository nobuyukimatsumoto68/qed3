// eig_wilson_lowmode_claude.cu
// _claude: WILSON LOW-MODE distribution scan (for the Hasenbusch ladder landscape,
// hasenbusch_massless_impl_plan_claude.md). For each thermalized config of a massless gsq ensemble it
// records lambda_min = smallest SINGULAR VALUE of D_W = sqrt( min eigenvalue of D_W^dag D_W ), the
// near-zero scale that drives the overlap force stiffness, plus lambda_max and the Zolotarev ratio
// lambda_min/lambda_max.
//
// NOTE (per NM): the 2-component Wilson operator in this setup has NO Hermitian H_W = gamma5 D_W (no
// gamma5/chirality in 3D), so this is a SINGULAR value of D_W (from D_W^dag D_W), not an eigenvalue of a
// Hermitian Wilson kernel. lambda_min/lambda_max are computed by OverlapWMass::compute_lambda_max (power
// iteration on D_W^dag D_W for the top, inverse power iteration for the bottom; overlap_wmass_claude.h)
// -- the SAME estimate production uses for the Zolotarev window, so the distribution is the one the
// sign-function / force actually sees.
//
// Usage: ./bin [gsq=8] [Nf=2] [stride=5] [kmax=0(all; <0 = FREE-ONLY, writes only the free .dat)]
//   (massless dir Nf<Nf>_gsq<gsq>...nt128L<LREF>/). kmax<0 records ONLY the free-field reference to
//   wilson_lowmode_free_L<LREF>_claude.dat and does NOT touch the config-scan .dat.
// Output: wilson_lowmode_Nf<Nf>_gsq<gsq>_L<LREF>_claude.dat  (one row per config, appended in scan order;
//   cols: k lambda_min lambda_max ratio -- the config id k is in the row, so the lowest-ev config is read
//   straight off the list). The header also records the FREE-field (U=0) lambda_min/lambda_max (a per-L
//   constant, computed before the loop). Build: tmp_eig_wilson_lowmode_local_claude.sh  (compiles L=1,2,4).

#include <typeinfo>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <filesystem>
#include <cstdint>
#include <complex>
#include <array>
#include <vector>
#include <map>
#include <memory>
#include <string>
#include <cmath>
#include <Eigen/Dense>

using Double = double;
using Idx = std::int32_t;
using Complex = std::complex<double>;
using Face = std::vector<Idx>;
using MS=Eigen::Matrix2cd;
using VD=Eigen::Vector2d;
using VE=Eigen::Vector3d;
using VC=Eigen::VectorXcd;

static constexpr int NS = 2;
static constexpr int DIM = 2;
static constexpr Complex I = Complex(0.0, 1.0);

#define InfoDelta

namespace Comp{
  constexpr bool is_compact=false;
  constexpr int NPARALLEL_DUPDATE=1;
  constexpr int NPARALLEL=NPARALLEL_DUPDATE;
  constexpr int NSTREAMS=NPARALLEL_DUPDATE;
  constexpr int NPARALLEL_GAUGE=16;
  constexpr int NPARALLEL_SORT=16;

#ifndef LREF
#define LREF 1
#endif
  constexpr int N_REFINE=LREF;
  constexpr int NS=2;
  constexpr int Nt=128;         // all massless gsq8 config dirs are nt128 (L1/L2/L4)

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;
  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-8;
}

// _fermilab: ABSOLUTE geometry paths (FNAL), mimicking hmc_hasenbusch_block_fermilab_claude.cu. The run CWD
// on FNAL is the /lustre2 output dir, so relative "../../geometry/..." breaks -- hardcode the absolute tree.
const std::string dir = "/project/affine/nmatsum/qed3/geometry/data/";
#include "/project/affine/nmatsum/qed3/geometry/geodesic.h"

#include "timer.h"
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

#include "sparse_matrix_claude.h"
#include "dirac_simp.h"
#include "dirac_ext.h"
#include "sparse_dirac_claude.h"
#include "matpoly_claude.h"
#define GRAD_L4
#include "includes/overlap_wmass_claude.h"

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(8);

  double gsq   = 8.0;
  int    Nf    = 2;
  int    stride= 5;
  int    kmax  = 0;   // 0 = all available
  if(argc>1) gsq    = atof(argv[1]);
  if(argc>2) Nf     = atoi(argv[2]);
  if(argc>3) stride = atoi(argv[3]);
  if(argc>4) kmax   = atoi(argv[4]);
  if(stride<1) stride = 1;

  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;
  using Fermion=OverlapWMass<WilsonDirac>;

  Base base(Comp::N_REFINE);
  const double M5 = -1.0, at = 0.2, nu0 = 1.0;
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);

  Gauge U(base);
  Fermion D(DW, Complex(0.0,0.0), 21, 0.001);   // massless overlap; only lambda_min/max are used here

  // ----- FREE-field reference (U = 0 = identity links): gsq/Nf-independent, per-L constant -----
  // U was just constructed = zero gauge = free field. Computed FIRST (before the scan .dat is opened) and
  // written to its OWN file, so a free-only invocation (kmax<0) never touches the scan .dat.
  D.unfreeze_lambda();
  D.update( U );
  const double lmin_free = D.lambda_min;
  const double lmax_free = D.lambda_max;
  std::cout << "# FREE (U=0): lambda_min=" << lmin_free << " lambda_max=" << lmax_free
            << " ratio=" << (lmin_free/lmax_free) << std::endl;
  const std::string freename = "wilson_lowmode_free_L"+std::to_string(Comp::N_REFINE)+"_claude.dat";
  std::ofstream fout( freename );
  fout << std::scientific << std::setprecision(10);
  fout << "# FREE-field (U=0) Wilson reference  L=" << Comp::N_REFINE << "  M5=" << M5
       << " at=" << at << " nu0=" << nu0 << std::endl;
  fout << "# lambda_min            lambda_max            ratio" << std::endl;
  fout << lmin_free << "   " << lmax_free << "   " << (lmin_free/lmax_free) << std::endl;
  fout.close();

  if(kmax < 0){   // FREE-ONLY: done, scan .dat untouched
    std::cout << "# free-only (kmax<0): wrote " << freename << "; scan .dat NOT touched" << std::endl;
    return 0;
  }

  // ----- massless config dir: two naming conventions in the tree --
  //   bare        Nf<Nf>_gsq<gsq>at<at>nu0<nu0>nt128L<L>/                        (gsq8 all L, gsq12 L1)
  //   explicit-0  Nf<Nf>_gsq<gsq>at<at>nu0<nu0>mRe0.000000mIm0.000000nt128L<L>/  (gsq12/16 L2/L4)
  const std::string head = "Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)
                         + "nu0"+std::to_string(nu0);
  const std::string tail = "nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
  std::string cfgdir = head + tail;
  if(!std::filesystem::exists(cfgdir)) cfgdir = head + "mRe0.000000mIm0.000000" + tail;
  // _fermilab: the redo dir3 appends a Hasenbusch "_hb<tag>" suffix (e.g. _hb1.000000 for L1/L2,
  // _hb0.400000-1.000000 for the 3-stage L3/L4) -- the bare/explicit-0 names above lack it. Glob for the
  // matching L<L>_hb* dir (massless: mRe0.000000mIm0.000000). Mirrors hmc_hasenbusch_block_fermilab dir3.
  if(!std::filesystem::exists(cfgdir)){
    const std::string pref = head + "mRe0.000000mIm0.000000nt"+std::to_string(Comp::Nt)
                           + "L"+std::to_string(Comp::N_REFINE)+"_hb";
    for(const auto& e : std::filesystem::directory_iterator(".")){
      if(!e.is_directory()) continue;
      const std::string nm = e.path().filename().string();
      if(nm.rfind("Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at",0)==0
         && nm.find("L"+std::to_string(Comp::N_REFINE)+"_hb")!=std::string::npos){ cfgdir = nm + "/"; break; }
    }
  }
  if(!std::filesystem::exists(cfgdir)){
    std::cout << "# ERROR: config dir not found (tried bare, mRe0.000000mIm0.000000, and L"
              << Comp::N_REFINE << "_hb* glob) for Nf" << Nf << " gsq" << gsq << std::endl;
    return 1;
  }

  const std::string outname = "wilson_lowmode_Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)
                            + "_L"+std::to_string(Comp::N_REFINE)+"_claude.dat";
  std::ofstream out( outname );
  out << std::scientific << std::setprecision(10);
  out << "# Wilson low-mode scan  dir=" << cfgdir << "  M5=" << M5 << " at=" << at << " nu0=" << nu0 << std::endl;
  out << "# lambda_min = smallest singular value of D_W = sqrt(min ev of D_W^dag D_W); ratio = lambda_min/lambda_max" << std::endl;
  out << "# FREE (U=0): lambda_min=" << lmin_free << "  lambda_max=" << lmax_free
      << "  ratio=" << (lmin_free/lmax_free) << std::endl;
  out << "# k   lambda_min            lambda_max            ratio" << std::endl;

  std::cout << "# ================ Wilson low-mode scan  L=" << Comp::N_REFINE
            << " gsq=" << gsq << " Nf=" << Nf << " (stride=" << stride << ") ================" << std::endl;
  std::cout << "# dir: " << cfgdir << "  -> " << outname << std::endl;

  int n_done = 0;
  double sum=0.0, sum2=0.0, mn=1.0e30, mx=0.0;
  Timer timer;

  for(int k=stride; ; k+=stride){
    if(kmax>0 && k>kmax) break;
    const std::string str_lat = cfgdir+"ckpoint_lat."+std::to_string(k);
    if(!std::filesystem::exists(str_lat)){
      // stop once we run past the last contiguous checkpoint (allow small gaps within one stride window)
      bool any_ahead = false;
      for(int j=k+1; j<=k+stride; j++){
        if(std::filesystem::exists(cfgdir+"ckpoint_lat."+std::to_string(j))){ any_ahead=true; break; }
      }
      if(!any_ahead) break;
      continue;
    }

    U.read( str_lat );
    D.unfreeze_lambda();
    D.update( U );                     // computes lambda_min, lambda_max

    const double lmin  = D.lambda_min;
    const double lmax  = D.lambda_max;
    const double ratio = lmin / lmax;
    out << k << "   " << lmin << "   " << lmax << "   " << ratio << std::endl;

    sum  += lmin;
    sum2 += lmin*lmin;
    mn = std::min(mn, lmin);
    mx = std::max(mx, lmin);
    n_done++;
    if(n_done % 50 == 0){
      out.flush();
      std::cout << "#   ..." << n_done << " configs (last k=" << k << ", lambda_min=" << lmin
                << ", " << timer.currentSeconds() << " s)" << std::endl;
    }
  }

  out.flush();

  const double mean = (n_done>0) ? sum/n_done : 0.0;
  const double var  = (n_done>0) ? std::max(0.0, sum2/n_done - mean*mean) : 0.0;
  std::cout << "# done: " << n_done << " configs in " << timer.currentSeconds() << " s" << std::endl;
  std::cout << "# lambda_min:  mean=" << mean << "  sd=" << std::sqrt(var)
            << "  min=" << mn << "  max=" << mx << std::endl;
  std::cout << "# wrote " << outname << std::endl;
  return 0;
}
