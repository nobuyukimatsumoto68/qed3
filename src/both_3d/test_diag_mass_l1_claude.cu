// test_diag_mass_l1_claude.cu
// _claude: L=1 equivalence check for the measure-weighted diagonal mass.
// At L=1 the icosahedral symmetry makes m_L = diag(m A_y/abar_s) UNIFORM, so the
// production OverlapWMass (diagonal m_L, physical m) must reduce to the frozen scalar
// obsolete::OverlapWMass with scalar mass c = m * mean_dual_area / mean_ell, applied to
// the same random vector. Tests mult, adj, DHD for m = 0, 0.1 (real), 0.1 i (complex).
// See mass_diag_l1_task_claude.md / mass_measure_factor_impl_plan_claude.md.
//
// Build: add a Makefile.fnal target (mirror the hmc_* targets); run on a GPU.
// Pass: all relative differences ~1e-13..1e-14.

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
#include <random>
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

  constexpr int N_REFINE=1;   // L = 1 (uniform measure -> diagonal reduces to scalar)
  constexpr int NS=2;
  constexpr int Nt=2;         // small; the mass check is t-independent

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;
  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-8;
}

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
#include "includes/overlap_wmass_obsolete_claude.h"   // obsolete::OverlapWMass (scalar mass)
#include "includes/overlap_wmass_claude.h"            // OverlapWMass (diagonal m_L)

// relative L2 difference ||a-b|| / max(||b||, tiny)
static double reldiff(const std::vector<Complex>& a, const std::vector<Complex>& b){
  double nd=0.0, nb=0.0;
  for(size_t i=0;i<a.size();i++){ nd += std::norm(a[i]-b[i]); nb += std::norm(b[i]); }
  return std::sqrt(nd) / std::max(std::sqrt(nb), 1.0e-300);
}

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);

  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();

  constexpr Idx N = Comp::N;
  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;

  Base base(Comp::N_REFINE);
  const double M5 = -1.0, at = 0.2, nu0 = 1.0;
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);
  Gauge U(base);                      // cold/identity gauge (mass reduction is gauge-independent)
  const int npole = 17;

  std::cout << "# L=1 check: N_SITES=" << Comp::N_SITES << " Nt=" << Comp::Nt << " N=" << N
            << " mean_dual_area=" << base.mean_dual_area << " mean_ell=" << base.mean_ell << std::endl;

  // fixed random source
  std::vector<Complex> xi(N);
  std::mt19937_64 gen(12345ull);
  std::normal_distribution<double> nd(0.0,1.0);
  for(Idx i=0;i<N;i++) xi[i] = Complex(nd(gen), nd(gen));

  // apply a device operator (void(CuC*,const CuC*)) to xi, return host result
  auto applyDev = [&](auto applyFn){
    CuC *d_in=nullptr, *d_out=nullptr;
    CUDA_CHECK(cudaMalloc(&d_in,  N*sizeof(CuC)));
    CUDA_CHECK(cudaMalloc(&d_out, N*sizeof(CuC)));
    CUDA_CHECK(cudaMemcpy(d_in, reinterpret_cast<const CuC*>(xi.data()), N*sizeof(CuC), cudaMemcpyHostToDevice));
    applyFn(d_out, d_in);
    CUDA_CHECK(cudaDeviceSynchronize());
    std::vector<Complex> out(N);
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(out.data()), d_out, N*sizeof(CuC), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaFree(d_in)); CUDA_CHECK(cudaFree(d_out));
    return out;
  };

  const std::array<Complex,3> masses = { Complex(0.0,0.0), Complex(0.1,0.0), Complex(0.0,0.1) };

  bool all_ok = true;
  for(const Complex m : masses){
    const Complex c = m * (base.mean_dual_area / base.mean_ell);   // uniform L=1 m_L value
    std::cout << "\n# === m = " << m << "  (scalar c = m*Abar/abar_s = " << c << ") ===" << std::endl;

    obsolete::OverlapWMass<WilsonDirac> Dold(DW, c, npole);   Dold.update(U);   // scalar reference
    OverlapWMass<WilsonDirac>           Dnew(DW, m, npole);   Dnew.update(U);   // diagonal production

    auto o_mult = applyDev([&](CuC* o,const CuC* i){ Dold.mult_deviceAsyncLaunch(o,i); });
    auto n_mult = applyDev([&](CuC* o,const CuC* i){ Dnew.mult_deviceAsyncLaunch(o,i); });
    auto o_adj  = applyDev([&](CuC* o,const CuC* i){ Dold.adj_deviceAsyncLaunch (o,i); });
    auto n_adj  = applyDev([&](CuC* o,const CuC* i){ Dnew.adj_deviceAsyncLaunch (o,i); });
    auto o_dhd  = applyDev([&](CuC* o,const CuC* i){ Dold.DHD_deviceAsyncLaunch (o,i); });
    auto n_dhd  = applyDev([&](CuC* o,const CuC* i){ Dnew.DHD_deviceAsyncLaunch (o,i); });

    const double r_mult = reldiff(n_mult, o_mult);
    const double r_adj  = reldiff(n_adj,  o_adj);
    const double r_dhd  = reldiff(n_dhd,  o_dhd);
    std::cout << "#   reldiff  mult = " << r_mult
              << "   adj = " << r_adj
              << "   DHD = " << r_dhd << std::endl;
    const double tol = 1.0e-10;
    const bool ok = (r_mult<tol && r_adj<tol && r_dhd<tol);
    std::cout << "#   -> " << (ok ? "PASS" : "FAIL") << std::endl;
    all_ok = all_ok && ok;
  }

  std::cout << "\n# ===== " << (all_ok ? "ALL PASS" : "SOME FAIL") << " =====" << std::endl;
  return all_ok ? 0 : 1;
}
