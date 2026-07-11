// test_hasenbusch_bilinear_l1_claude.cu
// _claude: L=1 force-vs-FD gate for the Hasenbusch external-bra force (C1 of
// hasenbusch_massless_impl_plan_claude.md). The ratio-frame cross term is
//   Term B = 2 Re <phi | K | eta>,   K = dD_ov/dU  (conserved current, mass-independent),
// and OverlapWMass::precalc_grad_bilinear_deviceAsyncLaunch_ms(U, phi, eta) + the usual
// grad_deviceAsyncLaunch(link, U, eta) must reproduce it. We check it against the CENTRAL
// finite difference of the scalar S_B(U) = 2 Re <phi | D_ov(U) | eta> with phi, eta FIXED
// (NOT re-solved: unlike a pseudo-fermion action there is no stationarity, so both the
// analytic and the FD keep phi, eta frozen). Via the GW relation D_ov^dag D_ov = D_ov +
// D_ov^dag the standard massless PF force is 2 Re <eta|K|eta>, i.e. the phi==eta special
// case of this test -- also checked as a sanity cross-reference.
//
// The KEY thing this gates: grad computes 2 Re <bra|K|ket> for a GENERIC bra (the massive
// path only ever used bra = (1+m_L)eta). If the two pole terms are not a clean z + z* for
// arbitrary bra, this FD FAILS and Term B would need an explicit bra<->ket symmetrization.
//
// Build (default reference grad; NO GRAD_L1/2/4): see tmp_claude.sh. Pass: |grad - fd| ~1e-5.

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

#ifndef LREF
#define LREF 1
#endif
  constexpr int N_REFINE=LREF;   // L; -DLREF=2 / -DLREF=4 to run the FD check at L=2,4 too
  constexpr int NS=2;
  constexpr int Nt=2;         // small; Term B is t-local in character

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;
  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-8;
}

const std::string dir = "../../geometry/data/";                              // local path
#include "../../geometry/geodesic.h"

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
// no GRAD_L1/2/4 -> default reference grad (add -DGRAD_L4 to also exercise the block grad).
#include "includes/overlap_wmass_claude.h"            // OverlapWMass (+ precalc_grad_bilinear_..._ms)
#include "pseudofermion_claude.h"                      // PseudoFermion (phi==eta cross-check)

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);

  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();

  constexpr Idx N = Comp::N;
  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;
  using Force=Gauge;

  Base base(Comp::N_REFINE);
  const double M5 = -1.0, at = 0.2, nu0 = 1.0;
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);

  using Rng = ParallelRngExt<Base, Comp::Nt>;
  Rng rng(base);
  Gauge U(base);
  U.gaussian(rng, 0.3);              // nontrivial, well-conditioned, reproducible
  const int npole = 17;

  std::cout << "# Hasenbusch bilinear-force L=" << Comp::N_REFINE << " FD gate: N=" << N
            << " Nt=" << Comp::Nt << " npole=" << npole << std::endl;

  // two DIFFERENT fixed Gaussian vectors: bra phi, ket eta (the essential bra != ket test)
  std::vector<Complex> phi_h(N), eta_h(N);
  { std::mt19937_64 gp(2718ull), ge(3141ull); std::normal_distribution<double> g01(0.0,1.0);
    for(Idx i=0;i<N;i++){ phi_h[i]=Complex(g01(gp),g01(gp)); eta_h[i]=Complex(g01(ge),g01(ge)); } }

  CuC *d_phi=nullptr, *d_eta=nullptr, *d_res=nullptr;
  CUDA_CHECK(cudaMalloc(&d_phi, N*sizeof(CuC)));
  CUDA_CHECK(cudaMalloc(&d_eta, N*sizeof(CuC)));
  CUDA_CHECK(cudaMalloc(&d_res, N*sizeof(CuC)));
  CUDA_CHECK(cudaMemcpy(d_phi, reinterpret_cast<const CuC*>(phi_h.data()), N*sizeof(CuC), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_eta, reinterpret_cast<const CuC*>(eta_h.data()), N*sizeof(CuC), cudaMemcpyHostToDevice));

  // massless operator: K = dD_ov/dU is mass-independent, and m_L=0 => mult applies bare D_ov.
  OverlapWMass<WilsonDirac> D(DW, Complex(0.0,0.0), npole);
  D.update(U);

  // --- analytic Term B force over ALL links: precalc_grad_bilinear(U, phi, eta) + compute(eta) ---
  D.precalc_grad_bilinear_deviceAsyncLaunch_ms(U, d_phi, d_eta);
  Force fB(base);
  fB.compute(U, d_eta, D);           // loops grad_deviceAsyncLaunch(link, U, eta) using the bilinear precalc

  // --- FD reference: S_B(U') = 2 Re <phi | D_ov(U') eta>, phi & eta FROZEN ---
  auto S_B = [&](const Gauge& Up)->double {
    D.update(Up);
    D.mult_deviceAsyncLaunch(d_res, d_eta);          // D_ov(Up) eta   (m_L=0)
    CUDA_CHECK(cudaDeviceSynchronize());
    std::vector<Complex> res(N);
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(res.data()), d_res, N*sizeof(CuC), cudaMemcpyDeviceToHost));
    Complex acc(0.0,0.0);
    for(Idx i=0;i<N;i++) acc += std::conj(phi_h[i])*res[i];
    return 2.0*acc.real();
  };

  const double eps = 1.0e-5;
  const int s = 0;
  const std::array<Idx,3> test_sp = {2, 7, 13};
  const std::array<Idx,2> test_tp = {3, 8};

  auto fd_sp = [&](Idx il){ Gauge UP(U); UP.sp(s,il)+=eps; double sp=S_B(UP);
                            Gauge UM(U); UM.sp(s,il)-=eps; double sm=S_B(UM); return (sp-sm)/(2.0*eps); };
  auto fd_tp = [&](Idx ix){ Gauge UP(U); UP.tp(s,ix)+=eps; double sp=S_B(UP);
                            Gauge UM(U); UM.tp(s,ix)-=eps; double sm=S_B(UM); return (sp-sm)/(2.0*eps); };

  std::cout << "\n# ===== Term B: analytic grad vs central-FD of 2Re<phi|D_ov|eta> (bra != ket) =====" << std::endl;
  bool ok = true;
  for(Idx il : test_sp){ double an=fB.sp(s,il), fd=fd_sp(il); double d=std::abs(an-fd);
    ok = ok && (d<1.0e-5);
    std::cout << "#   sp " << il << ": grad=" << an << " fd=" << fd << " |d|=" << d << std::endl; }
  for(Idx ix : test_tp){ double an=fB.tp(s,ix), fd=fd_tp(ix); double d=std::abs(an-fd);
    ok = ok && (d<1.0e-5);
    std::cout << "#   tp " << ix << ": grad=" << an << " fd=" << fd << " |d|=" << d << std::endl; }
  D.update(U);   // restore

  // --- sanity cross-check: phi == eta must reproduce the standard massless PF force ---
  // (standard: precalc_grad(eta) + compute -> 2Re<eta|K|eta>; bilinear with bra=ket=eta must match)
  {
    PseudoFermion<OverlapWMass<WilsonDirac>> pf(D);
    CUDA_CHECK(cudaMemcpy(pf.d_phi, d_eta, N*sizeof(CuC), cudaMemcpyDeviceToDevice)); // phi_pf := eta
    pf.update_eta();                                   // eta_pf = (D^dag D)^-1 eta   (a DIFFERENT vector)
    // Use the SAME device vector on both sides of the bilinear -> must equal the standard force on it.
    D.precalc_grad_bilinear_deviceAsyncLaunch_ms(U, pf.d_eta, pf.d_eta);
    Force fdiag_bil(base); fdiag_bil.compute(U, pf.d_eta, D);
    D.precalc_grad_deviceAsyncLaunch(U, pf.d_eta);
    Force fdiag_std(base); fdiag_std.compute(U, pf.d_eta, D);
    double maxd=0.0;
    for(int t=0;t<Comp::Nt;t++){
      for(Idx il=0; il<base.n_links; il++) maxd=std::max(maxd, std::abs(fdiag_bil.sp(t,il)-fdiag_std.sp(t,il)));
      for(Idx ix=0; ix<base.n_sites; ix++) maxd=std::max(maxd, std::abs(fdiag_bil.tp(t,ix)-fdiag_std.tp(t,ix)));
    }
    const bool okD = (maxd < 1.0e-9);
    ok = ok && okD;
    std::cout << "\n# phi==eta cross-check: max|bilinear - standard massless force| = " << maxd
              << "  -> " << (okD ? "PASS" : "FAIL") << std::endl;
    D.update(U);
  }

  CUDA_CHECK(cudaFree(d_phi)); CUDA_CHECK(cudaFree(d_eta)); CUDA_CHECK(cudaFree(d_res));
  std::cout << "\n# ===== Hasenbusch bilinear-force L=1 gate: " << (ok ? "PASS" : "FAIL")
            << " (grad-vs-FD tol ~1e-5; phi==eta tol ~1e-9) =====" << std::endl;
  return ok ? 0 : 1;
}
