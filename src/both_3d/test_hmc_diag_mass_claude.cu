// test_hmc_diag_mass_claude.cu
// _claude: FULL-HMC test for the measure-weighted diagonal mass m_L = diag(m A_y/abar_s) at L>2,
// real and imaginary m. Exercises the whole HMC stack (momentum refresh, pf heatbath, integrator,
// dH, Metropolis) with the converted operator + production force (GRAD_L4) -- complementary to the
// single-apply operator test test_diag_mass_l1_claude.cu.
// Plan: mass_hmc_test_impl_plan_claude.md.  Setup mirrors hmc_fermilab_wmass_L2_claude.cu.
//
// Build: -DLREF=2 / -DLREF=4 (lattice L). nvcc flags as tmp_claude.sh. Run on a GPU (handoff).
// Tests (per L x m): (1) reversibility, (2) dH ~ tau^2 scaling, (3) trajectory sniff <exp(-dH)>~1.

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
  constexpr int NPARALLEL_DUPDATE=1;   // omp threads / streams = 1 (single-stream test)
  constexpr int NPARALLEL=NPARALLEL_DUPDATE;
  constexpr int NSTREAMS=NPARALLEL_DUPDATE;
  constexpr int NPARALLEL_GAUGE=16;
  constexpr int NPARALLEL_SORT=16;

#ifndef LREF
#define LREF 2
#endif
  constexpr int N_REFINE=LREF;     // L; build -DLREF=2 / -DLREF=4
  constexpr int NS=2;
  constexpr int Nt=4;              // small Nt: m_L is t-uniform, cheap (production=128)

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;
  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-8;
}

const std::string dir = "../../geometry/data/";
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
#define GRAD_L4   // exercise the PRODUCTION force path in the HMC loop
#include "includes/overlap_wmass_claude.h"
#include "pseudofermion_claude.h"
#include "integrator.h"
#include "hmc.h"

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  // CLI: [gsq] [Nf]  (mass is looped internally: real + imaginary)
  double gsq = 8.0;
  if(argc>1) gsq = atof(argv[1]);
  int Nf = 2;
  if(argc>2) Nf = atoi(argv[2]);
  assert(Nf%2==0);

  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();

  using BaseLink = std::array<Idx,2>;
  constexpr Idx N = Comp::N;
  constexpr int Nt = Comp::Nt;
  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Force=GaugeExt<Base,Nt,Comp::is_compact>;
  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  using Action=U1WilsonExt<Base>;
  using Rng=ParallelRngExt<Base,Nt>;
  using Fermion=OverlapWMass<WilsonDirac>;
  using PF=PseudoFermion<Fermion>;

  Base base(Comp::N_REFINE);
  const double M5 = -1.0, at = 0.2, nu0 = 1.0;
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);
  Action SW(gsq, at, base);

  Rng rng(base);
  Gauge U(base);
  U.gaussian(rng, 0.3);            // nontrivial config (OPEN Q: thermalize first via run(no_reject))

  const int npole = (Comp::N_REFINE>=4) ? 13 : 17;   // match the wmass drivers
  const int n_pf  = Nf/2;

  // integrator + trajectory params (OPEN Q): mirror the wmass drivers
  const double tmax = 1.9;
  const int nsteps = (Comp::N_REFINE>=4) ? 12 : 9;
  const int nsteps_inner = 100;

  std::cout << "# L=" << Comp::N_REFINE << " Nt=" << Nt << " N=" << N
            << " gsq=" << gsq << " Nf=" << Nf << " npole=" << npole
            << " mean_dual_area=" << base.mean_dual_area << " mean_ell=" << base.mean_ell << std::endl;

  // OPEN Q (mass set): machinery test -> literal 0.1 real + 0.1 imaginary (+ massless baseline).
  const std::array<Complex,3> masses = { Complex(0.0,0.0), Complex(0.1,0.0), Complex(0.0,0.1) };

  bool all_ok = true;
  Gauge Uinit(U);   // each mass starts from the same gaussian config (tests 1-2 restore it; test 3 evolves)

  for(const Complex mass : masses){
    std::cout << "\n# =================== m = " << mass << " ===================" << std::endl;
    U = Uinit;
    Fermion D(DW, mass, npole);
    D.update(U);
    std::cout << "#   delta = " << D.Delta()
              << "  lambda_min/lambda_max = " << D.lambda_min/D.lambda_max << std::endl;

    std::vector<std::shared_ptr<PF>> pfs;
    for(int f=0; f<n_pf; f++) pfs.push_back( std::make_shared<PF>(D) );

    Force pi(base);
    MinimumNorm2 integrator(tmax, nsteps, nsteps_inner);
    HMC2 hmc(rng, &SW, &D, U, pi, pfs, &integrator);

    // ===== Test 1: reversibility =====
    // refresh pi + heatbath phi; integrate forward; flip pi; integrate back; must return to (U0,-pi0).
    // Limited by the inner-CG residual accumulated over 2*nsteps (~1e-7..1e-9); a wrong diagonal force
    // (e.g. dropped M_mass weight) breaks it grossly.
    {
      pi.gaussian(rng);
      for(auto& pf : pfs) pf->gen(rng);             // fixed phi for the whole round trip
      Gauge U0(U);
      Force pi0(pi);
      hmc.integrate();                              // forward tau
      pi = (-1.0)*pi;                               // flip momenta
      hmc.integrate();                              // backward tau
      double dU=0.0, dpi=0.0;
      for(int t=0;t<Nt;t++){
        for(Idx il=0; il<base.n_links; il++){
          dU  = std::max(dU,  std::abs(U.sp(t,il)  - U0.sp(t,il)));
          dpi = std::max(dpi, std::abs(pi.sp(t,il) + pi0.sp(t,il)));   // pi should be -pi0
        }
        for(Idx ix=0; ix<base.n_sites; ix++){
          dU  = std::max(dU,  std::abs(U.tp(t,ix)  - U0.tp(t,ix)));
          dpi = std::max(dpi, std::abs(pi.tp(t,ix) + pi0.tp(t,ix)));
        }
      }
      const bool ok = (dU < 1.0e-5 && dpi < 1.0e-5);
      std::cout << "#   [test 1 reversibility] max|U-U0|=" << dU << "  max|pi+pi0|=" << dpi
                << "  -> " << (ok?"PASS":"FAIL") << std::endl;
      all_ok = all_ok && ok;
      U = U0; D.update(U);
    }

    // ===== Test 2: dH ~ tau^2 scaling =====
    // same fixed (pi,phi); integrate from U0 with nsteps and 2*nsteps; MN2 is 2nd-order => |dH| ~ tau^2,
    // so doubling nsteps cuts |dH| ~4x. A broken force gives ratio ~1 or |dH| blowing up.
    // NB: the tau^2 law is ASYMPTOTIC (small tau). At the production tmax=1.9 the massive force's
    // higher-order terms spoil the clean 4x at coarse nsteps (the ratio came out ~1.3-2.4 for m!=0
    // -- NOT a force bug: L=1 force-vs-FD ~1e-8 for all m, reversibility ~1e-10). Use a SMALL
    // trajectory length tmax_sc so the leading tau^2 dominates and the ratio is ~4 for every m.
    {
      // FORCE FIXED (grad_diag_mass_force_bug_claude.md Sec 5'): the force is now exactly dS/dU, so dH
      // scales ~tau^2 at the PRODUCTION trajectory length -> ratio ~4 for every m (no more tmax workaround).
      const double tmax_sc = tmax;        // production traj length (was 0.2 diagnostic workaround)
      const int    ns_sc   = nsteps;
      // ONE shared start (U0, pi0, fixed phi); h0 computed ONCE so both integrations use the IDENTICAL
      // start Hamiltonian. Each j fully restores U=U0 (+ D.update), pi=pi0, eta@U0 before integrating.
      pi.gaussian(rng);
      for(auto& pf : pfs) pf->gen(rng);            // heatbath phi ONCE -> fixed for both ns
      Gauge U0(U);
      Force pi0(pi);
      for(auto& pf : pfs) pf->update_eta();         // eta at U0 with fixed phi
      double h0;
      {
        MinimumNorm2 integ0(tmax_sc, ns_sc, nsteps_inner);
        HMC2 h0op(rng, &SW, &D, U, pi, pfs, &integ0);
        h0 = h0op.H();                              // shared start Hamiltonian
      }
      double dHv[2];
      const int ns[2] = { ns_sc, 2*ns_sc };
      for(int j=0;j<2;j++){
        U = U0; D.update(U);                        // restore EXACT same config
        pi = pi0;
        for(auto& pf : pfs) pf->update_eta();        // eta@U0 (same phi) -> H() here == h0
        MinimumNorm2 integ(tmax_sc, ns[j], nsteps_inner);
        HMC2 h(rng, &SW, &D, U, pi, pfs, &integ);
        h.integrate();
        dHv[j] = h.H() - h0;
      }
      const double ratio = std::abs(dHv[0]) / std::max(std::abs(dHv[1]), 1.0e-300);
      const bool ok = (std::isfinite(ratio) && std::abs(dHv[0])<5.0 && ratio>2.5);
      std::cout << "#   [test 2 dH scaling] h0=" << h0
                << "  dH(ns="<<ns[0]<<")="<<dHv[0]
                << "  dH(ns="<<ns[1]<<")="<<dHv[1]
                << "  ratio="<<ratio<<" (expect ~4)  -> "<<(ok?"PASS":"FAIL")<<std::endl;
      all_ok = all_ok && ok;
      U = U0; D.update(U);
    }

    // ===== Test 3: trajectory sniff =====
    // K full hmc.run() trajectories; <exp(-dH)> should be ~1 (unbiased Metropolis). Loose statistical
    // gate (K small); a broken force shows up as runaway dH / 0% acceptance.
    {
      const int K = 8;
      double sum_exp = 0.0;
      int nacc = 0;
      std::cout << "#   [test 3 trajectory] per-traj dH / accept:" << std::endl;
      for(int k=0;k<K;k++){
        double r, dH;
        bool acc;
        hmc.run(r, dH, acc);
        sum_exp += std::exp(-dH);
        if(acc) nacc++;
        std::cout << "#     traj " << k << ": dH=" << dH << "  r=" << r << "  acc=" << acc << std::endl;
      }
      const double mean_exp = sum_exp / K;
      const bool ok = (std::isfinite(mean_exp) && mean_exp>0.5 && mean_exp<2.0 && nacc>0);
      std::cout << "#   [test 3] <exp(-dH)>=" << mean_exp << "  accept=" << nacc << "/" << K
                << "  -> " << (ok?"PASS":"FAIL") << std::endl;
      all_ok = all_ok && ok;
    }

    D.update(U);   // restore D(U) for the next mass
  }

  std::cout << "\n# ===== HMC diagonal-mass test: " << (all_ok ? "ALL PASS" : "SOME FAIL") << " =====" << std::endl;

  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
  return all_ok ? 0 : 1;
}
