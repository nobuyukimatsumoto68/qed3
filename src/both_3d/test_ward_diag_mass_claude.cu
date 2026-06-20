// test_ward_diag_mass_claude.cu
// _claude: OPERATOR-LEVEL Ward / current-conservation check for the measure-weighted diagonal mass
// D_m = D_ov + m_L (m_L = diag(m A_y/abar_s)). Two checks, looped over m in {0, 0.1, 0.1i}:
//   (A) gauge-commutator Ward identity:  [D_m, Theta] xi = sum_l dtheta_l K^l xi
//       (Theta = diag per-site phase; K = conserved current apply_k). Since m_L is diagonal and
//       gauge-singlet, [Theta, m_L]=0, so this must hold for ANY m exactly as for D_ov.
//       -> tests that the mass added to mult/adj is gauge-covariant (a gauge-dependent mass bug breaks it).
//   (B) current conservation:  sum_{z~w} tr(K^{wz} D_m^{-1}) = 0  (vector Ward identity, mass-independent
//       physically: the mass is vector-singlet). Uses the MASSIVE propagator D_m^{-1} -> tests the converted
//       D_m solve. Analytic (basis loop) at small lattice; stochastic Z2xZ2 otherwise.
//
// Ported from saved_scripts_claude/check_conserved_current_claude.cu (Overlap -> OverlapWMass + mass loop).
// Plan/refs: mass_measure_audit_handoff_claude.md (#5), check_conserved_current_impl_plan_claude.md.
// Build: -DLREF=1 (Nt=4, with analytic reference) / -DLREF=2. Runner: tmp_ward_claude.sh (GPU=1).

#include <typeinfo>
#include <cmath>
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
#include <functional>
#include <Eigen/Dense>

using Double  = double;
using Idx     = std::int32_t;
using Complex = std::complex<double>;
using Link    = std::array<Idx,2>;
using Face    = std::vector<Idx>;
using MS = Eigen::Matrix2cd;
using VD = Eigen::Vector2d;
using VE = Eigen::Vector3d;
using VC = Eigen::VectorXcd;

static constexpr int NS  = 2;
static constexpr int DIM = 2;
static constexpr Complex I = Complex(0.0, 1.0);

namespace Comp {
  constexpr bool is_compact = false;
  constexpr int NPARALLEL_DUPDATE = 1;
  constexpr int NPARALLEL         = NPARALLEL_DUPDATE;
  constexpr int NSTREAMS          = NPARALLEL_DUPDATE;
  constexpr int NPARALLEL_GAUGE   = NPARALLEL_DUPDATE;
  constexpr int NPARALLEL_SORT    = NPARALLEL_DUPDATE;

#ifndef LREF
#define LREF 1
#endif
  constexpr int N_REFINE = LREF;
  constexpr int NS       = 2;
  constexpr int Nt       = 4;          // small; analytic basis-loop reference at LREF=1, Nt<=4

  constexpr Idx N_SITES  = 10*N_REFINE*N_REFINE + 2;
  constexpr int N_LINKS  = 30*N_REFINE*N_REFINE;
  constexpr Idx Nx       = NS*N_SITES;
  constexpr Idx N        = Nx*Nt;

  const double TOL_INNER = 1.0e-9;
  const double TOL_OUTER = 1.0e-8;
}

const std::string dir = "../../geometry/data/";
#include "timer.h"
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

#include "valence_claude.h"
#include "gauge_ext.h"
#include "action_ext.h"

#include "sparse_matrix.h"
#include "dirac_base.h"
#include "dirac_simp.h"
#include "dirac_dual.h"
#include "dirac_ext.h"
#include "sparse_dirac_claude.h"
#include "matpoly_claude.h"

#include "overlap_wmass_claude.h"
#include "conserved_current_claude.h"

using Base        = S2Simp;
using WilsonDirac = DiracExt<Base, DiracS2Simp>;
using Gauge       = GaugeExt<Base, Comp::Nt, Comp::is_compact>;
using Fermion     = OverlapWMass<WilsonDirac>;

// (Theta v)[site, sigma] = h_theta[site] * v[site, sigma]; h_theta one real per spacetime site s*N_SITES+ix.
void multTheta(CuC* d_out, const CuC* d_in, const std::vector<double>& h_theta) {
  constexpr Idx N = Comp::N;
  std::vector<CuC> h(N);
  CUDA_CHECK(cudaMemcpy(h.data(), d_in, N*CD, D2H));
  for(int site = 0; site < (int)h_theta.size(); site++) {
    const double t = h_theta[site];
    h[2*site+0] = make_cuDoubleComplex(t*cuCreal(h[2*site+0]), t*cuCimag(h[2*site+0]));
    h[2*site+1] = make_cuDoubleComplex(t*cuCreal(h[2*site+1]), t*cuCimag(h[2*site+1]));
  }
  CUDA_CHECK(cudaMemcpy(d_out, h.data(), N*CD, H2D));
}

int main(int argc, char* argv[]) {
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  double gsq = 8.0;
  if(argc>1) gsq = atof(argv[1]);
  int Nf = 2;
  if(argc>2) Nf = atoi(argv[2]);
  int nhits = 32;
  if(argc>3) nhits = atoi(argv[3]);

  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();

  constexpr Idx N  = Comp::N;
  constexpr int Nt = Comp::Nt;

  Base base(Comp::N_REFINE);
  const double at = 0.2, M5 = -1.0, nu0 = 1.0;
  if(Nt != 1) assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);
  const int npole = 21;

  Gauge U(base);   // free-field (Ward identity holds config-by-config; free is simplest + has analytic ref)
  std::cout << "# Ward check  L=" << Comp::N_REFINE << " Nt=" << Nt << " N=" << N
            << " n_sites=" << base.n_sites << " n_links=" << base.n_links << " npole=" << npole << std::endl;

  // fixed random xi for the commutator check
  std::mt19937 rng_host(42);
  FermionVector h_xi;
  {
    std::normal_distribution<double> dist(0.0,1.0);
    for(int i=0;i<N;i++) h_xi.field[i] = Complex(dist(rng_host), dist(rng_host));
  }
  CuC* d_xi = nullptr;
  CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(h_xi.field), N*CD, H2D));

  // random per-site theta (shared across masses)
  std::vector<double> h_theta(Comp::N_SITES * Nt);
  {
    std::normal_distribution<double> rdist(0.0,1.0);
    for(auto& t : h_theta) t = rdist(rng_host);
  }

  const std::array<Complex,3> masses = { Complex(0.0,0.0), Complex(0.1,0.0), Complex(0.0,0.1) };
  bool all_ok = true;

  for(const Complex mass : masses) {
    std::cout << "\n# =================== m = " << mass << " ===================" << std::endl;
    Fermion D(DW, mass, npole);
    D.update(U);
    ConservedCurrent<Fermion,Gauge> kop(D);   // K is mass-independent (uses D_ov structure)
    MatPoly op_K; op_K.push_back(cplx(1.0), {&kop});

    // ===== (A) [D_m, Theta] xi = sum_l dtheta_l K^l xi =====
    {
      auto f_Dm = std::bind(&Fermion::mult_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);
      LinOpWrapper M_Dm(f_Dm);
      MatPoly op_Dm; op_Dm.push_back(cplx(1.0), {&M_Dm});

      CuC *d_th=nullptr,*d_chi1=nullptr,*d_chi2=nullptr,*d_tmp=nullptr;
      CUDA_CHECK(cudaMalloc(&d_th,N*CD)); CUDA_CHECK(cudaMalloc(&d_chi1,N*CD));
      CUDA_CHECK(cudaMalloc(&d_chi2,N*CD)); CUDA_CHECK(cudaMalloc(&d_tmp,N*CD));

      // chi1 = D_m(Theta xi) - Theta(D_m xi)
      multTheta(d_th, d_xi, h_theta);
      op_Dm.on_gpu<N>(d_chi1, d_th);   CUDA_CHECK(cudaDeviceSynchronize());
      op_Dm.on_gpu<N>(d_tmp, d_xi);    CUDA_CHECK(cudaDeviceSynchronize());
      multTheta(d_tmp, d_tmp, h_theta);
      Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_chi1, -1.0, d_tmp, d_chi1);
      CUDA_CHECK(cudaDeviceSynchronize());

      // chi2 = sum_l dtheta_l K^l xi  (spatial then temporal links)
      CUDA_CHECK(cudaMemset(d_chi2, 0, N*CD));
      for(int s=0;s<Nt;s++) {
        for(const auto& lk : base.links) {
          const double dt = h_theta[s*Comp::N_SITES + lk[0]] - h_theta[s*Comp::N_SITES + lk[1]];
          kop.apply_k(d_tmp, d_xi, U, std::pair<int,BaseLink>{s, lk});
          CUDA_CHECK(cudaDeviceSynchronize());
          Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_chi2, dt, d_tmp, d_chi2);
          CUDA_CHECK(cudaDeviceSynchronize());
        }
      }
      if(Nt > 1) {
        for(int s=0;s<Nt;s++) {
          const int s1 = (s+1)%Nt;
          for(int ix=0; ix<base.n_sites; ix++) {
            const double dt = h_theta[s*Comp::N_SITES + ix] - h_theta[s1*Comp::N_SITES + ix];
            kop.apply_k(d_tmp, d_xi, U, std::pair<int,Idx>{s, ix});
            CUDA_CHECK(cudaDeviceSynchronize());
            Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_chi2, dt, d_tmp, d_chi2);
            CUDA_CHECK(cudaDeviceSynchronize());
          }
        }
      }

      FermionVector h_chi1, h_chi2;
      CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(h_chi1.field), d_chi1, N*CD, D2H));
      CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(h_chi2.field), d_chi2, N*CD, D2H));
      double max_err=0.0, nrm=0.0;
      for(int i=0;i<N;i++){ max_err=std::max(max_err,std::abs(h_chi1.field[i]-h_chi2.field[i])); nrm=std::max(nrm,std::abs(h_chi1.field[i])); }
      const double rel = max_err/std::max(nrm,1.0e-300);
      const bool ok = (rel < 1.0e-6);
      std::cout << "#   (A) [D_m,Theta]xi=sum dtheta K^l xi : ||chi1-chi2||_inf=" << max_err
                << " (rel " << rel << ")  -> " << (ok?"PASS":"FAIL") << std::endl;
      all_ok = all_ok && ok;
      CUDA_CHECK(cudaFree(d_th)); CUDA_CHECK(cudaFree(d_chi1));
      CUDA_CHECK(cudaFree(d_chi2)); CUDA_CHECK(cudaFree(d_tmp));
    }

    // ===== (B) current conservation: sum_{z~w} tr(K^{wz} D_m^{-1}) = 0 at w=(s=0,ix=0) =====
    // J_{V,+}^{wz} = -Nf * (eta^dag K^{wz} phi),  D_m phi = eta (normal eqs D_m^dag D_m phi = D_m^dag eta).
    {
      constexpr int s_focus=0, ix_focus=0;
      auto f_DH  = std::bind(&Fermion::adj_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);
      auto f_DHD = std::bind(&Fermion::DHD_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);
      LinOpWrapper M_DH(f_DH), M_DHD(f_DHD);
      MatPoly op_DH;  op_DH.push_back(cplx(1.0), {&M_DH});
      MatPoly op_DHD; op_DHD.push_back(cplx(1.0), {&M_DHD});

      std::vector<Idx> nns_list(base.nns[ix_focus].begin(), base.nns[ix_focus].end());
      const int n_spatial = (int)nns_list.size();

      // analytic divergence (basis loop), small lattice only -> deterministic current-conservation gate
      if constexpr (Comp::N_REFINE == 1 && Comp::Nt <= 4) {
        Complex div_exact(0.0,0.0);
        FermionVector e_i, DH_e, phi_e, fv_Ke;
        std::vector<Complex> JV(n_spatial, Complex(0.0,0.0));
        Complex JV_t[2] = {Complex(0.0,0.0), Complex(0.0,0.0)};
        for(Idx k=0;k<N;k++) {
          e_i.set_pt_source(k/Comp::Nx, (k%Comp::Nx)/NS, k%NS);
          op_DH.from_cpu<N>(DH_e.field, e_i.field);
          op_DHD.solve<N>(phi_e.field, DH_e.field, Comp::TOL_OUTER);
          for(int j=0;j<n_spatial;j++) {
            kop.set(U, std::pair<int,BaseLink>{s_focus, BaseLink{ix_focus, nns_list[j]}}, false);
            op_K.from_cpu<N>(fv_Ke.field, phi_e.field);
            JV[j] += fv_Ke.field[k];
          }
          const std::pair<int,Idx> tl[2] = { {s_focus, ix_focus}, {(s_focus-1+Nt)%Nt, ix_focus} };
          const double tsign[2] = {1.0, -1.0};
          for(int t=0;t<2;t++) {
            kop.set(U, tl[t], false);
            op_K.from_cpu<N>(fv_Ke.field, phi_e.field);
            JV_t[t] += tsign[t]*fv_Ke.field[k];
          }
        }
        for(int j=0;j<n_spatial;j++) div_exact += -(double)Nf * JV[j];
        for(int t=0;t<2;t++)        div_exact += -(double)Nf * JV_t[t];
        const bool ok = (std::abs(div_exact) < 1.0e-6);
        std::cout << "#   (B) analytic sum_z J_V^{wz} (current conservation) = (" << div_exact.real()
                  << "," << div_exact.imag() << ")  -> " << (ok?"PASS":"FAIL") << std::endl;
        all_ok = all_ok && ok;
      } else {
        // stochastic divergence: should -> 0 within ~1/sqrt(nhits)
        using Rng5 = ParallelRngExt<Base, Comp::Nt>;
        Rng5 rng_st(base, 42);
        FermionVector eta, DH_eta, phi, fv_Kphi;
        Complex div_acc(0.0,0.0);
        std::vector<Complex> div_all;
        for(int h=0;h<nhits;h++) {
          eta.fill_z2_source(rng_st);
          op_DH.from_cpu<N>(DH_eta.field, eta.field);
          op_DHD.solve<N>(phi.field, DH_eta.field, Comp::TOL_OUTER);
          Complex div_h(0.0,0.0);
          for(int j=0;j<n_spatial;j++) {
            kop.set(U, std::pair<int,BaseLink>{s_focus, BaseLink{ix_focus, nns_list[j]}}, false);
            op_K.from_cpu<N>(fv_Kphi.field, phi.field);
            div_h += eta.dag(fv_Kphi);
          }
          const std::pair<int,Idx> tl[2] = { {s_focus, ix_focus}, {(s_focus-1+Nt)%Nt, ix_focus} };
          const double tsign[2] = {1.0,-1.0};
          for(int t=0;t<2;t++) {
            kop.set(U, tl[t], false);
            op_K.from_cpu<N>(fv_Kphi.field, phi.field);
            div_h += tsign[t]*eta.dag(fv_Kphi);
          }
          div_acc += div_h;
          div_all.push_back(div_h);
        }
        const Complex mean = -(double)Nf * div_acc / (double)nhits;
        double vr=0.0, vi=0.0;
        const Complex m0 = div_acc/(double)nhits;
        for(const Complex& x : div_all){ vr+=std::pow(x.real()-m0.real(),2); vi+=std::pow(x.imag()-m0.imag(),2); }
        const double n1=(double)(nhits-1)*nhits;
        const double er=(double)Nf*std::sqrt(vr/n1), ei=(double)Nf*std::sqrt(vi/n1);
        const bool ok = (std::abs(mean.real()) < 4.0*er+1.0e-8 && std::abs(mean.imag()) < 4.0*ei+1.0e-8);
        std::cout << "#   (B) stoch sum_z J_V^{wz} = (" << mean.real() << "+-" << er
                  << ", " << mean.imag() << "+-" << ei << ")  (nhits=" << nhits << ")  -> "
                  << (ok?"PASS (within 4 sigma of 0)":"CHECK") << std::endl;
        all_ok = all_ok && ok;
      }
    }
  }

  std::cout << "\n# ===== Ward diagonal-mass check: " << (all_ok ? "ALL PASS" : "SOME FAIL") << " =====" << std::endl;
  CUDA_CHECK(cudaFree(d_xi));
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
  return all_ok ? 0 : 1;
}
