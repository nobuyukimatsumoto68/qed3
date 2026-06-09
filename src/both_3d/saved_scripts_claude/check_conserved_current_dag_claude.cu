// Check program for the adjoint conserved-current kernel K^{wz\dag} (apply_k_dag).
// D2: adjoint-consistency  <eta, K phi> = <K^dag eta, phi>  per link.
// D3: commutator  [D_ov^dag, Theta] xi = -sum_l dtheta_l K^{l dag} xi.
// D4: tr(K^dag) (exact + stochastic) -> 0.
// D5: stochastic Ward identity (free-field + Gaussian U), analog of Chunks 5b/6.
// Mirrors check_conserved_current_claude.cu with apply_k -> apply_k_dag.

#include <typeinfo>
#include <cmath>
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
#include <random>
#include <Eigen/Dense>

using Double  = double;
using Idx     = std::int32_t;
using Complex = std::complex<double>;

using Link = std::array<Idx,2>;
using Face = std::vector<Idx>;

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

  constexpr int N_REFINE = 1;
  constexpr int NS       = 2;
  constexpr int Nt       = 4;

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

#include "sparse_dirac.h"
// #include "matpoly.h"
#include "matpoly_claude.h"

#include "overlap.h"
#include "conserved_current_claude.h"

#include <getopt.h>

// ---- type aliases (global scope, after all includes) ----
using Base        = S2Simp;
using WilsonDirac = DiracExt<Base, DiracS2Simp>;
using Gauge       = GaugeExt<Base, Comp::Nt, Comp::is_compact>;
using Fermion     = Overlap<WilsonDirac>;

// ---- CLI ----
void PrintHelp() {
  printf("Usage: ./a.out [options]\n");
  printf("  --gsq   <gsq>              Wilson coupling squared (default: 8.0)\n");
  printf("  --Nf    <Nf>               number of fermion flavors (default: 2)\n");
  printf("  --nu0   <nu0>              asymmetry parameter nu_0 (default: 1.0)\n");
  printf("  --gauge <trivial|gaussian|both>  gauge for adjoint-consistency check (default: both)\n");
  printf("  --nhits <nhits>            stochastic hits for tr(K^dag) / Ward identity (default: 16)\n");
  printf("  -h, --help                 show this help\n");
  exit(0);
}

void ParseArgs(int argc, char** argv,
               double& gsq, int& Nf, double& nu0,
               std::string& gauge_type, int& nhits) {
  const char* const short_opts = "gNnGHh";
  const option long_opts[] = {
    {"gsq",   required_argument, nullptr, 'g'},
    {"Nf",    required_argument, nullptr, 'N'},
    {"nu0",   required_argument, nullptr, 'n'},
    {"gauge", required_argument, nullptr, 'G'},
    {"nhits", required_argument, nullptr, 'H'},
    {"help",  no_argument,       nullptr, 'h'},
    {nullptr, no_argument,       nullptr,  0}
  };
  int option_index, opt;
  while((opt = getopt_long(argc, argv, short_opts, long_opts, &option_index)) != -1) {
    switch(opt) {
    case 'g': gsq        = std::stod(optarg); break;
    case 'N': Nf         = std::stoi(optarg); break;
    case 'n': nu0        = std::stod(optarg); break;
    case 'G': gauge_type = std::string(optarg); break;
    case 'H': nhits      = std::stoi(optarg); break;
    case 'h':
    default:  PrintHelp(); break;
    }
  }
}

// ---- helpers ----

// fill len complex entries at h_v with Gaussian random values
void fill_random(Complex* h_v, int len, std::mt19937& gen) {
  std::normal_distribution<double> dist(0.0, 1.0);
  for(int i = 0; i < len; i++) h_v[i] = Complex(dist(gen), dist(gen));
}

// multTheta: (Theta v)[site, sigma] = h_theta[site] * v[site, sigma]
// h_theta has one real scalar per spacetime site, ordered s*N_SITES + ix.
void multTheta(CuC* d_out, const CuC* d_in,
               const std::vector<double>& h_theta) {
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

// ---- main ----
int main(int argc, char* argv[]) {
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  double      gsq        = 8.0;
  int         Nf         = 2;
  double      nu0        = 1.0;
  std::string gauge_type = "both";
  int         nhits      = 16;
  ParseArgs(argc, argv, gsq, Nf, nu0, gauge_type, nhits);
  std::cout << "# gsq=" << gsq << "  Nf=" << Nf << "  nu0=" << nu0
            << "  gauge=" << gauge_type << "  nhits=" << nhits << std::endl;

  // --- GPU setup ---
  for(int i = 0; i < Comp::NSTREAMS; i++) d_MemorySets[i].allocate();
  int device;
  CUDA_CHECK(cudaGetDeviceCount(&device));
  cudaDeviceProp device_prop[device];
  cudaGetDeviceProperties(&device_prop[0], 0);
  std::cout << "# dev = " << device_prop[0].name << std::endl;
  CUDA_CHECK(cudaSetDevice(0));

  constexpr Idx N  = Comp::N;
  constexpr int Nt = Comp::Nt;

  // --- lattice, operators, gauge field ---
  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. n_sites=" << base.n_sites
            << "  n_links=" << base.n_links << std::endl;

  const double at = 0.2;
  const double M5 = -1.0;
  if(Nt != 1)
    assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);

  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);
  std::cout << "# DW set." << std::endl;

  Fermion D(DW, 21);
  std::cout << "# D set." << std::endl;

  Gauge U(base);   // free-field: all link phases = 0
  D.update(U);
  std::cout << "# D updated (free-field)." << std::endl;

  ConservedCurrent<Fermion,Gauge> kop(D);
  std::cout << "# ConservedCurrent set." << std::endl;
  std::cout << "# lambda_max = " << D.lambda_max << std::endl;

  // K applied through MatPoly::from_cpu (same path as the production measurement); configure
  // the link via kop.set(U, el, dag) just before each from_cpu.
  MatPoly op_K; op_K.push_back(cplx(1.0), {&kop});

  // --- device vectors ---
  CuC* d_phi    = nullptr;   // random phi (also reused as xi in the D3 commutator micro-test)
  CUDA_CHECK(cudaMalloc(&d_phi,   N*CD));

  // --- host vectors ---
  FermionVector h_phi, h_eta, h_Kphi, h_Kdeta;

  // draw random eta, phi; copy phi to device for the D3 micro-test
  std::mt19937 rng_host(42);
  fill_random(h_phi.field, N, rng_host);
  fill_random(h_eta.field, N, rng_host);
  CUDA_CHECK(cudaMemcpy(d_phi, reinterpret_cast<const CuC*>(h_phi.field), N*CD, H2D));
  std::cout << "# random eta, phi loaded." << std::endl;

  // ---- D2: adjoint-consistency  <eta, K phi> = <K^dag eta, phi>  per link ----
  // For each spatial and temporal link el: compare inner products of K phi and K^dag eta,
  // both applied via op_K.from_cpu (the production path). They must agree to solver tolerance.
  {
    auto run_adjoint_check = [&](const std::string& label) {
      std::cout << "\n# --- D2: adjoint-consistency check (" << label << ") ---" << std::endl;
      double max_err = 0.0; int n = 0;

      auto check_link = [&](auto el) {
        kop.set(U, el, /*dag=*/false);
        op_K.from_cpu<N>(h_Kphi.field, h_phi.field);    // K phi
        const Complex lhs = h_eta.dag(h_Kphi);          // <eta, K phi>

        kop.set(U, el, /*dag=*/true);
        op_K.from_cpu<N>(h_Kdeta.field, h_eta.field);   // K^dag eta
        const Complex rhs = h_Kdeta.dag(h_phi);         // <K^dag eta, phi>

        const double err = std::abs(lhs - rhs);
        if(err > max_err) max_err = err; n++;
      };

      // spatial links
      for(int s_lat = 0; s_lat < Nt; s_lat++)
        for(const auto& lk : base.links)
          check_link(std::pair<int,BaseLink>{s_lat, lk});
      // temporal links (skip for Nt=1)
      if(Nt > 1)
        for(int s_lat = 0; s_lat < Nt; s_lat++)
          for(int ix = 0; ix < base.n_sites; ix++)
            check_link(std::pair<int,Idx>{s_lat, ix});

      std::cout << "# n_links=" << n
                << "  max|<eta,K phi> - <K^dag eta,phi>|=" << max_err << std::endl;
      // bound set by inner-solver tolerance TOL_INNER=1e-9 compounded over the
      // R_m solves in apply_k / apply_k_dag; free-field is exact (~1e-15).
      assert(max_err < 1.0e-7);
      std::cout << "# PASS" << std::endl;
    };

    using Rng = ParallelRngExt<Base, Comp::Nt>;
    Rng rng(base, 42);
    const Gauge U_trivial = U;
    if(gauge_type == "trivial" || gauge_type == "both") run_adjoint_check("trivial");
    if(gauge_type == "gaussian" || gauge_type == "both") {
      U.gaussian(rng, 1.0); D.update(U); run_adjoint_check("gaussian");
      U = U_trivial; D.update(U);
    }
  }

  // ---- D3: commutator check  [D_ov^dag, Theta] xi = -sum_l dtheta_l K^{l dag} xi ----
  // Adjoint of the forward [D_ov,Theta] check: daggering sum_z K^{wz} = [D_ov, Theta_w]
  // gives sum_z K^{wz\dag} = [D_ov, Theta_w]^dag = -[D_ov^dag, Theta_w].
  // Runs on the free-field U restored at the end of D2.
  {
    auto f_DH = std::bind(&Fermion::adj_deviceAsyncLaunch, &D,
                          std::placeholders::_1, std::placeholders::_2);
    LinOpWrapper M_DH(f_DH);
    MatPoly op_DH;
    op_DH.push_back(cplx(1.0), {&M_DH});

    std::cout << "\n# --- D3: commutator [D_ov^dag,Theta]xi = -sum_l dtheta_l K^{l dag} xi ---" << std::endl;

    CuC* d_xi = d_phi; // reuse the random phi as xi (unmodified by D2)

    // sanity: Theta = c*I => [D_ov^dag, Theta] = 0 and all dtheta = 0
    {
      const double c = 2.5;
      std::vector<double> h_theta_c(Comp::N_SITES * Nt, c);
      CuC* d_c1 = nullptr; CuC* d_tmp = nullptr; CuC* d_th = nullptr;
      CUDA_CHECK(cudaMalloc(&d_c1,  N*CD));
      CUDA_CHECK(cudaMalloc(&d_tmp, N*CD));
      CUDA_CHECK(cudaMalloc(&d_th,  N*CD));
      multTheta(d_th, d_xi, h_theta_c);
      op_DH.on_gpu<N>(d_c1, d_th);     CUDA_CHECK(cudaDeviceSynchronize());
      op_DH.on_gpu<N>(d_tmp, d_xi);    CUDA_CHECK(cudaDeviceSynchronize());
      multTheta(d_tmp, d_tmp, h_theta_c);
      Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_c1, -1.0, d_tmp, d_c1);
      CUDA_CHECK(cudaDeviceSynchronize());
      FermionVector h1;
      CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(h1.field), d_c1, N*CD, D2H));
      double e1 = 0.0;
      for(int i = 0; i < N; i++) e1 = std::max(e1, std::abs(h1.field[i]));
      std::cout << "# Theta=c*I: ||chi1||_inf=" << e1 << std::endl;
      CUDA_CHECK(cudaFree(d_c1)); CUDA_CHECK(cudaFree(d_tmp)); CUDA_CHECK(cudaFree(d_th));
    }

    // random theta
    std::vector<double> h_theta(Comp::N_SITES * Nt);
    {
      std::normal_distribution<double> rdist(0.0, 1.0);
      for(auto& t : h_theta) t = rdist(rng_host);
    }

    CuC* d_Theta_xi = nullptr;
    CuC* d_chi1     = nullptr;
    CuC* d_chi2     = nullptr;
    CuC* d_tmp_th   = nullptr;
    CUDA_CHECK(cudaMalloc(&d_Theta_xi, N*CD));
    CUDA_CHECK(cudaMalloc(&d_chi1,     N*CD));
    CUDA_CHECK(cudaMalloc(&d_chi2,     N*CD));
    CUDA_CHECK(cudaMalloc(&d_tmp_th,   N*CD));

    // chi1 = [D_ov^dag, Theta] xi = D_ov^dag(Theta xi) - Theta(D_ov^dag xi)
    multTheta(d_Theta_xi, d_xi, h_theta);
    op_DH.on_gpu<N>(d_chi1, d_Theta_xi);   CUDA_CHECK(cudaDeviceSynchronize());
    op_DH.on_gpu<N>(d_tmp_th, d_xi);       CUDA_CHECK(cudaDeviceSynchronize());
    multTheta(d_tmp_th, d_tmp_th, h_theta);
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_chi1, -1.0, d_tmp_th, d_chi1);
    CUDA_CHECK(cudaDeviceSynchronize());

    // chi2 = - sum_l dtheta_l K^{l dag} xi  (note the overall minus sign)
    CUDA_CHECK(cudaMemset(d_chi2, 0, N*CD));
    // spatial links: canonical orientation lk={ix_w, ix_z}, dtheta = theta_w - theta_z
    for(int s = 0; s < Nt; s++) {
      for(const auto& lk : base.links) {
        const double dt = h_theta[s*Comp::N_SITES + lk[0]]
                        - h_theta[s*Comp::N_SITES + lk[1]];
        kop.apply_k_dag(d_tmp_th, d_xi, U, std::pair<int,BaseLink>{s, lk});
        CUDA_CHECK(cudaDeviceSynchronize());
        Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_chi2, -dt, d_tmp_th, d_chi2);
        CUDA_CHECK(cudaDeviceSynchronize());
      }
    }
    // temporal links: skip for Nt=1
    if(Nt > 1) {
      for(int s = 0; s < Nt; s++) {
        const int s1 = (s + 1) % Nt;
        for(int ix = 0; ix < base.n_sites; ix++) {
          const double dt = h_theta[s*Comp::N_SITES + ix]
                          - h_theta[s1*Comp::N_SITES + ix];
          kop.apply_k_dag(d_tmp_th, d_xi, U, std::pair<int,Idx>{s, ix});
          CUDA_CHECK(cudaDeviceSynchronize());
          Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_chi2, -dt, d_tmp_th, d_chi2);
          CUDA_CHECK(cudaDeviceSynchronize());
        }
      }
    }

    FermionVector h_chi1, h_chi2;
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(h_chi1.field), d_chi1, N*CD, D2H));
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(h_chi2.field), d_chi2, N*CD, D2H));
    double max_err = 0.0;
    for(int i = 0; i < N; i++)
      max_err = std::max(max_err, std::abs(h_chi1.field[i] - h_chi2.field[i]));
    std::cout << "# D_ov^dag commutator check: ||chi1 - chi2||_inf = " << max_err << std::endl;
    // bound set by the overlap pole solves inside D_ov^dag and apply_k_dag
    assert(max_err < 1.0e-6);
    std::cout << "# PASS" << std::endl;

    CUDA_CHECK(cudaFree(d_Theta_xi));
    CUDA_CHECK(cudaFree(d_chi1));
    CUDA_CHECK(cudaFree(d_chi2));
    CUDA_CHECK(cudaFree(d_tmp_th));
  }

  // ---- D4: tr(K^dag) checks (exact and stochastic) at w=(s=0,ix=0), per gauge background ----
  // tr(K^{wz\dag}) = conj(tr(K^{wz})) = 0 analytically; stochastic eta^dag K^dag eta -> 0.
  {
    using RngTrK = ParallelRngExt<Base, Comp::Nt>;
    RngTrK rng_gauge(base, 99);
    const Gauge U_free = U;

    auto run_trKdag = [&](const std::string& glabel) {
      constexpr int s_focus  = 0;
      constexpr int ix_focus = 0;

      std::vector<Idx> nns_list(base.nns[ix_focus].begin(), base.nns[ix_focus].end());
      const int n_spatial = static_cast<int>(nns_list.size());

      // analytic tr(K^dag) values, filled by the exact block and used by the stochastic summary
      std::vector<Complex> trK_exact(n_spatial, Complex(0.0, 0.0));
      Complex trK_exact_t[2] = {Complex(0.0,0.0), Complex(0.0,0.0)};

      // ---- Exact tr(K^dag): basis-vector loop (small lattice only) ----
      // tr(K^{wz\dag}) = sum_k e_k^dag K^{wz\dag} e_k; flat index k: s=k/Nx, ix=(k%Nx)/NS, spin=k%NS.
      if constexpr (Comp::N_REFINE == 1 && Comp::Nt <= 4) {
        std::cout << "\n# --- tr(K^dag) exact (" << glabel << ") at (s=0, ix=0) ---" << std::endl;

        FermionVector e_i, fv_Ke;

        auto trK_link = [&](auto el) -> Complex {
          kop.set(U, el, /*dag=*/true);
          Complex trK(0.0, 0.0);
          for(Idx k = 0; k < N; k++) {
            e_i.set_pt_source(k / Comp::Nx, (k % Comp::Nx) / NS, k % NS);
            op_K.from_cpu<N>(fv_Ke.field, e_i.field);   // K^dag e_k
            trK += fv_Ke.field[k];
          }
          return trK;
        };

        for(int j = 0; j < n_spatial; j++) {
          trK_exact[j] = trK_link(std::pair<int,BaseLink>{s_focus, BaseLink{ix_focus, nns_list[j]}});
          std::cout << "# trKdag  nn" << nns_list[j]
                    << " = (" << trK_exact[j].real() << "," << trK_exact[j].imag() << ")" << std::endl;
        }

        if(Nt > 1) {
          trK_exact_t[0] = trK_link(std::pair<int,Idx>{s_focus, ix_focus});
          trK_exact_t[1] = trK_link(std::pair<int,Idx>{(s_focus-1+Nt)%Nt, ix_focus});
          std::cout << "# trKdag  tlink_fwd"
                    << " = (" << trK_exact_t[0].real() << "," << trK_exact_t[0].imag() << ")" << std::endl;
          std::cout << "# trKdag  tlink_bwd"
                    << " = (" << trK_exact_t[1].real() << "," << trK_exact_t[1].imag() << ")" << std::endl;
        }
      }

      // ---- Stochastic tr(K^dag): Z2xZ2 noise ----
      // tr(K^{wz\dag}) ~ (1/nhits) sum_h eta_h^dag K^{wz\dag} eta_h; compared to exact value above.
      {
        std::cout << "\n# --- stochastic tr(K^dag) (" << glabel << ") at (s=0, ix=0) ---" << std::endl;

        RngTrK rng_trk(base, 42);

        FermionVector eta, fv_Keta;

        std::vector<Complex> acc(n_spatial, Complex(0.0, 0.0));
        std::vector<std::vector<Complex>> all_hits(n_spatial);
        // temporal links: index 0 = fwd (s_focus, ix_focus), 1 = bwd ((s_focus-1+Nt)%Nt, ix_focus)
        Complex acc_t[2] = {Complex(0.0,0.0), Complex(0.0,0.0)};
        std::vector<Complex> all_hits_t[2];

        for(int h = 0; h < nhits; h++) {
          eta.fill_z2_source(rng_trk);
          for(int j = 0; j < n_spatial; j++) {
            const BaseLink lk_j{ix_focus, nns_list[j]};
            kop.set(U, std::pair<int,BaseLink>{s_focus, lk_j}, /*dag=*/true);
            op_K.from_cpu<N>(fv_Keta.field, eta.field);   // K^dag eta
            const Complex c = eta.dag(fv_Keta);
            acc[j] += c;
            all_hits[j].push_back(c);
          }
          if(Nt > 1) {
            std::pair<int,Idx> tlinks[2] = {
              {s_focus, ix_focus},
              {(s_focus-1+Nt)%Nt, ix_focus}
            };
            for(int t = 0; t < 2; t++) {
              kop.set(U, tlinks[t], /*dag=*/true);
              op_K.from_cpu<N>(fv_Keta.field, eta.field);   // K^dag eta
              const Complex c = eta.dag(fv_Keta);
              acc_t[t] += c;
              all_hits_t[t].push_back(c);
            }
          }
          const double norm = 1.0 / (h+1);
          std::cout << "# hit " << std::setw(4) << (h+1);
          for(int j = 0; j < n_spatial; j++)
            std::cout << "  nn" << nns_list[j]
                      << "=(" << acc[j].real()*norm << "," << acc[j].imag()*norm << ")";
          if(Nt > 1)
            std::cout << "  tlink_fwd=(" << acc_t[0].real()*norm << "," << acc_t[0].imag()*norm << ")"
                      << "  tlink_bwd=(" << acc_t[1].real()*norm << "," << acc_t[1].imag()*norm << ")";
          std::cout << std::endl;
        }

        auto print_stderr = [&](const std::string& linklabel, Complex sum,
                                const std::vector<Complex>& samples, Complex exact) {
          const double mean_re = sum.real() / nhits;
          const double mean_im = sum.imag() / nhits;
          double vr = 0.0, vi = 0.0;
          for(int h = 0; h < nhits; h++) {
            vr += std::pow(samples[h].real() - mean_re, 2);
            vi += std::pow(samples[h].imag() - mean_im, 2);
          }
          const double n1 = static_cast<double>(nhits - 1) * nhits;
          std::cout << "#   " << linklabel
                    << "  stoch=(" << mean_re << "+-" << std::sqrt(vr/n1)
                    << ", "        << mean_im << "+-" << std::sqrt(vi/n1) << ")";
          if constexpr (Comp::N_REFINE == 1 && Comp::Nt <= 4)
            std::cout << "  exact=(" << exact.real() << "," << exact.imag() << ")";
          std::cout << std::endl;
        };

        std::cout << "# stoch trKdag mean +/- stderr  vs  exact:" << std::endl;
        for(int j = 0; j < n_spatial; j++)
          print_stderr("nn=" + std::to_string(nns_list[j]), acc[j], all_hits[j], trK_exact[j]);
        if(Nt > 1) {
          print_stderr("tlink_fwd", acc_t[0], all_hits_t[0], trK_exact_t[0]);
          print_stderr("tlink_bwd", acc_t[1], all_hits_t[1], trK_exact_t[1]);
        }
      }
    }; // end run_trKdag

    run_trKdag("free");

    U.gaussian(rng_gauge, 1.0);
    D.update(U);
    run_trKdag("gaussian");

    U = U_free;
    D.update(U);
  }

  // ---- D5: stochastic Ward identity for K^dag at w=(s=0,ix=0), free + Gaussian U ----
  // Adjoint of the forward Ward identity. Daggering sum_z K^{wz} = [D_ov, Theta_w] and using
  // cyclicity gives, config-by-config (both Re and Im):
  //   sum_{z~w} tr(K^{wz dag} D_ov^{-dag}) = -tr([D_ov^dag, Theta_w] D_ov^{-dag}) = 0.
  // Define J_dag^{wz} = -Nf tr(K^{wz dag} D_ov^{-dag}); the divergence sum_z J_dag^{wz} -> 0.
  // Stochastic: J_dag^{wz} ~ -Nf mean( eta^dag K^{wz dag} psi ), with psi = D^{-dag} eta solved
  // from the normal equations (D D^dag) psi = D eta.
  {
    using Rng5 = ParallelRngExt<Base, Comp::Nt>;
    Rng5 rng5_gauge(base, 99);
    const Gauge U_free5 = U;
    const int Nf5 = Nf;
    constexpr int s_focus = 0;
    constexpr int ix_focus = 0;

    auto f_D   = std::bind(&Fermion::mult_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);
    auto f_DDH = std::bind(&Fermion::DDH_deviceAsyncLaunch,  &D, std::placeholders::_1, std::placeholders::_2);
    LinOpWrapper M_D(f_D), M_DDH(f_DDH);
    MatPoly op_D;   op_D.push_back(cplx(1.0), {&M_D});
    MatPoly op_DDH; op_DDH.push_back(cplx(1.0), {&M_DDH});

    FermionVector eta, D_eta, psi, fv_Kpsi;

    const std::vector<Idx> nns5(base.nns[ix_focus].begin(), base.nns[ix_focus].end());
    const int n_sp5 = static_cast<int>(nns5.size());

    // temporal links incident to w: index 0 = fwd (sign +1), 1 = bwd (sign -1)
    const std::pair<int,Idx> tlinks5[2] = { {s_focus, ix_focus}, {(s_focus-1+Nt)%Nt, ix_focus} };
    const double tsign5[2] = {1.0, -1.0};

    auto run_ward = [&](const std::string& glabel) {
      std::cout << "\n# --- D5: Ward identity K^dag (" << glabel << ") at (s=0, ix=0) ---" << std::endl;

      // ---- analytic J_dag via basis-vector loop (small lattice only) ----
      // J_dag^{wz} = -Nf sum_k (K^{wz dag} psi_k)[k], with psi_k = D^{-dag} e_k.
      std::vector<Complex> JV_exact(n_sp5, Complex(0.0,0.0));
      Complex JV_exact_t[2] = {Complex(0.0,0.0), Complex(0.0,0.0)};
      Complex div_exact = Complex(0.0,0.0);
      if constexpr (Comp::N_REFINE == 1 && Comp::Nt <= 4) {
        std::cout << "# --- analytic J_dag via basis-vector loop ---" << std::endl;
        FermionVector e_i, D_e, psi_e, fv_Ke;
        for(Idx k = 0; k < N; k++) {
          e_i.set_pt_source(k / Comp::Nx, (k % Comp::Nx) / NS, k % NS);
          op_D.from_cpu<N>(D_e.field, e_i.field);                       // D e_k
          op_DDH.solve<N>(psi_e.field, D_e.field, Comp::TOL_OUTER);     // psi_e = D^{-dag} e_k
          for(int j = 0; j < n_sp5; j++) {
            kop.set(U, std::pair<int,BaseLink>{s_focus, BaseLink{ix_focus, nns5[j]}}, /*dag=*/true);
            op_K.from_cpu<N>(fv_Ke.field, psi_e.field);   // K^dag psi_k
            JV_exact[j] += fv_Ke.field[k];
          }
          if(Nt > 1) {
            for(int t = 0; t < 2; t++) {
              kop.set(U, tlinks5[t], /*dag=*/true);
              op_K.from_cpu<N>(fv_Ke.field, psi_e.field);   // K^dag psi_k
              JV_exact_t[t] += tsign5[t] * fv_Ke.field[k];
            }
          }
          if((k+1) % (N/4) == 0) std::cout << "# basis loop: k=" << (k+1) << "/" << N << std::endl;
        }
        for(int j = 0; j < n_sp5; j++) {
          JV_exact[j] *= -(double)Nf5;
          div_exact += JV_exact[j];
          std::cout << "# J_dag exact  nn=" << nns5[j]
                    << "  re=" << JV_exact[j].real() << "  im=" << JV_exact[j].imag() << std::endl;
        }
        if(Nt > 1) {
          for(int t = 0; t < 2; t++) JV_exact_t[t] *= -(double)Nf5;
          div_exact += JV_exact_t[0] + JV_exact_t[1];
          std::cout << "# J_dag exact  tlink_fwd  re=" << JV_exact_t[0].real() << "  im=" << JV_exact_t[0].imag() << std::endl;
          std::cout << "# J_dag exact  tlink_bwd  re=" << JV_exact_t[1].real() << "  im=" << JV_exact_t[1].imag() << std::endl;
        }
        std::cout << "# J_dag exact  SUM  re=" << div_exact.real() << "  im=" << div_exact.imag() << std::endl;
      }

      // ---- stochastic estimator ----
      Rng5 rng5_stoch(base, 42);
      std::vector<Complex> C_acc(n_sp5, Complex(0.0,0.0));
      std::vector<std::vector<Complex>> C_all(n_sp5);
      Complex C_t_acc[2] = {Complex(0.0,0.0), Complex(0.0,0.0)};
      std::vector<Complex> C_t_all[2];
      Complex div_acc = Complex(0.0,0.0);
      std::vector<Complex> div_all;

      for(int h = 0; h < nhits; h++) {
        eta.fill_z2_source(rng5_stoch);
        op_D.from_cpu<N>(D_eta.field, eta.field);                   // D eta
        op_DDH.solve<N>(psi.field, D_eta.field, Comp::TOL_OUTER);   // psi = D^{-dag} eta

        Complex div_h(0.0,0.0);
        for(int j = 0; j < n_sp5; j++) {
          kop.set(U, std::pair<int,BaseLink>{s_focus, BaseLink{ix_focus, nns5[j]}}, /*dag=*/true);
          op_K.from_cpu<N>(fv_Kpsi.field, psi.field);   // K^dag psi
          const Complex c = eta.dag(fv_Kpsi);
          C_acc[j] += c; C_all[j].push_back(c); div_h += c;
        }
        if(Nt > 1) {
          for(int t = 0; t < 2; t++) {
            kop.set(U, tlinks5[t], /*dag=*/true);
            op_K.from_cpu<N>(fv_Kpsi.field, psi.field);   // K^dag psi
            const Complex c = eta.dag(fv_Kpsi);
            C_t_acc[t] += tsign5[t] * c; C_t_all[t].push_back(tsign5[t] * c); div_h += tsign5[t] * c;
          }
        }
        div_acc += div_h; div_all.push_back(div_h);

        const double norm = 1.0 / (h+1);
        std::cout << "# hit " << std::setw(4) << (h+1)
                  << "  div_re=" << -((double)Nf5)*div_acc.real()*norm
                  << "  div_im=" << -((double)Nf5)*div_acc.imag()*norm << std::endl;
      }

      auto print_ward = [&](const std::string& label, const std::vector<Complex>& samples, Complex exact) {
        double mean_im = 0.0, mean_re = 0.0;
        for(const Complex& x : samples) { mean_im += x.imag(); mean_re += x.real(); }
        mean_im /= nhits; mean_re /= nhits;
        const double stoch_re = -((double)Nf5) * mean_re;
        const double stoch_im = -((double)Nf5) * mean_im;
        double var_im = 0.0, var_re = 0.0;
        for(const Complex& x : samples) {
          var_im += std::pow(x.imag() - mean_im, 2);
          var_re += std::pow(x.real() - mean_re, 2);
        }
        const double n1 = static_cast<double>(nhits - 1) * nhits;
        const double err_re = ((double)Nf5) * std::sqrt(var_re / n1);
        const double err_im = ((double)Nf5) * std::sqrt(var_im / n1);
        std::cout << "#   J_dag  " << label;
        if constexpr (Comp::N_REFINE == 1 && Comp::Nt <= 4)
          std::cout << "  Re_analy=" << exact.real();
        std::cout << "  Re_stoch=" << stoch_re << "  err_Re=" << err_re;
        if constexpr (Comp::N_REFINE == 1 && Comp::Nt <= 4)
          std::cout << "  sigma_Re=" << (stoch_re - exact.real()) / err_re;
        if constexpr (Comp::N_REFINE == 1 && Comp::Nt <= 4)
          std::cout << "  Im_analy=" << exact.imag();
        std::cout << "  Im_stoch=" << stoch_im << "  err_Im=" << err_im;
        if constexpr (Comp::N_REFINE == 1 && Comp::Nt <= 4)
          std::cout << "  sigma_Im=" << (stoch_im - exact.imag()) / err_im;
        std::cout << std::endl;
      };

      std::cout << "# per-link J_dag stoch vs exact (Nhits=" << nhits << ", Nf=" << Nf5 << ", " << glabel << "):" << std::endl;
      for(int j = 0; j < n_sp5; j++)
        print_ward("nn=" + std::to_string(nns5[j]), C_all[j], JV_exact[j]);
      if(Nt > 1) {
        print_ward("tlink_fwd", C_t_all[0], JV_exact_t[0]);
        print_ward("tlink_bwd", C_t_all[1], JV_exact_t[1]);
      }
      print_ward("SUM", div_all, div_exact);
    }; // end run_ward

    run_ward("free");

    U.gaussian(rng5_gauge, 2.0);
    D.update(U);
    run_ward("gaussian");

    U = U_free5;
    D.update(U);
  }

  // --- cleanup ---
  CUDA_CHECK(cudaFree(d_phi));
  for(int i = 0; i < Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();

  return 0;
}
