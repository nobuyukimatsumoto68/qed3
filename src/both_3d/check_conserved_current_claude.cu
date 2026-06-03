// Check program for the exactly conserved vector current kernel K^{wz} (eq. (3.34) of qed3int_v2-2.pdf).
// Chunks 1-3: boilerplate, construction, helper functions.
// Chunk 4: force-based cross-check Im[<xi,K xi>] = (lmax/2)*F^{wz}.
// Chunks 5-6 (free-field divergence, Ward identity) to be added.

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

// ---- Chunk 2: type aliases (global scope, after all includes) ----
using Base        = S2Simp;
using WilsonDirac = DiracExt<Base, DiracS2Simp>;
using Gauge       = GaugeExt<Base, Comp::Nt, Comp::is_compact>;
using Fermion     = Overlap<WilsonDirac>;

// ---- Chunk 1: CLI ----
void PrintHelp() {
  printf("Usage: ./a.out [options]\n");
  printf("  --gsq   <gsq>              Wilson coupling squared (default: 8.0)\n");
  printf("  --Nf    <Nf>               number of fermion flavors (default: 2)\n");
  printf("  --nu0   <nu0>              asymmetry parameter nu_0 (default: 1.0)\n");
  printf("  --gauge <trivial|gaussian|both>  gauge for Chunk 4 force check (default: both)\n");
  printf("  --nhits <nhits>            stochastic hits for Ward identity (default: 16)\n");
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

// ---- Chunk 3: helper functions ----

// use FermionVector::dag(b) for inner products; replaced inner_product.
// use cudaMemcpy(reinterpret_cast<CuC*>(fv.field), d_ptr, N*CD, D2H) for D2H; replaced device_to_host_vec.
/*
Complex inner_product(const std::vector<Complex>& eta, const std::vector<Complex>& v);
void device_to_host_vec(std::vector<Complex>& h_res, const CuC* d_res, int len);
*/

// accumulate scalar divergence: h_div[i_pos] += j, h_div[i_neg] -= j.
// j = 2*Re[<xi, K(el) xi>] is the oriented current through a link.
// i_pos/i_neg are flat site indices (s*n_sites + ix) for the two endpoints.
void accumulate_divergence(std::vector<double>& h_div, double j, int i_pos, int i_neg) {
  h_div[i_pos] += j;
  h_div[i_neg] -= j;
}

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

  // --- Chunk 2: lattice, operators, gauge field ---
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

  ConservedCurrent<Fermion> kop(D);
  std::cout << "# ConservedCurrent set." << std::endl;
  std::cout << "# lambda_max = " << D.lambda_max << std::endl;

  // --- device vectors ---
  CuC* d_xi   = nullptr;
  CuC* d_k_xi = nullptr;
  CUDA_CHECK(cudaMalloc(&d_xi,   N*CD));
  CUDA_CHECK(cudaMalloc(&d_k_xi, N*CD));

  // --- host vectors ---
  FermionVector h_xi, h_k_xi;

  // draw random \xi, copy to device
  std::mt19937 rng_host(42);
  fill_random(h_xi.field, N, rng_host);
  CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(h_xi.field), N*CD, H2D));
  std::cout << "# random xi loaded." << std::endl;

  // ---- Chunk 4: force-based cross-check ----
  {
    auto run_force_check = [&](const std::string& label) {
      std::cout << "\n# --- Chunk 4: force check (" << label << ") ---" << std::endl;
      D.precalc_grad_deviceAsyncLaunch(U, d_xi);
      double max_err = 0.0; int n = 0;
      // spatial links
      for(int s_lat = 0; s_lat < Nt; s_lat++) {
        for(const auto& lk : base.links) {
          const double force = D.grad_deviceAsyncLaunch(
              std::pair<int, BaseLink>{s_lat, lk}, U, d_xi);
          kop.apply_k_wz(d_k_xi, d_xi, U, s_lat, lk);
          CUDA_CHECK(cudaDeviceSynchronize());
          CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(h_k_xi.field), d_k_xi, N*CD, D2H));
          const Complex inner = h_xi.dag(h_k_xi);
          const double err = std::abs(std::imag(inner) + 0.5 * force);
          if(err > max_err) max_err = err; n++;
        }
      }
      // temporal links (skip for Nt=1: no temporal links)
      if(Nt > 1) {
        for(int s_lat = 0; s_lat < Nt; s_lat++) {
          for(int ix = 0; ix < base.n_sites; ix++) {
            const double force = D.grad_deviceAsyncLaunch(
                std::pair<int, Idx>{s_lat, ix}, U, d_xi);
            kop.apply_k(d_k_xi, d_xi, U, std::pair<int, Idx>{s_lat, ix});
            CUDA_CHECK(cudaDeviceSynchronize());
            CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(h_k_xi.field), d_k_xi, N*CD, D2H));
            const Complex inner = h_xi.dag(h_k_xi);
            const double err = std::abs(std::imag(inner) + 0.5 * force);
            if(err > max_err) max_err = err; n++;
          }
        }
      }
      std::cout << "# n_spatial=" << (Nt*(int)base.links.size())
                << "  n_temporal=" << (Nt*base.n_sites)
                << "  max|Im(<xi,K xi>) + 0.5*F|=" << max_err << std::endl;
      assert(max_err < 1.0e-4);
      std::cout << "# PASS" << std::endl;
    };
    using Rng = ParallelRngExt<Base, Comp::Nt>;
    Rng rng(base, 42);
    const Gauge U_trivial = U;
    if(gauge_type == "trivial" || gauge_type == "both") run_force_check("trivial");
    if(gauge_type == "gaussian" || gauge_type == "both") {
      U.gaussian(rng, 1.0); D.update(U); run_force_check("gaussian");
      U = U_trivial; D.update(U);
    }
  }

  // ---- Theta check: [X, Theta] xi = sum_l dtheta_l W^l xi, X = D_W / lambda_max ----
  // chi1 = X(Theta xi) - Theta(X xi)
  // chi2 = sum_l (theta_w - theta_z) W^l xi, over all spatial and temporal links
  // Expect ||chi1 - chi2||_inf ~ 1e-14
  {
    std::cout << "\n# --- Theta check: [X,Theta]xi = sum_l dtheta_l W^l xi ---" << std::endl;

    // sanity: Theta = c*I => [X,Theta]=0 and all dtheta=0; both chi1,chi2 must vanish
    {
      const double c = 2.5;
      std::vector<double> h_theta_c(Comp::N_SITES * Nt, c);
      CuC* d_sc1 = nullptr; CuC* d_sc2 = nullptr; CuC* d_stmp = nullptr;
      CUDA_CHECK(cudaMalloc(&d_sc1, N*CD));
      CUDA_CHECK(cudaMalloc(&d_sc2, N*CD));
      CUDA_CHECK(cudaMalloc(&d_stmp, N*CD));
      {
        MatPoly op_X; op_X.push_back(cplx(1.0/D.lambda_max), {&D.M_DW});
        CuC* d_sth = nullptr; CUDA_CHECK(cudaMalloc(&d_sth, N*CD));
        multTheta(d_sth, d_xi, h_theta_c);
        op_X.on_gpu<N>(d_sc1, d_sth); CUDA_CHECK(cudaDeviceSynchronize());
        op_X.on_gpu<N>(d_stmp, d_xi); CUDA_CHECK(cudaDeviceSynchronize());
        multTheta(d_stmp, d_stmp, h_theta_c);
        Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_sc1, -1.0, d_stmp, d_sc1);
        CUDA_CHECK(cudaDeviceSynchronize());
        CUDA_CHECK(cudaFree(d_sth));
      }
      CUDA_CHECK(cudaMemset(d_sc2, 0, N*CD));
      FermionVector h1, h2;
      CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(h1.field), d_sc1, N*CD, D2H));
      CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(h2.field), d_sc2, N*CD, D2H));
      double e1=0.0, e2=0.0;
      for(int i=0;i<N;i++) { e1=std::max(e1,std::abs(h1.field[i])); e2=std::max(e2,std::abs(h2.field[i])); }
      std::cout << "# Theta=c*I: ||chi1||_inf=" << e1 << "  ||chi2||_inf=" << e2 << std::endl;
      CUDA_CHECK(cudaFree(d_sc1)); CUDA_CHECK(cudaFree(d_sc2)); CUDA_CHECK(cudaFree(d_stmp));
    }

    // localized theta: 1 at ix=0 (all s), 0 elsewhere
    // temporal dtheta=0 for ix=0, so only spatial links touching ix=0 contribute to chi2
    {
      std::vector<double> h_thloc(Comp::N_SITES * Nt, 0.0);
      for(int s = 0; s < Nt; s++) h_thloc[s * Comp::N_SITES + 0] = 1.0;

      CuC* d_lth = nullptr; CuC* d_lc1 = nullptr;
      CuC* d_lc2 = nullptr; CuC* d_ltmp = nullptr;
      CUDA_CHECK(cudaMalloc(&d_lth,  N*CD));
      CUDA_CHECK(cudaMalloc(&d_lc1,  N*CD));
      CUDA_CHECK(cudaMalloc(&d_lc2,  N*CD));
      CUDA_CHECK(cudaMalloc(&d_ltmp, N*CD));

      {
        MatPoly op_X; op_X.push_back(cplx(1.0/D.lambda_max), {&D.M_DW});
        multTheta(d_lth, d_xi, h_thloc);
        op_X.on_gpu<N>(d_lc1, d_lth);   CUDA_CHECK(cudaDeviceSynchronize());
        op_X.on_gpu<N>(d_ltmp, d_xi);   CUDA_CHECK(cudaDeviceSynchronize());
        multTheta(d_ltmp, d_ltmp, h_thloc);
        Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_lc1, -1.0, d_ltmp, d_lc1);
        CUDA_CHECK(cudaDeviceSynchronize());
      }
      FermionVector h_lc1;
      CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(h_lc1.field), d_lc1, N*CD, D2H));
      std::cout << "# --- localized Theta (ix=0 all s): chi1=[X,Theta]xi ---" << std::endl;
      for(int i = 0; i < N; i++)
        if(std::abs(h_lc1.field[i]) > 1.0e-15)
          std::cout << "#  chi1[" << i << "] = " << h_lc1.field[i] << std::endl;

      CUDA_CHECK(cudaMemset(d_lc2, 0, N*CD));
      for(const auto& lk : base.links) {
        if(lk[0] != 0 && lk[1] != 0) continue;
        // dtheta = theta_{lk[0]} - theta_{lk[1]}: nonzero endpoint is ix=0 (theta=1), other is 0
        const double dt = (lk[0] == 0) ? 1.0 : -1.0;
        const Idx iy = (lk[0] == 0) ? lk[1] : lk[0];
        std::cout << "# --- nn iy=" << iy
                  << " lk={" << lk[0] << "," << lk[1] << "} dtheta=" << dt << " ---" << std::endl;
        for(int s = 0; s < Nt; s++) {
          COO<N> coo_W;
          kop.build_W(coo_W, U, std::pair<int,BaseLink>{s, lk});
          coo_W(d_ltmp, d_xi);
          CUDA_CHECK(cudaDeviceSynchronize());
          FermionVector h_wxi;
          CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(h_wxi.field), d_ltmp, N*CD, D2H));
          std::cout << "#  s=" << s << "  W^{lk}*xi (nonzero):" << std::endl;
          for(int i = 0; i < N; i++)
            if(std::abs(h_wxi.field[i]) > 1.0e-15)
              std::cout << "#    [" << i << "] = " << h_wxi.field[i]
                        << "  dtheta*W*xi=" << dt*h_wxi.field[i] << std::endl;
          Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_lc2, dt, d_ltmp, d_lc2);
          CUDA_CHECK(cudaDeviceSynchronize());
        }
      }
      FermionVector h_lc2;
      CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(h_lc2.field), d_lc2, N*CD, D2H));
      std::cout << "# --- chi2 = sum dtheta W^l xi (nonzero) ---" << std::endl;
      for(int i = 0; i < N; i++)
        if(std::abs(h_lc2.field[i]) > 1.0e-15)
          std::cout << "#  chi2[" << i << "] = " << h_lc2.field[i] << std::endl;
      double lerr = 0.0;
      for(int i = 0; i < N; i++) lerr = std::max(lerr, std::abs(h_lc1.field[i]-h_lc2.field[i]));
      std::cout << "# localized: ||chi1-chi2||_inf = " << lerr << std::endl;

      CUDA_CHECK(cudaFree(d_lth)); CUDA_CHECK(cudaFree(d_lc1));
      CUDA_CHECK(cudaFree(d_lc2)); CUDA_CHECK(cudaFree(d_ltmp));
    }

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

    {
      MatPoly op_X;
      op_X.push_back(cplx(1.0/D.lambda_max), {&D.M_DW});
      multTheta(d_Theta_xi, d_xi, h_theta);
      op_X.on_gpu<N>(d_chi1, d_Theta_xi);
      CUDA_CHECK(cudaDeviceSynchronize());
      op_X.on_gpu<N>(d_tmp_th, d_xi);
      CUDA_CHECK(cudaDeviceSynchronize());
      multTheta(d_tmp_th, d_tmp_th, h_theta);
      Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_chi1, -1.0, d_tmp_th, d_chi1);
      CUDA_CHECK(cudaDeviceSynchronize());
    }

    CUDA_CHECK(cudaMemset(d_chi2, 0, N*CD));
    // spatial links: canonical orientation lk={ix_w, ix_z}, dtheta = theta_w - theta_z
    for(int s = 0; s < Nt; s++) {
      for(const auto& lk : base.links) {
        const double dt = h_theta[s*Comp::N_SITES + lk[0]]
                        - h_theta[s*Comp::N_SITES + lk[1]];
        COO<N> coo_W;
        kop.build_W(coo_W, U, std::pair<int,BaseLink>{s, lk});
        coo_W(d_tmp_th, d_xi);
        CUDA_CHECK(cudaDeviceSynchronize());
        Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_chi2, dt, d_tmp_th, d_chi2);
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
          COO<N> coo_W;
          kop.build_W(coo_W, U, std::pair<int,Idx>{s, ix});
          coo_W(d_tmp_th, d_xi);
          CUDA_CHECK(cudaDeviceSynchronize());
          Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_chi2, dt, d_tmp_th, d_chi2);
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
    std::cout << "# Theta check: ||chi1 - chi2||_inf = " << max_err << std::endl;

    CUDA_CHECK(cudaFree(d_Theta_xi));
    CUDA_CHECK(cudaFree(d_chi1));
    CUDA_CHECK(cudaFree(d_chi2));
    CUDA_CHECK(cudaFree(d_tmp_th));
  }

  // ---- D_ov Theta check: [D_ov, Theta] xi = sum_l dtheta_l K^l xi (eq. 3.27 for overlap) ----
  {
    auto f_Dov = std::bind(&Fermion::mult_deviceAsyncLaunch, &D,
                           std::placeholders::_1, std::placeholders::_2);
    LinOpWrapper M_Dov(f_Dov);
    MatPoly op_Dov;
    op_Dov.push_back(cplx(1.0), {&M_Dov});

    std::cout << "\n# --- D_ov Theta check: [D_ov,Theta]xi = sum_l dtheta_l K^l xi ---" << std::endl;

    // sanity: Theta = c*I => [D_ov, Theta] = 0 and dtheta = 0
    {
      const double c = 2.5;
      std::vector<double> h_theta_c(Comp::N_SITES * Nt, c);
      CuC* d_sc1 = nullptr; CuC* d_sc2 = nullptr; CuC* d_stmp = nullptr; CuC* d_sth = nullptr;
      CUDA_CHECK(cudaMalloc(&d_sc1,  N*CD));
      CUDA_CHECK(cudaMalloc(&d_sc2,  N*CD));
      CUDA_CHECK(cudaMalloc(&d_stmp, N*CD));
      CUDA_CHECK(cudaMalloc(&d_sth,  N*CD));
      multTheta(d_sth, d_xi, h_theta_c);
      op_Dov.on_gpu<N>(d_sc1, d_sth);    CUDA_CHECK(cudaDeviceSynchronize());
      op_Dov.on_gpu<N>(d_stmp, d_xi);    CUDA_CHECK(cudaDeviceSynchronize());
      multTheta(d_stmp, d_stmp, h_theta_c);
      Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_sc1, -1.0, d_stmp, d_sc1);
      CUDA_CHECK(cudaDeviceSynchronize());
      CUDA_CHECK(cudaMemset(d_sc2, 0, N*CD));
      FermionVector h1, h2;
      CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(h1.field), d_sc1, N*CD, D2H));
      CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(h2.field), d_sc2, N*CD, D2H));
      double e1=0.0, e2=0.0;
      for(int i=0;i<N;i++) { e1=std::max(e1,std::abs(h1.field[i])); e2=std::max(e2,std::abs(h2.field[i])); }
      std::cout << "# Theta=c*I: ||chi1||_inf=" << e1 << "  ||chi2||_inf=" << e2 << std::endl;
      CUDA_CHECK(cudaFree(d_sc1)); CUDA_CHECK(cudaFree(d_sc2));
      CUDA_CHECK(cudaFree(d_stmp)); CUDA_CHECK(cudaFree(d_sth));
    }

    // localized: theta = 1 at ix=0 (all s), 0 elsewhere; temporal dtheta = 0
    {
      std::vector<double> h_thloc(Comp::N_SITES * Nt, 0.0);
      for(int s = 0; s < Nt; s++) h_thloc[s * Comp::N_SITES + 0] = 1.0;

      CuC* d_lth = nullptr; CuC* d_lc1 = nullptr;
      CuC* d_lc2 = nullptr; CuC* d_ltmp = nullptr;
      CUDA_CHECK(cudaMalloc(&d_lth,  N*CD));
      CUDA_CHECK(cudaMalloc(&d_lc1,  N*CD));
      CUDA_CHECK(cudaMalloc(&d_lc2,  N*CD));
      CUDA_CHECK(cudaMalloc(&d_ltmp, N*CD));

      multTheta(d_lth, d_xi, h_thloc);
      op_Dov.on_gpu<N>(d_lc1, d_lth);    CUDA_CHECK(cudaDeviceSynchronize());
      op_Dov.on_gpu<N>(d_ltmp, d_xi);    CUDA_CHECK(cudaDeviceSynchronize());
      multTheta(d_ltmp, d_ltmp, h_thloc);
      Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_lc1, -1.0, d_ltmp, d_lc1);
      CUDA_CHECK(cudaDeviceSynchronize());

      FermionVector h_lc1;
      CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(h_lc1.field), d_lc1, N*CD, D2H));
      std::cout << "# --- localized Theta (ix=0 all s): chi1=[D_ov,Theta]xi ---" << std::endl;
      for(int i = 0; i < N; i++)
        if(std::abs(h_lc1.field[i]) > 1.0e-15)
          std::cout << "#  chi1[" << i << "] = " << h_lc1.field[i] << std::endl;

      CUDA_CHECK(cudaMemset(d_lc2, 0, N*CD));
      for(const auto& lk : base.links) {
        if(lk[0] != 0 && lk[1] != 0) continue;
        const double dt = (lk[0] == 0) ? 1.0 : -1.0;
        for(int s = 0; s < Nt; s++) {
          kop.apply_k(d_ltmp, d_xi, U, std::pair<int,BaseLink>{s, lk});
          CUDA_CHECK(cudaDeviceSynchronize());
          Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_lc2, dt, d_ltmp, d_lc2);
          CUDA_CHECK(cudaDeviceSynchronize());
        }
      }
      FermionVector h_lc2;
      CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(h_lc2.field), d_lc2, N*CD, D2H));
      std::cout << "# --- chi2 = sum dtheta K^l xi (nonzero) ---" << std::endl;
      for(int i = 0; i < N; i++)
        if(std::abs(h_lc2.field[i]) > 1.0e-15)
          std::cout << "#  chi2[" << i << "] = " << h_lc2.field[i] << std::endl;
      double lerr = 0.0;
      for(int i = 0; i < N; i++) lerr = std::max(lerr, std::abs(h_lc1.field[i]-h_lc2.field[i]));
      std::cout << "# localized: ||chi1-chi2||_inf = " << lerr << std::endl;

      CUDA_CHECK(cudaFree(d_lth)); CUDA_CHECK(cudaFree(d_lc1));
      CUDA_CHECK(cudaFree(d_lc2)); CUDA_CHECK(cudaFree(d_ltmp));
    }

    // random theta
    std::vector<double> h_theta_ov(Comp::N_SITES * Nt);
    {
      std::normal_distribution<double> rdist(0.0, 1.0);
      for(auto& t : h_theta_ov) t = rdist(rng_host);
    }

    CuC* d_Theta_xi2 = nullptr;
    CuC* d_chi1_ov   = nullptr;
    CuC* d_chi2_ov   = nullptr;
    CuC* d_tmp_ov    = nullptr;
    CUDA_CHECK(cudaMalloc(&d_Theta_xi2, N*CD));
    CUDA_CHECK(cudaMalloc(&d_chi1_ov,   N*CD));
    CUDA_CHECK(cudaMalloc(&d_chi2_ov,   N*CD));
    CUDA_CHECK(cudaMalloc(&d_tmp_ov,    N*CD));

    multTheta(d_Theta_xi2, d_xi, h_theta_ov);
    op_Dov.on_gpu<N>(d_chi1_ov, d_Theta_xi2);
    CUDA_CHECK(cudaDeviceSynchronize());
    op_Dov.on_gpu<N>(d_tmp_ov, d_xi);
    CUDA_CHECK(cudaDeviceSynchronize());
    multTheta(d_tmp_ov, d_tmp_ov, h_theta_ov);
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_chi1_ov, -1.0, d_tmp_ov, d_chi1_ov);
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemset(d_chi2_ov, 0, N*CD));
    for(int s = 0; s < Nt; s++) {
      for(const auto& lk : base.links) {
        const double dt = h_theta_ov[s*Comp::N_SITES + lk[0]]
                        - h_theta_ov[s*Comp::N_SITES + lk[1]];
        kop.apply_k(d_tmp_ov, d_xi, U, std::pair<int,BaseLink>{s, lk});
        CUDA_CHECK(cudaDeviceSynchronize());
        Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_chi2_ov, dt, d_tmp_ov, d_chi2_ov);
        CUDA_CHECK(cudaDeviceSynchronize());
      }
    }
    if(Nt > 1) {
      for(int s = 0; s < Nt; s++) {
        const int s1 = (s + 1) % Nt;
        for(int ix = 0; ix < base.n_sites; ix++) {
          const double dt = h_theta_ov[s*Comp::N_SITES + ix]
                          - h_theta_ov[s1*Comp::N_SITES + ix];
          kop.apply_k(d_tmp_ov, d_xi, U, std::pair<int,Idx>{s, ix});
          CUDA_CHECK(cudaDeviceSynchronize());
          Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_chi2_ov, dt, d_tmp_ov, d_chi2_ov);
          CUDA_CHECK(cudaDeviceSynchronize());
        }
      }
    }

    FermionVector h_chi1_ov, h_chi2_ov;
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(h_chi1_ov.field), d_chi1_ov, N*CD, D2H));
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(h_chi2_ov.field), d_chi2_ov, N*CD, D2H));
    double max_err_ov = 0.0;
    for(int i = 0; i < N; i++)
      max_err_ov = std::max(max_err_ov, std::abs(h_chi1_ov.field[i] - h_chi2_ov.field[i]));
    std::cout << "# D_ov Theta check: ||chi1 - chi2||_inf = " << max_err_ov << std::endl;

    CUDA_CHECK(cudaFree(d_Theta_xi2));
    CUDA_CHECK(cudaFree(d_chi1_ov));
    CUDA_CHECK(cudaFree(d_chi2_ov));
    CUDA_CHECK(cudaFree(d_tmp_ov));
  }


  // ---- tr(K) checks (exact and stochastic) at w=(s=0,ix=0), looped over gauge backgrounds ----
  // tr(K^{wz}) = 0 analytically; stochastic estimator eta^dag K eta should also converge to 0.
  {
    using RngTrK = ParallelRngExt<Base, Comp::Nt>;
    RngTrK rng_gauge(base, 99);
    const Gauge U_free = U;

    auto run_trK = [&](const std::string& glabel) {
      constexpr int s_focus  = 0;
      constexpr int ix_focus = 0;

      std::vector<Idx> nns_list(base.nns[ix_focus].begin(), base.nns[ix_focus].end());
      const int n_spatial = static_cast<int>(nns_list.size());

      // analytic tr(K) values, filled by the exact block and used by the stochastic summary
      std::vector<Complex> trK_exact(n_spatial, Complex(0.0, 0.0));
      Complex trK_exact_t[2] = {Complex(0.0,0.0), Complex(0.0,0.0)};

      // ---- Exact tr(K): basis-vector loop (small lattice only) ----
      // tr(K^{wz}) = sum_k e_k^dag K^{wz} e_k; flat index k: s=k/Nx, ix=(k%Nx)/NS, spin=k%NS.
      if constexpr (Comp::N_REFINE == 1 && Comp::Nt <= 4) {
        std::cout << "\n# --- tr(K) exact (" << glabel << ") at (s=0, ix=0) ---" << std::endl;

        FermionVector e_i, fv_Ke;
        CuC* d_ei = nullptr;
        CUDA_CHECK(cudaMalloc(&d_ei, N*CD));

        auto trK_link = [&](auto el) -> Complex {
          Complex trK(0.0, 0.0);
          for(Idx k = 0; k < N; k++) {
            e_i.set_pt_source(k / Comp::Nx, (k % Comp::Nx) / NS, k % NS);
            CUDA_CHECK(cudaMemcpy(d_ei, reinterpret_cast<const CuC*>(e_i.field), N*CD, H2D));
            kop.apply_k(d_k_xi, d_ei, U, el);
            CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(fv_Ke.field), d_k_xi, N*CD, D2H));
            trK += fv_Ke.field[k];
          }
          return trK;
        };

        for(int j = 0; j < n_spatial; j++) {
          trK_exact[j] = trK_link(std::pair<int,BaseLink>{s_focus, BaseLink{ix_focus, nns_list[j]}});
          std::cout << "# trK  nn" << nns_list[j]
                    << " = (" << trK_exact[j].real() << "," << trK_exact[j].imag() << ")" << std::endl;
        }

        if(Nt > 1) {
          trK_exact_t[0] = trK_link(std::pair<int,Idx>{s_focus, ix_focus});
          trK_exact_t[1] = trK_link(std::pair<int,Idx>{(s_focus-1+Nt)%Nt, ix_focus});
          std::cout << "# trK  tlink_fwd"
                    << " = (" << trK_exact_t[0].real() << "," << trK_exact_t[0].imag() << ")" << std::endl;
          std::cout << "# trK  tlink_bwd"
                    << " = (" << trK_exact_t[1].real() << "," << trK_exact_t[1].imag() << ")" << std::endl;
        }

        CUDA_CHECK(cudaFree(d_ei));
      }

      // ---- Stochastic tr(K): Z2xZ2 noise ----
      // tr(K^{wz}) ~ (1/nhits) sum_h eta_h^dag K^{wz} eta_h; compared to exact value above.
      {
        std::cout << "\n# --- stochastic tr(K) (" << glabel << ") at (s=0, ix=0) ---" << std::endl;

        RngTrK rng_trk(base, 42);

        FermionVector eta, fv_Keta;
        CuC* d_eta = nullptr;
        CUDA_CHECK(cudaMalloc(&d_eta, N*CD));

        std::vector<Complex> acc(n_spatial, Complex(0.0, 0.0));
        std::vector<std::vector<Complex>> all_hits(n_spatial);
        // temporal links: index 0 = fwd (s_focus, ix_focus), 1 = bwd ((s_focus-1+Nt)%Nt, ix_focus)
        Complex acc_t[2] = {Complex(0.0,0.0), Complex(0.0,0.0)};
        std::vector<Complex> all_hits_t[2];

        for(int h = 0; h < nhits; h++) {
          eta.fill_z2_source(rng_trk);
          CUDA_CHECK(cudaMemcpy(d_eta, reinterpret_cast<const CuC*>(eta.field), N*CD, H2D));
          for(int j = 0; j < n_spatial; j++) {
            const BaseLink lk_j{ix_focus, nns_list[j]};
            kop.apply_k(d_k_xi, d_eta, U, std::pair<int,BaseLink>{s_focus, lk_j});
            CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(fv_Keta.field), d_k_xi, N*CD, D2H));
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
              kop.apply_k(d_k_xi, d_eta, U, tlinks[t]);
              CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(fv_Keta.field), d_k_xi, N*CD, D2H));
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

        std::cout << "# stoch trK mean +/- stderr  vs  exact:" << std::endl;
        for(int j = 0; j < n_spatial; j++)
          print_stderr("nn=" + std::to_string(nns_list[j]), acc[j], all_hits[j], trK_exact[j]);
        if(Nt > 1) {
          print_stderr("tlink_fwd", acc_t[0], all_hits_t[0], trK_exact_t[0]);
          print_stderr("tlink_bwd", acc_t[1], all_hits_t[1], trK_exact_t[1]);
        }

        CUDA_CHECK(cudaFree(d_eta));
      }
    }; // end run_trK

    run_trK("free");

    U.gaussian(rng_gauge, 1.0);
    D.update(U);
    run_trK("gaussian");

    U = U_free;
    D.update(U);
  }


  // ---- Chunk 5b: stochastic Ward identity at w=(s=0,ix=0), free-field U ----
  // <J_V^{wz}> = -i Nf Im[ eta^dag K^{wz} phi ], D phi = eta, Z2xZ2 noise.
  // Ward identity: sum_{z~w} <J_V^{wz}> = 0 (eq. 3.32); purely imaginary, sum is 0.
  // Code reports coefficient of -i: +Nf * mean( Im[eta^dag K phi] ).
  // Forward solve D phi=eta via normal eqs: D^dag D phi = D^dag eta.
  {
    std::cout << "\n# --- Chunk 5b: stochastic Ward identity at (s=0, ix=0) ---" << std::endl;

    const int Nhits   = nhits;
    const int Nf_val  = 2;
    constexpr int s_focus  = 0;
    constexpr int ix_focus = 0;

    std::cout << "# base.nns[0]:";
    for(Idx iy : base.nns[0]) std::cout << " " << iy;
    std::cout << std::endl;

    using Rng5 = ParallelRngExt<Base, Comp::Nt>;
    Rng5 rng_stoch(base, 42);

    auto f_DH  = std::bind(&Fermion::adj_deviceAsyncLaunch,  &D, std::placeholders::_1, std::placeholders::_2);
    auto f_DHD = std::bind(&Fermion::DHD_deviceAsyncLaunch,  &D, std::placeholders::_1, std::placeholders::_2);
    LinOpWrapper M_DH(f_DH), M_DHD(f_DHD);
    MatPoly op_DH;  op_DH.push_back(cplx(1.0), {&M_DH});
    MatPoly op_DHD; op_DHD.push_back(cplx(1.0), {&M_DHD});

    FermionVector eta, DH_eta, phi;
    CuC* d_phi_dev = nullptr;
    CUDA_CHECK(cudaMalloc(&d_phi_dev, N*CD));

    std::vector<Idx> nns_list(base.nns[ix_focus].begin(), base.nns[ix_focus].end());
    const int n_spatial = static_cast<int>(nns_list.size());

    // ---- Analytic <J_V^{wz}>: basis-vector loop (small lattice only) ----
    // J_V^{wz} = -i Nf Im[ tr(K^{wz} D^{-1}) ] = -i Nf sum_k Im( (K^{wz} phi_k)[k] ), D phi_k = e_k.
    // JV_exact stores the coefficient of -i: +Nf * sum_k Im( (K phi_k)[k] ).
    // One solve per basis vector; all links share the same set of phi_k.
    std::vector<Complex> JV_exact(n_spatial, Complex(0.0,0.0));
    Complex JV_exact_t[2] = {Complex(0.0,0.0), Complex(0.0,0.0)};
    Complex div_exact = Complex(0.0,0.0);
    if constexpr (Comp::N_REFINE == 1 && Comp::Nt <= 4) {
      std::cout << "\n# --- analytic <J_V> via basis-vector loop ---" << std::endl;
      FermionVector e_i, DH_e, phi_e, fv_Ke;
      for(Idx k = 0; k < N; k++) {
        e_i.set_pt_source(k / Comp::Nx, (k % Comp::Nx) / NS, k % NS);
        op_DH.from_cpu<N>(DH_e.field, e_i.field);
        op_DHD.solve<N>(phi_e.field, DH_e.field, Comp::TOL_OUTER);
        CUDA_CHECK(cudaMemcpy(d_phi_dev, reinterpret_cast<const CuC*>(phi_e.field), N*CD, H2D));
        for(int j = 0; j < n_spatial; j++) {
          const BaseLink lk_j{ix_focus, nns_list[j]};
          kop.apply_k(d_k_xi, d_phi_dev, U, std::pair<int,BaseLink>{s_focus, lk_j});
          CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(fv_Ke.field), d_k_xi, N*CD, D2H));
          JV_exact[j] += fv_Ke.field[k];
        }
        if(Nt > 1) {
          const std::pair<int,Idx> tlinks[2] = {
            {s_focus, ix_focus},
            {(s_focus-1+Nt)%Nt, ix_focus}
          };
          const double tsign[2] = {1.0, -1.0};
          for(int t = 0; t < 2; t++) {
            kop.apply_k(d_k_xi, d_phi_dev, U, tlinks[t]);
            CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(fv_Ke.field), d_k_xi, N*CD, D2H));
            JV_exact_t[t] += tsign[t] * fv_Ke.field[k];
          }
        }
        if((k+1) % (N/4) == 0)
          std::cout << "# basis loop: k=" << (k+1) << "/" << N << std::endl;
      }
      for(int j = 0; j < n_spatial; j++) {
        JV_exact[j] *= -(double)Nf_val;
        div_exact += JV_exact[j];
        std::cout << "# J_{V,+} exact  nn=" << nns_list[j]
                  << "  re=" << JV_exact[j].real() << "  im=" << JV_exact[j].imag() << std::endl;
      }
      if(Nt > 1) {
        for(int t = 0; t < 2; t++) JV_exact_t[t] *= -(double)Nf_val;
        div_exact += JV_exact_t[0] + JV_exact_t[1];
        std::cout << "# J_{V,+} exact  tlink_fwd"
                  << "  re=" << JV_exact_t[0].real() << "  im=" << JV_exact_t[0].imag() << std::endl;
        std::cout << "# J_{V,+} exact  tlink_bwd"
                  << "  re=" << JV_exact_t[1].real() << "  im=" << JV_exact_t[1].imag() << std::endl;
      }
      std::cout << "# J_{V,+} exact  SUM"
                << "  re=" << div_exact.real() << "  im=" << div_exact.imag() << std::endl;
    }

    std::vector<Complex> C_acc(n_spatial, Complex(0.0,0.0));
    std::vector<std::vector<Complex>> C_all(n_spatial);
    // temporal: index 0=fwd (sign +1), 1=bwd (sign -1)
    Complex C_t_acc[2] = {Complex(0.0,0.0), Complex(0.0,0.0)};
    std::vector<Complex> C_t_all[2];
    Complex div_acc = Complex(0.0,0.0);
    std::vector<Complex> div_all;

    FermionVector fv_Kphi;

    for(int h = 0; h < Nhits; h++) {
      eta.fill_z2_source(rng_stoch);
      op_DH.from_cpu<N>(DH_eta.field, eta.field);
      op_DHD.solve<N>(phi.field, DH_eta.field, Comp::TOL_OUTER);
      CUDA_CHECK(cudaMemcpy(d_phi_dev, reinterpret_cast<const CuC*>(phi.field), N*CD, H2D));

      Complex div_h = Complex(0.0,0.0);
      for(int j = 0; j < n_spatial; j++) {
        const BaseLink lk_j{ix_focus, nns_list[j]};
        kop.apply_k(d_k_xi, d_phi_dev, U, std::pair<int,BaseLink>{s_focus, lk_j});
        CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(fv_Kphi.field), d_k_xi, N*CD, D2H));
        const Complex c = eta.dag(fv_Kphi);
        C_acc[j] += c;
        C_all[j].push_back(c);
        div_h += c;
      }
      if(Nt > 1) {
        const std::pair<int,Idx> tlinks[2] = {
          {s_focus, ix_focus},
          {(s_focus-1+Nt)%Nt, ix_focus}
        };
        const double tsign[2] = {1.0, -1.0};
        for(int t = 0; t < 2; t++) {
          kop.apply_k(d_k_xi, d_phi_dev, U, tlinks[t]);
          CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(fv_Kphi.field), d_k_xi, N*CD, D2H));
          const Complex c = eta.dag(fv_Kphi);
          C_t_acc[t] += tsign[t] * c;
          C_t_all[t].push_back(tsign[t] * c);
          div_h += tsign[t] * c;
        }
      }
      div_acc += div_h;
      div_all.push_back(div_h);

      const double norm = 1.0 / (h+1);
      std::cout << "# hit " << std::setw(4) << (h+1)
                << "  div_im=" << ((double)Nf_val)*div_acc.imag()*norm
                << "  div_re=" << ((double)Nf_val)*div_acc.real()*norm;
      for(int j = 0; j < n_spatial; j++)
        std::cout << "  nn" << nns_list[j]
                  << "_im=" << C_acc[j].imag()*norm
                  << "_re=" << C_acc[j].real()*norm;
      if(Nt > 1)
        std::cout << "  tfwd_im=" << C_t_acc[0].imag()*norm << "_re=" << C_t_acc[0].real()*norm
                  << "  tbwd_im=" << C_t_acc[1].imag()*norm << "_re=" << C_t_acc[1].real()*norm;
      std::cout << std::endl;
    }

    auto print_ward = [&](const std::string& label, const std::vector<Complex>& samples,
                          Complex exact) {
      double mean_im = 0.0, mean_re = 0.0;
      for(const Complex& x : samples) { mean_im += x.imag(); mean_re += x.real(); }
      mean_im /= Nhits; mean_re /= Nhits;
      const double stoch_re = -((double)Nf_val) * mean_re;
      const double stoch_im = -((double)Nf_val) * mean_im;
      double var_im = 0.0, var_re = 0.0;
      for(const Complex& x : samples) {
        var_im += std::pow(x.imag() - mean_im, 2);
        var_re += std::pow(x.real() - mean_re, 2);
      }
      const double n1 = static_cast<double>(Nhits - 1) * Nhits;
      const double err_re = ((double)Nf_val) * std::sqrt(var_re / n1);
      const double err_im = ((double)Nf_val) * std::sqrt(var_im / n1);
      std::cout << "#   J_{V,+}  " << label;
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

    std::cout << "# per-link <J_V> stoch vs exact (Nhits=" << Nhits << ", Nf=" << Nf_val << "):" << std::endl;
    for(int j = 0; j < n_spatial; j++)
      print_ward("nn=" + std::to_string(nns_list[j]), C_all[j], JV_exact[j]);
    if(Nt > 1) {
      print_ward("tlink_fwd", C_t_all[0], JV_exact_t[0]);
      print_ward("tlink_bwd", C_t_all[1], JV_exact_t[1]);
    }
    print_ward("SUM", div_all, div_exact);

    CUDA_CHECK(cudaFree(d_phi_dev));
  }

  // ---- Chunk 6: stochastic Ward identity on Gaussian gauge config ----
  // Same estimator as Chunk 5b; gauge background is a single Gaussian-random U.
  {
    std::cout << "\n# --- Chunk 6: Ward identity on Gaussian gauge config ---" << std::endl;

    using Rng6 = ParallelRngExt<Base, Comp::Nt>;
    Rng6 rng6_gauge(base, 99);
    U.gaussian(rng6_gauge, 2.0);
    D.update(U);
    std::cout << "# Gaussian U generated (width=2.0), D updated." << std::endl;

    const int Nhits6       = nhits;
    constexpr int s_focus6 = 0;
    constexpr int ix_focus6 = 0;

    auto f_DH6  = std::bind(&Fermion::adj_deviceAsyncLaunch,  &D, std::placeholders::_1, std::placeholders::_2);
    auto f_DHD6 = std::bind(&Fermion::DHD_deviceAsyncLaunch,  &D, std::placeholders::_1, std::placeholders::_2);
    LinOpWrapper M_DH6(f_DH6), M_DHD6(f_DHD6);
    MatPoly op_DH6;  op_DH6.push_back(cplx(1.0), {&M_DH6});
    MatPoly op_DHD6; op_DHD6.push_back(cplx(1.0), {&M_DHD6});

    Rng6 rng6_stoch(base, 42);
    FermionVector eta6, DH_eta6, phi6, fv_Kphi6;
    CuC* d_phi6 = nullptr;
    CUDA_CHECK(cudaMalloc(&d_phi6, N*CD));

    const std::vector<Idx> nns6(base.nns[ix_focus6].begin(), base.nns[ix_focus6].end());
    const int n_sp6 = static_cast<int>(nns6.size());

    // ---- Analytic <J_V^{wz}> on Gaussian U: basis-vector loop (small lattice only) ----
    std::vector<Complex> JV_exact6(n_sp6, Complex(0.0,0.0));
    Complex JV_exact6_t[2] = {Complex(0.0,0.0), Complex(0.0,0.0)};
    Complex div_exact6 = Complex(0.0,0.0);
    if constexpr (Comp::N_REFINE == 1 && Comp::Nt <= 4) {
      std::cout << "\n# --- analytic <J_V> (Gaussian U) via basis-vector loop ---" << std::endl;
      FermionVector e_i6, DH_e6, phi_e6, fv_Ke6;
      for(Idx k = 0; k < N; k++) {
        e_i6.set_pt_source(k / Comp::Nx, (k % Comp::Nx) / NS, k % NS);
        op_DH6.from_cpu<N>(DH_e6.field, e_i6.field);
        op_DHD6.solve<N>(phi_e6.field, DH_e6.field, Comp::TOL_OUTER);
        CUDA_CHECK(cudaMemcpy(d_phi6, reinterpret_cast<const CuC*>(phi_e6.field), N*CD, H2D));
        for(int j = 0; j < n_sp6; j++) {
          const BaseLink lk_j{ix_focus6, nns6[j]};
          kop.apply_k(d_k_xi, d_phi6, U, std::pair<int,BaseLink>{s_focus6, lk_j});
          CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(fv_Ke6.field), d_k_xi, N*CD, D2H));
          JV_exact6[j] += fv_Ke6.field[k];
        }
        if(Nt > 1) {
          const std::pair<int,Idx> tlinks6_e[2] = {
            {s_focus6, ix_focus6},
            {(s_focus6-1+Nt)%Nt, ix_focus6}
          };
          const double tsign6_e[2] = {1.0, -1.0};
          for(int t = 0; t < 2; t++) {
            kop.apply_k(d_k_xi, d_phi6, U, tlinks6_e[t]);
            CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(fv_Ke6.field), d_k_xi, N*CD, D2H));
            JV_exact6_t[t] += tsign6_e[t] * fv_Ke6.field[k];
          }
        }
        if((k+1) % (N/4) == 0)
          std::cout << "# basis loop: k=" << (k+1) << "/" << N << std::endl;
      }
      for(int j = 0; j < n_sp6; j++) {
        JV_exact6[j] *= -(double)Nf;
        div_exact6 += JV_exact6[j];
        std::cout << "# J_{V,+} exact  nn=" << nns6[j]
                  << "  re=" << JV_exact6[j].real() << "  im=" << JV_exact6[j].imag() << std::endl;
      }
      if(Nt > 1) {
        for(int t = 0; t < 2; t++) JV_exact6_t[t] *= -(double)Nf;
        div_exact6 += JV_exact6_t[0] + JV_exact6_t[1];
        std::cout << "# J_{V,+} exact  tlink_fwd"
                  << "  re=" << JV_exact6_t[0].real() << "  im=" << JV_exact6_t[0].imag() << std::endl;
        std::cout << "# J_{V,+} exact  tlink_bwd"
                  << "  re=" << JV_exact6_t[1].real() << "  im=" << JV_exact6_t[1].imag() << std::endl;
      }
      std::cout << "# J_{V,+} exact  SUM"
                << "  re=" << div_exact6.real() << "  im=" << div_exact6.imag() << std::endl;
    }

    std::vector<Complex> C_acc6(n_sp6, Complex(0.0,0.0));
    std::vector<std::vector<Complex>> C_all6(n_sp6);
    Complex C_t_acc6[2] = {Complex(0.0,0.0), Complex(0.0,0.0)};
    std::vector<Complex> C_t_all6[2];
    Complex div_acc6 = Complex(0.0,0.0);
    std::vector<Complex> div_all6;

    for(int h = 0; h < Nhits6; h++) {
      eta6.fill_z2_source(rng6_stoch);
      op_DH6.from_cpu<N>(DH_eta6.field, eta6.field);
      op_DHD6.solve<N>(phi6.field, DH_eta6.field, Comp::TOL_OUTER);
      CUDA_CHECK(cudaMemcpy(d_phi6, reinterpret_cast<const CuC*>(phi6.field), N*CD, H2D));

      Complex div_h6 = Complex(0.0,0.0);
      for(int j = 0; j < n_sp6; j++) {
        const BaseLink lk_j{ix_focus6, nns6[j]};
        kop.apply_k(d_k_xi, d_phi6, U, std::pair<int,BaseLink>{s_focus6, lk_j});
        CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(fv_Kphi6.field), d_k_xi, N*CD, D2H));
        const Complex c = eta6.dag(fv_Kphi6);
        C_acc6[j] += c; C_all6[j].push_back(c);
        div_h6 += c;
      }
      if(Nt > 1) {
        const std::pair<int,Idx> tlinks6[2] = {
          {s_focus6, ix_focus6},
          {(s_focus6-1+Nt)%Nt, ix_focus6}
        };
        const double tsign6[2] = {1.0, -1.0};
        for(int t = 0; t < 2; t++) {
          kop.apply_k(d_k_xi, d_phi6, U, tlinks6[t]);
          CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(fv_Kphi6.field), d_k_xi, N*CD, D2H));
          const Complex c = eta6.dag(fv_Kphi6);
          C_t_acc6[t] += tsign6[t] * c; C_t_all6[t].push_back(tsign6[t] * c);
          div_h6 += tsign6[t] * c;
        }
      }
      div_acc6 += div_h6;
      div_all6.push_back(div_h6);

      const double norm6 = 1.0 / (h+1);
      std::cout << "# hit " << std::setw(4) << (h+1)
                << "  div_im=" << ((double)Nf)*div_acc6.imag()*norm6
                << "  div_re=" << ((double)Nf)*div_acc6.real()*norm6;
      for(int j = 0; j < n_sp6; j++)
        std::cout << "  nn" << nns6[j]
                  << "_im=" << C_acc6[j].imag()*norm6
                  << "_re=" << C_acc6[j].real()*norm6;
      if(Nt > 1)
        std::cout << "  tfwd_im=" << C_t_acc6[0].imag()*norm6 << "_re=" << C_t_acc6[0].real()*norm6
                  << "  tbwd_im=" << C_t_acc6[1].imag()*norm6 << "_re=" << C_t_acc6[1].real()*norm6;
      std::cout << std::endl;
    }

    auto print_ward6 = [&](const std::string& label, const std::vector<Complex>& samples,
                           Complex exact) {
      double mean_im = 0.0, mean_re = 0.0;
      for(const Complex& x : samples) { mean_im += x.imag(); mean_re += x.real(); }
      mean_im /= Nhits6; mean_re /= Nhits6;
      const double stoch_re = -((double)Nf) * mean_re;
      const double stoch_im = -((double)Nf) * mean_im;
      double var_im = 0.0, var_re = 0.0;
      for(const Complex& x : samples) {
        var_im += std::pow(x.imag() - mean_im, 2);
        var_re += std::pow(x.real() - mean_re, 2);
      }
      const double n1 = static_cast<double>(Nhits6 - 1) * Nhits6;
      const double err_re = ((double)Nf) * std::sqrt(var_re / n1);
      const double err_im = ((double)Nf) * std::sqrt(var_im / n1);
      std::cout << "#   J_{V,+}  " << label;
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

    std::cout << "# per-link <J_V> (Nhits=" << Nhits6 << ", Nf=" << Nf << ", Gaussian U):" << std::endl;
    for(int j = 0; j < n_sp6; j++)
      print_ward6("nn=" + std::to_string(nns6[j]), C_all6[j], JV_exact6[j]);
    if(Nt > 1) {
      print_ward6("tlink_fwd", C_t_all6[0], JV_exact6_t[0]);
      print_ward6("tlink_bwd", C_t_all6[1], JV_exact6_t[1]);
    }
    print_ward6("SUM", div_all6, div_exact6);

    CUDA_CHECK(cudaFree(d_phi6));
  }

  // --- cleanup ---
  CUDA_CHECK(cudaFree(d_xi));
  CUDA_CHECK(cudaFree(d_k_xi));
  for(int i = 0; i < Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();

  return 0;
}
