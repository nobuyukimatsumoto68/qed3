// jj_kbuild_exact_claude.cu
// -----------------------------------------------------------------------------
// STAGE 1b of the exact current-correlator validation (jj_exact_freefield_impl_plan_claude.md):
// build & save the DENSE projected conserved current K^P(t) (N x N), per gauge config.
//
// K is MASS-INDEPENDENT, so it lives in its OWN (mass-free) directory and is built ONCE per config/L,
// then reused by jj_contract_exact_claude.cu for every mass.  This program holds all the overlap /
// multishift GPU work for the current side; the contraction stage is then pure linear algebra.
//
//   K^P(t) e_j = sum_a w^P_a K(a,t) e_j      (P=tp: a over sites, w_tp; P=sp: a over links, w_sp)
// built column-by-column with the SAME ConservedCurrent kop / op_K + weights as the stochastic
// estimator (apply_k_ms multishift).  free field => only t=0 needed (the contraction time-shifts it);
// interacting => all t=0..Nt-1 (no translation invariance).
//
// L (= N_REFINE) COMPILE-TIME: -DN_REFINE_CLI=<L> (default 1).  CLI: --ens-dir (omit => free, t=0
// only), --ninter, --gpu.   Output: data_<free|ensbase>_Kdense_L<L>/K.<k>.h5  (datasets <proj>_t<t>
// row-major N^2; atomic write).   No mass argument (K is mass-independent).
// -----------------------------------------------------------------------------

#include <typeinfo>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <highfive/H5File.hpp>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <filesystem>
#include <memory>
#include <chrono>
#include <cstdint>
#include <complex>
#include <array>
#include <vector>
#include <map>
#include <string>
#include <ctime>
#include <getopt.h>
#include <Eigen/Dense>

using Double = double;
using Idx = std::int32_t;
using Complex = std::complex<double>;
using Link = std::array<Idx,2>;
using Face = std::vector<Idx>;

using MS = Eigen::Matrix2cd;
using VD = Eigen::Vector2d;
using VE = Eigen::Vector3d;
using VC = Eigen::VectorXcd;

static constexpr int NS = 2;
static constexpr int DIM = 2;
static constexpr Complex I = Complex(0.0, 1.0);

namespace Comp{
  constexpr bool is_compact = false;
  constexpr int NPARALLEL_DUPDATE = 1;
  constexpr int NPARALLEL = NPARALLEL_DUPDATE;
  constexpr int NSTREAMS = 4;
  constexpr int NPARALLEL_GAUGE = NPARALLEL_DUPDATE;
  constexpr int NPARALLEL_SORT = NPARALLEL_DUPDATE;
#ifndef N_REFINE_CLI
#define N_REFINE_CLI 1
#endif
  constexpr int N_REFINE = N_REFINE_CLI;
  constexpr int NS = 2;
#ifndef NT_CLI
#define NT_CLI 128                     // temporal extent; override: nvcc -DNT_CLI=16 ... (small-Nt checks)
#endif
  constexpr int Nt = NT_CLI;
  constexpr Idx N_SITES = 10*N_REFINE*N_REFINE + 2;
  constexpr int N_LINKS = 30*N_REFINE*N_REFINE;
  constexpr Idx Nx = NS*N_SITES;
  constexpr Idx N = Nx*Nt;
  const double TOL_INNER = 1.0e-9;
  const double TOL_OUTER = 1.0e-5;
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
#include "matpoly_claude.h"
#include "dirac_pf.h"
#include "overlap_wmass_claude.h"
#include "conserved_current_claude.h"

using BaseLink = std::array<Idx,2>;
using BaseFace = std::vector<Idx>;


// ---- CLI ----
struct Args {
  double nu0 = 1.0, nu1 = -1.0;
  std::string ens_dir;
  int ninter = 10, gpu = 0;
  int tp_only_t = -1;       // >=0: build ONLY proj=tp at ONLY this slice t (cheap single-slice check)
};

void PrintHelp(){
  printf("Usage: ./a.out [options]   (build dense K^P(t), mass-INDEPENDENT)\n");
  printf("  --ens-dir <path>     sea config dir; OMIT => free field (U=1, build t=0 only)\n");
  printf("  --ninter <N>         ensemble stride ckpoint_lat.k (default 10)\n");
  printf("  --nu0 <x> --nu1 <y>  (default 1.0 ; nu1 defaults to nu0)\n");
  printf("  --tp-only-t <t>      build ONLY tp at ONLY slice t (single-slice cross-check)\n");
  printf("  --gpu <id>           CUDA device (default 0; via CUDA_VISIBLE_DEVICES)\n");
  printf("  -h, --help\n");
}

void ParseArgs(int argc, char** argv, Args& a){
  static struct option lo[] = {
    {"nu0",       required_argument, 0, 'n'},
    {"nu1",       required_argument, 0, 'm'},
    {"ens-dir",   required_argument, 0, 'e'},
    {"ninter",    required_argument, 0, 'I'},
    {"tp-only-t", required_argument, 0, 'T'},
    {"gpu",       required_argument, 0, 'G'},
    {"help",      no_argument,       0, 'h'},
    {0,0,0,0}
  };
  int opt, idx;
  while((opt = getopt_long(argc, argv, "n:m:e:I:T:G:h", lo, &idx)) != -1){
    switch(opt){
      case 'n': a.nu0    = std::stod(optarg); break;
      case 'm': a.nu1    = std::stod(optarg); break;
      case 'e': a.ens_dir = optarg;           break;
      case 'I': a.ninter = std::stoi(optarg); break;
      case 'T': a.tp_only_t = std::stoi(optarg); break;
      case 'G': a.gpu    = std::stoi(optarg); break;
      case 'h': default: PrintHelp(); std::exit(0);
    }
  }
}


int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  Args a;
  ParseArgs(argc, argv, a);
  if(a.nu1 < 0.0) a.nu1 = a.nu0;
  (void)a.gpu;
  CUDA_CHECK(cudaSetDevice(0));
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  std::cout << "# dev = " << prop.name << std::endl;

  constexpr Idx N  = Comp::N;
  constexpr int Nt = Comp::Nt;
  const bool free_field = a.ens_dir.empty();
  std::cout << "# K-build: N=" << N
            << (free_field ? "  [free: t=0 only]" : "  [interacting: all t]") << std::endl;

  using Base = S2Simp;
  using WilsonDirac = DiracExt<Base, DiracS2Simp>;
  using Gauge = GaugeExt<Base, Comp::Nt, Comp::is_compact>;
  using Fermion = OverlapWMass<WilsonDirac>;

  Base base(Comp::N_REFINE);
  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();
  const double r = 1.0, M5 = -1.0, at = 0.2;
  WilsonDirac DW(base, 0.0, r, M5, at, a.nu1);
  Gauge U(base);
  Fermion D(DW, Complex(0.0), 21);                 // massless D_ov; K is mass-independent
  ConservedCurrent<Fermion,Gauge> kop(D);
  MatPoly op_K;
  op_K.push_back(cplx(1.0), {&kop});

  // ---- projection weights (Eq. 4.32 tp ; Eq. 4.29 sp) ----
  const int n_sites = (int)base.n_sites, n_links = (int)base.n_links;
  std::vector<double> w_tp(n_sites), w_sp(n_links);
  for(int n=0; n<n_sites; n++){ const double kt = DW.kappa_t[n];     w_tp[n] = base.dual_areas[n]/(kt*kt); }
  for(int il=0; il<n_links; il++){ const double ks = DW.bd.kappa[il]; w_sp[il] = base.link_volume[il]/(ks*ks); }

  // ---- output dir: MASS-FREE (K is mass-independent => shared across masses) ----
  std::string ens_base = a.ens_dir;
  if(!ens_base.empty() && ens_base.back()=='/') ens_base.pop_back();
  { const auto s = ens_base.find_last_of('/'); if(s!=std::string::npos) ens_base = ens_base.substr(s+1); }
  const std::string kdir = "data_" + (free_field ? std::string("free") : ens_base)
                         + "_Kdense_L" + std::to_string(Comp::N_REFINE) + "/";
  std::filesystem::create_directories(kdir);
  std::cout << "# kdir = " << kdir << std::endl;

  std::vector<Complex> ej(N), out(N), K0;
  Timer timer;

  // ---- config loop ----
  const int k_lo   = 0;
  const int k_step = free_field ? 1 : a.ninter;
  const int k_hi   = free_field ? 1 : 1000000;
  for(int k=k_lo; k<k_hi; k+=k_step){
    if(!free_field){
      const std::string lat = a.ens_dir + "ckpoint_lat." + std::to_string(k);
      if(!std::filesystem::exists(lat)){ if(k==0) continue; else break; }
      U.read(lat);
    }
    const std::string h5path = kdir + "K." + std::to_string(k) + ".h5";
    if(std::filesystem::exists(h5path)){
      bool c=false;
      try { HighFive::File f(h5path, HighFive::File::ReadOnly); c = f.exist("complete"); } catch(...) {}
      if(c){ std::cout << "# skip k=" << k << " (complete)\n"; if(free_field) break; else continue; }
    }
    D.update(U);
    std::cout << "# k=" << k << (free_field ? " (free)" : "") << "  building K^{tp,sp}(t) ...\n";

    // ATOMIC output: write to .tmp, rename after the 'complete' sentinel + close.
    const std::string h5tmp = h5path + ".tmp";
    auto h5p = std::make_unique<HighFive::File>(h5tmp,
                 HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
    HighFive::File& h5 = *h5p;
    h5.createDataSet("N",    std::vector<int>{(int)N});
    h5.createDataSet("Nt",   std::vector<int>{Nt});
    h5.createDataSet("free", std::vector<int>{free_field ? 1 : 0});

    const int tmax = free_field ? 1 : Nt;          // free: t=0 only; interacting: all t
    // single-slice cross-check: build ONLY tp at ONLY slice a.tp_only_t
    std::vector<std::string> projs = {std::string("tp"), std::string("sp")};
    if(a.tp_only_t >= 0) projs = {std::string("tp")};
    for(const std::string proj : projs){
      const int nins = (proj=="tp") ? n_sites : n_links;
      const int t_lo = (a.tp_only_t >= 0) ? a.tp_only_t   : 0;
      const int t_hi = (a.tp_only_t >= 0) ? a.tp_only_t+1 : tmax;
      for(int t=t_lo; t<t_hi; t++){
        // build dense K^proj(t) column-by-column: K0[i*N+j] = (K^proj(t))_{i,j}
        K0.assign((size_t)N*N, Complex(0,0));
        for(Idx j=0; j<N; j++){
          std::fill(ej.begin(), ej.end(), Complex(0,0));
          ej[j] = Complex(1,0);
          for(int ins=0; ins<nins; ins++){
            if(proj=="tp") kop.set_temporal(U, t, (Idx)ins, /*dag=*/false);
            else           kop.set_spatial (U, t, base.links[ins], /*dag=*/false);
            op_K.from_cpu<N>(out.data(), ej.data());            // out = K(ins,t) e_j
            const double w = (proj=="tp") ? w_tp[ins] : w_sp[ins];
            for(Idx i=0; i<N; i++) K0[(size_t)i*N + j] += w*out[i];
          }
        }
        std::vector<double> re(K0.size()), im(K0.size());
        for(size_t i=0; i<K0.size(); i++){ re[i] = K0[i].real(); im[i] = K0[i].imag(); }
        const std::string key = proj + "_t" + std::to_string(t);
        h5.createDataSet(key + "/real", re);
        h5.createDataSet(key + "/imag", im);
        std::cout << "#   wrote " << key << " [" << timer.currentSeconds() << " s]\n";
      }
    }

    h5.createDataSet("complete", std::vector<int>{1});
    h5p.reset();
    std::filesystem::rename(h5tmp, h5path);
    std::cout << "# wrote " << h5path << "\n";
    if(free_field) break;
  }

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();
  return 0;
}
