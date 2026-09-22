// distill_peram_claude.cu
// EXACT (deterministic) distillation for the sigma\sigma - F^2 (0++) mixing program (Chester-Pufu,
// arXiv:1603.05582).  Method: exact distillation (Peardon 0905.2160); NO stochastic estimators.
// Plan: distillation_impl_plan_claude.md ; method note: distillation_for_two_meson_claude.md.
//
// CHUNK 1 (this stage): the DISTILLATION BASIS.  Per (config, timeslice t) build the covariant TIMESLICE
// WILSON operator D_{W,2}(t) (= the 2D S^2 Wilson-Dirac DiracS2Simp at that timeslice's gauge field; it
// carries the spatial covariant hop kappa e^{i\theta} AND the spin connection \Omega), form the Hermitian
// PSD normal operator D_{W,2}^dag D_{W,2}, and diagonalize it densely (Eigen SelfAdjointEigenSolver).  The
// lowest N_v eigenvectors w_k (plain-orthonormal spinor modes on the 2 N_s space) are the basis V(t).
// Spin is FOLDED INTO the mode index (no separate spin index downstream).  At N_v = 2 N_s the basis is the
// FULL space => V V^dag = 1 (exact local operator; L1: N_v = 24).
//
// Basis operator = timeslice WILSON (cheap, ultralocal), NOT overlap (overlap restricted to a slice is
// time-nonlocal / costly).  The perambulator (chunk 2) STILL inverts the massless OVERLAP D_ov^{-1}; only
// the smearing basis is Wilson low modes.  Convention (D) of the plan; M5 = -1 matches the overlap kernel.
//
// The stored COO value IS the true matrix element (sparse_dirac_claude.h:184 coo2csr does plain assignment;
// the "-tmp" comments in DiracS2Simp::coo_format are stale), so the dense assembly below reproduces exactly
// the operator the code applies.
//
// CHUNK 2 (this stage too): the PERAMBULATOR.  For each source timeslice t in a window [tsrc0, tsrc0+twin)
// and each mode l, embed w_l(t) as a full-lattice spinor source (nonzero only at t), forward-solve the
// MASSLESS OVERLAP psi_l = D_ov^{-1} w_l (op_DH then op_Dsq.solve; D_ov normal), and contract per sink t' in
// the window: tau_kl(t',t) = V(t')^dag psi_l|_{t'}.  Also tau'_kl(t',t) = V(t')^dag (1-D_ov^dag) psi_l|_{t'}
// (FS furnished leg; ONE extra D_ov^dag apply, no second solve -- plan key-simpl #4).  Spin folded into k,l.
//   Window (plan): source = sink = [tsrc0, tsrc0+twin); twin=32 brackets the L1 scalar plateau dt[8,24].
//   Validation: T2a (reconstruct V(t')tau V(t)^dag vs a DIRECT D_ov^{-1} unit-source solve = gold L1 exact),
//   T2b ((IV.17): backward bar_tau = delta - tau via one D_ov^{-dag} solve), T2c (tau' consistency).
//
// Output: data_<ens>/distill_Nv<Nv>/peram.<k>.h5  with /meta, /evals (Nt,Nv), /V (Nt,Nv,2*N_s)
// (spinor component j = NS*x + s; reshape offline to (N_s,2)), and /peram/tau, /peram/tau_gw (twin,twin,Nv,Nv).
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
#include <random>
#include <Eigen/Dense>

using Double = double;
using Idx = std::int32_t;
using Complex = std::complex<double>;

using Link = std::array<Idx,2>;
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

  constexpr int NPARALLEL_DUPDATE=1;
  constexpr int NPARALLEL=NPARALLEL_DUPDATE;
  constexpr int NSTREAMS=4;
  constexpr int NPARALLEL_GAUGE=NPARALLEL_DUPDATE;
  constexpr int NPARALLEL_SORT=NPARALLEL_DUPDATE;

#ifndef N_REFINE_CLI
#define N_REFINE_CLI 1
#endif
  constexpr int N_REFINE=N_REFINE_CLI;
  constexpr int NS=2;
  constexpr int Nt=128;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;

  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-5;
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

#include "overlap_wmass_claude.h"        // (chunk 2) complex-mass overlap, massless at mass=0

//------------------------------------------
#include <getopt.h>

// extract a_t from an ensemble dir name (unique 'at' immediately followed by a digit); -1 if not found.
// (a_t is NOT used for the basis operator, which is purely spatial; kept for the chunk-2 overlap kernel.)
static double at_from_ensdir(const std::string& s){
  for(std::size_t p = 0; p + 2 < s.size(); ++p){
    const bool is_at = (s[p] == 'a' && s[p+1] == 't');
    const bool digit_after = (s[p+2] >= '0' && s[p+2] <= '9');
    if(is_at && digit_after) return std::stod(s.substr(p+2));
  }
  return -1.0;
}

// Fixed-timeslice gauge adaptor: DiracS2Simp::coo_format calls u(Link{ix,iy}) and expects the covariant
// U(1) phase \theta_{xy} (orientation-signed) on the spatial link.  U.sp(t,link) already folds map2sign, so
// this returns exactly the phase the 3D DiracExt uses for its spatial hop at timeslice t.
template<typename Gauge>
struct SliceGauge{
  const Gauge& U;
  int t;
  SliceGauge(const Gauge& U_, const int t_) : U(U_), t(t_) {}
  double operator()(const Link& ell) const { return U.sp(t, ell); }
};

// Assemble the dense 2 N_s x 2 N_s timeslice Wilson operator D_{W,2}(t) from the DiracS2Simp COO.  The COO
// (is,js,v_coo) IS the true matrix (plain assignment in coo2csr); duplicate (i,j) entries accumulate.
template<typename Dirac2D, typename Gauge>
static void assemble_Dslice(Eigen::MatrixXcd& M, const Dirac2D& D2,
                            const std::vector<Idx>& is, const std::vector<Idx>& js,
                            const Gauge& U, const int t){
  const Idx n = NS*Comp::N_SITES;
  M.setZero(n, n);
  SliceGauge<Gauge> ug(U, t);
  std::vector<Complex> v_coo(is.size(), Complex(0.0,0.0));
  D2.coo_format(v_coo, ug);
  for(std::size_t e=0; e<is.size(); e++){
    M(is[e], js[e]) += v_coo[e];
  }
}

// Embed mode l of the timeslice-t basis V_t (column l, a 2 N_s spinor) as a full-lattice source: nonzero
// ONLY on timeslice t (all sites, both spins), zero elsewhere.  src(t,x,s) = V_t(NS*x + s, l).
static void embed_mode_source(FermionVector& src, const Eigen::MatrixXcd& V_t, const int l, const int t){
  memset(src.field, 0, Comp::N*CD);
  for(Idx x=0; x<Comp::N_SITES; x++){
    for(int s=0; s<NS; s++){
      src(t, x, s) = V_t(NS*x+s, l);
    }
  }
}

// Extract the 2 N_s spinor subvector of a full-lattice field at timeslice tp: p[NS*x + s] = psi(tp,x,s).
static void extract_slice(Eigen::VectorXcd& p, const FermionVector& psi, const int tp){
  for(Idx x=0; x<Comp::N_SITES; x++){
    for(int s=0; s<NS; s++){
      p(NS*x+s) = psi(tp, x, s);
    }
  }
}

void PrintHelp(){
  printf("distill_peram (chunk 1): distillation BASIS V(t) = low modes of the timeslice Wilson normal op\n");
  printf("  D_{W,2}^dag D_{W,2} (spinor modes, spin connection included) for the sigma\\sigma - F^2 program.\n");
  printf("  --gsq <x>        Wilson coupling squared (ensemble id; default 8.0)\n");
  printf("  --Nf <n>         number of fermion flavors (ensemble id; default 2)\n");
  printf("  --nu0 <x>        sea quark asymmetry (ensemble id; default 1.0)\n");
  printf("  --nu1 <x>        valence Wilson-Dirac asymmetry for the overlap kernel (default nu0)\n");
  printf("  --at <x>         valence temporal spacing a_t (overlap kernel; default: derive from ens-dir; free 0.2)\n");
  printf("  --nv <n>         number of basis modes N_v (default 2*N_SITES = full = exact local)\n");
  printf("  --ens-dir <path> sea config dir; OMIT => free field (U=1)\n");
  printf("  --stride <s>     ensemble config stride (default 10)\n");
  printf("  --kmin <a> --kmax <b>  config range [a,b) (default 0..1e6)\n");
  printf("  --tsrc0 <t>      base source timeslice (default 0)\n");
  printf("  --twin <W>       source+sink window length; peram is (W,W,Nv,Nv) over [tsrc0,tsrc0+W) (default 32)\n");
  printf("  --tol <x>        outer CG tolerance for the overlap solves (default 1e-5; massless 1e-8 ~100x slower)\n");
  printf("  --gauge-check    (T1c) apply a random U(1) gauge transform on the first config and verify\n");
  printf("                   the eigenvalues are unchanged (covariance of the Laplacian).\n");
  printf("  -h, --help\n");
  exit(0);
}

void ParseArgs(int argc, char* argv[],
               double& gsq, int& Nf, double& nu0, double& nu1, double& at_cli,
               int& nv, std::string& ens_dir, int& stride,
               int& kmin, int& kmax, int& tsrc0, int& twin, double& tol, bool& gauge_check){
  static struct option long_opts[] = {
    {"gsq",     required_argument, nullptr, 'g'},
    {"Nf",      required_argument, nullptr, 'N'},
    {"nu0",     required_argument, nullptr, 'n'},
    {"nu1",     required_argument, nullptr, 'm'},
    {"at",      required_argument, nullptr, 'A'},
    {"nv",      required_argument, nullptr, 'v'},
    {"ens-dir", required_argument, nullptr, 'e'},
    {"stride",  required_argument, nullptr, 'I'},
    {"kmin",    required_argument, nullptr, 'a'},
    {"kmax",    required_argument, nullptr, 'b'},
    {"tsrc0",   required_argument, nullptr, 'S'},
    {"twin",    required_argument, nullptr, 'W'},
    {"tol",     required_argument, nullptr, 'O'},
    {"gauge-check", no_argument,   nullptr, 'C'},
    {"help",    no_argument,       nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
  };
  int opt, idx;
  while((opt = getopt_long(argc, argv, "g:N:n:m:A:v:e:I:a:b:S:W:O:Ch", long_opts, &idx)) != -1){
    switch(opt){
    case 'g': gsq     = std::stod(optarg); break;
    case 'N': Nf      = std::stoi(optarg); break;
    case 'n': nu0     = std::stod(optarg); break;
    case 'm': nu1     = std::stod(optarg); break;
    case 'A': at_cli  = std::stod(optarg); break;
    case 'v': nv      = std::stoi(optarg); break;
    case 'e': ens_dir = optarg; break;
    case 'I': stride  = std::stoi(optarg); break;
    case 'a': kmin    = std::stoi(optarg); break;
    case 'b': kmax    = std::stoi(optarg); break;
    case 'S': tsrc0   = std::stoi(optarg); break;
    case 'W': twin    = std::stoi(optarg); break;
    case 'O': tol     = std::stod(optarg); break;
    case 'C': gauge_check = true; break;
    case 'h':
    case '?':
    default:  PrintHelp(); break;
    }
  }
}
//------------------------------------------

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  double gsq=8.0;  int Nf=2;  double nu0=1.0;  double nu1=-1.0;
  int nv=-1;                    // -1 => full 2*N_SITES (exact local)
  std::string ens_dir="";
  int stride=10;
  int kmin=0;
  int kmax=1000000;
  int tsrc0=0;
  int twin=32;                  // source+sink window (brackets the L1 scalar plateau dt[8,24])
  double at_cli=-1.0;           // -1 => derive a_t from ens-dir (free field: 0.2)
  double tol=1.0e-5;            // outer CG tolerance.  MASSLESS overlap: near-zero modes make CG to 1e-8 ~100x
                               // slower (11 h/cfg observed) for NO physics gain -- at 1e-5 the reconstruction
                               // is already ~1e-6 (T2a=7e-7, correlated errors cancel) << gauge noise.  Tighten
                               // only via --tol if ever needed.
  bool gauge_check=false;

  ParseArgs(argc, argv, gsq, Nf, nu0, nu1, at_cli, nv, ens_dir, stride, kmin, kmax, tsrc0, twin, tol, gauge_check);
  if(nu1 < 0.0) nu1 = nu0;

  const bool free_field = ens_dir.empty();
  const Idx n2 = NS*Comp::N_SITES;               // spinor-space dimension per timeslice (= 2 N_s)
  if(nv < 0) nv = n2;
  assert(nv >= 1 && nv <= n2 && "N_v out of range [1, 2*N_SITES]");

  const double M5 = -1.0;                         // matches the overlap kernel D_W(M5=-1)
  double at = at_cli;
  if(at < 0.0) at = free_field ? 0.2 : at_from_ensdir(ens_dir);
  assert(at > 0.0 && "could not resolve a_t; pass --at explicitly");
  assert(twin >= 1 && tsrc0 >= 0 && tsrc0 + twin <= Comp::Nt && "require 1<=twin, 0<=tsrc0, tsrc0+twin<=Nt");

  std::cout << "# distill_peram (chunk1+2): gsq="<<gsq<<" Nf="<<Nf<<" nu0="<<nu0<<" nu1="<<nu1<<" at="<<at
            << " nv="<<nv<<" / 2N_s="<<n2<<" (N_SITES="<<Comp::N_SITES<<", L="<<Comp::N_REFINE<<")"
            << " ens_dir="<<(free_field?std::string("<free-field U=1>"):ens_dir)
            << " M5="<<M5<<" tsrc0="<<tsrc0<<" twin="<<twin << std::endl;

  // device setup (chunk 1 math is CPU/Eigen; this matches the other drivers + readies chunk-2 overlap solves)
  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();
  {
    int device;
    CUDA_CHECK(cudaGetDeviceCount(&device));
    cudaDeviceProp device_prop[device];
    cudaGetDeviceProperties(&device_prop[0], 0);
    std::cout << "# dev = " << device_prop[0].name << std::endl;
    CUDA_CHECK(cudaSetDevice(0));
  }

  // basis operator: 2D S^2 Wilson-Dirac (m=0, r=1, M5).  Reads omega/alpha spin-structure files.
  using Base=S2Simp;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;
  Base base(Comp::N_REFINE);
  std::cout << "# lattice set." << std::endl;
  DiracS2Simp D2(base, 0.0, 1.0, M5);
  std::cout << "# 2D Wilson-Dirac (basis operator) set." << std::endl;

  // COO structure is config/timeslice INDEPENDENT: build once.
  std::vector<Idx> is, js;
  D2.coo_structure(is, js);
  std::cout << "# COO structure built (nnz="<<is.size()<<")." << std::endl;

  Gauge U(base);

  // ---- overlap kernel + operators for the perambulator (chunk 2): massless D_ov, forward solves only ----
  // D_ov^{-1} b = op_Dsq^{-1}(op_DH b)  (D_ov normal);  D_ov^{-dag} b = op_D(op_Dsq^{-1} b)  (T2b backward).
  constexpr Idx N = Comp::N;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Fermion=OverlapWMass<WilsonDirac>;
  if(Comp::Nt!=1) assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu1);
  Fermion D(DW, Complex(0.0), 11);
  std::cout << "# overlap operator set: D_ov (massless, M5="<<M5<<")." << std::endl;

  auto f_D   = std::bind(&Fermion::mult_deviceAsyncLaunch_ms, &D, std::placeholders::_1, std::placeholders::_2);
  auto f_DH  = std::bind(&Fermion::adj_deviceAsyncLaunch_ms,  &D, std::placeholders::_1, std::placeholders::_2);
  auto f_Dsq = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms,  &D, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_D(f_D), M_DH(f_DH), M_Dsq(f_Dsq);
  MatPoly op_D;   op_D.push_back(cplx(1.0), {&M_D});
  MatPoly op_DH;  op_DH.push_back(cplx(1.0), {&M_DH});
  MatPoly op_Dsq; op_Dsq.push_back(cplx(1.0), {&M_Dsq});
  // op_oneMinusDdag : v -> (1 - D_ov^dag) v  (FS furnished leg for tau').
  MatPoly op_oneMinusDdag;
  op_oneMinusDdag.push_back(cplx( 1.0), {});
  op_oneMinusDdag.push_back(cplx(-1.0), {&M_DH});

  // perambulator work vectors (host-pinned FermionVectors) reused across (source t, mode l)
  FermionVector src, tmp, psi, dpsi, uni, psi_d, dpsi_d, chi;

  // ---- output dir: data_<ESNID>/distill_Nv<Nv>/peram.<k>.h5 ----
  std::string ens_base = ens_dir;
  if(!ens_base.empty() && ens_base.back()=='/') ens_base.pop_back();
  { const auto slash = ens_base.find_last_of('/'); if(slash!=std::string::npos) ens_base = ens_base.substr(slash+1); }
  const std::string esnid = (free_field ? std::string("free") : ens_base);
  const std::string dir_out = "data_"+esnid+"/distill_Nv"+std::to_string(nv)+"/";
  std::filesystem::create_directories(dir_out);
  std::cout << "# dir_out = " << dir_out << std::endl;

  const int Nt = Comp::Nt;

  const int k_ckpoint = free_field ? 1 : stride;
  const int k_lo      = free_field ? 0 : kmin;
  const int k_hi      = free_field ? 1 : kmax;

  for(int k = k_lo; k < k_hi; k += k_ckpoint){
    std::string str_lat;
    if(!free_field){
      str_lat = ens_dir + "ckpoint_lat." + std::to_string(k);
      if(!std::filesystem::exists(str_lat)){ if(k==0) continue; else break; }
    }
    const std::string h5path = dir_out + "peram." + std::to_string(k) + ".h5";
    // skip if the perambulator is already present.
    if(std::filesystem::exists(h5path)){
      bool has_peram=false;
      try {
        HighFive::File f(h5path, HighFive::File::ReadOnly);
        if(f.exist("peram")) has_peram = f.getGroup("peram").exist("tau");
      } catch(...) {}
      if(has_peram){ std::cout<<"# skip k="<<k<<" (peram present)"<<std::endl; continue; }
    }
    if(!free_field) U.read(str_lat);
    D.update(U);

    const auto t_k0 = std::chrono::steady_clock::now();

    // evals[t][k], V[t] as (n2 x nv) column matrix.  Stored flat for h5 below.
    std::vector<std::vector<double>> evals(Nt, std::vector<double>(nv, 0.0));
    std::vector<Eigen::MatrixXcd> Vt(Nt);

    // validation accumulators (per config)
    double min_eval = 1.0e300;
    double max_ortho = 0.0;      // T1a: max |V^dag V - I|
    double max_compl = 0.0;      // T1d: max |V V^dag - I| (only meaningful at nv==n2)

    // BASIS_SYM (compile-time, -DBASIS_SYM=1): distillation-basis operator whose low modes form V(t).
    //   0 (default) = D_W^H D_W (original, ONE-SIDED normal op).  Its eigenbasis does NOT carry the theory's
    //     symmetry (sigma3/parity maps D_W^H D_W eigenmodes to D_W D_W^H eigenmodes -- a DIFFERENT operator),
    //     so a TRUNCATED basis (Nv < 2 N_s) spans a NON-symmetry-invariant subspace and breaks symmetry-
    //     protected zeros -- e.g. the single meson m_PS leaks into sigma^2 at L2 Nv=24-of-84 (dominates the P+
    //     GEVP ground; verified free: complete Nv84 -> two-meson, truncated Nv24 -> m_PS).
    //   1 = the SYMMETRIZED operator D_W^H D_W + D_W D_W^H, which IS invariant under the symmetry that swaps the
    //     two normal ops, so its low-mode truncation spans a symmetry-invariant subspace and the protection
    //     survives truncation.  (Identical to option 0 at the COMPLETE basis and whenever D_W is normal; it only
    //     differs, and only helps, when D_W is non-normal -- which the Wilson term makes true.)
    //   Rationale + the truncation-leak diagnosis: fs_gw_collapse_v_agent_two_meson_claude.md.
#ifndef BASIS_SYM
#define BASIS_SYM 0
#endif
    Eigen::MatrixXcd M(n2, n2), MhM(n2, n2);
    for(int t=0; t<Nt; t++){
      assemble_Dslice(M, D2, is, js, U, t);
#if BASIS_SYM
      MhM.noalias() = M.adjoint() * M + M * M.adjoint();   // symmetrized basis (symmetry-respecting -> truncation-safe)
#else
      MhM.noalias() = M.adjoint() * M;             // Hermitian PSD; low modes = small singular vals of D2 (original one-sided)
#endif
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(MhM);   // eigenvalues ASCENDING
      const Eigen::VectorXd&  lam = es.eigenvalues();
      const Eigen::MatrixXcd& Q   = es.eigenvectors();

      Vt[t] = Q.leftCols(nv);                       // lowest nv modes
      for(int a=0; a<nv; a++) evals[t][a] = lam(a);
      if(lam(0) < min_eval) min_eval = lam(0);

      // T1a orthonormality residual (columns are orthonormal by construction; confirm numerically)
      Eigen::MatrixXcd G = Vt[t].adjoint() * Vt[t];
      for(int a=0; a<nv; a++) G(a,a) -= 1.0;
      const double o = G.cwiseAbs().maxCoeff();
      if(o > max_ortho) max_ortho = o;

      // T1d completeness (exact local): only at full rank nv==n2
      if(nv == n2){
        Eigen::MatrixXcd P = Vt[t] * Vt[t].adjoint();
        for(int a=0; a<n2; a++) P(a,a) -= 1.0;
        const double c = P.cwiseAbs().maxCoeff();
        if(c > max_compl) max_compl = c;
      }
    } // t

    // ---- (T1c) gauge covariance on this config (first config only, if requested) ----
    if(gauge_check && (free_field || k==k_lo)){
      Gauge Ug(U);
      ParallelRngExt<Base,Comp::Nt> rngc(base, 4321);
      Ug.random_gauge_trsf(rngc, 1.0);
      double max_dlam = 0.0;
      Eigen::MatrixXcd Mg(n2,n2), MhMg(n2,n2);
      for(int t=0; t<Nt; t++){
        assemble_Dslice(Mg, D2, is, js, Ug, t);
        MhMg.noalias() = Mg.adjoint() * Mg;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> esg(MhMg);
        const Eigen::VectorXd& lg = esg.eigenvalues();
        for(int a=0; a<nv; a++){
          const double d = std::abs(lg(a) - evals[t][a]);
          if(d > max_dlam) max_dlam = d;
        }
      }
      std::cout << "#   [T1c gauge covariance] max |lambda(U^g) - lambda(U)| = " << max_dlam
                << "  (should be ~1e-12)" << std::endl;
    }

    // ---- validation summary (per config) ----
    std::cout << "#   [T1a] min eigenvalue = " << min_eval
              << " (>=0 ?)   max|V^dag V - 1| = " << max_ortho << std::endl;
    if(nv == n2){
      std::cout << "#   [T1d] N_v=2N_s exact local: max|V V^dag - 1| = " << max_compl
                << "  (should be ~1e-13)" << std::endl;
    }
    // shell structure (T1b): print the lowest eigenvalues at t=0 to eyeball spinor-harmonic degeneracies.
    {
      std::cout << "#   [T1b shells] lowest evals (t=0):";
      const int nshow = std::min(nv, 16);
      for(int a=0; a<nshow; a++) std::cout << " " << evals[0][a];
      std::cout << std::endl;
    }

    // ================= CHUNK 2: perambulator on the window [tsrc0, tsrc0+twin) =================
    // Tau[a'][a](k,l) = tau_{kl}(t'=tsrc0+a', t=tsrc0+a) = w_k(t')^dag D_ov^{-1} w_l(t).
    // Taup = tau' (FS furnished leg): w_k(t')^dag (1 - D_ov^dag) D_ov^{-1} w_l(t).
    std::vector<std::vector<Eigen::MatrixXcd>>
      Tau (twin, std::vector<Eigen::MatrixXcd>(twin, Eigen::MatrixXcd::Zero(nv,nv))),
      Taup(twin, std::vector<Eigen::MatrixXcd>(twin, Eigen::MatrixXcd::Zero(nv,nv)));
    {
      Eigen::VectorXcd p(n2), pd(n2);
      for(int a=0; a<twin; a++){
        const int t = tsrc0 + a;
        for(int l=0; l<nv; l++){
          embed_mode_source(src, Vt[t], l, t);
          op_DH.from_cpu<N>(tmp.field, src.field);                     // D_ov^dag w_l
          op_Dsq.solve<N>(psi.field, tmp.field, tol);      // psi = D_ov^{-1} w_l(t)
          op_oneMinusDdag.from_cpu<N>(dpsi.field, psi.field);          // (1 - D_ov^dag) psi
          for(int ap=0; ap<twin; ap++){
            const int tp = tsrc0 + ap;
            extract_slice(p,  psi,  tp);
            extract_slice(pd, dpsi, tp);
            Tau [ap][a].col(l) = Vt[tp].adjoint() * p;
            Taup[ap][a].col(l) = Vt[tp].adjoint() * pd;
          }
        }
      }
    }

    // (T2a) reconstruct V(t') tau(t',tsrc0) [V(tsrc0)^dag e_{j0}] vs a DIRECT D_ov^{-1} unit-source solve;
    // (T2c) same for tau' vs (1-D_ov^dag) D_ov^{-1}.  Exact at nv=2N_s (-> ~solve tol).
    {
      const int t0 = tsrc0;
      const int j0 = 0;                              // unit source at (tsrc0, x=0, s=0)
      memset(uni.field, 0, Comp::N*CD);
      uni(t0, 0, 0) = Complex(1.0, 0.0);
      op_DH.from_cpu<N>(tmp.field, uni.field);
      op_Dsq.solve<N>(psi_d.field, tmp.field, tol);        // psi_d = D_ov^{-1} e_{j0}(t0)
      op_oneMinusDdag.from_cpu<N>(dpsi_d.field, psi_d.field);          // (1 - D_ov^dag) psi_d
      Eigen::VectorXcd w0dag = Vt[t0].row(j0).adjoint();              // V(t0)^dag e_{j0}  (length nv)
      Eigen::VectorXcd p(n2), pd(n2);
      double max_a=0.0, max_c=0.0;
      for(int ap=0; ap<twin; ap++){
        const int tp = tsrc0 + ap;
        Eigen::VectorXcd recon  = Vt[tp] * (Tau [ap][0] * w0dag);
        Eigen::VectorXcd reconc = Vt[tp] * (Taup[ap][0] * w0dag);
        extract_slice(p,  psi_d,  tp);
        extract_slice(pd, dpsi_d, tp);
        const double da = (recon  - p ).cwiseAbs().maxCoeff();
        const double dc = (reconc - pd).cwiseAbs().maxCoeff();
        if(da>max_a) max_a=da;
        if(dc>max_c) max_c=dc;
      }
      std::cout << "#   [T2a] max|V tau V^dag - D_ov^{-1}| (unit src) = " << max_a
                << "  (exact at nv="<<n2<<" -> ~solve tol)" << std::endl;
      std::cout << "#   [T2c] max|V tau' V^dag - (1-D_ov^dag)D_ov^{-1}| = " << max_c << std::endl;
    }

    // (T2b) (IV.17): backward peram bar_tau(t',t0) = w_k(t')^dag D_ov^{-dag} w_l(t0) should equal
    // delta_{t't0} I - tau(t',t0).  D_ov^{-dag} b = op_D(op_Dsq^{-1} b).  First config only (nv extra solves).
    if(free_field || k==k_lo){
      const int t0 = tsrc0;
      std::vector<Eigen::MatrixXcd> barTau(twin, Eigen::MatrixXcd::Zero(nv,nv));
      Eigen::VectorXcd p(n2);
      for(int l=0; l<nv; l++){
        embed_mode_source(src, Vt[t0], l, t0);
        op_Dsq.solve<N>(tmp.field, src.field, tol);        // (D D^dag)^{-1} w_l
        op_D.from_cpu<N>(chi.field, tmp.field);                       // chi = D_ov^{-dag} w_l
        for(int ap=0; ap<twin; ap++){
          const int tp = tsrc0 + ap;
          extract_slice(p, chi, tp);
          barTau[ap].col(l) = Vt[tp].adjoint() * p;
        }
      }
      double max_b=0.0;
      for(int ap=0; ap<twin; ap++){
        Eigen::MatrixXcd resid = barTau[ap] + Tau[ap][0];             // bar_tau + tau
        if(ap==0) for(int kk=0; kk<nv; kk++) resid(kk,kk) -= 1.0;     // - delta_{t't0} I
        const double b = resid.cwiseAbs().maxCoeff();
        if(b>max_b) max_b=b;
      }
      std::cout << "#   [T2b] (IV.17) max|bar_tau - (delta - tau)| = " << max_b
                << "  (should be ~solve tol)" << std::endl;
    }

    // ---- write h5: /meta, /evals (Nt,Nv), /V (Nt,Nv,2*N_s), /peram/{tau,tau_gw} (twin,twin,Nv,Nv) ----
    const std::string h5tmp = h5path + ".tmp";
    {
      std::unique_ptr<HighFive::File> h5p = std::make_unique<HighFive::File>(h5tmp,
              HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
      HighFive::File& h5 = *h5p;
      h5.createDataSet("meta/Nv",       std::vector<int>{nv});
      h5.createDataSet("meta/Nt",       std::vector<int>{Nt});
      h5.createDataSet("meta/Ns",       std::vector<int>{(int)Comp::N_SITES});
      h5.createDataSet("meta/Nspin",    std::vector<int>{NS});
      h5.createDataSet("meta/config_k", std::vector<int>{k});
      h5.createDataSet("meta/M5",       std::vector<double>{M5});
      h5.createDataSet("meta/at",       std::vector<double>{at});
      h5.createDataSet("meta/nu0",      std::vector<double>{nu0});
      h5.createDataSet("meta/nu1",      std::vector<double>{nu1});
      h5.createDataSet("meta/overlap_tol", std::vector<double>{tol});
      h5.createDataSet("meta/tsrc0",    std::vector<int>{tsrc0});
      h5.createDataSet("meta/twin",     std::vector<int>{twin});
      h5.createDataSet("meta/method",   std::string("exact_distillation_wilson_basis"));
      h5.createDataSet("meta/ensemble", esnid);

      // evals: (Nt, Nv)
      std::vector<std::vector<double>> ev(Nt, std::vector<double>(nv));
      for(int t=0;t<Nt;t++) for(int a=0;a<nv;a++) ev[t][a]=evals[t][a];
      h5.createDataSet("evals", ev);

      // V: (Nt, Nv, 2*N_s) split real/imag.  V[t][a][j] = eigvec component j (= NS*x + s) of mode a.
      // Spinor component index j runs 0..2N_s-1; reshape offline to (N_s, 2) as x=j/2, s=j%2.
      std::vector<std::vector<std::vector<double>>>
        Vre(Nt, std::vector<std::vector<double>>(nv, std::vector<double>(n2, 0.0))),
        Vim = Vre;
      for(int t=0;t<Nt;t++){
        for(int a=0;a<nv;a++){
          for(Idx j=0;j<n2;j++){
            const Complex c = Vt[t](j, a);
            Vre[t][a][j] = c.real();
            Vim[t][a][j] = c.imag();
          }
        }
      }
      h5.createDataSet("V/real", Vre);
      h5.createDataSet("V/imag", Vim);

      // peram: (twin, twin, Nv, Nv) split re/im, row-major [a_snk, a_src, k, l].  Flat + DataSpace so the
      // 4D shape survives for h5py (nested-vector depth kept modest).  t' = tsrc0 + a_snk, t = tsrc0 + a_src.
      const std::size_t nT = (std::size_t)twin*twin*nv*nv;
      std::vector<double> Tre(nT), Tim(nT), Tpre(nT), Tpim(nT);
      std::size_t idx=0;
      for(int ap=0; ap<twin; ap++){
        for(int a=0; a<twin; a++){
          for(int kk=0; kk<nv; kk++){
            for(int ll=0; ll<nv; ll++){
              const Complex c  = Tau [ap][a](kk,ll);
              const Complex cp = Taup[ap][a](kk,ll);
              Tre[idx]  = c.real();
              Tim[idx]  = c.imag();
              Tpre[idx] = cp.real();
              Tpim[idx] = cp.imag();
              idx++;
            }
          }
        }
      }
      const std::vector<std::size_t> pdims = { (std::size_t)twin, (std::size_t)twin, (std::size_t)nv, (std::size_t)nv };
      h5.createDataSet<double>("peram/tau/real",    HighFive::DataSpace(pdims)).write_raw(Tre.data());
      h5.createDataSet<double>("peram/tau/imag",    HighFive::DataSpace(pdims)).write_raw(Tim.data());
      h5.createDataSet<double>("peram/tau_gw/real", HighFive::DataSpace(pdims)).write_raw(Tpre.data());
      h5.createDataSet<double>("peram/tau_gw/imag", HighFive::DataSpace(pdims)).write_raw(Tpim.data());
    }
    std::filesystem::rename(h5tmp, h5path);
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now()-t_k0).count();
    std::cout << "# k="<<k<<(free_field?" (free field)":"")<<" basis+peram done ["<<secs<<" s] -> "<<h5path<<std::endl;
  } // k

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();
  return 0;
}
