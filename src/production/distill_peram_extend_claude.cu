// distill_peram_extend_claude.cu
// VARIANT of distill_peram_mrhs_claude.cu (original untouched): INCREMENTAL basis extension N_o -> N_n without redoing
// the stored D_ov solves (NM 2026-11-03: L3 24->32, L4 48->64).  Plan: distill_peram_extend_impl_plan_claude.md;
// derivation + check: final/analysis_axial/peram_extend_nv_note_claude.md.  Blocks of tau_ij(t',t) (old=0..N_o-1, new=N_o..N_n-1):
//   old-old : READ from the stored file (no work);
//   any-new : NEW solves for the new sources only (block width NSTACK = N_n - N_o), projected on all N_n modes;
//   new-old : identity tau_ij(t',t) = delta_{t't} delta_ij - conj(tau_ji(t,t')) from (D_ov^{-1})^dag = 1 - D_ov^{-1}
//             (D_ov = 1 + V, V unitary; paper-21 Eq. A.25 / IV.16), using the new-source solves (i != j: no delta).
// Old modes = the STORED V (exact consistency with the stored tau); new modes = fresh BASIS_SYM low modes N_o..N_n-1,
// re-orthonormalized against the stored ones.  Per-config checks: stored evals reproduced, projector onto the old
// span reproduced, eigenvalue gaps at the cuts; first config: identity block vs a DIRECT old-source solve (window 0, a=0).
// tau_gw is NOT written (NM 2026-11-03: extended perams are for the axial sp only).  meta/tau_gw_present = 0.
// ------------------------------------------------------------------------------------------------------------------
// (original header of distill_peram_mrhs_claude.cu follows)
// distill_peram_mrhs_claude.cu
// VARIANT of distill_peram_claude.cu (do NOT touch the original: a job runs off distill_peram_L1_claude.o).
// Two changes vs the original, plan: distill_peram_mrhs_impl_plan_claude.md.
//   (1) MRHS-BLOCK the N_v mode solves.  The per-mode single-RHS op_Dsq.solve (D_ov^dag D_ov)^{-1} loop is
//       replaced by ONE BlockedMat::solve_sq_from_cpu over an NSTACK-wide host block (mrhs block CG, bit-
//       identical per column -- validated C6a-d in blocked_mat_claude.h).  NSTACK is COMPILE-TIME = 2 N_s
//       (= Comp::Nx, the full basis; 24 @L1).  If --nv < 2 N_s, the extra block columns carry a zero RHS and
//       freeze immediately in the block CG (b_norm_sq ~ 0), so a full-width block is always safe.
//   (2) SECOND SOURCE WINDOW: --tsrc-list / --nsrc select several base source timeslices (nsrc=2 -> {0,Nt/2}).
//       V(t) is all-t (built once/config); only the solves+contraction repeat per window.  h5 gains a leading
//       nsrc axis on tau/tau_gw (nsrc=1 keeps the OLD layout for backward-compatible readers).
// UNCERTAINTIES (resolve if a test flags them): NSTACK is fixed at the FULL basis 2 N_s = Comp::Nx (not the
//   runtime --nv) so the block width is a compile-time constant; solve_sq_from_cpu(x_host,b_host,tol) is the
//   host-block entry (device I/O hidden) and is in-place safe (x may alias b), per blocked_mat_claude.h:448.
//
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
// Output: data_<ens>/distill_Nv<Nv>_v2/peram.<k>.h5  (SEPARATE "_v2" dir -- never overwrites the existing
//         distill_Nv<Nv>/ production perams)  with /meta, /evals (Nt,Nv), /V (Nt,Nv,2*N_s)
// (spinor component j = NS*x + s; reshape offline to (N_s,2)), and /peram/tau, /peram/tau_gw of shape
// (nsrc,twin,twin,Nv,Nv) when nsrc>1 (leading window axis; meta/tsrc_list gives the base t per window), or
// (twin,twin,Nv,Nv) when nsrc==1 (OLD layout, backward-compatible with existing readers).
#include <typeinfo>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
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
#include "blocked_mat_claude.h"          // BlockedMat<N,NSTACK,Op> mrhs block solves (solve_sq_from_cpu); C7

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

// FUSED multi-source (--fused-sources): embed mode l at EVERY window's source timeslice ts[j] from its own
// basis Vts[j], summed into ONE full-lattice source.  By linearity D_ov^{-1} of this sum = sum of the
// per-source solutions; well-separated timeslices ($N_t/2$ apart) do not interfere within each window (the
// far-source tail is $e^{-m_\text{gap}\Delta}$-suppressed), so one solve yields all windows' perambulators.
static void embed_mode_source_fused(FermionVector& src, const std::vector<const Eigen::MatrixXcd*>& Vts,
                                    const std::vector<int>& ts, const int l){
  memset(src.field, 0, Comp::N*CD);
  for(std::size_t j=0; j<ts.size(); j++){
    const int t = ts[j];
    const Eigen::MatrixXcd& V_t = *Vts[j];
    for(Idx x=0; x<Comp::N_SITES; x++){
      for(int s=0; s<NS; s++){
        src(t, x, s) = V_t(NS*x+s, l);
      }
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
  printf("  --out-suffix <s> output dir = distill_Nv<nv><s>/ (default \"_v2\"); \"\" merges into distill_Nv<nv>/\n");
  printf("  --stride <s>     ensemble config stride (default 10)\n");
  printf("  --kmin <a> --kmax <b>  config range [a,b) (default 0..1e6)\n");
  printf("  --tsrc0 <t>      base source timeslice (default 0)\n");
  printf("  --nsrc <n>       number of source windows; n>1 -> base t = tsrc0 + i*(Nt/n), i=0..n-1 (default 1)\n");
  printf("  --tsrc-list a,b  explicit comma-separated base source timeslices (overrides --nsrc and --tsrc0)\n");
  printf("  --twin <W>       source+sink window length; peram is (W,W,Nv,Nv) per window over [t,t+W) (default 32)\n");
  printf("  --tol <x>        outer CG tolerance for the overlap solves (default 1e-5; massless 1e-8 ~100x slower)\n");
  printf("  --gauge-check    (T1c) apply a random U(1) gauge transform on the first config and verify\n");
  printf("                   the eigenvalues are unchanged (covariance of the Laplacian).\n");
  printf("  --fused-sources  nsrc>1: solve ONE combined source (sum over windows) per offset and extract every\n");
  printf("                   window from it (well-separated -> no interference); solves = twin, not nsrc*twin.\n");
  printf("  --nv-old <n>     EXTEND: stored basis size N_o (required); --nv = new size N_n > N_o\n");
  printf("  --old-suffix <s> EXTEND: stored dir = distill_Nv<N_o><s>/ (default \"_sym\")\n");
  printf("  -h, --help\n");
  exit(0);
}

void ParseArgs(int argc, char* argv[],
               double& gsq, int& Nf, double& nu0, double& nu1, double& at_cli,
               int& nv, std::string& ens_dir, int& stride,
               int& kmin, int& kmax, int& tsrc0, int& twin, double& tol, bool& gauge_check,
               int& nsrc, std::string& tsrc_list_str, std::string& out_suffix, bool& fused,
               int& nv_old, std::string& old_suffix){
  static struct option long_opts[] = {
    {"out-suffix", required_argument, nullptr, 'X'},   // output dir = distill_Nv<nv><suffix>/ (default "_v2"); "" merges into distill_Nv<nv>/
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
    {"nsrc",    required_argument, nullptr, 'R'},
    {"tsrc-list", required_argument, nullptr, 'L'},
    {"twin",    required_argument, nullptr, 'W'},
    {"tol",     required_argument, nullptr, 'O'},
    {"gauge-check", no_argument,   nullptr, 'C'},
    {"fused-sources", no_argument, nullptr, 'F'},
    {"nv-old",  required_argument, nullptr, 'o'},
    {"old-suffix", required_argument, nullptr, 'P'},
    {"help",    no_argument,       nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
  };
  int opt, idx;
  while((opt = getopt_long(argc, argv, "X:g:N:n:m:A:v:e:I:a:b:S:R:L:W:O:CFo:P:h", long_opts, &idx)) != -1){
    switch(opt){
    case 'X': out_suffix = optarg; break;
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
    case 'R': nsrc    = std::stoi(optarg); break;
    case 'L': tsrc_list_str = optarg; break;
    case 'W': twin    = std::stoi(optarg); break;
    case 'O': tol     = std::stod(optarg); break;
    case 'C': gauge_check = true; break;
    case 'F': fused = true; break;
    case 'o': nv_old  = std::stoi(optarg); break;
    case 'P': old_suffix = optarg; break;
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
  int nsrc=1;                   // number of source windows (multi-source); 1 => original single-window layout
  std::string tsrc_list_str=""; // explicit comma list of base source timeslices (overrides nsrc/tsrc0)
  std::string out_suffix="_sym_ext"; // EXTEND: output dir distill_Nv<N_n>_sym_ext/ (never the stored dirs)
  int nv_old=-1;                // EXTEND: stored basis size N_o (required)
  std::string old_suffix="_sym"; // EXTEND: stored dir distill_Nv<N_o><old_suffix>/
  bool fused=false;             // --fused-sources: solve ONE combined source per offset, extract all windows

  ParseArgs(argc, argv, gsq, Nf, nu0, nu1, at_cli, nv, ens_dir, stride, kmin, kmax, tsrc0, twin, tol, gauge_check,
            nsrc, tsrc_list_str, out_suffix, fused, nv_old, old_suffix);
#if !BASIS_SYM
  static_assert(BASIS_SYM, "extension driver: the stored perams use the symmetrized basis -> compile with -DBASIS_SYM=1");
#endif
  if(nu1 < 0.0) nu1 = nu0;

  // ---- build the list of base source timeslices (multi-source change 2) ----
  // Priority: explicit --tsrc-list ; else --nsrc>1 spreads evenly (tsrc0 + i*Nt/nsrc) ; else single {tsrc0}.
  std::vector<int> tsrc_list;
  if(!tsrc_list_str.empty()){
    std::string tok;
    std::stringstream ss(tsrc_list_str);
    while(std::getline(ss, tok, ',')){
      if(!tok.empty()) tsrc_list.push_back(std::stoi(tok));
    }
    assert(!tsrc_list.empty() && "--tsrc-list parsed to empty");
  }
  else if(nsrc > 1){
    const int dt = Comp::Nt / nsrc;                  // window spacing (Nt/2=64 for nsrc=2)
    for(int i=0; i<nsrc; i++) tsrc_list.push_back(tsrc0 + i*dt);
  }
  else{
    tsrc_list.push_back(tsrc0);
  }
  nsrc = (int)tsrc_list.size();

  const bool free_field = ens_dir.empty();
  const Idx n2 = NS*Comp::N_SITES;               // spinor-space dimension per timeslice (= 2 N_s)
  if(nv < 0) nv = n2;
  assert(nv >= 1 && nv <= n2 && "N_v out of range [1, 2*N_SITES]");
  assert(nv_old >= 1 && nv_old < nv && "EXTEND: require 1 <= --nv-old < --nv");
  const int nnew = nv - nv_old;                  // number of NEW source modes (= block solve width)
  assert(out_suffix != old_suffix && "EXTEND: output suffix must differ from the stored dir suffix (no overwrite)");

  const double M5 = -1.0;                         // matches the overlap kernel D_W(M5=-1)
  double at = at_cli;
  if(at < 0.0) at = free_field ? 0.2 : at_from_ensdir(ens_dir);
  assert(at > 0.0 && "could not resolve a_t; pass --at explicitly");
  assert(twin >= 1 && "require twin>=1");
  for(int s=0; s<nsrc; s++){
    assert(tsrc_list[s] >= 0 && tsrc_list[s] + twin <= Comp::Nt && "require 0<=tsrc, tsrc+twin<=Nt for every window");
  }

  std::cout << "# distill_peram_mrhs (chunk1+2): gsq="<<gsq<<" Nf="<<Nf<<" nu0="<<nu0<<" nu1="<<nu1<<" at="<<at
            << " nv="<<nv<<" / 2N_s="<<n2<<" (N_SITES="<<Comp::N_SITES<<", L="<<Comp::N_REFINE<<")"
            << " ens_dir="<<(free_field?std::string("<free-field U=1>"):ens_dir)
            << " M5="<<M5<<" twin="<<twin<<" nsrc="<<nsrc<<" tsrc_list={";
  for(int s=0; s<nsrc; s++) std::cout << (s?",":"") << tsrc_list[s];
  std::cout << "}" << std::endl;

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

  // ---- mrhs (change 1): ONE BlockedMat over the FULL basis width, created ONCE (heavy ctor: allocs the
  // block scratch/pool), reused every (config, source window, source timeslice).  NSTACK is COMPILE-TIME =
  // 2 N_s (= Comp::Nx, the complete basis; 24 @L1) so the block width is a constant; --nv < 2N_s just leaves
  // the trailing columns with a zero RHS (they freeze immediately in the block CG).  solve_sq_from_cpu does
  // (D_ov^dag D_ov)^{-1} per column, bit-identical to op_Dsq.solve<N> (validated C6a-d in blocked_mat).
  // NSTACK = block-solve width.  Default = Comp::Nx (full 2 N_s basis).  Override with -DNSTACK_CLI=<n> to make
  // the block matvec only n-wide -- REQUIRED for a real speedup when TRUNCATING the basis (--nv n < 2N_s): with
  // the default NSTACK the frozen trailing columns still ride the n-wide-vs-full matvec, so --nv alone saves the
  // eigensolve/contraction but NOT the dominant block solve.  Require nv <= NSTACK.
#ifndef NSTACK_CLI
  constexpr int NSTACK = Comp::Nx;                 // = NS * N_SITES = full distillation basis width
#else
  constexpr int NSTACK = NSTACK_CLI;               // truncated block width (e.g. 24 to match L1's basis size)
#endif
  assert(nnew <= NSTACK && "EXTEND: require N_n - N_o <= NSTACK (compile -DNSTACK_CLI=<N_n-N_o>)");
  BlockedMat<N, NSTACK, Fermion> blk_Dsq(D);       // holds a const ref to D: per-config coeffs read at solve time
  std::vector<Complex> hblk((std::size_t)N * NSTACK);   // host RHS/solution block staging (device I/O hidden)
  std::cout << "# mrhs block engine set: NSTACK="<<NSTACK<<" (2N_s="<<Comp::Nx<<"), one solve_sq_from_cpu per source t."
            << std::endl;

  // perambulator work vectors (host-pinned FermionVectors) reused across (source t, mode l)
  FermionVector src, tmp, psi, dpsi, uni, psi_d, dpsi_d, chi;

  // ---- output dir: data_<ESNID>/distill_Nv<Nv>_v2/peram.<k>.h5  (SEPARATE "_v2" dir: this mrhs/multi-source
  //      variant NEVER writes into the existing distill_Nv<Nv>/ -- no risk of overwriting production perams) ----
  std::string ens_base = ens_dir;
  if(!ens_base.empty() && ens_base.back()=='/') ens_base.pop_back();
  { const auto slash = ens_base.find_last_of('/'); if(slash!=std::string::npos) ens_base = ens_base.substr(slash+1); }
  const std::string esnid = (free_field ? std::string("free") : ens_base);
  // output dir suffix (default "_v2"); pass --out-suffix "" to MERGE into the production distill_Nv<nv>/ (the
  // nsrc==1 h5 is byte-identical to _v1, and the skip-if-present guard protects existing perams -> no overwrite).
  const std::string dir_out = "data_"+esnid+"/distill_Nv"+std::to_string(nv)+out_suffix+"/";
  std::filesystem::create_directories(dir_out);
  std::cout << "# dir_out = " << dir_out << std::endl;
  const std::string dir_old = "data_"+esnid+"/distill_Nv"+std::to_string(nv_old)+old_suffix+"/";
  std::cout << "# EXTEND: N_o="<<nv_old<<" -> N_n="<<nv<<" ("<<nnew<<" new source modes), stored dir = "<<dir_old << std::endl;

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
    // ---- EXTEND: the stored N_o perambulator for this config (skip LOUDLY if missing) ----
    const std::string h5old = dir_old + "peram." + std::to_string(k) + ".h5";
    if(!std::filesystem::exists(h5old)){
      std::cout << "# SKIP k="<<k<<": stored peram missing ("<<h5old<<")" << std::endl;
      continue;
    }
    std::vector<double> Vo_re((std::size_t)Comp::Nt*nv_old*n2), Vo_im(Vo_re.size());
    std::vector<double> ev_old((std::size_t)Comp::Nt*nv_old);
    std::vector<double> To_re((std::size_t)nsrc*twin*twin*nv_old*nv_old), To_im(To_re.size());
    {
      HighFive::File fo(h5old, HighFive::File::ReadOnly);
      std::vector<int> nvo, two, nso, tlo;
      std::vector<double> tolo;
      fo.getDataSet("meta/Nv").read(nvo);
      fo.getDataSet("meta/twin").read(two);
      fo.getDataSet("meta/nsrc").read(nso);
      fo.getDataSet("meta/tsrc_list").read(tlo);
      fo.getDataSet("meta/overlap_tol").read(tolo);
      assert(nvo[0] == nv_old && "EXTEND: stored meta/Nv != --nv-old");
      assert(two[0] == twin && nso[0] == nsrc && "EXTEND: stored twin/nsrc differ from this run");
      for(int s=0; s<nsrc; s++) assert(tlo[s] == tsrc_list[s] && "EXTEND: stored tsrc_list differs");
      assert(std::abs(tolo[0] - tol) <= 1e-12*tol && "EXTEND: --tol must equal the stored meta/overlap_tol");
      fo.getDataSet("V/real").read_raw(Vo_re.data());
      fo.getDataSet("V/imag").read_raw(Vo_im.data());
      fo.getDataSet("evals").read_raw(ev_old.data());
      fo.getDataSet("peram/tau/real").read_raw(To_re.data());
      fo.getDataSet("peram/tau/imag").read_raw(To_im.data());
    }
    if(!free_field) U.read(str_lat);
    D.update(U);

    const auto t_k0 = std::chrono::steady_clock::now();

    // evals[t][k], V[t] as (n2 x nv) column matrix.  Stored flat for h5 below.
    std::vector<std::vector<double>> evals(Nt, std::vector<double>(nv, 0.0));
    std::vector<Eigen::MatrixXcd> Vt(Nt);
    std::vector<double> lam_cut(Nt, 0.0);      // EXTEND: lambda_{N_n} (first mode beyond the new cut), if N_n < 2N_s

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
      if(nv < n2) lam_cut[t] = lam(nv);
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
#if BASIS_SYM
        MhMg.noalias() = Mg.adjoint() * Mg + Mg * Mg.adjoint();   // symmetrized (match the V(t) operator)
#else
        MhMg.noalias() = Mg.adjoint() * Mg;
#endif
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

    // ================= EXTEND: basis = [stored old modes | fresh new modes] =================
    // Stored V (t, a, j) row-major -> Vs(t) (n2 x N_o).  Checks (per config): stored evals reproduced; projector onto
    // the old span reproduced (phase/rotation-invariant); gaps at the cuts.  Any failure -> NO output for this k.
    double max_dev_eval = 0.0, max_dev_proj = 0.0, max_new_overlap = 0.0, min_gap_old = 1.0e300, min_gap_new = 1.0e300;
    {
      for(int t=0; t<Nt; t++){
        Eigen::MatrixXcd Vs(n2, nv_old);
        for(int a=0; a<nv_old; a++){
          for(Idx j=0; j<n2; j++){
            const std::size_t q = ((std::size_t)t*nv_old + a)*n2 + j;
            Vs(j, a) = Complex(Vo_re[q], Vo_im[q]);
          }
          max_dev_eval = std::max(max_dev_eval, std::abs(ev_old[(std::size_t)t*nv_old + a] - evals[t][a]));
        }
        const Eigen::MatrixXcd Qo = Vt[t].leftCols(nv_old);
        const Eigen::MatrixXcd dP = Vs * Vs.adjoint() - Qo * Qo.adjoint();
        max_dev_proj = std::max(max_dev_proj, dP.cwiseAbs().maxCoeff());
        // new modes: remove any stored-span component, then re-orthonormalize (thin Householder QR)
        Eigen::MatrixXcd Qn = Vt[t].rightCols(nnew);
        max_new_overlap = std::max(max_new_overlap, (Vs.adjoint() * Qn).cwiseAbs().maxCoeff());
        Qn -= Vs * (Vs.adjoint() * Qn);
        Eigen::HouseholderQR<Eigen::MatrixXcd> qr(Qn);
        Eigen::MatrixXcd Qthin = qr.householderQ() * Eigen::MatrixXcd::Identity(n2, nnew);
        // keep each new column's phase/orientation close to the fresh eigenvector (QR may flip signs)
        for(int b=0; b<nnew; b++){
          const Complex ph = Qthin.col(b).dot(Qn.col(b));
          if(std::abs(ph) > 0.0) Qthin.col(b) *= ph / std::abs(ph);
        }
        Vt[t].leftCols(nv_old) = Vs;
        Vt[t].rightCols(nnew) = Qthin;
        min_gap_old = std::min(min_gap_old, evals[t][nv_old] - evals[t][nv_old-1]);
      }
      // gap at the NEW cut: lambda_{N_n} - lambda_{N_n - 1} (only defined for N_n < 2N_s)
      if(nv < n2){
        for(int t=0; t<Nt; t++) min_gap_new = std::min(min_gap_new, lam_cut[t] - evals[t][nv-1]);
      }
      // orthonormality of the combined basis
      double max_ortho2 = 0.0;
      for(int t=0; t<Nt; t++){
        Eigen::MatrixXcd G = Vt[t].adjoint() * Vt[t];
        for(int a=0; a<nv; a++) G(a,a) -= 1.0;
        max_ortho2 = std::max(max_ortho2, G.cwiseAbs().maxCoeff());
      }
      std::cout << "#   [EXT basis] max|eval_stored - eval_fresh| = " << max_dev_eval
                << " , max|P_old(stored) - P_old(fresh)| = " << max_dev_proj
                << " , max|V_old^dag V_new(fresh)| = " << max_new_overlap
                << " , combined max|V^dag V - 1| = " << max_ortho2 << std::endl;
      std::cout << "#   [EXT gaps] min_t (lambda_{N_o+1} - lambda_{N_o}) = " << min_gap_old
                << " , min_t (lambda_{N_n+1} - lambda_{N_n}) = " << min_gap_new << std::endl;
      if(max_dev_eval > 1.0e-9 || max_dev_proj > 1.0e-8 || max_ortho2 > 1.0e-10){
        std::cout << "# FAILED k="<<k<<": stored basis NOT reproduced (see [EXT basis]); NO output written" << std::endl;
        continue;
      }
    }

    // ================= EXTEND CHUNK 2: new-source solves + block assembly =================
    // Tnew[s][ap][a] (N_n x nnew) = V(t'=src0_s+ap)^dag D_ov^{-1} w_{N_o+b}(t=src0_s+a), b = 0..nnew-1.
    std::vector<std::vector<std::vector<Eigen::MatrixXcd>>> Tau_all(nsrc), Tnew(nsrc);
    for(int s=0; s<nsrc; s++){
      Tnew[s].assign(twin, std::vector<Eigen::MatrixXcd>(twin, Eigen::MatrixXcd::Zero(nv, nnew)));
    }
    const auto t_new0 = std::chrono::steady_clock::now();
    {
      Eigen::VectorXcd p(n2);
      if(fused && nsrc>1){
        // FUSED (as the stored run): ONE combined source per offset a, summed over the windows
        std::vector<const Eigen::MatrixXcd*> Vsrc(nsrc);
        std::vector<int> tsr(nsrc);
        for(int a=0; a<twin; a++){
          for(int l=nnew; l<NSTACK; l++){
            std::fill(hblk.begin()+(std::size_t)l*N, hblk.begin()+(std::size_t)(l+1)*N, Complex(0.0,0.0));
          }
          for(int bb=0; bb<nnew; bb++){
            for(int s=0; s<nsrc; s++){
              tsr[s]  = tsrc_list[s] + a;
              Vsrc[s] = &Vt[tsr[s]];
            }
            embed_mode_source_fused(src, Vsrc, tsr, nv_old + bb);         // sum_s w_{N_o+b}(src0_s + a)
            op_DH.from_cpu<N>(hblk.data() + (std::size_t)bb*N, src.field);
          }
          blk_Dsq.solve_sq_from_cpu(hblk.data(), hblk.data(), tol);       // col b = D_ov^{-1} (fused source)
          for(int bb=0; bb<nnew; bb++){
            for(Idx i=0; i<N; i++) psi.field[i] = hblk[(std::size_t)bb*N + i];
            for(int s=0; s<nsrc; s++){
              for(int ap=0; ap<twin; ap++){
                const int tp = tsrc_list[s] + ap;
                extract_slice(p, psi, tp);
                Tnew[s][ap][a].col(bb) = Vt[tp].adjoint() * p;
              }
            }
          }
        }
      }
      else{
        for(int s=0; s<nsrc; s++){
          const int src0 = tsrc_list[s];
          for(int a=0; a<twin; a++){
            const int t = src0 + a;
            for(int l=nnew; l<NSTACK; l++){
              std::fill(hblk.begin()+(std::size_t)l*N, hblk.begin()+(std::size_t)(l+1)*N, Complex(0.0,0.0));
            }
            for(int bb=0; bb<nnew; bb++){
              embed_mode_source(src, Vt[t], nv_old + bb, t);
              op_DH.from_cpu<N>(hblk.data() + (std::size_t)bb*N, src.field);
            }
            blk_Dsq.solve_sq_from_cpu(hblk.data(), hblk.data(), tol);
            for(int bb=0; bb<nnew; bb++){
              for(Idx i=0; i<N; i++) psi.field[i] = hblk[(std::size_t)bb*N + i];
              for(int ap=0; ap<twin; ap++){
                const int tp = src0 + ap;
                extract_slice(p, psi, tp);
                Tnew[s][ap][a].col(bb) = Vt[tp].adjoint() * p;
              }
            }
          }
        }
      }
    }
    const double secs_new = std::chrono::duration<double>(std::chrono::steady_clock::now()-t_new0).count();
    std::cout << "#   [EXT peram] new-source solves ("<<nnew<<" modes x "<<twin<<" offsets"<<(fused&&nsrc>1?", fused":"")
              <<") done ["<<secs_new<<" s]" << std::endl;

    // ---- assemble tau (N_n x N_n): old-old stored ; any-new from Tnew ; new-old by the identity ----
    double max_id_newnew = 0.0, max_id_oldold = 0.0;
    for(int s=0; s<nsrc; s++){
      Tau_all[s].assign(twin, std::vector<Eigen::MatrixXcd>(twin, Eigen::MatrixXcd::Zero(nv, nv)));
      for(int ap=0; ap<twin; ap++){
        for(int a=0; a<twin; a++){
          Eigen::MatrixXcd& T = Tau_all[s][ap][a];
          for(int i=0; i<nv_old; i++){
            for(int j=0; j<nv_old; j++){
              const std::size_t q = ((((std::size_t)s*twin + ap)*twin + a)*nv_old + i)*nv_old + j;
              T(i, j) = Complex(To_re[q], To_im[q]);
            }
          }
          T.rightCols(nnew) = Tnew[s][ap][a];
          for(int bb=0; bb<nnew; bb++){
            for(int j=0; j<nv_old; j++){
              T(nv_old + bb, j) = -std::conj(Tnew[s][a][ap](j, bb));       // tau_ij(t',t) = -conj(tau_ji(t,t')), i != j
            }
          }
        }
      }
    }
    // identity residuals on the blocks that are known BOTH ways (diagnostic of solver tolerance, every config)
    for(int s=0; s<nsrc; s++){
      for(int ap=0; ap<twin; ap++){
        for(int a=0; a<twin; a++){
          const Eigen::MatrixXcd& T  = Tau_all[s][ap][a];
          const Eigen::MatrixXcd& Tr = Tau_all[s][a][ap];
          Eigen::MatrixXcd R = T + Tr.adjoint();
          if(ap == a) for(int i=0; i<nv; i++) R(i, i) -= 1.0;
          max_id_newnew = std::max(max_id_newnew, R.bottomRightCorner(nnew, nnew).cwiseAbs().maxCoeff());
          max_id_oldold = std::max(max_id_oldold, R.topLeftCorner(nv_old, nv_old).cwiseAbs().maxCoeff());
        }
      }
    }
    std::cout << "#   [EXT identity] max|tau(t',t) + tau(t,t')^dag - delta| : new-new (solved) = " << max_id_newnew
              << " , old-old (stored) = " << max_id_oldold << "  (~ solver tol " << tol << " scale)" << std::endl;

    // first config: DIRECT check of the identity-filled new-old block (window 0, source offset a=0): N_o clean
    // single-source solves of the OLD modes at t = src0, projected on the NEW sink modes.
    if(k == k_lo){
      const int src0 = tsrc_list[0];
      Eigen::VectorXcd p(n2);
      double max_direct = 0.0, max_direct_rel = 0.0;
      for(int j=0; j<nv_old; j++){
        embed_mode_source(src, Vt[src0], j, src0);
        op_DH.from_cpu<N>(tmp.field, src.field);
        op_Dsq.solve<N>(psi.field, tmp.field, tol);                       // psi = D_ov^{-1} w_j(src0)
        for(int ap=0; ap<twin; ap++){
          const int tp = src0 + ap;
          extract_slice(p, psi, tp);
          const Eigen::VectorXcd direct = Vt[tp].rightCols(nnew).adjoint() * p;   // tau_{new, j}(t', src0)
          const Eigen::VectorXcd fill   = Tau_all[0][ap][0].col(j).tail(nnew);
          const double d = (direct - fill).cwiseAbs().maxCoeff();
          const double sc = std::max(direct.cwiseAbs().maxCoeff(), 1.0e-300);
          max_direct = std::max(max_direct, d);
          max_direct_rel = std::max(max_direct_rel, d / sc);
        }
      }
      std::cout << "#   [EXT direct] new-old block, identity vs direct solve (window 0, a=0): max abs = " << max_direct
                << " , max rel = " << max_direct_rel << std::endl;
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
      h5.createDataSet("meta/tsrc0",    std::vector<int>{tsrc_list[0]});   // base of window 0 (backward compat)
      h5.createDataSet("meta/twin",     std::vector<int>{twin});
      h5.createDataSet("meta/nsrc",     std::vector<int>{nsrc});           // number of source windows
      h5.createDataSet("meta/tsrc_list", tsrc_list);                       // base source timeslice per window
      h5.createDataSet("meta/method",   std::string("exact_distillation_wilson_basis_extended"));
      h5.createDataSet("meta/nv_old",   std::vector<int>{nv_old});
      h5.createDataSet("meta/extended_from", dir_old);
      h5.createDataSet("meta/identity_block", std::vector<int>{1});        // new-old block from (D_ov^{-1})^dag = 1 - D_ov^{-1}
      h5.createDataSet("meta/tau_gw_present", std::vector<int>{0});
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

      // peram: split re/im, row-major.  MULTI-SOURCE (change 2): with nsrc>1 the array carries a LEADING
      // window axis -> shape (nsrc, twin, twin, Nv, Nv), index [s, a_snk, a_src, k, l], with base source t of
      // window s given by meta/tsrc_list[s] (t' = tsrc_list[s] + a_snk, t = tsrc_list[s] + a_src).  With
      // nsrc==1 the OLD layout (twin, twin, Nv, Nv) is written verbatim so existing readers keep working.
      const std::size_t nWindow = (std::size_t)twin*twin*nv*nv;   // entries per window
      const std::size_t nT = (std::size_t)nsrc*nWindow;
      std::vector<double> Tre(nT), Tim(nT);
      std::size_t idx=0;
      for(int s=0; s<nsrc; s++){
        const std::vector<std::vector<Eigen::MatrixXcd>>& Tau  = Tau_all [s];
        for(int ap=0; ap<twin; ap++){
          for(int a=0; a<twin; a++){
            for(int kk=0; kk<nv; kk++){
              for(int ll=0; ll<nv; ll++){
                const Complex c  = Tau [ap][a](kk,ll);
                Tre[idx]  = c.real();
                Tim[idx]  = c.imag();
                idx++;
              }
            }
          }
        }
      }
      std::vector<std::size_t> pdims;
      if(nsrc == 1){
        pdims = { (std::size_t)twin, (std::size_t)twin, (std::size_t)nv, (std::size_t)nv };            // OLD layout
      }
      else{
        pdims = { (std::size_t)nsrc, (std::size_t)twin, (std::size_t)twin, (std::size_t)nv, (std::size_t)nv }; // + nsrc axis
      }
      h5.createDataSet<double>("peram/tau/real",    HighFive::DataSpace(pdims)).write_raw(Tre.data());
      h5.createDataSet<double>("peram/tau/imag",    HighFive::DataSpace(pdims)).write_raw(Tim.data());
      // tau_gw NOT written (extended perams are for the axial sp only; meta/tau_gw_present = 0)
    }
    std::filesystem::rename(h5tmp, h5path);
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now()-t_k0).count();
    std::cout << "# k="<<k<<(free_field?" (free field)":"")<<" basis+peram done (nsrc="<<nsrc<<") ["<<secs
              <<" s] -> "<<h5path<<std::endl;
  } // k

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();
  return 0;
}
