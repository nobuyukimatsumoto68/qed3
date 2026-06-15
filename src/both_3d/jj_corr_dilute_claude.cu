// jj_corr_dilute_claude.cu
// -----------------------------------------------------------------------------
// DILUTED + MASTER-FIELD variant of jj_corr_block_t_claude.cu.  In ONE pass sharing the diluted sink leg
// phi' = D_m^{-1} eta, it computes BOTH the exact conserved current K AND the local (ultralocal e.sigma)
// current:
//   exact K : connected tp + sp, disc, axial.   (ylm REMOVED -> local-ylm is a SEPARATE one-end-trick file.)
//   local   : connected s1,s2,s3 (vector) + LOCAL AXIAL s1,s2,s3 (Apm; t0-leg D_m^{-dag}) + local disc J_a(t).
// exact and local differ only in the INSERTION (K vs point-sigma); the local sink is a cheap local
// sigma-multiply on the SAME phi' (no extra solve), local adds only its own source legs.
//
// Variance reduction:
//   (1) spin x even/odd-time DILUTION = 4 patterns (eta.time_spin_dilution, t_block = Nt/2 => interval 2):
//       diluted estimator = sum over the 4 disjoint patterns (unbiased); removes the dominant source/sink
//       time cross-noise + spin cross-terms.  Ref: Foley/Juge/O'Cais/Peardon/Ryan/Skullerud, hep-lat/0505023.
//   (2) SOURCE-ORIGIN SUPERPOSITION = master-field: insert at t0 in {0, Nt/2} SIMULTANEOUSLY (the source is
//       their SUM, ONE solve per insertion); the measured correlator is C(t) = G(t) + G(t - Nt/2), FOLDED in
//       analysis (wrap-around accepted).  Both for exact and local.  Ref: Luscher arXiv:1707.09758 (orig);
//       Francis/Fritzsch/Luscher/Rago arXiv:1911.04533; our framing arXiv:2301.08696.
//
// Output: data_<ESNID>/corr_dil_nt0<N>_nhits<H>/corr.<k>.h<h>.h5 (one file per hit; origins superposed, the
//   origin set t0s stored for the fold; the RNG seed string stored as rng_seed).  Production
//   jj_corr_block_t_claude.cu is left untouched (A/B).  RNG seeded from a string via std::seed_seq
//   (seed_from_string), so a run is reproducible from the recorded rng_seed.
//
// SINK lever (inherited): t-BLOCKED K apply (ConservedCurrentBlockT::apply_k[_dag]_block_t) computes
// K(.,t) phi' for ALL t in one pass (~3x, the C6f lever).  Solver: overlap multishift
// (*_deviceAsyncLaunch_ms); kernel = multishift apply_k_ms via ConservedCurrent::operator().
// Conventions (v3): valence mass --mass-re/--mass-im (D_m = D_ov + m); --ens-dir omit => free field U = 1;
// parity (m purely imaginary) dagger leg uses \tilde D_{m_P} = D_ov + m_P/(1-m_P).
// Plan: jj_corr_dilute_impl_plan_claude.md ; estimator algebra: jj_local_stoch_estimator_design_claude.md.

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
  constexpr int NSTREAMS=4; // NPARALLEL_DUPDATE;
  constexpr int NPARALLEL_GAUGE=NPARALLEL_DUPDATE;
  constexpr int NPARALLEL_SORT=NPARALLEL_DUPDATE;

  constexpr int N_REFINE=1;
  constexpr int NS=2;
  constexpr int Nt=128;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;

  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;   // inner Zolotarev pole solves
  const double TOL_OUTER=1.0e-5;   // outer CG; above the ~1e-15 machine-precision residual floor
                                   // for the small-norm sink RHS (plenty for the correlator)
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

// #include "sparse_dirac.h"
#include "sparse_dirac_claude.h"   // O(len) bucketing CSR build (was O(N*len)); -DCSR_VERIFY to check
// #include "matpoly.h"
#include "matpoly_claude.h"

#include "overlap_wmass_claude.h"        // complex-mass overlap (massless at mass=0)
#include "blocked_mat_claude.h"      // C7: BlockedMat<N,NSTACK,Op> (mrhs block solves, host-block solve_sq_from_cpu)
#include "conserved_current_claude.h"   // ConservedCurrent: apply_k / apply_k_dag (multishift apply_k_ms)
#include "conserved_current_block_claude.h"  // C6f-c: ConservedCurrentBlockT::apply_k[_dag]_block_t (t-blocked sink)

//------------------------------------------
#include <getopt.h>

// Stable, portable string -> RNG seed (std::seed_seq mixes the chars deterministically; NOT std::hash,
// which is implementation-defined).  Lets a run be reproduced from the rng_seed string stored in the .h5.
static int seed_from_string(const std::string& s){
  std::seed_seq seq(s.begin(), s.end());
  std::uint32_t w;
  seq.generate(&w, &w + 1);          // one 32-bit word
  return static_cast<int>(w);
}

void PrintHelp(){
  printf("Usage: ./a.out [options]   (UNIFIED: computes disc + conn tp/sp/ylm, vector + axial)\n");
  printf("  --gsq <gsq>          Wilson coupling squared (ensemble id; default: 8.0)\n");
  printf("  --Nf <Nf>            number of fermion flavors (ensemble id; default: 2)\n");
  printf("  --nu0 <nu0>          sea quark asymmetry (ensemble id; default: 1.0)\n");
  printf("  --nu1 <nu1>          valence Wilson-Dirac asymmetry (operator; default: nu0)\n");
  printf("  --mass-re <x>        valence mass Re (default: 0.0)\n");
  printf("  --mass-im <y>        valence mass Im (default: 0.0)\n");
  printf("  --ens-dir <path>     sea config directory; OMIT => free field (U=1) check\n");
  printf("  --nhits <n>          stochastic hits (default: 1)\n");
  printf("  --n-t0 <N>           number of source-time origins t0=b*(Nt/N), b=0..N-1 (default: 2)\n");
  printf("  --stride <N>         ensemble config stride: measure ckpoint_lat.k for k=kmin,kmin+N,... (default: 10)\n");
  printf("  --kmin <k>           first ckpoint index, inclusive (default: 0)\n");
  printf("  --kmax <k>           one past last ckpoint index, EXCLUSIVE: kmin <= k < kmax (default: 1000000)\n");
  printf("  --spin-dilution      turn ON spin dilution (NS=2 spin classes); default OFF (1 spin class)\n");
  printf("  --time-dilution <td> number of even/interval time-dilution classes (default: 2); require (Nt/2)%%td==0\n");
  printf("                       => total dilution patterns = (spin-dilution? 2:1) x td.  spin-only = --spin-dilution --time-dilution 1\n");
  printf("  -h, --help           show this help\n");
  exit(0);
}

void ParseArgs(int argc, char* argv[],
               double& gsq, int& Nf, double& nu0, double& nu1,
               double& mass_re, double& mass_im,
               std::string& ens_dir, int& nhits, int& n_t0, int& stride,
               int& kmin, int& kmax, bool& spin_dilution, int& time_dilution){
  static struct option long_opts[] = {
    {"gsq",     required_argument, nullptr, 'g'},
    {"Nf",      required_argument, nullptr, 'N'},
    {"nu0",     required_argument, nullptr, 'n'},
    {"nu1",     required_argument, nullptr, 'm'},
    {"mass-re", required_argument, nullptr, 'r'},
    {"mass-im", required_argument, nullptr, 'i'},
    {"ens-dir", required_argument, nullptr, 'e'},
    {"nhits",   required_argument, nullptr, 'H'},
    {"n-t0",    required_argument, nullptr, 'T'},
    {"stride",  required_argument, nullptr, 'I'},
    {"kmin",    required_argument, nullptr, 'a'},
    {"kmax",    required_argument, nullptr, 'b'},
    {"spin-dilution", no_argument,       nullptr, 's'},
    {"time-dilution", required_argument, nullptr, 'd'},
    {"help",    no_argument,       nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
  };
  int opt, idx;
  while((opt = getopt_long(argc, argv, "g:N:n:m:r:i:e:H:T:I:a:b:sd:h", long_opts, &idx)) != -1){
    switch(opt){
    case 'g': gsq     = std::stod(optarg); break;
    case 'N': Nf      = std::stoi(optarg); break;
    case 'n': nu0     = std::stod(optarg); break;
    case 'm': nu1     = std::stod(optarg); break;
    case 'r': mass_re = std::stod(optarg); break;
    case 'i': mass_im = std::stod(optarg); break;
    case 'e': ens_dir = optarg; break;
    case 'H': nhits   = std::stoi(optarg); break;
    case 'T': n_t0    = std::stoi(optarg); break;
    case 'I': stride  = std::stoi(optarg); break;
    case 'a': kmin    = std::stoi(optarg); break;
    case 'b': kmax    = std::stoi(optarg); break;
    case 's': spin_dilution = true; break;
    case 'd': time_dilution = std::stoi(optarg); break;
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

  double gsq=8.0;  int Nf=2;  double nu0=1.0;  double nu1=-1.0;   // nu1<0 => default to nu0
  double mass_re=0.0, mass_im=0.0;
  std::string ens_dir="";     // empty => free-field mode
  int nhits=1;
  int n_t0=2;   // number of source-time origins (Sec. 3.7)
  int stride=10;   // ensemble config stride (k = kmin, kmin+stride, ...)
  int kmin=0;          // first ckpoint index, inclusive
  int kmax=1000000;    // one past last, exclusive: kmin <= k < kmax (existence-break ends it earlier)
  bool spin_dilution=false;  // --spin-dilution: NS=2 spin classes (else 1, no spin dilution)
  int  time_dilution=2;      // --time-dilution td: td even/interval time classes; require (Nt/2)%td==0

  ParseArgs(argc, argv, gsq, Nf, nu0, nu1, mass_re, mass_im, ens_dir, nhits, n_t0, stride, kmin, kmax,
            spin_dilution, time_dilution);
  assert(time_dilution >= 1 && (Comp::Nt/2) % time_dilution == 0
         && "--time-dilution td must satisfy (Nt/2) % td == 0");
  const int n_spin_classes = spin_dilution ? NS : 1;     // dilution patterns = n_spin_classes * time_dilution
  const int t_block_dil = Comp::Nt / time_dilution;      // time_spin/time_dilution interval = Nt/t_block = td
  if(nu1 < 0.0) nu1 = nu0;    // valence asymmetry defaults to the sea value nu0 (knob retained)

  const Complex valence_mass(mass_re, mass_im);
  const bool free_field = ens_dir.empty();

  // parity case: purely imaginary valence mass -> dagger leg uses \tilde D_{m_P}
  const bool parity = (std::abs(mass_im) > 1.0e-15) && (std::abs(mass_re) <= 1.0e-15);
  // flavor case: purely real valence mass -> axial uses massless D_ov legs (m_F formulas)
  const bool flavor = (std::abs(mass_re) > 1.0e-15) && (std::abs(mass_im) <= 1.0e-15);
  (void)flavor;  // DEAD FLAG: m_F now uses the D_m path everywhere (same as massless); kept only for the
                 // commented A/B reference lines in the axial block.  No live branch reads `flavor`.

  std::cout << "# gsq="<<gsq<<" Nf="<<Nf<<" nu0="<<nu0<<" nu1="<<nu1
            << " mass=("<<mass_re<<","<<mass_im<<")"
            << " ens_dir="<<(free_field?std::string("<free-field U=1>"):ens_dir)
            << " nhits="<<nhits << std::endl;

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();

  int device;
  CUDA_CHECK(cudaGetDeviceCount(&device));
  cudaDeviceProp device_prop[device];
  cudaGetDeviceProperties(&device_prop[0], 0);
  std::cout << "# dev = " << device_prop[0].name << std::endl;
  CUDA_CHECK(cudaSetDevice(0));

  constexpr Idx N  = Comp::N;
  constexpr int Nt = Comp::Nt;

  // source-time origins (Sec. 3.7): n_t0 evenly-spaced t0 = b*(Nt/n_t0); FULL dt-loop dt=0..Nt-1;
  // store every (hit, t0) raw (no in-program averaging); average + jackknife downstream.
  assert(n_t0 >= 1 && Nt % n_t0 == 0 && "Nt must be divisible by n_t0");
  const int t0_spacing = Nt / n_t0;
  std::vector<int> t0s(n_t0);
  for(int b=0; b<n_t0; b++) t0s[b] = b*t0_spacing;
  std::cout << "# n_t0=" << n_t0 << " source origins (t0=b*"<<t0_spacing<<"), full dt; one file/config" << std::endl;

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  using Rng=ParallelRngExt<Base,Nt>;
  using Fermion=OverlapWMass<WilsonDirac>;

  // ---- operators (grouped as in hmc_w_mass_claude.cu) -----------------------------
  Base base(Comp::N_REFINE);
  std::cout << "# lattice set." << std::endl;

  const double r  = 1.0;
  const double M5 = -1.0;
  const double at = 0.2;
  if(Nt!=1) assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu1);
  std::cout << "# DW set." << std::endl;

  Gauge U(base);                       // free field (U=1) unless a config is read below
  Rng rng(base, 1234);

  // Overlap operators:
  //   D    = D_ov              (massless; axial GW factors + flavor-axial legs)
  //   Dm   = D_ov + m          (Eq. 3.59; vector (++) both legs; tp/sp)
  //   Dtil = D_ov + m/(1-m)    (\tilde D_{m_P}, Eq. 3.62; parity (-) dagger leg)
  Fermion D   (DW, Complex(0.0), 11);                              // npole=11 (match the HMC; was 21)
  Fermion Dm  (DW, valence_mass, 11);
  Fermion Dtil(DW, valence_mass / (Complex(1.0) - valence_mass), 11);
  std::cout << "# overlap operators set: D_ov, D_m, tilde D_{m_P} (M5="<<M5<<")." << std::endl;

  ConservedCurrent<Fermion,Gauge> kop(Dm);   // K is mass-independent; multishift apply_k_ms via operator()
  MatPoly op_K; op_K.push_back(cplx(1.0), {&kop});   // apply K via op_K.from_cpu (source K^dag(.,t0) applies)
  // C6f-c: t-blocked sink kernel.  One wrapper (full pool ~189 MB + KBlockScratch ~19 MB), created once;
  // apply_k[_dag]_block_t computes K(t,fixed) sinkvec for ALL t in one pass (the SINK lever).
  ConservedCurrentBlockT<Fermion,Gauge,Comp::Nt> kblock(kop);

  // Uniform operator set, multishift (_ms) entry points (~4x).  For each overlap X:
  //   mult = X,  H = X^dag,  sq = X^dag X = X X^dag (D_ov+m normal -> fused DDH=DHD).
  //   X^{-1} b  = op_XH (RHS X^dag b) + op_Xsq (CG);  X^{-dag} b = op_X (RHS X b) + op_Xsq.
  // massless D_ov: op_oneMinusD = (1 - D_ov) GW factor (apply); op_DH/op_Dsq for the flavor-axial solve.
  auto f_D   = std::bind(&Fermion::mult_deviceAsyncLaunch_ms, &D, std::placeholders::_1, std::placeholders::_2);
  auto f_DH  = std::bind(&Fermion::adj_deviceAsyncLaunch_ms,  &D, std::placeholders::_1, std::placeholders::_2);
  auto f_Dsq = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms,  &D, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_D(f_D), M_DH(f_DH), M_Dsq(f_Dsq);
  MatPoly op_oneMinusD;
  op_oneMinusD.push_back(cplx( 1.0), {});       // identity term (empty product)
  op_oneMinusD.push_back(cplx(-1.0), {&M_D});   // - D_ov   => op_oneMinusD : v -> (1 - D_ov) v
  MatPoly op_DH;  op_DH.push_back(cplx(1.0), {&M_DH});
  MatPoly op_Dsq; op_Dsq.push_back(cplx(1.0), {&M_Dsq});

  auto f_Dm   = std::bind(&Fermion::mult_deviceAsyncLaunch_ms, &Dm, std::placeholders::_1, std::placeholders::_2);
  auto f_DmH  = std::bind(&Fermion::adj_deviceAsyncLaunch_ms,  &Dm, std::placeholders::_1, std::placeholders::_2);
  auto f_Dmsq = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms,  &Dm, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_Dm(f_Dm), M_DmH(f_DmH), M_Dmsq(f_Dmsq);
  MatPoly op_Dm;   op_Dm.push_back(cplx(1.0), {&M_Dm});
  MatPoly op_DmH;  op_DmH.push_back(cplx(1.0), {&M_DmH});
  MatPoly op_Dmsq; op_Dmsq.push_back(cplx(1.0), {&M_Dmsq});

  auto f_tilDm   = std::bind(&Fermion::mult_deviceAsyncLaunch_ms, &Dtil, std::placeholders::_1, std::placeholders::_2);
  auto f_tilDmH  = std::bind(&Fermion::adj_deviceAsyncLaunch_ms,  &Dtil, std::placeholders::_1, std::placeholders::_2);
  auto f_tilDmsq = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms,  &Dtil, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_tilDm(f_tilDm), M_tilDmH(f_tilDmH), M_tilDmsq(f_tilDmsq);
  MatPoly op_tilDm;   op_tilDm.push_back(cplx(1.0), {&M_tilDm});
  MatPoly op_tilDmH;  op_tilDmH.push_back(cplx(1.0), {&M_tilDmH});
  MatPoly op_tilDmsq; op_tilDmsq.push_back(cplx(1.0), {&M_tilDmsq});

  // ---- mrhs (C6e) block-apply ops: DDH(=op_*sq) on an NSTACK-wide block, used by
  // MatPoly::solve_block_cg in the connected SOURCE loops. tp width = N_SITES, sp width = N_LINKS.
  // These REPLACE the per-site/per-link single op_*sq.solve (left commented in the loops); the op_*sq
  // objects still serve the single phi'/ylm/tilphi solves.
  constexpr int NSTACK_TP = Comp::N_SITES;        // tp batch width
  constexpr int NSTACK_SP = Comp::N_LINKS;        // sp batch width
  // mrhs block engines (C7): one BlockedMat per (operator, width). Own their scratch (RAII),
  // hold a const ref to the fermion (per-config coeffs read at solve time), expose solve_sq_from_cpu
  // (host-block, hides device I/O). Created once; reused across hits/loops.
  BlockedMat<N,NSTACK_TP,Fermion> blk_tp_D(D),  blk_tp_Dm(Dm),  blk_tp_Dtil(Dtil);
  BlockedMat<N,NSTACK_SP,Fermion> blk_sp_D(D),  blk_sp_Dm(Dm),  blk_sp_Dtil(Dtil);

  // ---- geometry weights ------------------------------------------------------------------------
  const double inv4pi = 1.0/(4.0*std::acos(-1.0));
  const int n_sites = static_cast<int>(base.n_sites);
  const int n_links = static_cast<int>(base.links.size());

  // ---- mrhs (C6e) block scratch: batch the connected SOURCE solves over n_sites (tp) / n_links (sp).
  // The per-site/per-link op_*sq.solve loops become ONE block solve (MatPoly::solve_block_cg over
  // DDH_..._ms_block<NSTACK>). NSTACK is compile-time (Comp:: constexpr); cap = max(N_SITES, N_LINKS);
  // buffers in the PARENT-stream set. D/Dm/Dtil all n=11 -> npole = size-1 identical.
  assert(n_sites==NSTACK_TP && n_links==NSTACK_SP);   // NSTACK_TP/SP defined with the block-apply ops above
  const int nstack_cap = (NSTACK_TP > NSTACK_SP ? NSTACK_TP : NSTACK_SP);
  std::vector<Complex> hblk((size_t)Comp::N * nstack_cap);   // host RHS/solution block staging (in-place; device I/O is inside BlockedMat)

  // temporal projection (Eq. 4.32): w_tp[n] = A_n / kappa_t(n)^2   (kappa^2)
  std::vector<double> w_tp(base.n_sites);
  for(Idx n=0; n<base.n_sites; n++){ const double kt=DW.kappa_t[n]; w_tp[n]=base.dual_areas[n]/(kt*kt); }

  // spatial projection (Eq. 4.29): w_sp[il] = A_{nn'} / kappa^{(0)2}_{nn'} = link_volume[il]/kappa[il]^2
  std::vector<double> w_sp(base.n_links);
  for(Idx il=0; il<base.n_links; il++){ const double ks=DW.bd.kappa[il]; w_sp[il]=base.link_volume[il]/(ks*ks); }

  // LOCAL (ultralocal) site current weight: the bare \sigma carries no \kappa -> w_site = dual_areas
  // (Eq. 4.28/4.31; jj_local_deter_claude.cu).  NO 1/kappa^2.
  std::vector<double> w_site(base.n_sites);
  for(Idx n=0; n<base.n_sites; n++) w_site[n]=base.dual_areas[n];

  // ylm m-summed tower (Sec. 3.6): keep l_1=l_2=l, sum m_1,m_2 over -l..l.  ONE solve per l via the
  // COMBINED real weight W_ell[l][n] = (A_n / kappa_t(n)) * sum_{m=-l}^{l} Y_lm(n^)  (kappa^1, Eq. 4.36).
  constexpr int L_MAX = 2;
  // YLM STRIPPED (Chunk 1b): the Y_lm tower is a LOCAL-current observable computed in a SEPARATE
  // one-end-trick file (jj_local_ylm_stoch).  Setting N_ELL=0 makes EVERY ylm loop/array/output inert
  // here (the ylm code is left in place, just never executed).  Was: N_ELL = L_MAX + 1.
  constexpr int N_ELL = 0;
  std::vector<int> ls(N_ELL); for(int l=0;l<N_ELL;l++) ls[l]=l;
  std::vector<std::vector<double>> W_ell(N_ELL, std::vector<double>(n_sites, 0.0));
  for(int l=0; l<N_ELL; l++)
    for(int n=0; n<n_sites; n++){
      const double kt = DW.kappa_t[n];
      double s = 0.0;
      for(int m=-l; m<=l; m++) s += Ylm_real(l, m, base.sites[n]);   // sum over m at this current-site
      W_ell[l][n] = base.dual_areas[n] * s / kt;
    }

  // ---- output: data_<ESNID>/corr_nt0<N>_nhits<H>/corr.<config>.h5  (connected + folded disc) -------
  std::string ens_base = ens_dir;
  if(!ens_base.empty() && ens_base.back()=='/') ens_base.pop_back();
  { const auto slash = ens_base.find_last_of('/'); if(slash!=std::string::npos) ens_base = ens_base.substr(slash+1); }
  const std::string esnid = (free_field ? std::string("free") : ens_base)
                          + "_vmRe"+std::to_string(mass_re)+"vmIm"+std::to_string(mass_im);
  // dilution descriptor in the dir name so different schemes (spin on/off, td) do NOT collide.
  const std::string dir_out = "data_"+esnid+"/corr_dil_nt0"+std::to_string(n_t0)+"_nhits"+std::to_string(nhits)
                            + "_s"+std::to_string(spin_dilution?1:0)+"td"+std::to_string(time_dilution)+"/";
  std::filesystem::create_directories(dir_out);
  std::cout << "# dir_out = " << dir_out
            << "  (n_sites="<<n_sites<<", n_links="<<n_links<<", n_ell="<<N_ELL<<")" << std::endl;

  // ---- host buffers ----------------------------------------------------------------------------
  //   eta = shared Z_2 source; phi = phi'_++ = D_m^{-1} eta; phimm = phi'_-- = tilde D_{m_P}^{-dag} eta
  //   (parity sink leg); kphi = K(.,t) phi' (sink apply, shared tp+ylm); rho = K^dag-applied source;
  //   tmp = preconditioned CG RHS; srcL[l]/PhiL[l] = ylm m-summed source/sink towers.
  // CONNECTED-only file (disc = standalone jj_disc_claude.cu; see plan Sec. 3.8 solve-sharing note).
  // Connected source legs held across the shared sink pass (one per (insertion,t0)); std::vector is safe
  // now that FermionVector has a move ctor (valence_claude.h).  Indexers ITP/IYL flatten (insertion,b).
  FermionVector eta, phi, phimm, kphi, rho, tmp, chi, tilphi;   // tilphi = tilde D_{m_P}^{-1} eta (parity disc)
  FermionVector eta_A;   // area-weighted source eta_A(t,n,i) = dual_areas[n] eta(t,n,i) (condensate trace)
  std::array<FermionVector, N_ELL> srcL, PhiL;
  // SUPERPOSITION (master-field): one source leg PER INSERTION, summed over the n_t0 origins
  //   psi_tp[n] = D_m^{-dag} ( sum_b K^dag(n,t0_b) eta ).  Origins are SUPERPOSED in the source (one solve
  // per insertion), so the per-(n,b) dimension is gone; ITP/ISP ignore b (call sites unchanged).
  std::vector<FermionVector> psi_tp(n_sites);          // tp source leg  psi_tp[n]
  std::vector<FermionVector> psi_yl(N_ELL);            // (ylm stripped; N_ELL=0)
  std::vector<FermionVector> psi_sp(n_links);          // sp source leg  psi_sp[a]
  // Reuse A: cache the K^dag-applied (summed-origin) sources rho for the axial passes (1-D_ov)*rho.
  std::vector<FermionVector> rho_tp(n_sites);          // = sum_b K^dag(n ,t0_b) eta  (sites)
  std::vector<FermionVector> rho_sp(n_links);          // = sum_b K^dag(lk,t0_b) eta  (links)
  auto ITP = [&](int n,int b){ (void)b; return n; };   // origins SUPERPOSED -> b ignored
  auto IYL = [&](int l,int b){ (void)b; return l; };
  auto ISP = [&](int a,int b){ (void)b; return a; };

  // C6f-c sink block-t buffers:  d_sinkvec = the sink vector (phi'/phimm/chi/tilphi) uploaded once per
  // pass;  d_kphi_block = K(t,fixed) sinkvec for all t (N*Nt device, apply_k_block_t output);  kblk =
  // its host copy (D2H once per insertion);  PhiLt[l*Nt+t] = ylm m-summed sink tower per timeslice (was
  // PhiL[l]; must hold all t across the outer-n loop).  IPL flattens (l,t).
  CuC *d_sinkvec=nullptr, *d_kphi_block=nullptr;
  CUDA_CHECK(cudaMalloc(&d_sinkvec,    (size_t)N*CD));
  CUDA_CHECK(cudaMalloc(&d_kphi_block, (size_t)N*Nt*CD));
  std::vector<Complex> kblk((size_t)N*Nt);
  std::vector<FermionVector> PhiLt(N_ELL*Nt);     // ylm sink tower per (l,t)
  auto IPL = [&](int l,int t){ return l*Nt + t; };
  // copy column t of the K-block (host kblk) into kphi.field, so the existing per-(n,t) accumulation
  // lines (eta.dag(kphi), psi.dag(kphi), PhiLt += w*kphi) run UNCHANGED.
  auto kcol = [&](int t){ memcpy(kphi.field, kblk.data()+(size_t)t*N, (size_t)N*sizeof(Complex)); };

  // helpers: write a length-Nt complex correlator (1/4pi folded) under <key>/{real,imag};
  // write_corr_conj writes the elementwise conjugate (massless / m_F: Vmm = conj(Vpp)).
  auto write_corr = [&](HighFive::File& h5, const std::string& key, const std::vector<Complex>& C){
    std::vector<double> re(Nt), im(Nt);
    for(int t=0;t<Nt;t++){ const Complex g=inv4pi*C[t]; re[t]=g.real(); im[t]=g.imag(); }
    h5.createDataSet(key+"/real", re);  h5.createDataSet(key+"/imag", im);
  };
  auto write_corr_conj = [&](HighFive::File& h5, const std::string& key, const std::vector<Complex>& C){
    std::vector<double> re(Nt), im(Nt);
    for(int t=0;t<Nt;t++){ const Complex g=inv4pi*C[t]; re[t]=g.real(); im[t]=-g.imag(); }
    h5.createDataSet(key+"/real", re);  h5.createDataSet(key+"/imag", im);
  };
  // disconnected single-current traces J(t) are RAW (no 1/4pi fold), matching jj_disc_claude.cu.
  auto write_vec = [&](HighFive::File& h5, const std::string& key, const std::vector<Complex>& C){
    std::vector<double> re(Nt), im(Nt);
    for(int t=0;t<Nt;t++){ re[t]=C[t].real(); im[t]=C[t].imag(); }
    h5.createDataSet(key+"/real", re);  h5.createDataSet(key+"/imag", im);
  };

  // free field: single deterministic config (k=0), U=1.  ensemble: loop ckpoint_lat.k in ens_dir
  // with config stride (--stride; default 10).
  const int k_ckpoint = free_field ? 1 : stride;
  const int k_lo      = free_field ? 0 : kmin;   // inclusive
  const int k_hi      = free_field ? 1 : kmax;   // exclusive: kmin <= k < kmax

  for(int k = k_lo; k < k_hi; k += k_ckpoint){
    if(!free_field){
      const std::string str_lat = ens_dir + "ckpoint_lat." + std::to_string(k);
      if(!std::filesystem::exists(str_lat)){ if(k==0) continue; else break; }
      U.read(str_lat);
    }
    // ONE FILE PER HIT (so an interrupted run keeps the finished hits): corr.<k>.h<h>.h5, each atomic
    // (.tmp + rename) with its own "complete" sentinel; data under h0/.  Resume = skip completed hits.
    D.update(U);  Dm.update(U);  Dtil.update(U);
    std::cout << "# k="<<k<<(free_field?" (free field)":"")
              << "  lambda_min/max="<<Dm.lambda_min<<"/"<<Dm.lambda_max<<std::endl;

    for(int h=0; h<nhits; h++){
      const std::string h5path_h = dir_out + "corr." + std::to_string(k) + ".h" + std::to_string(h) + ".h5";
      if(std::filesystem::exists(h5path_h)){
        bool c=false; try { HighFive::File f(h5path_h,HighFive::File::ReadOnly); c=f.exist("complete"); } catch(...) {}
        if(c){ std::cout<<"# skip k="<<k<<" hit "<<h<<" (complete)"<<std::endl; continue; }
      }
      // RNG seed from a STRING (reproducible; stored as rng_seed).  Seed per (config, hit); the 4 dilution
      // patterns advance the same stream.  std::seed_seq via seed_from_string (no rng.h change).
      const std::string seed_str = esnid + "_k" + std::to_string(k) + "_h" + std::to_string(h);
      rng.reseed(seed_from_string(seed_str));
      const auto t_hit0 = std::chrono::steady_clock::now();
      auto elapsed = [&](){ return std::chrono::duration<double>(std::chrono::steady_clock::now()-t_hit0).count(); };
      std::cout << "# k="<<k<<" hit "<<(h+1)<<"/"<<nhits<<"  (n_t0="<<n_t0<<", n_sites="<<n_sites
                << ", n_links="<<n_links<<", seed='"<<seed_str<<"')" << std::endl;
      const std::string hp = "h0/";                  // one hit per file => always h0
      const std::string h5tmp = h5path_h + ".tmp";
      auto h5p = std::make_unique<HighFive::File>(h5tmp,
                   HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
      HighFive::File& h5 = *h5p;
      h5.createDataSet("t0s",   t0s);
      h5.createDataSet("n_t0",  std::vector<int>{n_t0});
      h5.createDataSet("nhits", std::vector<int>{nhits});
      h5.createDataSet("hit",   std::vector<int>{h});
      h5.createDataSet("rng_seed", seed_str);
      h5.createDataSet("spin_dilution", std::vector<int>{spin_dilution?1:0});
      h5.createDataSet("time_dilution", std::vector<int>{time_dilution});
      h5.createDataSet("n_patterns",    std::vector<int>{n_spin_classes*time_dilution});

      // HIT-SCOPE accumulators, SUMMED over the n_spin_classes*time_dilution dilution patterns (and over insertions/origins):
      //   conn (superposed, absolute t): Ctp/Csp (vector), Atp/Asp (axial);  disc (RAW): Jtp/Jsp.
      std::vector<Complex> Ctp(Nt,Complex(0,0)), Csp(Nt,Complex(0,0));
      std::vector<Complex> Atp(Nt,Complex(0,0)), Asp(Nt,Complex(0,0));
      std::vector<Complex> Jtp(Nt,Complex(0,0)), Jsp(Nt,Complex(0,0));
      // LOCAL current accumulators: conn Cs[c-1], disc Js[c-1] for channels c=1,2,3 (summed over patterns).
      std::vector<std::vector<Complex>> Cs(3, std::vector<Complex>(Nt,Complex(0,0)));
      std::vector<std::vector<Complex>> Js(3, std::vector<Complex>(Nt,Complex(0,0)));
      // LOCAL AXIAL conn accumulator CsA[c-1] (same channels; t0-leg D_m^{-dag} -> single complex Apm).
      std::vector<std::vector<Complex>> CsA(3, std::vector<Complex>(Nt,Complex(0,0)));
      // PARITY (m_P) "-" channels: SUPERPOSED + hit-scope accumulated, written ONCE after the loop (mirror
      // the (++) frame; audit I1/I2).  Ctp_mm/Csp_mm = exact-K V-- (Eq 3.65, NO (1+m_P) factor: kernel-
      // cancelled); JtpT_til/JspT_til = disc tilde traces; Cs_mm[c] = local vector V-- ((1+m_P)^{-2} at output).
      std::vector<Complex> Ctp_mm(Nt,Complex(0,0)), Csp_mm(Nt,Complex(0,0));
      std::vector<Complex> JtpT_til(Nt,Complex(0,0)), JspT_til(Nt,Complex(0,0));
      std::vector<std::vector<Complex>> Cs_mm(3, std::vector<Complex>(Nt,Complex(0,0)));
      // CONDENSATE traces (Eqs. 1.23 PS, 1.55 FS): single-propagator traces, reuse the dilute solves.
      // The forward leg xi -> eta^dag is D_m^{-1} in ALL mass cases (Eqs. 3.50, 3.57, 3.60); the
      // backward leg eta -> xi^dag is its dagger for massless/m_F (3.51, 3.58) but the tilde-D propagator
      // (1+m_P)^{-1} tilde D_{m_P}^{-dag} for m_P (3.61).  The spacetime average uses the site-area weights
      // A_n=dual_areas, so every trace is AREA-WEIGHTED, tr[A M]=sum_{n,t} A_n M((n,t),(n,t)), estimated with
      // the area-weighted source eta_A = A eta (A real diagonal -> eta_A^dag = eta^dag A).  Per pattern:
      //   acc_etadag_phi          = sum_d eta_A^dag phi'_d              (phi' = D_m^{-1} eta)  -> tr[A D_m^{-1}]
      //   acc_etadag_1mD_phi      = sum_d eta_A^dag (1-D_ov) phi'_d     (massless/m_F FS leg)
      //   acc_etadag_phimm        = sum_d eta_A^dag phimm_d            (phimm = tilde D_{m_P}^{-dag} eta; m_P only)
      //   acc_etadag_1mDdag_phimm = sum_d eta_A^dag (1-D_ov^dag) phimm_d (m_P FS leg)
      // (full area-weighted spacetime-summed traces; patterns partition the source -> sum is the FULL trace.)
      Complex acc_etadag_phi(0,0);
      Complex acc_etadag_1mD_phi(0,0);
      Complex acc_etadag_phimm(0,0);
      Complex acc_etadag_1mDdag_phimm(0,0);

      // ===== DILUTION: n_spin_classes (spin) x time_dilution (interval time classes) patterns.
      // --spin-dilution => NS spin classes (single spin component each); else 1 class (both spins).
      // --time-dilution td => td classes t_s=0..td-1, slices t == t_s (mod td) (t_block=Nt/td, interval=td).
      // Each pattern is its own diluted source; the estimator is the SUM over the disjoint patterns. =====
      for(int t_s=0; t_s<time_dilution; t_s++) for(int sc=0; sc<n_spin_classes; sc++){
        if(spin_dilution) eta.time_spin_dilution(rng, t_s, t_block_dil, sc);  // single spin component sc
        else              eta.time_dilution(rng, t_s, t_block_dil);           // both spins (no spin dilution)
        std::cout << "#   [dilution pattern t_s="<<t_s<<"/"<<time_dilution
                  << (spin_dilution ? (" spin="+std::to_string(sc)) : std::string(" (both spins)")) <<"]" << std::endl;

        // shared forward leg phi' = D_m^{-1} eta (reused by tp/sp/local/disc as the sink leg K(.,t)phi').
        std::cout << "#   phi' = D_m^{-1} eta : solving ..." << std::flush;
        op_DmH.from_cpu<N>(tmp.field, eta.field);
        op_Dmsq.solve<N>(phi.field, tmp.field, Comp::TOL_OUTER);
        std::cout << " done ["<<elapsed()<<" s]" << std::endl;

        // CONDENSATE: spacetime average uses the site-area weights A_n=dual_areas, so the trace is the
        // AREA-WEIGHTED tr[A M] = sum_{n,t} A_n <sigma(n,t)>, NOT the uniform trace.  Build the area-weighted
        // source eta_A(t,n,i) = A_n eta(t,n,i) (A diagonal, broadcast over t and spin); then eta_A.dag(M eta)
        // estimates tr[A M] (A real -> eta_A^dag = (A eta)^dag = eta^dag A).  tmp free here (held D_m^dag eta).
        for(int tt=0; tt<Comp::Nt; tt++)
          for(int n=0; n<n_sites; n++){
            const double a = w_site[n];
            for(int i=0; i<NS; i++) eta_A(tt,(Idx)n,i) = a*eta(tt,(Idx)n,i);
          }
        acc_etadag_phi += eta_A.dag(phi);                  // tr[A D_m^{-1}]  (etadag_xi leg, ALL mass cases)
        // FS backward leg for massless/m_F (xi^dag eta = (D_m^{-1})^dag; reuse phi' via trace cyclicity,
        // hence the conj at output).  For m_P the leg uses tilde D (phimm) -- handled in the parity block.
        if(!parity){
          op_oneMinusD.from_cpu<N>(tmp.field, phi.field);   // tmp = (1 - D_ov) phi'
          acc_etadag_1mD_phi += eta_A.dag(tmp);             // tr[A (1-D_ov) D_m^{-1}]
        }

      // ============ CONNECTED VECTOR -- temporal tp + ylm (shared phi' + shared K(n,t)phi' pass) =======
      // (++) source legs at all origins (the per-site K^dag(n,t0)eta apply is shared by tp and ylm):
      //   tp:  psi_tp[ITP(n,b)] = D_m^{-dag} K^dag(n,t0_b) eta
      //   ylm: psi_yl[IYL(l,b)] = D_m^{-dag} sum_n W_ell[l][n] K^dag(n,t0_b) eta   (m-summed tower)
      // SUPERPOSED source (master-field): rho_tp[n] = sum_b K^dag(n,t0_b) eta (summed over origins), then
      // ONE mrhs block solve psi_tp = D_m^{-dag} rho_tp.  (ylm stripped: N_ELL=0.)
      std::cout << "#   [vec tp ++] source solves ("<<n_sites<<" sites, "<<n_t0<<" origins SUPERPOSED) ..." << std::endl;
      for(int n=0;n<n_sites;n++){
        memset(rho_tp[n].field, 0, Comp::N*CD);
        for(int b=0;b<n_t0;b++){
          kop.set_temporal(U, t0s[b], (Idx)n, /*dag=*/true);
          op_K.from_cpu<N>(rho.field, eta.field);                    // rho = K^dag(n,t0_b) eta
          for(Idx i=0;i<N;i++) rho_tp[n].field[i]+=rho.field[i];     // SUM over origins
        }
        op_Dm.from_cpu<N>(hblk.data() + (size_t)n*N, rho_tp[n].field);   // RHS block col n = D_m (summed K^dag eta)
      }
      blk_tp_Dm.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);  // psi block = D_m^{-dag} (summed K^dag eta)
      for(int n=0;n<n_sites;n++) for(Idx i=0;i<N;i++) psi_tp[n].field[i]=hblk[(size_t)n*N+i];
      std::cout << "#     tp source done ["<<elapsed()<<" s]" << std::endl;
      // (++) sink pass: kphi = K(n,t) phi' (block-t) ONCE per (n,t); SUPERPOSED conn C(t)=G(t)+G(t-Nt/2)
      // accumulated at ABSOLUTE t (psi_tp is the summed-origin source); disc Jtp rides this apply.
      std::cout << "#   [vec tp ++] sink pass (block-t: "<<n_sites<<" applies x "<<Nt<<" t) ..." << std::flush;
      {
        CUDA_CHECK(cudaMemcpy(d_sinkvec, reinterpret_cast<CuC*>(phi.field), N*CD, H2D));
        for(int n=0;n<n_sites;n++){
          kblock.apply_k_block_t(d_kphi_block, d_sinkvec, U, (Idx)n);     // K(n,t) phi' for ALL t
          CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(kblk.data()), d_kphi_block, (size_t)N*Nt*CD, D2H));
          for(int t=0;t<Nt;t++){
            kcol(t);                                                     // kphi = column t = K(n,t) phi'
            Jtp[t] += w_tp[n]*eta.dag(kphi);                             // disc tp (rides this apply; summed over patterns)
            Ctp[t] += w_tp[n]*psi_tp[n].dag(kphi);                       // conn tp (superposed; absolute t; summed over patterns)
          }
        }
      }
      std::cout << " done ["<<elapsed()<<" s]" << std::endl;

      // (--) parity: independent tilde legs (operator-adjoint mirror of (++)).  Writes Vmm only.
      if(parity){
        std::cout << "#   [vec tp --] tilde source+sink (superposed) ..." << std::flush;
        op_tilDm.from_cpu<N>(tmp.field, eta.field);                  // tmp = tilde eta
        op_tilDmsq.solve<N>(phimm.field, tmp.field, Comp::TOL_OUTER);// phimm = tilde D_{m_P}^{-dag} eta

        // CONDENSATE m_P backward leg (Eq. 3.61): xi^dag eta = (1+m_P)^{-1} tilde D_{m_P}^{-dag}.
        // Reuse phimm = tilde D_{m_P}^{-dag} eta (just solved); (1-D_ov^dag) phimm = phimm - D_ov^dag phimm.
        // (1+m_P)^{-1} factor applied once at output.  tmp free here (held tilde eta).
        // area-weighted (eta_A built in the condensate block above, same pattern): tr[A (...) tilde D^{-dag}]
        acc_etadag_phimm += eta_A.dag(phimm);
        op_DH.from_cpu<N>(tmp.field, phimm.field);                       // tmp = D_ov^dag phimm
        acc_etadag_1mDdag_phimm += eta_A.dag(phimm) - eta_A.dag(tmp);    // eta^dag A (1-D_ov^dag) phimm

        // exact V-- (Eq. 3.65): SUMMED-origin source rho = sum_b K(n,t0_b) eta (dag=false); psi = tilde D^{-1} rho
        // (op_tilDmH + blk_tp_Dtil); sink K^dag(n,t) phimm -> Ctp_mm[t] (absolute t, superposed, summed over
        // patterns; mirrors the (++) frame).  rho/tmp scratch only -> rho_tp preserved for the axial.  NO factor
        // (exact-K kernel correction cancels (1+m_P)^{-1}).  Written ONCE after the dilution loop.
        for(int n=0;n<n_sites;n++){
          memset(rho.field, 0, Comp::N*CD);
          for(int b=0;b<n_t0;b++){
            kop.set_temporal(U, t0s[b], (Idx)n, /*dag=*/false);
            op_K.from_cpu<N>(tmp.field, eta.field);                  // tmp = K(n,t0_b) eta
            for(Idx i=0;i<N;i++) rho.field[i]+=tmp.field[i];         // SUM over origins
          }
          op_tilDmH.from_cpu<N>(hblk.data() + (size_t)n*N, rho.field);  // RHS block col n = tilde^dag (summed K eta)
        }
        blk_tp_Dtil.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);  // psi = tilde D^{-1}(summed K eta)
        for(int n=0;n<n_sites;n++) for(Idx i=0;i<N;i++) psi_tp[n].field[i]=hblk[(size_t)n*N+i];
        CUDA_CHECK(cudaMemcpy(d_sinkvec, reinterpret_cast<CuC*>(phimm.field), N*CD, H2D));
        for(int n=0;n<n_sites;n++){
          kblock.apply_k_dag_block_t(d_kphi_block, d_sinkvec, U, (Idx)n);  // K^dag(n,t) phimm for ALL t
          CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(kblk.data()), d_kphi_block, (size_t)N*Nt*CD, D2H));
          // disc V-- RIDES this sink (cyclicity: tr[tilde D^{-dag} K^dag] = tr[K^dag tilde D^{-dag}] =
          //   eta^dag K^dag(n,t) phimm; phimm = tilde D^{-dag} eta is this sink leg) -> NO separate disc sweep.
          for(int t=0;t<Nt;t++){ kcol(t); Ctp_mm[t] += w_tp[n]*psi_tp[n].dag(kphi); JtpT_til[t] += w_tp[n]*eta.dag(kphi); }  // conn V-- + disc -- (rides)
        }
        std::cout << " done ["<<elapsed()<<" s]" << std::endl;
      }

      // ============ CONNECTED VECTOR -- spatial sp (own K(l,t)phi' pass over links) =================
      // SUPERPOSED source: rho_sp[a] = sum_b K^dag(lk,t0_b) eta; ONE mrhs block solve psi_sp = D_m^{-dag} rho_sp.
      std::cout << "#   [vec sp ++] source ("<<n_links<<" links, "<<n_t0<<" origins SUPERPOSED) + sink ..." << std::flush;
      for(int a=0;a<n_links;a++){
        const BaseLink lk = base.links[a];
        memset(rho_sp[a].field, 0, Comp::N*CD);
        for(int b=0;b<n_t0;b++){
          kop.set_spatial(U, t0s[b], lk, /*dag=*/true);
          op_K.from_cpu<N>(rho.field, eta.field);                    // rho = K^dag(lk,t0_b) eta
          for(Idx i=0;i<N;i++) rho_sp[a].field[i]+=rho.field[i];     // SUM over origins
        }
        op_Dm.from_cpu<N>(hblk.data() + (size_t)a*N, rho_sp[a].field);  // RHS block col a = D_m (summed K^dag eta)
      }
      blk_sp_Dm.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);  // psi block = D_m^{-dag} (summed K^dag eta)
      for(int a=0;a<n_links;a++) for(Idx i=0;i<N;i++) psi_sp[a].field[i]=hblk[(size_t)a*N+i];
      // (++) sink pass: kphi = K(lk,t) phi' once per (a,t); SUPERPOSED conn at ABSOLUTE t.
      {
        CUDA_CHECK(cudaMemcpy(d_sinkvec, reinterpret_cast<CuC*>(phi.field), N*CD, H2D));
        for(int a=0;a<n_links;a++){
          const BaseLink lk = base.links[a];
          const Idx il = base.map2il.at(lk);
          kblock.apply_k_block_t(d_kphi_block, d_sinkvec, U, lk);    // K(lk,t) phi' for ALL t
          CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(kblk.data()), d_kphi_block, (size_t)N*Nt*CD, D2H));
          for(int t=0;t<Nt;t++){
            kcol(t);                                                // kphi = column t = K(lk,t) phi'
            Jsp[t] += w_sp[il]*eta.dag(kphi);                       // disc sp (rides this apply; summed over patterns)
            Csp[t] += w_sp[il]*psi_sp[a].dag(kphi);                 // conn sp (superposed; absolute t; summed over patterns)
          }
        }
      }
      std::cout << " done ["<<elapsed()<<" s]" << std::endl;
      // (--) parity sp: tilde mirror.  phimm = tilde^{-dag} eta reused from the temporal (--) block above.
      if(parity){
        std::cout << "#   [vec sp --] tilde source+sink (superposed) ..." << std::flush;
        // exact V-- sp (Eq. 3.65): SUMMED-origin source rho = sum_b K(lk,t0_b) eta; psi = tilde D^{-1} rho;
        // sink K^dag(lk,t) phimm -> Csp_mm[t] (absolute t, superposed).  rho/tmp scratch (rho_sp preserved).
        for(int a=0;a<n_links;a++){
          const BaseLink lk = base.links[a];
          memset(rho.field, 0, Comp::N*CD);
          for(int b=0;b<n_t0;b++){
            kop.set_spatial(U, t0s[b], lk, /*dag=*/false);
            op_K.from_cpu<N>(tmp.field, eta.field);                  // tmp = K(lk,t0_b) eta
            for(Idx i=0;i<N;i++) rho.field[i]+=tmp.field[i];         // SUM over origins
          }
          op_tilDmH.from_cpu<N>(hblk.data() + (size_t)a*N, rho.field);  // RHS block col a = tilde^dag (summed K eta)
        }
        blk_sp_Dtil.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);  // psi = tilde D^{-1}(summed K eta)
        for(int a=0;a<n_links;a++) for(Idx i=0;i<N;i++) psi_sp[a].field[i]=hblk[(size_t)a*N+i];
        CUDA_CHECK(cudaMemcpy(d_sinkvec, reinterpret_cast<CuC*>(phimm.field), N*CD, H2D));
        for(int a=0;a<n_links;a++){
          const BaseLink lk = base.links[a];
          const Idx il = base.map2il.at(lk);
          kblock.apply_k_dag_block_t(d_kphi_block, d_sinkvec, U, lk);  // K^dag(lk,t) phimm for ALL t
          CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(kblk.data()), d_kphi_block, (size_t)N*Nt*CD, D2H));
          for(int t=0;t<Nt;t++){ kcol(t); Csp_mm[t] += w_sp[il]*psi_sp[a].dag(kphi); JspT_til[t] += w_sp[il]*eta.dag(kphi); }  // conn V-- + disc -- (rides)
        }
        std::cout << " done ["<<elapsed()<<" s]" << std::endl;
      }

      // DISCONNECTED Jtp/Jsp rode the (++) sink applies above (RAW; summed over the dilution patterns).
      // WRITTEN ONCE per hit AFTER the dilution loop (below) -- writing here would collide across patterns.
      // parity disc -- (JtpT_til/JspT_til) NOW RIDES the connected (--) sinks above (cyclicity:
      //   tr[tilde D^{-dag} K^dag] = tr[K^dag tilde D^{-dag}] = eta^dag K^dag(n,t) phimm), so the standalone
      //   tilphi solve + 2 K-apply sweeps below are REDUNDANT (~20% of the m_P cost) -- commented out.
      //   OLD standalone version kept for A/B (uses the other ordering K(n,t) tilphi; same trace, more cost):
      // if(parity){
      //   std::cout << "#   [disc --] tilde trace (tilphi solve + sink) ..." << std::flush;
      //   op_tilDmH.from_cpu<N>(tmp.field, eta.field);
      //   op_tilDmsq.solve<N>(tilphi.field, tmp.field, Comp::TOL_OUTER);   // tilphi = tilde D_{m_P}^{-1} eta
      //   CUDA_CHECK(cudaMemcpy(d_sinkvec, reinterpret_cast<CuC*>(tilphi.field), N*CD, H2D));
      //   for(int n=0;n<n_sites;n++){
      //     kblock.apply_k_block_t(d_kphi_block, d_sinkvec, U, (Idx)n);  // K(n,t) tilphi for ALL t
      //     CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(kblk.data()), d_kphi_block, (size_t)N*Nt*CD, D2H));
      //     for(int t=0;t<Nt;t++){ kcol(t); JtpT_til[t] += w_tp[n]*kphi.dag(eta); }  // (K tilphi)^dag eta
      //   }
      //   for(int a=0;a<n_links;a++){
      //     const BaseLink lk = base.links[a]; const Idx il = base.map2il.at(lk);
      //     kblock.apply_k_block_t(d_kphi_block, d_sinkvec, U, lk);      // K(lk,t) tilphi for ALL t
      //     CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(kblk.data()), d_kphi_block, (size_t)N*Nt*CD, D2H));
      //     for(int t=0;t<Nt;t++){ kcol(t); JspT_til[t] += w_sp[il]*kphi.dag(eta); }
      //   }
      //   std::cout << " done ["<<elapsed()<<" s]" << std::endl;
      // }

      // ============ LOCAL (ultralocal) current s1,s2,s3 (shares phi'; local sigma-sink, no extra K solve) ==
      // SUPERPOSED localized source chi[c][n] = D_m^{-dag} ( sum_b sigma_c(n,t0_b) eta ); local sink
      // sigphi_c = sigma_c phi'; conn Cs_c(t) += w_site[n] chi[c][n](t,n)^dag sigphi_c(t,n) (absolute t,
      // superposed); disc Js_c(t) += w_site[n] eta(t,n)^dag sigphi_c(t,n).  Diagonal-in-site (the localized
      // source isolates site n in expectation).  MUST run before axial (which overwrites phi/chi/rho/psi_tp).
      // (jj_local_stoch_estimator_design_claude.md.)
      {
        std::cout << "#   [local s1,s2,s3] ("<<n_sites<<" sites, superposed) ..." << std::flush;
        for(int c=1;c<=3;c++){
          // localized summed-origin source per site n, then ONE mrhs block solve -> psi_tp[n] = chi[c][n].
          for(int n=0;n<n_sites;n++){
            memset(rho.field, 0, Comp::N*CD);
            for(int b=0;b<n_t0;b++) for(int s2=0;s2<NS;s2++) rho(t0s[b], (Idx)n, s2) = eta(t0s[b], (Idx)n, s2);
            rho.mult_sigma(c);                                       // sigma_c on the (n,t0_b) supports
            op_Dm.from_cpu<N>(hblk.data()+(size_t)n*N, rho.field);   // RHS block col n = D_m (sum_b sigma_c eta)
          }
          blk_tp_Dm.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);  // chi[c][n] = D_m^{-dag}(...)
          for(int n=0;n<n_sites;n++) for(Idx i=0;i<N;i++) psi_tp[n].field[i]=hblk[(size_t)n*N+i];
          chi = phi; chi.mult_sigma(c);                              // local sink sigphi_c = sigma_c phi'
          for(int n=0;n<n_sites;n++) for(int t=0;t<Nt;t++){
            Complex cdot(0,0), ddot(0,0);
            for(int s2=0;s2<NS;s2++){
              cdot += std::conj(psi_tp[n](t,(Idx)n,s2)) * chi(t,(Idx)n,s2);   // conn: chi[c][n]^dag sigphi_c
              ddot += std::conj(eta(t,(Idx)n,s2))       * chi(t,(Idx)n,s2);   // disc: eta^dag sigphi_c
            }
            Cs[c-1][t] += w_site[n]*cdot;
            Js[c-1][t] += w_site[n]*ddot;
          }
          if(parity){   // local V-- (m_P): BOTH legs tilde D^{-dag} (Eq. 3.65 local) -> (1+m_P)^{-2} at output.
            // psi_mm = tilde D^{-1}(sum_b sigma_c(n,t0_b)eta); sink sigma_c phimm (phimm = tilde D^{-dag} eta).
            for(int n=0;n<n_sites;n++){
              memset(rho.field, 0, Comp::N*CD);
              for(int b=0;b<n_t0;b++) for(int s2=0;s2<NS;s2++) rho(t0s[b], (Idx)n, s2) = eta(t0s[b], (Idx)n, s2);
              rho.mult_sigma(c);
              op_tilDmH.from_cpu<N>(hblk.data()+(size_t)n*N, rho.field);  // RHS = tilde^dag (sum_b sigma_c eta)
            }
            blk_tp_Dtil.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);  // psi_mm[n] = tilde D^{-1}(...)
            for(int n=0;n<n_sites;n++) for(Idx i=0;i<N;i++) psi_tp[n].field[i]=hblk[(size_t)n*N+i];
            chi = phimm; chi.mult_sigma(c);                            // local sink sigma_c phimm
            for(int n=0;n<n_sites;n++) for(int t=0;t<Nt;t++){
              Complex cdot(0,0);
              for(int s2=0;s2<NS;s2++) cdot += std::conj(psi_tp[n](t,(Idx)n,s2)) * chi(t,(Idx)n,s2);
              Cs_mm[c-1][t] += w_site[n]*cdot;
            }
          }
        }
        std::cout << " done ["<<elapsed()<<" s]" << std::endl;
      }

      // ============ LOCAL AXIAL s1,s2,s3 (rip off (1-D_ov)): tr[sigma(t0) [t0-leg] sigma(t) D_m^{-1}] ====
      // Local sibling of the connected axial ("replace K -> sigma, rip off (1-D_ov)"): SAME as the local
      // VECTOR above but the t0-leg propagator is the BACKWARD one (eta xi^dag).  Mass cases (Eqs. 3.50/3.51,
      // 3.57/3.58, 3.60/3.61): massless/m_F -> t0-leg = D_m^{-dag} (psi_A = D_m^{-1} sigma_c eta via op_DmH);
      // m_P -> t0-leg = (1+m_P)^{-1} tilde D_{m_P}^{-dag} (psi_A = tilde D^{-1} sigma_c eta via op_tilDmH; the
      // (1+m_P)^{-1} factor is applied at output -- the BARE local current has NO conserved-kernel correction
      // that cancels it, unlike the exact-K axial Eq. 3.67).  The local sink chi=sigma_c phi' (phi'=D_m^{-1}eta
      // = forward leg, STILL VALID here -- the connected axial overwrites phi BELOW) is SHARED.  Channel "Apm".
      {
        std::cout << "#   [local axial s1,s2,s3] ("<<n_sites<<" sites, superposed) ..." << std::flush;
        for(int c=1;c<=3;c++){
          for(int n=0;n<n_sites;n++){
            memset(rho.field, 0, Comp::N*CD);
            for(int b=0;b<n_t0;b++) for(int s2=0;s2<NS;s2++) rho(t0s[b], (Idx)n, s2) = eta(t0s[b], (Idx)n, s2);
            rho.mult_sigma(c);                                       // sigma_c on the (n,t0_b) supports
            if(parity) op_tilDmH.from_cpu<N>(hblk.data()+(size_t)n*N, rho.field);  // m_P: tilde D t0-leg (Eq. 3.61)
            else       op_DmH.from_cpu<N>(hblk.data()+(size_t)n*N, rho.field);     // D_m t0-leg
          }
          if(parity) blk_tp_Dtil.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER); // psi_A = tilde D^{-1}(sigma eta)
          else       blk_tp_Dm.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);   // psi_A = D_m^{-1}(sigma eta)
          for(int n=0;n<n_sites;n++) for(Idx i=0;i<N;i++) psi_tp[n].field[i]=hblk[(size_t)n*N+i];
          chi = phi; chi.mult_sigma(c);                              // local sink sigphi_c = sigma_c phi' (shared)
          for(int n=0;n<n_sites;n++) for(int t=0;t<Nt;t++){
            Complex cdot(0,0);
            for(int s2=0;s2<NS;s2++) cdot += std::conj(psi_tp[n](t,(Idx)n,s2)) * chi(t,(Idx)n,s2);
            CsA[c-1][t] += w_site[n]*cdot;
          }
        }
        std::cout << " done ["<<elapsed()<<" s]" << std::endl;
      }

      // ============ CONNECTED AXIAL -- C_{A+-} (tp + sp; own GW chi=(1-D_ov)phi' and K^dag sink) ======
      // Valence legs (Sec. 1.1): flavor m_F -> D_m (m=m_F) legs (GW (1-D_ov) stays MASSLESS); parity m_P
      //   -> D_m forward + tilde D sink; massless -> D_m (= D_ov).  (m_F leg changed from massless D_ov.)
      // Only C_{A+-} (Apm) is computed; C_{A-+} = reflection dt->Nt-dt (Eq. 3.57) is reconstructed downstream.
      // psi_tp/psi_sp are REUSED (the vector results above are already written to h5); phi is overwritten
      // with the axial forward leg.
      {
        // forward leg phi'_A = D_m^{-1} eta (m_F: m=m_F; m_P/massless: D_m);  chi = (1 - D_ov) phi'_A (massless GW)
        std::cout << "#   [axial] forward leg + chi=(1-D_ov)phi' ..." << std::flush;
        // if(flavor){ op_DH.from_cpu<N>(tmp.field, eta.field);  op_Dsq.solve<N>(phi.field, tmp.field, Comp::TOL_OUTER); }  // OLD: massless D_ov fwd for m_F
        // else      { op_DmH.from_cpu<N>(tmp.field, eta.field); op_Dmsq.solve<N>(phi.field, tmp.field, Comp::TOL_OUTER); }
        op_DmH.from_cpu<N>(tmp.field, eta.field); op_Dmsq.solve<N>(phi.field, tmp.field, Comp::TOL_OUTER);  // D_m forward (ALL cases)
        op_oneMinusD.from_cpu<N>(chi.field, phi.field);            // chi = (1 - D_ov) phi'_A
        std::cout << " done ["<<elapsed()<<" s]" << std::endl;

        // --- axial tp ---  source: psi_tp[n] = X_sink^{-1} (1 - D_ov) rho_tp[n]   (rho_tp = summed-origin
        // K^dag eta, REUSED from vec ++; superposed, ylm stripped).
        std::cout << "#   [axial tp] source ("<<n_sites<<" sites, superposed) + sink ..." << std::flush;
        for(int n=0;n<n_sites;n++){
          op_oneMinusD.from_cpu<N>(rho.field, rho_tp[n].field);          // rho = (1 - D_ov) (summed K^dag eta)
          // if(flavor)    op_DH.from_cpu<N>(hblk.data()+(size_t)n*N, rho.field);   // OLD: massless D_ov^dag sink for m_F
          if(parity)  op_tilDmH.from_cpu<N>(hblk.data()+(size_t)n*N, rho.field);
          else        op_DmH.from_cpu<N>(hblk.data()+(size_t)n*N, rho.field);       // m_F + massless: D_m^dag
        }
        // if(flavor)   blk_tp_D.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);   // OLD: massless block
        if(parity) blk_tp_Dtil.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);
        else       blk_tp_Dm.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);       // m_F + massless: D_m
        for(int n=0;n<n_sites;n++) for(Idx i=0;i<N;i++) psi_tp[n].field[i]=hblk[(size_t)n*N+i];
        {
          CUDA_CHECK(cudaMemcpy(d_sinkvec, reinterpret_cast<CuC*>(chi.field), N*CD, H2D));
          for(int n=0;n<n_sites;n++){
            kblock.apply_k_dag_block_t(d_kphi_block, d_sinkvec, U, (Idx)n);  // K^dag(n,t) chi for ALL t
            CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(kblk.data()), d_kphi_block, (size_t)N*Nt*CD, D2H));
            for(int t=0;t<Nt;t++){ kcol(t); Atp[t] += w_tp[n]*psi_tp[n].dag(kphi); }   // superposed; absolute t; summed over patterns
          }
        }
        std::cout << " done ["<<elapsed()<<" s]" << std::endl;

        // --- axial sp ---  source: psi_sp[a] = X_sink^{-1} (1 - D_ov) rho_sp[a]   (summed-origin, reused).
        std::cout << "#   [axial sp] source ("<<n_links<<" links, superposed) + sink ..." << std::flush;
        for(int a=0;a<n_links;a++){
          op_oneMinusD.from_cpu<N>(rho.field, rho_sp[a].field);          // rho = (1 - D_ov) (summed K^dag eta)
          // if(flavor)    op_DH.from_cpu<N>(hblk.data()+(size_t)a*N, rho.field);   // OLD: massless D_ov^dag sink for m_F
          if(parity)  op_tilDmH.from_cpu<N>(hblk.data()+(size_t)a*N, rho.field);
          else        op_DmH.from_cpu<N>(hblk.data()+(size_t)a*N, rho.field);       // m_F + massless: D_m^dag
        }
        // if(flavor)   blk_sp_D.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);   // OLD: massless block
        if(parity) blk_sp_Dtil.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);
        else       blk_sp_Dm.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);       // m_F + massless: D_m
        for(int a=0;a<n_links;a++) for(Idx i=0;i<N;i++) psi_sp[a].field[i]=hblk[(size_t)a*N+i];
        {
          CUDA_CHECK(cudaMemcpy(d_sinkvec, reinterpret_cast<CuC*>(chi.field), N*CD, H2D));
          for(int a=0;a<n_links;a++){
            const BaseLink lk = base.links[a];
            const Idx il = base.map2il.at(lk);
            kblock.apply_k_dag_block_t(d_kphi_block, d_sinkvec, U, lk);  // K^dag(lk,t) chi for ALL t
            CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(kblk.data()), d_kphi_block, (size_t)N*Nt*CD, D2H));
            for(int t=0;t<Nt;t++){ kcol(t); Asp[t] += w_sp[il]*psi_sp[a].dag(kphi); }   // superposed; absolute t; summed over patterns
          }
        }
        std::cout << " done ["<<elapsed()<<" s]" << std::endl;
      }
      } // ===== end 4-pattern dilution loop =====

      // ----- write the 4-pattern-summed correlators (one file per hit; superposed origins -> fold in analysis) -----
      write_corr(h5, hp+"tp/Vpp", Ctp);  if(!parity) write_corr_conj(h5, hp+"tp/Vmm", Ctp);
      write_corr(h5, hp+"sp/Vpp", Csp);  if(!parity) write_corr_conj(h5, hp+"sp/Vmm", Csp);
      if(parity){   // m_P: exact-K V-- (Eq. 3.65; tilde D, NO (1+m_P) factor -- kernel-cancelled) + tilde disc.
        write_corr(h5, hp+"tp/Vmm", Ctp_mm);    // SAME superposed frame as Vpp (audit I1/I2 fix)
        write_corr(h5, hp+"sp/Vmm", Csp_mm);
        write_vec (h5, hp+"disc/tp/Jtil", JtpT_til);
        write_vec (h5, hp+"disc/sp/Jtil", JspT_til);
        // A/B (audit I1/I2): these REPLACE the OLD per-origin writes that ran INSIDE the dilution loop:
        //   write_corr(h5, hp+"t0_<b>/tp/Vmm", Ctp[b]); write_corr(h5, hp+"t0_<b>/sp/Vmm", Csp[b]);
        //   write_vec (h5, hp+"disc/{tp,sp}/Jtil", JtpT/JspT);   // crashed on pattern 2 (duplicate key);
        // now superposed (h0/-level, no t0_<b> nesting) + hit-scope-accumulated + written once here.
      }
      write_corr(h5, hp+"axial/tp/Apm", Atp);
      write_corr(h5, hp+"axial/sp/Apm", Asp);
      write_vec (h5, hp+"disc/tp/J", Jtp);
      write_vec (h5, hp+"disc/sp/J", Jsp);
      for(int c=1;c<=3;c++){                                 // LOCAL s1,s2,s3 (conn Vpp/Vmm + disc J)
        write_corr(h5, hp+"s"+std::to_string(c)+"/Vpp", Cs[c-1]);
        if(!parity) write_corr_conj(h5, hp+"s"+std::to_string(c)+"/Vmm", Cs[c-1]);
        else {   // m_P local V-- : BOTH legs tilde -> (1+m_P)^{-2} (bare current, no kernel cancellation)
          const Complex inv1pmP2 = Complex(1.0,0.0)/((Complex(1.0,0.0)+valence_mass)*(Complex(1.0,0.0)+valence_mass));
          for(auto& z : Cs_mm[c-1]) z *= inv1pmP2;
          write_corr(h5, hp+"s"+std::to_string(c)+"/Vmm", Cs_mm[c-1]);   // NEW (audit I3): was silently dropped for m_P
        }
        write_vec (h5, hp+"disc/s"+std::to_string(c)+"/J", Js[c-1]);
        // m_P: bare local axial carries the (1+m_P)^{-1} factor of the backward propagator (Eq. 3.61); it is
        // NOT cancelled here (no conserved-kernel correction, unlike the exact-K axial Eq. 3.67).
        if(parity){ const Complex inv1pmP = Complex(1.0,0.0)/(Complex(1.0,0.0)+valence_mass);
                    for(auto& z : CsA[c-1]) z *= inv1pmP; }
        write_corr(h5, hp+"axial/s"+std::to_string(c)+"/Apm", CsA[c-1]);   // LOCAL AXIAL s_c (single complex channel)
      }

      // ----- CONDENSATES (Eqs. 1.23 PS, 1.55 FS): three bilinear traces (full spacetime-summed). -----
      // Mass-pattern-dependent propagators (Eqs. 3.50/3.51 massless, 3.57/3.58 m_F, 3.60/3.61 m_P):
      //   etadag_xi        = <eta^dag xi>             = -tr D_m^{-1}                 = -acc_etadag_phi  (ALL cases)
      //   xidag_eta        = <xi^dag eta>:
      //        massless/m_F : = conj(etadag_xi)                       (backward leg = (D_m^{-1})^dag)
      //        m_P          : = -(1+m_P)^{-1} acc_etadag_phimm        (backward leg = (1+m_P)^{-1} tilde D^{-dag})
      //   xidag_1mDdag_eta = <xi^dag (1-D_ov^dag) eta>:
      //        massless/m_F : = -conj(acc_etadag_1mD_phi)             (phi'-reuse via trace cyclicity)
      //        m_P          : = -(1+m_P)^{-1} acc_etadag_1mDdag_phimm
      // Stored = AREA-WEIGHTED spacetime-summed traces (site weight A_n=dual_areas folded in via eta_A).
      // Analysis: spacetime-AVERAGE by /(Nt * sum_n dual_areas[n])  (= Nt*4pi on the unit sphere); then
      //   sigma_PS = etadag_xi + xidag_eta;   sigma_FS = etadag_xi - xidag_1mDdag_eta.
      const Complex etadag_xi = -acc_etadag_phi;
      Complex xidag_eta, xidag_1mDdag_eta;
      if(parity){
        const Complex inv1pmP = Complex(1.0,0.0) / (Complex(1.0,0.0) + valence_mass);  // 1/(1+m_P)
        xidag_eta        = -inv1pmP * acc_etadag_phimm;
        xidag_1mDdag_eta = -inv1pmP * acc_etadag_1mDdag_phimm;
      } else {
        xidag_eta        =  std::conj(etadag_xi);
        xidag_1mDdag_eta = -std::conj(acc_etadag_1mD_phi);
      }
      write_vec(h5, hp+"condensate/etadag_xi",        std::vector<Complex>{ etadag_xi });
      write_vec(h5, hp+"condensate/xidag_eta",        std::vector<Complex>{ xidag_eta });
      write_vec(h5, hp+"condensate/xidag_1mDdag_eta", std::vector<Complex>{ xidag_1mDdag_eta });

      h5.createDataSet("complete", std::vector<int>{1});   // per-hit sentinel (written LAST)
      h5p.reset();                                          // close before the atomic rename
      std::filesystem::rename(h5tmp, h5path_h);            // atomic publish: corr.<k>.h<h>.h5
      const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now()-t_hit0).count();
      std::cout << "#   hit "<<(h+1)<<" done (exact tp+sp+axial + local s1,s2,s3 + local axial; diluted, origins SUPERPOSED) ["<<secs<<" s] -> "<<h5path_h << std::endl;
    } // hits
  } // k

  CUDA_CHECK(cudaFree(d_sinkvec)); CUDA_CHECK(cudaFree(d_kphi_block));  // C6f-c sink block-t buffers
  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();  // BlockedMat engines free their own scratch (dtor)
  return 0;
}
