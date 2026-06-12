// jj_corr_dilute_claude.cu
// -----------------------------------------------------------------------------
// DILUTED + MASTER-FIELD variant of jj_corr_block_t_claude.cu.  In ONE pass sharing the diluted sink leg
// phi' = D_m^{-1} eta, it computes BOTH the exact conserved current K AND the local (ultralocal e.sigma)
// current:
//   exact K : connected tp + sp, disc, axial.   (ylm REMOVED -> local-ylm is a SEPARATE one-end-trick file.)
//   local   : connected sp/tp diagonal channels s1,s2,s3, + local disc J_a(t).  weight w_site = dual_areas.
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
  printf("  --ninter <N>         ensemble config stride: measure ckpoint_lat.k for k=kmin,kmin+N,... (default: 10)\n");
  printf("  --kmin <k>           first ckpoint index, inclusive (default: 0)\n");
  printf("  --kmax <k>           one past last ckpoint index, EXCLUSIVE: kmin <= k < kmax (default: 1000000)\n");
  printf("  -h, --help           show this help\n");
  exit(0);
}

void ParseArgs(int argc, char* argv[],
               double& gsq, int& Nf, double& nu0, double& nu1,
               double& mass_re, double& mass_im,
               std::string& ens_dir, int& nhits, int& n_t0, int& ninter,
               int& kmin, int& kmax){
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
    {"ninter",  required_argument, nullptr, 'I'},
    {"kmin",    required_argument, nullptr, 'a'},
    {"kmax",    required_argument, nullptr, 'b'},
    {"help",    no_argument,       nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
  };
  int opt, idx;
  while((opt = getopt_long(argc, argv, "g:N:n:m:r:i:e:H:T:I:a:b:h", long_opts, &idx)) != -1){
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
    case 'I': ninter  = std::stoi(optarg); break;
    case 'a': kmin    = std::stoi(optarg); break;
    case 'b': kmax    = std::stoi(optarg); break;
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
  int ninter=10;   // ensemble config stride (k = kmin, kmin+ninter, ...)
  int kmin=0;          // first ckpoint index, inclusive
  int kmax=1000000;    // one past last, exclusive: kmin <= k < kmax (existence-break ends it earlier)

  ParseArgs(argc, argv, gsq, Nf, nu0, nu1, mass_re, mass_im, ens_dir, nhits, n_t0, ninter, kmin, kmax);
  if(nu1 < 0.0) nu1 = nu0;    // valence asymmetry defaults to the sea value nu0 (knob retained)

  const Complex valence_mass(mass_re, mass_im);
  const bool free_field = ens_dir.empty();

  // parity case: purely imaginary valence mass -> dagger leg uses \tilde D_{m_P}
  const bool parity = (std::abs(mass_im) > 1.0e-15) && (std::abs(mass_re) <= 1.0e-15);
  // flavor case: purely real valence mass -> axial uses massless D_ov legs (m_F formulas)
  const bool flavor = (std::abs(mass_re) > 1.0e-15) && (std::abs(mass_im) <= 1.0e-15);

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
  //   Dm   = D_ov + m          (Eq. 3.60; vector (++) both legs; tp/sp/ylm)
  //   Dtil = D_ov + m/(1-m)    (\tilde D_{m_P}, Eq. 3.63; parity (-) dagger leg)
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
  const std::string dir_out = "data_"+esnid+"/corr_dil_nt0"+std::to_string(n_t0)+"_nhits"+std::to_string(nhits)+"/";
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
  // with stride ninter (--ninter; default 10).
  const int k_ckpoint = free_field ? 1 : ninter;
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

      // HIT-SCOPE accumulators, SUMMED over the 4 dilution patterns (and over insertions/origins):
      //   conn (superposed, absolute t): Ctp/Csp (vector), Atp/Asp (axial);  disc (RAW): Jtp/Jsp.
      std::vector<Complex> Ctp(Nt,Complex(0,0)), Csp(Nt,Complex(0,0));
      std::vector<Complex> Atp(Nt,Complex(0,0)), Asp(Nt,Complex(0,0));
      std::vector<Complex> Jtp(Nt,Complex(0,0)), Jsp(Nt,Complex(0,0));
      // LOCAL current accumulators: conn Cs[c-1], disc Js[c-1] for channels c=1,2,3 (summed over patterns).
      std::vector<std::vector<Complex>> Cs(3, std::vector<Complex>(Nt,Complex(0,0)));
      std::vector<std::vector<Complex>> Js(3, std::vector<Complex>(Nt,Complex(0,0)));

      // ===== spin x even/odd-time DILUTION: 4 patterns (eta.time_spin_dilution, t_block=Nt/2 => interval 2).
      // Each pattern is its own diluted source; the estimator is the SUM over the 4 disjoint patterns. =====
      for(int t_s=0; t_s<2; t_s++) for(int spin=0; spin<NS; spin++){
        eta.time_spin_dilution(rng, t_s, Comp::Nt/2, spin);   // even/odd timeslices x spin component
        std::cout << "#   [dilution pattern t_s="<<t_s<<" spin="<<spin<<"]" << std::endl;

        // shared forward leg phi' = D_m^{-1} eta (reused by tp/sp/local/disc as the sink leg K(.,t)phi').
        std::cout << "#   phi' = D_m^{-1} eta : solving ..." << std::flush;
        op_DmH.from_cpu<N>(tmp.field, eta.field);
        op_Dmsq.solve<N>(phi.field, tmp.field, Comp::TOL_OUTER);
        std::cout << " done ["<<elapsed()<<" s]" << std::endl;

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
        std::cout << "#   [vec tp+ylm --] tilde source+sink ..." << std::flush;
        op_tilDm.from_cpu<N>(tmp.field, eta.field);                  // tmp = tilde eta
        op_tilDmsq.solve<N>(phimm.field, tmp.field, Comp::TOL_OUTER);// phimm = tilde D_{m_P}^{-dag} eta
        for(int b=0;b<n_t0;b++){
          const int t0=t0s[b];
          for(int l=0;l<N_ELL;l++) memset(srcL[l].field, 0, Comp::N*CD);
          for(int n=0;n<n_sites;n++){    // mrhs (C6e): build RHS block {tilde^dag K(n,t0)eta} then one block solve
            kop.set_temporal(U, t0, (Idx)n, /*dag=*/false);
            op_K.from_cpu<N>(rho.field, eta.field);                  // rho = K(n,t0) eta
            op_tilDmH.from_cpu<N>(hblk.data() + (size_t)n*N, rho.field);  // RHS block col n = tilde^dag rho
            // [pre-mrhs single solve, replaced by the block solve below]
            // op_tilDmH.from_cpu<N>(tmp.field, rho.field);             // tmp = tilde^dag rho
            // op_tilDmsq.solve<N>(psi_tp[ITP(n,b)].field, tmp.field, Comp::TOL_OUTER);  // tilde^{-1} K eta
            for(int l=0;l<N_ELL;l++){ const double w=W_ell[l][n]; for(Idx i=0;i<N;i++) srcL[l].field[i]+=w*rho.field[i]; }
          }
          blk_tp_Dtil.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);  // tilde^{-1} K eta (in-place host)
          for(int n=0;n<n_sites;n++) for(Idx i=0;i<N;i++) psi_tp[ITP(n,b)].field[i]=hblk[(size_t)n*N+i];
          for(int l=0;l<N_ELL;l++){
            op_tilDmH.from_cpu<N>(tmp.field, srcL[l].field);
            op_tilDmsq.solve<N>(psi_yl[IYL(l,b)].field, tmp.field, Comp::TOL_OUTER);
          }
        }
        std::vector<std::vector<Complex>> Ctp(n_t0, std::vector<Complex>(Nt, Complex(0,0)));
        std::vector<std::vector<Complex>> Gyl(N_ELL*n_t0, std::vector<Complex>(Nt, Complex(0,0)));
        for(int l=0;l<N_ELL;l++) for(int t=0;t<Nt;t++) memset(PhiLt[IPL(l,t)].field, 0, Comp::N*CD);
        CUDA_CHECK(cudaMemcpy(d_sinkvec, reinterpret_cast<CuC*>(phimm.field), N*CD, H2D));
        for(int n=0;n<n_sites;n++){
          kblock.apply_k_dag_block_t(d_kphi_block, d_sinkvec, U, (Idx)n);  // K^dag(n,t) phimm for ALL t
          CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(kblk.data()), d_kphi_block, (size_t)N*Nt*CD, D2H));
          for(int t=0;t<Nt;t++){
            kcol(t);                                                       // kphi = column t = K^dag(n,t) phimm
            for(int b=0;b<n_t0;b++){ const int dt=(t - t0s[b] + Nt)%Nt; Ctp[b][dt] += w_tp[n]*psi_tp[ITP(n,b)].dag(kphi); }
            for(int l=0;l<N_ELL;l++){ const double w=W_ell[l][n]; FermionVector& P=PhiLt[IPL(l,t)]; for(Idx i=0;i<N;i++) P.field[i]+=w*kphi.field[i]; }
          }
        }
        for(int t=0;t<Nt;t++){
          for(int b=0;b<n_t0;b++){ const int dt=(t - t0s[b] + Nt)%Nt; for(int l=0;l<N_ELL;l++) Gyl[IYL(l,b)][dt] += psi_yl[IYL(l,b)].dag(PhiLt[IPL(l,t)]); }
        }
        for(int b=0;b<n_t0;b++){
          const std::string kp = hp + "t0_" + std::to_string(b) + "/";
          write_corr(h5, kp+"tp/Vmm", Ctp[b]);
          for(int l=0;l<N_ELL;l++) write_corr(h5, kp+"ylm/Vmm/l"+std::to_string(l), Gyl[IYL(l,b)]);
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
        std::cout << "#   [vec sp --] tilde source+sink ..." << std::flush;
        for(int b=0;b<n_t0;b++){
          const int t0=t0s[b];
          for(int a=0;a<n_links;a++){    // mrhs (C6e): build RHS block {tilde^dag K(lk,t0)eta}
            const BaseLink lk = base.links[a];
            kop.set_spatial(U, t0, lk, /*dag=*/false);
            op_K.from_cpu<N>(rho.field, eta.field);                  // rho = K(lk,t0) eta
            op_tilDmH.from_cpu<N>(hblk.data() + (size_t)a*N, rho.field);  // RHS block col a = tilde^dag rho
            // [pre-mrhs single solve, replaced by the block solve below]
            // op_tilDmH.from_cpu<N>(tmp.field, rho.field);             // tmp = tilde^dag rho
            // op_tilDmsq.solve<N>(psi_sp[ISP(a,b)].field, tmp.field, Comp::TOL_OUTER);  // tilde^{-1} K eta
          }
          blk_sp_Dtil.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);  // tilde^{-1} K eta (in-place host)
          for(int a=0;a<n_links;a++) for(Idx i=0;i<N;i++) psi_sp[ISP(a,b)].field[i]=hblk[(size_t)a*N+i];
        }
        std::vector<std::vector<Complex>> Csp(n_t0, std::vector<Complex>(Nt, Complex(0,0)));
        CUDA_CHECK(cudaMemcpy(d_sinkvec, reinterpret_cast<CuC*>(phimm.field), N*CD, H2D));
        for(int a=0;a<n_links;a++){
          const BaseLink lk = base.links[a];
          const Idx il = base.map2il.at(lk);
          kblock.apply_k_dag_block_t(d_kphi_block, d_sinkvec, U, lk);  // K^dag(lk,t) phimm for ALL t
          CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(kblk.data()), d_kphi_block, (size_t)N*Nt*CD, D2H));
          for(int t=0;t<Nt;t++){
            kcol(t);                                                  // kphi = column t = K^dag(lk,t) phimm
            for(int b=0;b<n_t0;b++){ const int dt=(t - t0s[b] + Nt)%Nt; Csp[b][dt] += w_sp[il]*psi_sp[ISP(a,b)].dag(kphi); }
          }
        }
        for(int b=0;b<n_t0;b++){
          const std::string kp = hp + "t0_" + std::to_string(b) + "/";
          write_corr(h5, kp+"sp/Vmm", Csp[b]);
        }
        std::cout << " done ["<<elapsed()<<" s]" << std::endl;
      }

      // DISCONNECTED Jtp/Jsp rode the (++) sink applies above (RAW; summed over the dilution patterns).
      // WRITTEN ONCE per hit AFTER the dilution loop (below) -- writing here would collide across patterns.
      // parity: dagger-leg tilde trace \tilde T(a,t) = (K(a,t) tilphi)^dag eta, tilphi = tilde D_{m_P}^{-1} eta.
      // Cannot ride the connected parity sink (that applies K^dag phimm) -> own forward solve + K applies.
      if(parity){
        std::cout << "#   [disc --] tilde trace (tilphi solve + sink) ..." << std::flush;
        op_tilDmH.from_cpu<N>(tmp.field, eta.field);
        op_tilDmsq.solve<N>(tilphi.field, tmp.field, Comp::TOL_OUTER);   // tilphi = tilde D_{m_P}^{-1} eta
        std::vector<Complex> JtpT(Nt, Complex(0,0)), JspT(Nt, Complex(0,0));
        std::vector<std::vector<Complex>> JylT(N_ELL, std::vector<Complex>(Nt, Complex(0,0)));
        CUDA_CHECK(cudaMemcpy(d_sinkvec, reinterpret_cast<CuC*>(tilphi.field), N*CD, H2D));
        for(int n=0;n<n_sites;n++){
          kblock.apply_k_block_t(d_kphi_block, d_sinkvec, U, (Idx)n);  // K(n,t) tilphi for ALL t
          CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(kblk.data()), d_kphi_block, (size_t)N*Nt*CD, D2H));
          for(int t=0;t<Nt;t++){
            kcol(t);                                                  // kphi = column t = K(n,t) tilphi
            const Complex TtilNT = kphi.dag(eta);                     // \tilde T(n,t) = (K tilphi)^dag eta
            JtpT[t] += w_tp[n]*TtilNT;
            for(int l=0;l<N_ELL;l++) JylT[l][t] += W_ell[l][n]*TtilNT;
          }
        }
        for(int a=0;a<n_links;a++){
          const BaseLink lk = base.links[a];
          const Idx il = base.map2il.at(lk);
          kblock.apply_k_block_t(d_kphi_block, d_sinkvec, U, lk);      // K(lk,t) tilphi for ALL t
          CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(kblk.data()), d_kphi_block, (size_t)N*Nt*CD, D2H));
          for(int t=0;t<Nt;t++){
            kcol(t);                                                  // kphi = column t = K(lk,t) tilphi
            JspT[t] += w_sp[il]*kphi.dag(eta);
          }
        }
        write_vec(h5, hp+"disc/tp/Jtil", JtpT);
        write_vec(h5, hp+"disc/sp/Jtil", JspT);
        for(int l=0;l<N_ELL;l++) write_vec(h5, hp+"disc/ylm/l"+std::to_string(l)+"/Jtil", JylT[l]);
        std::cout << " done ["<<elapsed()<<" s]" << std::endl;
      }

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
        }
        std::cout << " done ["<<elapsed()<<" s]" << std::endl;
      }

      // ============ CONNECTED AXIAL -- C_{A+-} (tp + sp; own GW chi=(1-D_ov)phi' and K^dag sink) ======
      // Valence legs (Sec. 1.1): flavor m_F -> massless D_ov both legs; parity m_P -> sink tilde; else D_m.
      // Only C_{A+-} (Apm) is computed; C_{A-+} = reflection dt->Nt-dt (Eq. 3.57) is reconstructed downstream.
      // psi_tp/psi_sp are REUSED (the vector results above are already written to h5); phi is overwritten
      // with the axial forward leg.
      {
        // forward leg phi'_A = X^{-1} eta  (X = D_ov if flavor, else D_m);  chi = (1 - D_ov) phi'_A
        std::cout << "#   [axial] forward leg + chi=(1-D_ov)phi' ..." << std::flush;
        if(flavor){ op_DH.from_cpu<N>(tmp.field, eta.field);  op_Dsq.solve<N>(phi.field, tmp.field, Comp::TOL_OUTER); }
        else      { op_DmH.from_cpu<N>(tmp.field, eta.field); op_Dmsq.solve<N>(phi.field, tmp.field, Comp::TOL_OUTER); }
        op_oneMinusD.from_cpu<N>(chi.field, phi.field);            // chi = (1 - D_ov) phi'_A
        std::cout << " done ["<<elapsed()<<" s]" << std::endl;

        // --- axial tp ---  source: psi_tp[n] = X_sink^{-1} (1 - D_ov) rho_tp[n]   (rho_tp = summed-origin
        // K^dag eta, REUSED from vec ++; superposed, ylm stripped).
        std::cout << "#   [axial tp] source ("<<n_sites<<" sites, superposed) + sink ..." << std::flush;
        for(int n=0;n<n_sites;n++){
          op_oneMinusD.from_cpu<N>(rho.field, rho_tp[n].field);          // rho = (1 - D_ov) (summed K^dag eta)
          if(flavor)       op_DH.from_cpu<N>(hblk.data()+(size_t)n*N, rho.field);    // RHS block col n = X^dag rho
          else if(parity)  op_tilDmH.from_cpu<N>(hblk.data()+(size_t)n*N, rho.field);
          else             op_DmH.from_cpu<N>(hblk.data()+(size_t)n*N, rho.field);
        }
        if(flavor)      blk_tp_D.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);
        else if(parity) blk_tp_Dtil.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);
        else            blk_tp_Dm.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);
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
          if(flavor)       op_DH.from_cpu<N>(hblk.data()+(size_t)a*N, rho.field);    // RHS block col a = X^dag rho
          else if(parity)  op_tilDmH.from_cpu<N>(hblk.data()+(size_t)a*N, rho.field);
          else             op_DmH.from_cpu<N>(hblk.data()+(size_t)a*N, rho.field);
        }
        if(flavor)      blk_sp_D.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);
        else if(parity) blk_sp_Dtil.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);
        else            blk_sp_Dm.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);
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
      write_corr(h5, hp+"axial/tp/Apm", Atp);
      write_corr(h5, hp+"axial/sp/Apm", Asp);
      write_vec (h5, hp+"disc/tp/J", Jtp);
      write_vec (h5, hp+"disc/sp/J", Jsp);
      for(int c=1;c<=3;c++){                                 // LOCAL s1,s2,s3 (conn Vpp/Vmm + disc J)
        write_corr(h5, hp+"s"+std::to_string(c)+"/Vpp", Cs[c-1]);
        if(!parity) write_corr_conj(h5, hp+"s"+std::to_string(c)+"/Vmm", Cs[c-1]);
        write_vec (h5, hp+"disc/s"+std::to_string(c)+"/J", Js[c-1]);
      }

      h5.createDataSet("complete", std::vector<int>{1});   // per-hit sentinel (written LAST)
      h5p.reset();                                          // close before the atomic rename
      std::filesystem::rename(h5tmp, h5path_h);            // atomic publish: corr.<k>.h<h>.h5
      const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now()-t_hit0).count();
      std::cout << "#   hit "<<(h+1)<<" done (exact tp+sp+axial + local s1,s2,s3; 4-pattern diluted, origins SUPERPOSED) ["<<secs<<" s] -> "<<h5path_h << std::endl;
    } // hits
  } // k

  CUDA_CHECK(cudaFree(d_sinkvec)); CUDA_CHECK(cudaFree(d_kphi_block));  // C6f-c sink block-t buffers
  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();  // BlockedMat engines free their own scratch (dtor)
  return 0;
}
