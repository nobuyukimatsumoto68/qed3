// jj_corr_block_t_claude.cu  (C6f-c: t-BLOCKED sink variant of jj_corr_mrhs_claude.cu)
// Same physics/output as jj_corr_mrhs; the SINK passes are restructured outer-n + block-over-t: each
// pass `for t: for n/a: op_K.from_cpu(kphi,sinkvec)` becomes `for n/a: kblock.apply_k[_dag]_block_t(
// d_kphi_block, sinkvec, U, fixed); D2H; for t: <same accumulation on column t>`. Sink K applies drop
// from Nt*n per pass to n block applies (~3x apply, the C6f lever). Numerically equals jj_corr_mrhs to
// ~1e-15 (block K's Term B solve_shift_block) -> validate with h5diff -d 1e-10 (NOT bit-exact). The
// mrhs SOURCE-solve batching (blk_*.solve_sq_from_cpu) is KEPT unchanged.
//
// UNIFIED current-current correlator program (plan: conserved_current_correlators_impl_plan_v3_claude.md
// Sec. 3.8).  ONE program computes EVERYTHING per config -- connected AND disconnected -- sharing
// phi'=D_m^{-1} eta and the K phi' applies between the two pieces.  Output: data_<ESNID>/corr_nt0<N>_nhits<H>/.
//   - disconnected (vector): single-time traces J(t) = sum_a w_a eta^dag K(a,t) phi'  (raw; tp/sp/ylm + parity Jtil).
//   - connected tp/sp (vector + axial): source leg psi_a(t0)=D^{-dag} K^dag(a,t0) eta, sink = K(a,t) phi'.
//   - connected ylm (vector): m-summed tower G_l(t), ONE solve per l (combined weight W_l, Sec. 3.6).
//
// Unification lever (Sec. 3.8): the disc trace and the connected SINK both evaluate K(a,t) phi' with the
// SAME phi'.  The (++) tp+ylm sink loop feeds disc-tp + conn-tp + conn-ylm; the (++) sp sink loop feeds
// disc-sp + conn-sp -- so the non-parity disc costs NO extra K applies.  Only the connected SOURCE legs
// cost solves; the parity disc tilde trace (\tilde D_{m_P}^{-1} eta) needs its own solve + sink pass.
// The standalone jj_disc_claude.cu is kept as the cheap disc-only path (many configs, no connected solves).
//
// Conventions (v3): valence mass via --mass-re/--mass-im (D_m = D_ov + m); ensemble via --ens-dir
// (omit => free field, U = 1); kernel K is always the massless-form Noether kernel (mass-independent).
// Parity ( m purely imaginary ) dagger leg uses \tilde D_{m_P} = D_ov + m_P/(1-m_P); massless / m_F uses
// D_m itself.  All overlap applies/solves go through the multi-shift (_ms) entry points.
//
// Solver: overlap multishift (*_deviceAsyncLaunch_ms) + the kernel is the multishift apply_k_ms via
// ConservedCurrent::operator() (inherited through op_K.from_cpu).  --n-t0 sets the source origins.

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

#include "sparse_dirac.h"
// #include "matpoly.h"
#include "matpoly_claude.h"

#include "overlap_wmass_claude.h"        // complex-mass overlap (massless at mass=0)
#include "blocked_mat_claude.h"      // C7: BlockedMat<N,NSTACK,Op> (mrhs block solves, host-block solve_sq_from_cpu)
#include "conserved_current_claude.h"   // ConservedCurrent: apply_k / apply_k_dag (multishift apply_k_ms)
#include "conserved_current_block_claude.h"  // C6f-c: ConservedCurrentBlockT::apply_k[_dag]_block_t (t-blocked sink)

//------------------------------------------
#include <getopt.h>

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

  // ylm m-summed tower (Sec. 3.6): keep l_1=l_2=l, sum m_1,m_2 over -l..l.  ONE solve per l via the
  // COMBINED real weight W_ell[l][n] = (A_n / kappa_t(n)) * sum_{m=-l}^{l} Y_lm(n^)  (kappa^1, Eq. 4.36).
  constexpr int L_MAX = 2;
  constexpr int N_ELL = L_MAX + 1;                 // l = 0,1,2  (l=0 ~ 0 by charge conservation, a check)
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
  const std::string dir_out = "data_"+esnid+"/corr_nt0"+std::to_string(n_t0)+"_nhits"+std::to_string(nhits)+"/";
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
  std::vector<FermionVector> psi_tp(n_sites * n_t0);   // tp source legs  psi_tp[ITP(n,b)]
  std::vector<FermionVector> psi_yl(N_ELL   * n_t0);   // ylm source towers psi_yl[IYL(l,b)]
  std::vector<FermionVector> psi_sp(n_links * n_t0);   // sp source legs  psi_sp[ISP(a,b)]
  // Reuse A: cache the K^dag-applied sources rho = K^dag(insertion,t0) eta produced by the vector (++)
  // source passes (always run, dag=true).  The axial source passes then reuse them as (1-D_ov)*rho_*
  // instead of re-applying K^dag, saving n_t0*(n_sites+n_links) kernel applies per hit.
  std::vector<FermionVector> rho_tp(n_sites * n_t0);   // = K^dag(n ,t0_b) eta  (sites)
  std::vector<FermionVector> rho_sp(n_links * n_t0);   // = K^dag(lk,t0_b) eta  (links)
  auto ITP = [&](int n,int b){ return n*n_t0 + b; };
  auto IYL = [&](int l,int b){ return l*n_t0 + b; };
  auto ISP = [&](int a,int b){ return a*n_t0 + b; };

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
    // ONE file per config.  Resume: skip ONLY if the "complete" sentinel (written LAST, after every
    // dataset) is present -- an interrupted write lacks it -> recompute.  Read-only open.
    const std::string h5path = dir_out + "corr." + std::to_string(k) + ".h5";
    if(std::filesystem::exists(h5path)){
      bool complete = false;
      try { HighFive::File h5c(h5path, HighFive::File::ReadOnly); complete = h5c.exist("complete"); }
      catch(...) {}
      if(complete){ std::cout<<"# skip k="<<k<<" (complete)"<<std::endl; continue; }
      std::cout<<"# k="<<k<<" exists but INCOMPLETE -> recompute"<<std::endl;
    }
    D.update(U);  Dm.update(U);  Dtil.update(U);
    std::cout << "# k="<<k<<(free_field?" (free field)":"")
              << "  lambda_min/max="<<Dm.lambda_min<<"/"<<Dm.lambda_max<<std::endl;

    // ATOMIC WRITE: write to "<h5path>.tmp", then rename() to the final name after the 'complete'
    // sentinel + close.  rename(2) is atomic on POSIX, so readers / the resume check above never see a
    // partial file; an interrupted run leaves only a stale ".tmp" (ignored by the *.h5 globs).
    const std::string h5tmp = h5path + ".tmp";
    auto h5p = std::make_unique<HighFive::File>(h5tmp,
                 HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
    HighFive::File& h5 = *h5p;
    h5.createDataSet("t0s",   t0s);
    h5.createDataSet("n_t0",  std::vector<int>{n_t0});
    h5.createDataSet("nhits", std::vector<int>{nhits});
    h5.createDataSet("ls",    ls);

    for(int h=0; h<nhits; h++){
      const auto t_hit0 = std::chrono::steady_clock::now();
      auto elapsed = [&](){ return std::chrono::duration<double>(std::chrono::steady_clock::now()-t_hit0).count(); };
      std::cout << "# k="<<k<<" hit "<<(h+1)<<"/"<<nhits<<"  (n_t0="<<n_t0<<", n_sites="<<n_sites
                << ", n_links="<<n_links<<", n_ell="<<N_ELL<<")" << std::endl;
      eta.fill_z2_source(rng);
      const std::string hp = "h" + std::to_string(h) + "/";   // key prefix /h{h}/

      // shared forward leg phi' = D_m^{-1} eta (op_DmH RHS-former + op_Dmsq CG); reused by ALL connected
      // projections (tp/sp/ylm) as the sink leg K(.,t)phi'.
      std::cout << "#   phi' = D_m^{-1} eta : solving ..." << std::flush;
      op_DmH.from_cpu<N>(tmp.field, eta.field);
      op_Dmsq.solve<N>(phi.field, tmp.field, Comp::TOL_OUTER);
      std::cout << " done ["<<elapsed()<<" s]" << std::endl;

      // disconnected single-current traces J(t) = sum_a w_a eta^dag K(a,t) phi' (RAW); these RIDE the
      // connected (++) sink passes below (the K(.,t)phi' applies are shared) -- zero extra K applies in
      // the non-parity case.  Origin-independent: accumulated over the whole hit, written once.
      std::vector<Complex> Jtp(Nt, Complex(0,0)), Jsp(Nt, Complex(0,0));
      std::vector<std::vector<Complex>> Jyl(N_ELL, std::vector<Complex>(Nt, Complex(0,0)));

      // ============ CONNECTED VECTOR -- temporal tp + ylm (shared phi' + shared K(n,t)phi' pass) =======
      // (++) source legs at all origins (the per-site K^dag(n,t0)eta apply is shared by tp and ylm):
      //   tp:  psi_tp[ITP(n,b)] = D_m^{-dag} K^dag(n,t0_b) eta
      //   ylm: psi_yl[IYL(l,b)] = D_m^{-dag} sum_n W_ell[l][n] K^dag(n,t0_b) eta   (m-summed tower)
      std::cout << "#   [vec tp+ylm ++] source solves ("<<n_t0<<" t0 x ("<<n_sites<<" tp + "<<N_ELL<<" ylm)) ..." << std::endl;
      for(int b=0;b<n_t0;b++){
        const int t0=t0s[b];
        for(int l=0;l<N_ELL;l++) memset(srcL[l].field, 0, Comp::N*CD);
        // mrhs (C6e): build the n_sites RHS block {D_m K^dag(n,t0) eta} then ONE block solve
        // (op_Dmsq.solve_block_cg over DDH_..._ms_block) instead of n_sites single op_Dmsq solves.
        for(int n=0;n<n_sites;n++){
          kop.set_temporal(U, t0, (Idx)n, /*dag=*/true);
          op_K.from_cpu<N>(rho_tp[ITP(n,b)].field, eta.field);       // rho = K^dag(n,t0) eta (CACHED for axial reuse)
          op_Dm.from_cpu<N>(hblk.data() + (size_t)n*N, rho_tp[ITP(n,b)].field);   // RHS block col n = D_m rho
          // [pre-mrhs single solve, replaced by the block solve below]
          // op_Dm.from_cpu<N>(tmp.field, rho_tp[ITP(n,b)].field);                // tmp = D_m rho
          // op_Dmsq.solve<N>(psi_tp[ITP(n,b)].field, tmp.field, Comp::TOL_OUTER);  // = D_m^{-dag} K^dag eta
          for(int l=0;l<N_ELL;l++){ const double w=W_ell[l][n]; for(Idx i=0;i<N;i++) srcL[l].field[i]+=w*rho_tp[ITP(n,b)].field[i]; }
        }
        blk_tp_Dm.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);  // psi block = D_m^{-dag} K^dag eta (in-place host)
        for(int n=0;n<n_sites;n++) for(Idx i=0;i<N;i++) psi_tp[ITP(n,b)].field[i]=hblk[(size_t)n*N+i];
        for(int l=0;l<N_ELL;l++){
          op_Dm.from_cpu<N>(tmp.field, srcL[l].field);               // tmp = D_m srcL[l]
          op_Dmsq.solve<N>(psi_yl[IYL(l,b)].field, tmp.field, Comp::TOL_OUTER);
        }
        std::cout << "#     t0="<<t0<<" ("<<(b+1)<<"/"<<n_t0<<") source done ["<<elapsed()<<" s]" << std::endl;
      }
      // (++) shared sink pass: kphi = K(n,t) phi' ONCE per (n,t); feeds tp (per (n,b)) + ylm (accumulate)
      std::cout << "#   [vec tp+ylm ++] sink pass (block-t: "<<n_sites<<" applies x "<<Nt<<" t) ..." << std::flush;
      {
        std::vector<std::vector<Complex>> Ctp(n_t0, std::vector<Complex>(Nt, Complex(0,0)));
        std::vector<std::vector<Complex>> Gyl(N_ELL*n_t0, std::vector<Complex>(Nt, Complex(0,0)));
        for(int l=0;l<N_ELL;l++) for(int t=0;t<Nt;t++) memset(PhiLt[IPL(l,t)].field, 0, Comp::N*CD);
        CUDA_CHECK(cudaMemcpy(d_sinkvec, reinterpret_cast<CuC*>(phi.field), N*CD, H2D));
        for(int n=0;n<n_sites;n++){
          kblock.apply_k_block_t(d_kphi_block, d_sinkvec, U, (Idx)n);     // K(n,t) phi' for ALL t
          CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(kblk.data()), d_kphi_block, (size_t)N*Nt*CD, D2H));
          for(int t=0;t<Nt;t++){
            kcol(t);                                                     // kphi = column t = K(n,t) phi'
            Jtp[t] += w_tp[n]*eta.dag(kphi);                             // disc tp: T(n,t)=eta^dag kphi (rides this apply)
            for(int b=0;b<n_t0;b++){ const int dt=(t - t0s[b] + Nt)%Nt; Ctp[b][dt] += w_tp[n]*psi_tp[ITP(n,b)].dag(kphi); }
            for(int l=0;l<N_ELL;l++){ const double w=W_ell[l][n]; FermionVector& P=PhiLt[IPL(l,t)]; for(Idx i=0;i<N;i++) P.field[i]+=w*kphi.field[i]; }
          }
        }
        for(int t=0;t<Nt;t++){
          for(int l=0;l<N_ELL;l++) Jyl[l][t] += eta.dag(PhiLt[IPL(l,t)]);   // disc ylm: J_l(t)=eta^dag PhiL[l](t)
          for(int b=0;b<n_t0;b++){ const int dt=(t - t0s[b] + Nt)%Nt; for(int l=0;l<N_ELL;l++) Gyl[IYL(l,b)][dt] += psi_yl[IYL(l,b)].dag(PhiLt[IPL(l,t)]); }
        }
        for(int b=0;b<n_t0;b++){
          const std::string kp = hp + "t0_" + std::to_string(b) + "/";
          write_corr(h5, kp+"tp/Vpp", Ctp[b]);
          if(!parity) write_corr_conj(h5, kp+"tp/Vmm", Ctp[b]);      // massless/m_F: Vmm = conj(Vpp)
          for(int l=0;l<N_ELL;l++){
            write_corr(h5, kp+"ylm/Vpp/l"+std::to_string(l), Gyl[IYL(l,b)]);
            if(!parity) write_corr_conj(h5, kp+"ylm/Vmm/l"+std::to_string(l), Gyl[IYL(l,b)]);
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
      // (++) source legs (insertion-diagonal over links): psi_sp[ISP(a,b)] = D_m^{-dag} K^dag(lk,t0_b) eta
      std::cout << "#   [vec sp ++] source solves ("<<n_t0<<" t0 x "<<n_links<<" links) + sink pass ..." << std::flush;
      for(int b=0;b<n_t0;b++){
        const int t0=t0s[b];
        for(int a=0;a<n_links;a++){    // mrhs (C6e): build RHS block {D_m K^dag(lk,t0)eta} then one block solve
          const BaseLink lk = base.links[a];
          kop.set_spatial(U, t0, lk, /*dag=*/true);
          op_K.from_cpu<N>(rho_sp[ISP(a,b)].field, eta.field);       // rho = K^dag(lk,t0) eta (CACHED for axial reuse)
          op_Dm.from_cpu<N>(hblk.data() + (size_t)a*N, rho_sp[ISP(a,b)].field);  // RHS block col a = D_m rho
          // [pre-mrhs single solve, replaced by the block solve below]
          // op_Dm.from_cpu<N>(tmp.field, rho_sp[ISP(a,b)].field);      // tmp = D_m rho
          // op_Dmsq.solve<N>(psi_sp[ISP(a,b)].field, tmp.field, Comp::TOL_OUTER);  // = D_m^{-dag} K^dag eta
        }
        blk_sp_Dm.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);  // psi block = D_m^{-dag} K^dag eta (in-place host)
        for(int a=0;a<n_links;a++) for(Idx i=0;i<N;i++) psi_sp[ISP(a,b)].field[i]=hblk[(size_t)a*N+i];
      }
      // (++) sink pass: kphi = K(lk,t) phi' once per (a,t); pair into Csp[b][dt]
      {
        std::vector<std::vector<Complex>> Csp(n_t0, std::vector<Complex>(Nt, Complex(0,0)));
        CUDA_CHECK(cudaMemcpy(d_sinkvec, reinterpret_cast<CuC*>(phi.field), N*CD, H2D));
        for(int a=0;a<n_links;a++){
          const BaseLink lk = base.links[a];
          const Idx il = base.map2il.at(lk);
          kblock.apply_k_block_t(d_kphi_block, d_sinkvec, U, lk);    // K(lk,t) phi' for ALL t
          CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(kblk.data()), d_kphi_block, (size_t)N*Nt*CD, D2H));
          for(int t=0;t<Nt;t++){
            kcol(t);                                                // kphi = column t = K(lk,t) phi'
            Jsp[t] += w_sp[il]*eta.dag(kphi);                       // disc sp: rides this apply
            for(int b=0;b<n_t0;b++){ const int dt=(t - t0s[b] + Nt)%Nt; Csp[b][dt] += w_sp[il]*psi_sp[ISP(a,b)].dag(kphi); }
          }
        }
        for(int b=0;b<n_t0;b++){
          const std::string kp = hp + "t0_" + std::to_string(b) + "/";
          write_corr(h5, kp+"sp/Vpp", Csp[b]);
          if(!parity) write_corr_conj(h5, kp+"sp/Vmm", Csp[b]);      // massless/m_F: Vmm = conj(Vpp)
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

      // ============ DISCONNECTED single-current traces (folded in; rode the (++) sink applies) ===========
      // Jtp/Jsp/Jyl were accumulated inside the (++) tp+ylm and sp sink passes above at zero extra K-applies;
      // origin-independent, written once per hit.  RAW (no 1/4pi), matching jj_disc_claude.cu.
      write_vec(h5, hp+"disc/tp/J", Jtp);
      write_vec(h5, hp+"disc/sp/J", Jsp);
      for(int l=0;l<N_ELL;l++) write_vec(h5, hp+"disc/ylm/l"+std::to_string(l)+"/J", Jyl[l]);
      std::cout << "#   [disc] J(t) written (rode the (++) sink passes)" << std::endl;
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

        // --- axial tp ---  source: psi_tp[ITP(n,b)] = X_sink^{-1} (1 - D_ov) K^dag(n,t0) eta
        std::cout << "#   [axial tp] source solves ("<<n_t0<<" t0 x "<<n_sites<<") + sink pass ..." << std::flush;
        for(int b=0;b<n_t0;b++){
          for(int n=0;n<n_sites;n++){    // mrhs (C6e): build RHS block {X^dag (1-D_ov) K^dag eta}
            op_oneMinusD.from_cpu<N>(rho.field, rho_tp[ITP(n,b)].field);  // rho = (1 - D_ov) K^dag(n,t0) eta  (K^dag reused from vec ++)
            if(flavor)       op_DH.from_cpu<N>(hblk.data()+(size_t)n*N, rho.field);      // RHS block col n = X^dag rho
            else if(parity)  op_tilDmH.from_cpu<N>(hblk.data()+(size_t)n*N, rho.field);
            else             op_DmH.from_cpu<N>(hblk.data()+(size_t)n*N, rho.field);
            // [pre-mrhs single solves, replaced by the block solve below]
            // if(flavor){      op_DH.from_cpu<N>(tmp.field, rho.field);     op_Dsq.solve<N>(psi_tp[ITP(n,b)].field, tmp.field, Comp::TOL_OUTER); }
            // else if(parity){ op_tilDmH.from_cpu<N>(tmp.field, rho.field); op_tilDmsq.solve<N>(psi_tp[ITP(n,b)].field, tmp.field, Comp::TOL_OUTER); }
            // else {           op_DmH.from_cpu<N>(tmp.field, rho.field);    op_Dmsq.solve<N>(psi_tp[ITP(n,b)].field, tmp.field, Comp::TOL_OUTER); }
          }
          if(flavor)      blk_tp_D.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);
          else if(parity) blk_tp_Dtil.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);
          else            blk_tp_Dm.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);
          for(int n=0;n<n_sites;n++) for(Idx i=0;i<N;i++) psi_tp[ITP(n,b)].field[i]=hblk[(size_t)n*N+i];
          // AXIAL YLM tower source: psi_yl[IYL(l,b)] = X_sink^{-1} (1-D_ov) srcL[l],
          //   srcL[l] = sum_n W_ell[l][n] rho_tp[ITP(n,b)]  (rho_tp = K^dag(n,t0) eta, cached from vec ++).
          //   (1-D_ov) pulls out of the m-sum; single solves (N_ELL per b, small).  Reuses srcL/psi_yl.
          for(int l=0;l<N_ELL;l++){
            memset(srcL[l].field, 0, Comp::N*CD);
            for(int n=0;n<n_sites;n++){ const double w=W_ell[l][n]; for(Idx i=0;i<N;i++) srcL[l].field[i]+=w*rho_tp[ITP(n,b)].field[i]; }
            op_oneMinusD.from_cpu<N>(rho.field, srcL[l].field);                 // rho = (1-D_ov) srcL[l]
            if(flavor)      { op_DH.from_cpu<N>(tmp.field, rho.field);    op_Dsq.solve<N>(psi_yl[IYL(l,b)].field, tmp.field, Comp::TOL_OUTER); }
            else if(parity) { op_tilDmH.from_cpu<N>(tmp.field, rho.field); op_tilDmsq.solve<N>(psi_yl[IYL(l,b)].field, tmp.field, Comp::TOL_OUTER); }
            else            { op_DmH.from_cpu<N>(tmp.field, rho.field);   op_Dmsq.solve<N>(psi_yl[IYL(l,b)].field, tmp.field, Comp::TOL_OUTER); }
          }
        }
        {
          std::vector<std::vector<Complex>> Atp(n_t0, std::vector<Complex>(Nt, Complex(0,0)));
          std::vector<std::vector<Complex>> Ayl(N_ELL*n_t0, std::vector<Complex>(Nt, Complex(0,0)));  // axial ylm tower
          for(int l=0;l<N_ELL;l++) for(int t=0;t<Nt;t++) memset(PhiLt[IPL(l,t)].field, 0, Comp::N*CD);  // zero ylm sink tower
          CUDA_CHECK(cudaMemcpy(d_sinkvec, reinterpret_cast<CuC*>(chi.field), N*CD, H2D));
          for(int n=0;n<n_sites;n++){
            kblock.apply_k_dag_block_t(d_kphi_block, d_sinkvec, U, (Idx)n);  // K^dag(n,t) chi for ALL t
            CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(kblk.data()), d_kphi_block, (size_t)N*Nt*CD, D2H));
            for(int t=0;t<Nt;t++){
              kcol(t);                                                      // kphi = column t = K^dag(n,t) chi
              for(int b=0;b<n_t0;b++){ const int dt=(t - t0s[b] + Nt)%Nt; Atp[b][dt] += w_tp[n]*psi_tp[ITP(n,b)].dag(kphi); }
              // AXIAL YLM sink tower: rides this loop at ZERO extra K applies (PhiLt[l](t) = sum_n W_ell K^dag(n,t) chi)
              for(int l=0;l<N_ELL;l++){ const double w=W_ell[l][n]; FermionVector& P=PhiLt[IPL(l,t)]; for(Idx i=0;i<N;i++) P.field[i]+=w*kphi.field[i]; }
            }
          }
          for(int t=0;t<Nt;t++) for(int b=0;b<n_t0;b++){ const int dt=(t - t0s[b] + Nt)%Nt;
            for(int l=0;l<N_ELL;l++) Ayl[IYL(l,b)][dt] += psi_yl[IYL(l,b)].dag(PhiLt[IPL(l,t)]); }
          for(int b=0;b<n_t0;b++){ const std::string kp=hp+"t0_"+std::to_string(b)+"/";
            write_corr(h5, kp+"axial/tp/Apm", Atp[b]);
            for(int l=0;l<N_ELL;l++) write_corr(h5, kp+"axial/ylm/Apm/l"+std::to_string(l), Ayl[IYL(l,b)]); }  // NEW: axial ylm tower
        }
        std::cout << " done ["<<elapsed()<<" s]" << std::endl;

        // --- axial sp ---  source: psi_sp[ISP(a,b)] = X_sink^{-1} (1 - D_ov) K^dag(lk,t0) eta
        std::cout << "#   [axial sp] source solves ("<<n_t0<<" t0 x "<<n_links<<") + sink pass ..." << std::flush;
        for(int b=0;b<n_t0;b++){
          for(int a=0;a<n_links;a++){    // mrhs (C6e): build RHS block {X^dag (1-D_ov) K^dag eta}
            op_oneMinusD.from_cpu<N>(rho.field, rho_sp[ISP(a,b)].field);  // rho = (1 - D_ov) K^dag(lk,t0) eta  (K^dag reused from vec ++)
            if(flavor)       op_DH.from_cpu<N>(hblk.data()+(size_t)a*N, rho.field);      // RHS block col a = X^dag rho
            else if(parity)  op_tilDmH.from_cpu<N>(hblk.data()+(size_t)a*N, rho.field);
            else             op_DmH.from_cpu<N>(hblk.data()+(size_t)a*N, rho.field);
            // [pre-mrhs single solves, replaced by the block solve below]
            // if(flavor){      op_DH.from_cpu<N>(tmp.field, rho.field);     op_Dsq.solve<N>(psi_sp[ISP(a,b)].field, tmp.field, Comp::TOL_OUTER); }
            // else if(parity){ op_tilDmH.from_cpu<N>(tmp.field, rho.field); op_tilDmsq.solve<N>(psi_sp[ISP(a,b)].field, tmp.field, Comp::TOL_OUTER); }
            // else {           op_DmH.from_cpu<N>(tmp.field, rho.field);    op_Dmsq.solve<N>(psi_sp[ISP(a,b)].field, tmp.field, Comp::TOL_OUTER); }
          }
          if(flavor)      blk_sp_D.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);
          else if(parity) blk_sp_Dtil.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);
          else            blk_sp_Dm.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);
          for(int a=0;a<n_links;a++) for(Idx i=0;i<N;i++) psi_sp[ISP(a,b)].field[i]=hblk[(size_t)a*N+i];
        }
        {
          std::vector<std::vector<Complex>> Asp(n_t0, std::vector<Complex>(Nt, Complex(0,0)));
          CUDA_CHECK(cudaMemcpy(d_sinkvec, reinterpret_cast<CuC*>(chi.field), N*CD, H2D));
          for(int a=0;a<n_links;a++){
            const BaseLink lk = base.links[a];
            const Idx il = base.map2il.at(lk);
            kblock.apply_k_dag_block_t(d_kphi_block, d_sinkvec, U, lk);  // K^dag(lk,t) chi for ALL t
            CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(kblk.data()), d_kphi_block, (size_t)N*Nt*CD, D2H));
            for(int t=0;t<Nt;t++){
              kcol(t);                                                  // kphi = column t = K^dag(lk,t) chi
              for(int b=0;b<n_t0;b++){ const int dt=(t - t0s[b] + Nt)%Nt; Asp[b][dt] += w_sp[il]*psi_sp[ISP(a,b)].dag(kphi); }
            }
          }
          for(int b=0;b<n_t0;b++){ const std::string kp=hp+"t0_"+std::to_string(b)+"/"; write_corr(h5, kp+"axial/sp/Apm", Asp[b]); }
        }
        std::cout << " done ["<<elapsed()<<" s]" << std::endl;
      }

      const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now()-t_hit0).count();
      std::cout << "#   hit "<<(h+1)<<" done (tp+sp+ylm vector + axial) ["<<secs<<" s]" << std::endl;
    } // hits

    h5.createDataSet("complete", std::vector<int>{1});   // sentinel: ALL datasets present (written LAST)
    h5p.reset();                                          // CLOSE the .tmp file before publishing
    std::filesystem::rename(h5tmp, h5path);               // atomic publish: now visible as corr.<k>.h5
    std::cout << "# wrote " << h5path << std::endl;
  } // k

  CUDA_CHECK(cudaFree(d_sinkvec)); CUDA_CHECK(cudaFree(d_kphi_block));  // C6f-c sink block-t buffers
  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();  // BlockedMat engines free their own scratch (dtor)
  return 0;
}
