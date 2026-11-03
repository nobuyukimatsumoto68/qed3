// jj_local_ylm_scalar_disc_stoch_claude.cu
// STOCHASTIC DISCONNECTED scalar one-point loop (identity vertex, l=0) for the flavor-singlet
// sigma_FS two-point function.  MASSLESS ONLY.  FORWARD SOLVE ONLY.
//
// COPY of jj_local_ylm_disc_stoch_claude.cu (production vector disc driver: at_from_ensdir, per-hit
// atomic write, "complete" gating, early config skip).  Changes vs the vector driver:
//   (1) loop = scalar identity vertex, l=0 only:
//         J_1(t) = tr[ W_00(t) D_ov^{-1} ],  W_00(t) = sum_n A_n Y_00(n^) 1_spin (n,t).
//       The vector currents s1/s2/s3 and the sigma_PS condensate are NOT computed here.
//   (2) NO (1 - D_ov) furnishing and NO D^{-dag} (backward) solve.  At m=0 Ginsparg-Wilson gives
//         D_ov^{-1} + D_ov^{-dag} = 1   =>   (1 - D_ov^dag) D_ov^{-dag} = - D_ov^{-1},
//       so the furnished backward leg is minus the forward leg and
//         J_FS = 2 J_1 = tr W + 2 i Im J_1 ,   J_PS = 2 Re J_1 = tr W  (config independent).
//       The analysis uses 2 i Im J_1.  (Derivation: fs_furnishing_derivation_claude.md, qed3-5c;
//       plan: final/analysis_scalar/scalar_fs_disc_impl_plan_claude.md.)
//   (3) the NSTACK hits of one config are solved TOGETHER as one mrhs block (BlockedMat::solve_sq),
//       NSTACK = -DNB_CLI (default 2).  -DNB_CLI=1 gives the unblocked reference (one column).
//       Each hit keeps its own RNG stream, seed string (esnid,k,h), and its own output file.
//   (4) output dir = data_<ESNID>/corr_ylm_disc_scalar_tb<tb>/  (the vector disc files are untouched).
//
// Estimator (per dilution pattern, per hit): with the volume source eta and phi = D_ov^{-1} eta,
//   J_1(t) += sum_n A_n Y_00(n^) conj(eta(t,n,spin)) phi(t,n,spin)     (RAW; summed over patterns).
// The spin loop reads one spin row per pattern; summing both patterns gives the full spin trace.
// NOTE: for a stochastic hit Re J_1 equals (1/2) tr W only up to noise; (1/2) tr W is stored as
// h0/disc/ylm/s0/half_trW for the statistical check.
//
// References: Ginsparg-Wilson, Phys. Rev. D 25, 2649 (1982); M. Luscher, hep-lat/9802011 (exact
// lattice chiral symmetry); J. Foley et al., hep-lat/0505023 (dilution of stochastic sources);
// B. Jegerlehner, hep-lat/9612014 (multishift solver inside BlockedMat).
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

  // constexpr int N_REFINE=1;   // L1 default; now compile-time via -DN_REFINE_CLI
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

// mrhs block width = number of hits solved together (compile-time; -DNB_CLI=1 => unblocked reference)
#ifndef NB_CLI
#define NB_CLI 2
#endif

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

#include "overlap_wmass_claude.h"        // complex-mass overlap (massless at mass=0)
#include "blocked_mat_claude.h"          // BlockedMat<N,NSTACK,Op> mrhs block solves (solve_sq_from_cpu)

//------------------------------------------
#include <getopt.h>

static int seed_from_string(const std::string& s){
  std::seed_seq seq(s.begin(), s.end());
  std::uint32_t w;
  seq.generate(&w, &w + 1);
  return static_cast<int>(w);
}

// extract a_t from an ensemble dir name.  The a_t token is the unique 'at' immediately
// followed by a digit, e.g. "...gsq2.000000at0.100000nu0..." -> 0.1 .  Returns -1.0 if not found.
// Used so the valence temporal spacing ALWAYS tracks the sea ensemble (no silent 0.2 mismatch).
static double at_from_ensdir(const std::string& s){
  for(std::size_t p = 0; p + 2 < s.size(); ++p){
    const bool is_at = (s[p] == 'a' && s[p+1] == 't');
    const bool digit_after = (s[p+2] >= '0' && s[p+2] <= '9');
    if(is_at && digit_after) return std::stod(s.substr(p+2));
  }
  return -1.0;
}

// true if the per-hit output file exists and carries the "complete" flag
static bool hit_is_complete(const std::string& h5path){
  if(!std::filesystem::exists(h5path)) return false;
  bool c = false;
  try {
    HighFive::File f(h5path, HighFive::File::ReadOnly);
    c = f.exist("complete");
  } catch(...) {}
  return c;
}

void PrintHelp(){
  printf("jj_local_ylm_scalar_disc_stoch: stochastic DISCONNECTED scalar l=0 loop J_1 = tr[W_00 D_ov^{-1}].\n");
  printf("  MASSLESS only; forward solve only; NB_CLI hits per config solved as one mrhs block.\n");
  printf("  --gsq <x>            Wilson coupling squared (ensemble id; default 8.0)\n");
  printf("  --Nf <n>             number of fermion flavors (ensemble id; default 2)\n");
  printf("  --nu0 <x>            sea quark asymmetry (ensemble id; default 1.0)\n");
  printf("  --nu1 <x>            valence Wilson-Dirac asymmetry (operator; default nu0)\n");
  printf("  --at <x>             valence temporal spacing a_t (default: derive from ens-dir; free field 0.2)\n");
  printf("  --ens-dir <path>     sea config dir; OMIT => free field (U=1)\n");
  printf("  --nhits <n>          stochastic hits (default NB_CLI; must be a multiple of NB_CLI)\n");
  printf("  --disc-tblock <tb>   timeslices per time-dilution class (default 2; interval=Nt/tb classes)\n");
  printf("  --stride <s>         ensemble config stride (default 10)\n");
  printf("  --kmin <a> --kmax <b> config range [a,b) (default 0..1e6)\n");
  printf("  --out-suffix <s>     appended to the output dir name (smoke/variant runs; default empty)\n");
  printf("  -h, --help\n");
  exit(0);
}

void ParseArgs(int argc, char* argv[],
               double& gsq, int& Nf, double& nu0, double& nu1,
               std::string& ens_dir, int& nhits, int& disc_tblock, int& stride,
               int& kmin, int& kmax, double& at_cli, std::string& out_suffix){
  static struct option long_opts[] = {
    {"gsq",     required_argument, nullptr, 'g'},
    {"Nf",      required_argument, nullptr, 'N'},
    {"nu0",     required_argument, nullptr, 'n'},
    {"nu1",     required_argument, nullptr, 'm'},
    {"at",      required_argument, nullptr, 'A'},
    {"ens-dir", required_argument, nullptr, 'e'},
    {"nhits",   required_argument, nullptr, 'H'},
    {"disc-tblock", required_argument, nullptr, 'T'},
    {"stride",  required_argument, nullptr, 'I'},
    {"kmin",    required_argument, nullptr, 'a'},
    {"kmax",    required_argument, nullptr, 'b'},
    {"out-suffix", required_argument, nullptr, 'S'},
    {"help",    no_argument,       nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
  };
  int opt, idx;
  while((opt = getopt_long(argc, argv, "g:N:n:m:A:e:H:T:I:a:b:S:h", long_opts, &idx)) != -1){
    switch(opt){
    case 'g': gsq     = std::stod(optarg); break;
    case 'N': Nf      = std::stoi(optarg); break;
    case 'n': nu0     = std::stod(optarg); break;
    case 'm': nu1     = std::stod(optarg); break;
    case 'A': at_cli  = std::stod(optarg); break;
    case 'e': ens_dir = optarg; break;
    case 'H': nhits   = std::stoi(optarg); break;
    case 'T': disc_tblock = std::stoi(optarg); break;
    case 'I': stride  = std::stoi(optarg); break;
    case 'a': kmin    = std::stoi(optarg); break;
    case 'b': kmax    = std::stoi(optarg); break;
    case 'S': out_suffix = optarg; break;
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

  constexpr int NSTACK = NB_CLI;      // hits per block solve

  double gsq=8.0;
  int Nf=2;
  double nu0=1.0;
  double nu1=-1.0;
  // MASSLESS ONLY: the forward-only form (no furnishing, no backward solve) relies on GW at m=0.
  const double mass_re=0.0;
  const double mass_im=0.0;
  std::string ens_dir="";
  int nhits=NSTACK;
  int disc_tblock=2;
  int stride=10;
  int kmin=0;
  int kmax=1000000;
  double at_cli=-1.0;       // -1 => derive a_t from ens-dir (free field: 0.2)
  std::string out_suffix="";

  ParseArgs(argc, argv, gsq, Nf, nu0, nu1, ens_dir, nhits, disc_tblock, stride, kmin, kmax, at_cli, out_suffix);
  if(nu1 < 0.0) nu1 = nu0;

  const Complex valence_mass(mass_re, mass_im);
  const bool free_field = ens_dir.empty();

  // resolve valence temporal spacing a_t: --at overrides; else derive from the ensemble dir name
  // (so a_t ALWAYS tracks the sea ensemble); free field defaults to 0.2.
  double at = at_cli;
  if(at < 0.0) at = free_field ? 0.2 : at_from_ensdir(ens_dir);
  assert(at > 0.0 && "could not resolve a_t from ens-dir; pass --at explicitly");
  assert(disc_tblock >= 1 && Comp::Nt % disc_tblock == 0 && "Nt % disc_tblock must be 0");
  assert(nhits >= NSTACK && nhits % NSTACK == 0 && "nhits must be a multiple of NB_CLI (block width)");

  std::cout << "# gsq="<<gsq<<" Nf="<<Nf<<" nu0="<<nu0<<" nu1="<<nu1<<" at="<<at
            << " mass=("<<mass_re<<","<<mass_im<<")"
            << " ens_dir="<<(free_field?std::string("<free-field U=1>"):ens_dir)
            << " nhits="<<nhits<<" disc_tblock="<<disc_tblock
            << " NSTACK="<<NSTACK << std::endl;

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();

  int device;
  CUDA_CHECK(cudaGetDeviceCount(&device));
  cudaDeviceProp device_prop[device];
  cudaGetDeviceProperties(&device_prop[0], 0);
  std::cout << "# dev = " << device_prop[0].name << std::endl;
  CUDA_CHECK(cudaSetDevice(0));

  constexpr Idx N  = Comp::N;
  constexpr int Nt = Comp::Nt;

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  using Rng=ParallelRngExt<Base,Nt>;
  using Fermion=OverlapWMass<WilsonDirac>;

  Base base(Comp::N_REFINE);
  std::cout << "# lattice set." << std::endl;

  const double M5 = -1.0;
  if(Nt!=1) assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu1);
  std::cout << "# DW set." << std::endl;

  Gauge U(base);

  // one RNG per block column (= per hit), each reseeded from its own (esnid,k,h) string
  // OLD (single stream, one hit at a time):
  //   Rng rng(base, 1234);
  std::vector<std::unique_ptr<Rng>> rng_col(NSTACK);
  for(int c=0; c<NSTACK; c++) rng_col[c] = std::make_unique<Rng>(base, 1234);

  Fermion Dm(DW, valence_mass, 11);
  std::cout << "# overlap operator set: D_m (M5="<<M5<<")." << std::endl;

  // phi = D_ov^{-1} eta = (D^dag D)^{-1} (D^dag eta).  Forward only.
  auto f_DmH  = std::bind(&Fermion::adj_deviceAsyncLaunch_ms,  &Dm, std::placeholders::_1, std::placeholders::_2);
  // OLD (single-RHS outer CG operator; replaced by the block engine below):
  //   auto f_Dmsq = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms,  &Dm, std::placeholders::_1, std::placeholders::_2);
  //   LinOpWrapper M_DmH(f_DmH), M_Dmsq(f_Dmsq);
  //   MatPoly op_Dmsq; op_Dmsq.push_back(cplx(1.0), {&M_Dmsq});
  LinOpWrapper M_DmH(f_DmH);
  MatPoly op_DmH;
  op_DmH.push_back(cplx(1.0), {&M_DmH});

  // mrhs block engine, created ONCE (heavy ctor); holds a const ref to Dm, so the per-config
  // coefficients set by Dm.update(U) are read at solve time.
  BlockedMat<N, NSTACK, Fermion> blk_Dmsq(Dm);
  std::vector<Complex> hblk((std::size_t)N * NSTACK);      // host RHS/solution block
  std::cout << "# mrhs block engine set: NSTACK="<<NSTACK<<" (hits per block solve)." << std::endl;

  const int n_sites = static_cast<int>(base.n_sites);
  const int interval = Nt / disc_tblock;   // number of time-dilution classes

  // (1/2) tr W_00 on one timeslice = sum_n A_n Y_00(n^)  (spin trace 2 x contact 1/2).
  // Computed with the SAME mult_Ylm_real that builds the vertex, so the normalization matches J_1.
  double half_trW = 0.0;
  {
    FermionVector one;
    for(Idx i=0; i<N; i++) one.field[i] = Complex(0.0, 0.0);
    for(int n=0; n<n_sites; n++) one(0, (Idx)n, 0) = Complex(1.0, 0.0);
    one.mult_Ylm_real(0, 0, base);
    for(int n=0; n<n_sites; n++) half_trW += one(0, (Idx)n, 0).real();
  }
  std::cout << "# half_trW = (1/2) tr W_00 = " << half_trW << std::endl;

  // ---- output: data_<ESNID>/corr_ylm_disc_scalar_tb<tb>/corr.<k>.h<h>.h5
  // Separate from the vector disc dir corr_ylm_disc_tb<tb>/ (never opened here).
  // Hits are indexed per file (h<h>) and gated individually on "complete"; the RNG seed is
  // (esnid, k, h), independent of nhits and of the block width.
  std::string ens_base = ens_dir;
  if(!ens_base.empty() && ens_base.back()=='/') ens_base.pop_back();
  {
    const auto slash = ens_base.find_last_of('/');
    if(slash!=std::string::npos) ens_base = ens_base.substr(slash+1);
  }
  const std::string esnid = (free_field ? std::string("free") : ens_base)
                          + "_vmRe"+std::to_string(mass_re)+"vmIm"+std::to_string(mass_im);
  // OLD (vector disc dir):
  //   const std::string dir_out = "data_"+esnid+"/corr_ylm_disc_tb"+std::to_string(disc_tblock)+"/";
  const std::string dir_out = "data_"+esnid+"/corr_ylm_disc_scalar_tb"+std::to_string(disc_tblock)
                            + out_suffix + "/";
  std::filesystem::create_directories(dir_out);
  std::cout << "# dir_out = " << dir_out << std::endl;

  // disc J(t) is RAW (no inv4pi fold), matching the vector disc driver.
  auto write_vec = [&](HighFive::File& h5, const std::string& key, const std::vector<Complex>& C){
    std::vector<double> re(C.size()), im(C.size());
    for(size_t t=0;t<C.size();t++){ re[t]=C[t].real(); im[t]=C[t].imag(); }
    h5.createDataSet(key+"/real", re);  h5.createDataSet(key+"/imag", im);
  };

  std::vector<std::unique_ptr<FermionVector>> eta_col(NSTACK);
  for(int c=0; c<NSTACK; c++) eta_col[c] = std::make_unique<FermionVector>();
  FermionVector phi, Gamma_phi;

  const int k_ckpoint = free_field ? 1 : stride;
  const int k_lo      = free_field ? 0 : kmin;
  const int k_hi      = free_field ? 1 : kmax;

  for(int k = k_lo; k < k_hi; k += k_ckpoint){
    // Checkpoint presence FIRST -- the break on a missing ckpoint_lat must still fire before any
    // skip logic, or a strided worker would run past the end of the stream.
    std::string str_lat;
    if(!free_field){
      str_lat = ens_dir + "ckpoint_lat." + std::to_string(k);
      if(!std::filesystem::exists(str_lat)){ if(k==0) continue; else break; }
    }

    // EARLY SKIP: if EVERY hit of this config is already "complete", skip it WITHOUT U.read +
    // Dm.update (Dm.update recomputes lambda_min/lambda_max, the expensive part).
    // SAFE: the RNG is reseeded per hit from the string (esnid,k,h).
    {
      bool all_done = true;
      for(int h=0; h<nhits; h++){
        const std::string h5p = dir_out + "corr." + std::to_string(k) + ".h" + std::to_string(h) + ".h5";
        if(!hit_is_complete(h5p)){
          all_done = false;
          break;
        }
      }
      if(all_done){
        std::cout << "# skip k=" << k << " (all " << nhits << " hits done; no U.read/update)" << std::endl;
        continue;
      }
    }

    if(!free_field) U.read(str_lat);
    Dm.update(U);
    std::cout << "# k="<<k<<(free_field?" (free field)":"")
              << "  lambda_min/max="<<Dm.lambda_min<<"/"<<Dm.lambda_max<<std::endl;

    // OLD (one hit at a time):
    //   for(int h=0; h<nhits; h++){
    // NEW: hits in groups of NSTACK; group g holds hits h = h0 .. h0+NSTACK-1 as block columns.
    for(int h0=0; h0<nhits; h0+=NSTACK){
      std::vector<std::string> h5path_col(NSTACK);
      std::vector<std::string> seed_col(NSTACK);
      std::vector<char> done_col(NSTACK, 0);
      bool group_done = true;
      for(int c=0; c<NSTACK; c++){
        const int h = h0 + c;
        h5path_col[c] = dir_out + "corr." + std::to_string(k) + ".h" + std::to_string(h) + ".h5";
        seed_col[c]   = esnid + "_k" + std::to_string(k) + "_h" + std::to_string(h);
        done_col[c]   = hit_is_complete(h5path_col[c]) ? 1 : 0;
        if(!done_col[c]) group_done = false;
      }
      if(group_done){
        std::cout<<"# skip k="<<k<<" hits "<<h0<<".."<<(h0+NSTACK-1)<<" (complete)"<<std::endl;
        continue;
      }
      for(int c=0; c<NSTACK; c++){
        rng_col[c]->reseed(seed_from_string(seed_col[c]));
      }
      const auto t_grp0 = std::chrono::steady_clock::now();
      std::cout << "# k="<<k<<" hits "<<h0<<".."<<(h0+NSTACK-1)<<" of "<<nhits
                <<"  (disc_tblock="<<disc_tblock<<", "<<interval<<" time x "<<NS
                <<" spin classes, NSTACK="<<NSTACK<<")" << std::endl;

      // accumulators JS[c][t]: scalar l=0 one-point loop of hit column c (SUMMED over patterns)
      std::vector<std::vector<Complex>> JS(NSTACK, std::vector<Complex>(Nt, Complex(0,0)));

      // ===== TIME+SPIN dilution sweep =====
      for(int t_s=0; t_s<interval; t_s++){
        for(int spin=0; spin<NS; spin++){
          // block RHS: column c = D^dag eta_c  (a finished hit gets a zero column, which freezes)
          for(int c=0; c<NSTACK; c++){
            Complex* col = hblk.data() + (std::size_t)c*N;
            if(done_col[c]){
              std::fill(col, col + N, Complex(0.0, 0.0));
              continue;
            }
            eta_col[c]->time_spin_dilution(*rng_col[c], t_s, disc_tblock, spin);   // volume Z2
            op_DmH.from_cpu<N>(col, eta_col[c]->field);
          }

          // phi_c = D_ov^{-1} eta_c for all columns at once
          // OLD (single RHS):
          //   op_DmH.from_cpu<N>(tmp.field, eta.field);
          //   op_Dmsq.solve<N>(phi.field, tmp.field, Comp::TOL_OUTER);
          blk_Dmsq.solve_sq_from_cpu(hblk.data(), hblk.data(), Comp::TOL_OUTER);

          // OLD (vector disc: condensate sigma_PS + currents s1,s2,s3, l<=3) -- not computed here:
          //   acc_etadag_xi += eta_A.dag(phi);
          //   Gamma_phi = phi; Gamma_phi.mult_sigma(a); Gamma_phi.mult_Ylm_real(l, m, base);
          //   eta.accumulate_loop_raw(J[a-1][l][m+l], Gamma_phi, t_s, disc_tblock, spin);
          // NEW: scalar identity vertex, l=0:  W = A_n Y_00 phi;  JS += sum_n conj(eta) W.
          for(int c=0; c<NSTACK; c++){
            if(done_col[c]) continue;
            const Complex* col = hblk.data() + (std::size_t)c*N;
            for(Idx i=0; i<N; i++) phi.field[i] = col[i];
            Gamma_phi = phi;
            Gamma_phi.mult_Ylm_real(0, 0, base);   // folds A_n Y_00
            eta_col[c]->accumulate_loop_raw(JS[c], Gamma_phi, t_s, disc_tblock, spin);
          }
        }
      }
      const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now()-t_grp0).count();

      // ---- write: one file per hit (atomic: .tmp then rename) ----
      for(int c=0; c<NSTACK; c++){
        if(done_col[c]) continue;
        const int h = h0 + c;
        const std::string hp = "h0/";
        const std::string h5tmp = h5path_col[c] + ".tmp";
        auto h5p = std::make_unique<HighFive::File>(h5tmp,
                     HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
        HighFive::File& h5 = *h5p;
        h5.createDataSet("nhits", std::vector<int>{nhits});
        h5.createDataSet("hit",   std::vector<int>{h});
        h5.createDataSet("rng_seed", seed_col[c]);
        h5.createDataSet("disc_tblock", std::vector<int>{disc_tblock});
        h5.createDataSet("nstack", std::vector<int>{NSTACK});
        h5.createDataSet("a_t", std::vector<double>{at});
        write_vec(h5, hp+"disc/ylm/s0/l0/m0/J", JS[c]);
        h5.createDataSet(hp+"disc/ylm/s0/half_trW", std::vector<double>{half_trW});
        h5.createDataSet("complete", std::vector<int>{1});
        h5p.reset();
        std::filesystem::rename(h5tmp, h5path_col[c]);
        std::cout << "#   hit "<<h<<" written -> "<<h5path_col[c] << std::endl;
      }
      std::cout << "#   hits "<<h0<<".."<<(h0+NSTACK-1)<<" done (scalar l=0 loop J_1) ["<<secs
                <<" s total, "<<(secs/NSTACK)<<" s per hit]" << std::endl;
    } // hit groups
  } // k

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();
  return 0;
}
