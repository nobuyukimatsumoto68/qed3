// jj_scalar_condensate_eo_stoch_claude.cu   (stripped from jj_local_ylm_scalar_disc_stoch_claude.cu)
// STANDALONE driver for the SCALAR CONDENSATES <sigma_PS>, <sigma_FS> (a=0, l=0, spacetime-summed one-point
// densities = SSB order parameters).  Split OUT of the (dropped) ylm disc driver.  Dilution = ONLY even/odd
// timeslice + spin (nothing finer).  Massless / real m_F; m_P (parity) OUT OF SCOPE.
// Plan: scalar_ylm_corr_impl_plan_claude.md (2026-07-13 section).  Note: condensate_contact_massive_claude.md.
//
// Operator defs (current, unchanged): sigma_PS = eta^dag xi + xi^dag eta (naive);
//                                     sigma_FS = eta^dag xi - xi^dag (1-D_ov^dag) eta (furnished).
//
// Estimator (area-weighted A_n = dual_areas, spacetime-summed; eta = volume Z2, <eta eta^dag> = 1).  TWO
// building blocks, each needing its OWN solve per dilution pattern (phi = D_m^{-1} eta forward;
// psibar = D_m^{-dag} eta backward -- the sigma_FS term carries (1-D_ov^dag) on the BACKWARD leg, so it is
// NOT reusable from phi; see condensate_contact_massive_claude.md Eq (301,303)):
//   B_1 = tr[A D_m^{-1}]                 <- sum_d (A eta)^dag phi                         (phi = D_m^{-1} eta)
//   B_2 = tr[A (1-D_ov^dag) D_m^{-dag}]  <- sum_d (A eta)^dag w,  w = (1-D_ov^dag) psibar (psibar = D_m^{-dag} eta)
// Store the RAW loops (loop sign -, matching the note's naming):
//   etadag_xi        = -sum_patterns (A eta)^dag phi   (= -B_1)
//   xidag_1mDdag_eta = -sum_patterns (A eta)^dag w     (= -B_2)
// Analysis: <sigma_PS> = etadag_xi + conj(etadag_xi) = 2 Re(etadag_xi) ;
//           <sigma_FS> = etadag_xi - xidag_1mDdag_eta   (condensate_contact_massive_claude.md Eq (303)).
// NOTE the FS term uses (1-D_ov^dag) D_m^{-dag} (NOT the forward (1-D_ov) D_m^{-1}, which coincides only at L1
// where dual_areas/m_L are uniform; they differ at L2/L4 by operator ordering).
//
// Dilution = e/o + spin (4 patterns/hit): time_spin_dilution(rng, t_s, t_block=Nt/2, spin) -> interval=2 ->
//   t_s=0 EVEN {0,2,4,...}, t_s=1 ODD {1,3,5,...}; x spin in {0,1}.  Sum all 4 -> unbiased full spin+time
//   trace (site stochastic).  2 legs x 4 patterns = 8 CG solves/config/hit.  Per-hit atomic + "complete"-gated.
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

//------------------------------------------
#include <getopt.h>

static int seed_from_string(const std::string& s){
  std::seed_seq seq(s.begin(), s.end());
  std::uint32_t w;
  seq.generate(&w, &w + 1);
  return static_cast<int>(w);
}

void PrintHelp(){
  printf("jj_scalar_condensate_eo_stoch: stochastic SCALAR condensates <sigma_PS>,<sigma_FS> (e/o + spin dilution).\n");
  printf("  --gsq <x>            Wilson coupling squared (ensemble id; default 8.0)\n");
  printf("  --Nf <n>             number of fermion flavors (ensemble id; default 2)\n");
  printf("  --nu0 <x>            sea quark asymmetry (ensemble id; default 1.0)\n");
  printf("  --nu1 <x>            valence Wilson-Dirac asymmetry (operator; default nu0)\n");
  printf("  --mass-re <x>        valence mass Re (default 0.0; m_F is real)\n");
  printf("  --mass-im <y>        valence mass Im (default 0.0; m_P parity NOT supported here)\n");
  printf("  --ens-dir <path>     sea config dir; OMIT => free field (U=1)\n");
  printf("  --nhits <n>          stochastic hits (default 1)\n");
  printf("  --stride <s>         ensemble config stride (default 10)\n");
  printf("  --kmin <a> --kmax <b> config range [a,b) (default 0..1e6)\n");
  printf("  -h, --help\n");
  exit(0);
}

void ParseArgs(int argc, char* argv[],
               double& gsq, int& Nf, double& nu0, double& nu1,
               double& mass_re, double& mass_im,
               std::string& ens_dir, int& nhits, int& stride,
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
    {"stride",  required_argument, nullptr, 'I'},
    {"kmin",    required_argument, nullptr, 'a'},
    {"kmax",    required_argument, nullptr, 'b'},
    {"help",    no_argument,       nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
  };
  int opt, idx;
  while((opt = getopt_long(argc, argv, "g:N:n:m:r:i:e:H:I:a:b:h", long_opts, &idx)) != -1){
    switch(opt){
    case 'g': gsq     = std::stod(optarg); break;
    case 'N': Nf      = std::stoi(optarg); break;
    case 'n': nu0     = std::stod(optarg); break;
    case 'm': nu1     = std::stod(optarg); break;
    case 'r': mass_re = std::stod(optarg); break;
    case 'i': mass_im = std::stod(optarg); break;
    case 'e': ens_dir = optarg; break;
    case 'H': nhits   = std::stoi(optarg); break;
    case 'I': stride  = std::stoi(optarg); break;
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

  double gsq=8.0;  int Nf=2;  double nu0=1.0;  double nu1=-1.0;
  double mass_re=0.0, mass_im=0.0;
  std::string ens_dir="";
  int nhits=1;
  int stride=10;
  int kmin=0;
  int kmax=1000000;

  ParseArgs(argc, argv, gsq, Nf, nu0, nu1, mass_re, mass_im, ens_dir, nhits, stride, kmin, kmax);
  if(nu1 < 0.0) nu1 = nu0;

  const Complex valence_mass(mass_re, mass_im);
  const bool free_field = ens_dir.empty();
  const bool parity = (std::abs(mass_im) > 1.0e-15) && (std::abs(mass_re) <= 1.0e-15);
  assert(!parity && "m_P needs the tilde-D backward leg -- OUT OF SCOPE (massless/m_F only)");
  assert(Comp::Nt % 2 == 0 && "e/o dilution needs Nt even");

  std::cout << "# gsq="<<gsq<<" Nf="<<Nf<<" nu0="<<nu0<<" nu1="<<nu1
            << " mass=("<<mass_re<<","<<mass_im<<")"
            << " ens_dir="<<(free_field?std::string("<free-field U=1>"):ens_dir)
            << " nhits="<<nhits<<" (dilution=e/o timeslice x spin)" << std::endl;

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
  const double at = 0.2;
  if(Nt!=1) assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu1);
  std::cout << "# DW set." << std::endl;

  Gauge U(base);
  Rng rng(base, 1234);

  // D_m = D_ov + m_L (m_L = measure-weighted diagonal mass, NOT scalar).  D = D_ov (massless) furnishes the
  // GW factor (1-D_ov^dag) on the sigma_FS building block B_2 = tr[A (1-D_ov^dag) D_m^{-dag}].
  Fermion Dm(DW, valence_mass, 11);
  Fermion D (DW, Complex(0.0), 11);
  std::cout << "# overlap operators set: D_m, D_ov (M5="<<M5<<")." << std::endl;

  // phi    = D_m^{-1}   eta = op_Dmsq^{-1}(op_DmH eta)   (sigma_PS leg).
  // psibar = D_m^{-dag} eta = op_Dm(op_Dmsq^{-1} eta)    (sigma_FS leg; D_m normal so DDH=DHD).
  auto f_DmH  = std::bind(&Fermion::adj_deviceAsyncLaunch_ms,  &Dm, std::placeholders::_1, std::placeholders::_2);
  auto f_Dmsq = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms,  &Dm, std::placeholders::_1, std::placeholders::_2);
  auto f_Dm   = std::bind(&Fermion::mult_deviceAsyncLaunch_ms, &Dm, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_DmH(f_DmH), M_Dmsq(f_Dmsq), M_Dm(f_Dm);
  MatPoly op_DmH;  op_DmH.push_back(cplx(1.0), {&M_DmH});
  MatPoly op_Dmsq; op_Dmsq.push_back(cplx(1.0), {&M_Dmsq});
  MatPoly op_Dm;   op_Dm.push_back(cplx(1.0), {&M_Dm});     // D_m forward (for D_m^{-dag} = op_Dm . op_Dmsq^{-1})

  // massless D_ov^dag apply: op_oneMinusDdag : v -> (1 - D_ov^dag) v  (the sigma_FS building block carries the
  // GW factor (1-D_ov^dag) on the BACKWARD leg D_m^{-dag}: tr[A (1-D_ov^dag) D_m^{-dag}]).
  auto f_DH = std::bind(&Fermion::adj_deviceAsyncLaunch_ms, &D, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_DH(f_DH);
  MatPoly op_oneMinusDdag;
  op_oneMinusDdag.push_back(cplx( 1.0), {});       // identity
  op_oneMinusDdag.push_back(cplx(-1.0), {&M_DH});  // - D_ov^dag

  const int n_sites = static_cast<int>(base.n_sites);
  // e/o time dilution: interval = 2 classes (even/odd), t_block = Nt/2 timeslices per class.
  const int t_block  = Nt / 2;
  const int interval = 2;

  // LOCAL site area weight A_n = dual_areas (bare sigma, no kappa).
  std::vector<double> w_site(base.n_sites);
  for(Idx n=0; n<base.n_sites; n++) w_site[n]=base.dual_areas[n];

  // ---- output: data_<ESNID>/corr_condensate_eo_nhits<H>/corr.<k>.h<h>.h5
  std::string ens_base = ens_dir;
  if(!ens_base.empty() && ens_base.back()=='/') ens_base.pop_back();
  { const auto slash = ens_base.find_last_of('/'); if(slash!=std::string::npos) ens_base = ens_base.substr(slash+1); }
  const std::string esnid = (free_field ? std::string("free") : ens_base)
                          + "_vmRe"+std::to_string(mass_re)+"vmIm"+std::to_string(mass_im);
  const std::string dir_out = "data_"+esnid+"/corr_condensate_eo_nhits"+std::to_string(nhits)+"/";
  std::filesystem::create_directories(dir_out);
  std::cout << "# dir_out = " << dir_out << std::endl;

  // condensate scalars are RAW (no inv4pi fold), matching disc_claude.cu / the dilute disc.
  auto write_vec = [&](HighFive::File& h5, const std::string& key, const std::vector<Complex>& C){
    std::vector<double> re(C.size()), im(C.size());
    for(size_t t=0;t<C.size();t++){ re[t]=C[t].real(); im[t]=C[t].imag(); }
    h5.createDataSet(key+"/real", re);  h5.createDataSet(key+"/imag", im);
  };

  FermionVector eta, eta_A, tmp, phi, psibar, w;   // psibar = D_m^{-dag} eta ; w = (1-D_ov^dag) psibar (FS leg)

  const int k_ckpoint = free_field ? 1 : stride;
  const int k_lo      = free_field ? 0 : kmin;
  const int k_hi      = free_field ? 1 : kmax;

  for(int k = k_lo; k < k_hi; k += k_ckpoint){
    std::string str_lat;
    if(!free_field){
      str_lat = ens_dir + "ckpoint_lat." + std::to_string(k);
      if(!std::filesystem::exists(str_lat)){ if(k==0) continue; else break; }
    }
    // CHEAP pre-skip (BEFORE any construction): skip the whole config if every hit is already done -- only
    // the output .h5 is needed, NOT U.read or the overlap Dm.update/D.update (expensive lambda_min/max).
    // "done" = complete AND has the NEW FS key xidag_1mDdag_eta, so OLD-driver files (etadag_1mD only) are
    // recomputed automatically (no rm needed).
    {
      bool all_done = true;
      for(int h=0; h<nhits; h++){
        const std::string h5p = dir_out + "corr." + std::to_string(k) + ".h" + std::to_string(h) + ".h5";
        bool done_h = false;
        if(std::filesystem::exists(h5p)){
          try {
            HighFive::File f(h5p, HighFive::File::ReadOnly);
            // non-throwing key check (getDataSet throws + HDF5 spams stderr when absent): navigate groups.
            if(f.exist("complete") && f.exist("h0")){
              auto g0 = f.getGroup("h0");
              if(g0.exist("condensate")) done_h = g0.getGroup("condensate").exist("xidag_1mDdag_eta");
            }
          } catch(...) {}
        }
        if(!done_h){ all_done = false; break; }
      }
      if(all_done){ std::cout<<"# skip k="<<k<<" (all "<<nhits<<" hits done; no U.read/update)"<<std::endl; continue; }
    }
    if(!free_field) U.read(str_lat);
    Dm.update(U);
    D.update(U);
    std::cout << "# k="<<k<<(free_field?" (free field)":"")
              << "  lambda_min/max="<<Dm.lambda_min<<"/"<<Dm.lambda_max<<std::endl;

    for(int h=0; h<nhits; h++){
      const std::string h5path_h = dir_out + "corr." + std::to_string(k) + ".h" + std::to_string(h) + ".h5";
      if(std::filesystem::exists(h5path_h)){
        bool done_h=false;
        try {
          HighFive::File f(h5path_h, HighFive::File::ReadOnly);
          if(f.exist("complete") && f.exist("h0")){
            auto g0 = f.getGroup("h0");
            if(g0.exist("condensate")) done_h = g0.getGroup("condensate").exist("xidag_1mDdag_eta");
          }
        } catch(...) {}
        if(done_h){ std::cout<<"# skip k="<<k<<" hit "<<h<<" (done: has xidag_1mDdag_eta)"<<std::endl; continue; }
      }
      const std::string seed_str = esnid + "_k" + std::to_string(k) + "_h" + std::to_string(h);
      rng.reseed(seed_from_string(seed_str));
      const auto t_hit0 = std::chrono::steady_clock::now();
      std::cout << "# k="<<k<<" hit "<<(h+1)<<"/"<<nhits<<"  (e/o x spin = "<<interval<<" time x "<<NS
                <<" spin classes, seed='"<<seed_str<<"')" << std::endl;

      // condensate building blocks (SUMMED over the 4 dilution patterns): acc_B1 = B_1 = tr[A D_m^{-1}] ;
      // acc_B2 = B_2 = tr[A (1-D_ov^dag) D_m^{-dag}].
      Complex acc_B1(0,0);
      Complex acc_B2(0,0);

      // ===== e/o TIME x SPIN dilution sweep =====
      for(int t_s=0; t_s<interval; t_s++){
        for(int spin=0; spin<NS; spin++){
          eta.time_spin_dilution(rng, t_s, t_block, spin);   // volume Z2 on this (e/o class, spin)
          // phi    = D_m^{-1}   eta = op_Dmsq^{-1}(op_DmH eta)          (sigma_PS leg)
          // psibar = D_m^{-dag} eta = op_Dm(op_Dmsq^{-1} eta) ; w = (1-D_ov^dag) psibar   (sigma_FS leg)
          op_DmH.from_cpu<N>(tmp.field, eta.field);
          op_Dmsq.solve<N>(phi.field, tmp.field, Comp::TOL_OUTER);
          op_Dmsq.solve<N>(tmp.field, eta.field, Comp::TOL_OUTER);
          op_Dm.from_cpu<N>(psibar.field, tmp.field);
          op_oneMinusDdag.from_cpu<N>(w.field, psibar.field);

          // eta_A = A_n eta (area weight; zero off the diluted support since eta is zero there).
          for(int tt=0; tt<Nt; tt++)
            for(int n=0; n<n_sites; n++){
              const double aw = w_site[n];
              for(int i=0;i<NS;i++) eta_A(tt,(Idx)n,i) = aw*eta(tt,(Idx)n,i);
            }
          acc_B1 += eta_A.dag(phi);   // -> tr[A D_m^{-1}]                 (etadag_xi        = -B1)
          acc_B2 += eta_A.dag(w);     // -> tr[A (1-D_ov^dag) D_m^{-dag}]  (xidag_1mDdag_eta = -B2)
        }
      }

      // ---- write (fresh only; per-hit atomic via .tmp + rename) ----
      const std::string hp = "h0/";
      const std::string h5tmp = h5path_h + ".tmp";
      std::unique_ptr<HighFive::File> h5p = std::make_unique<HighFive::File>(h5tmp,
              HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
      HighFive::File& h5 = *h5p;
      h5.createDataSet("nhits",    std::vector<int>{nhits});
      h5.createDataSet("hit",      std::vector<int>{h});
      h5.createDataSet("rng_seed", seed_str);
      h5.createDataSet("dilution", std::string("eo_spin"));
      // loop sign - (matching condensate_contact_massive_claude.md Eq (301,303)).  Analysis:
      //   <sigma_PS> = etadag_xi + conj(etadag_xi) = 2 Re(etadag_xi)
      //   <sigma_FS> = etadag_xi - xidag_1mDdag_eta
      const Complex etadag_xi        = -acc_B1;
      const Complex xidag_1mDdag_eta = -acc_B2;
      write_vec(h5, hp+"condensate/etadag_xi",        std::vector<Complex>{ etadag_xi        });
      write_vec(h5, hp+"condensate/xidag_1mDdag_eta", std::vector<Complex>{ xidag_1mDdag_eta });
      h5.createDataSet("complete", std::vector<int>{1});
      h5p.reset();
      std::filesystem::rename(h5tmp, h5path_h);
      const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now()-t_hit0).count();
      std::cout << "#   hit "<<(h+1)<<" done (condensate <sigma_PS>,<sigma_FS>, e/o+spin) ["<<secs<<" s] -> "
                << h5path_h << std::endl;
    } // hits
  } // k

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();
  return 0;
}
