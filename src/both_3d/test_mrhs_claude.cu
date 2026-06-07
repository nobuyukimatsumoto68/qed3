// test_mrhs_claude.cu
// End-to-end A/B for the multi-shift inner-pole optimization (C4): solves
// op_Dmsq = D_m^dag D_m two ways -- the reference DDH (per-pole Zolotarev loop)
// vs DDH_ms (multi-shift CG, solve_multishift) -- on the same RHS. Confirms the
// solutions match to ~tol (residual + x_ms-vs-x_ref) and times both for the
// realistic end-to-end speedup.
//
// The operator construction is lifted verbatim from jj_conn_tpproj_claude.cu so
// the validated CG path is the SAME one the production code drives:
//   op_Dmsq = D_m^dag D_m,  and the jj forward solve  phi' = D_m^{-1} eta  is
//   reproduced as  tmp = D_m^dag eta  (op_DmH);  x = op_Dmsq^{-1} tmp.
// Free field only (U = 1): the solver math is gauge-independent, so no ensemble
// is needed for an exactness/benchmark test.
//
// Usage:
//   ./test_mrhs_claude.o [--mass-re x] [--mass-im y] [--tol t] [--nrhs n] [--reps r]
// Compile: handled by the both_3d Makefile (auto-picks every *.cu).

#include <typeinfo>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
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

  constexpr int N_REFINE=1;
  constexpr int NS=2;
  constexpr int Nt=128;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;

  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;   // inner Zolotarev pole solves
  const double TOL_OUTER=1.0e-5;   // outer CG default (overridable via --tol)
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

#include "overlap_wmass_claude.h"
#include "conserved_current_claude.h"

//------------------------------------------
#include <getopt.h>

void PrintHelp(){
  printf("Usage: ./test_mrhs_claude.o [options]\n");
  printf("  --mass-re <x>   valence mass Re (default: 0.0)\n");
  printf("  --mass-im <y>   valence mass Im (default: 0.0)\n");
  printf("  --tol <t>       outer CG tolerance (default: %.1e)\n", Comp::TOL_OUTER);
  printf("  --nrhs <n>      number of (looped, single-RHS) solves to validate/time (default: 1)\n");
  printf("  --reps <r>      benchmark repetitions of the full nrhs batch (default: 10)\n");
  printf("  -h, --help      show this help\n");
  exit(0);
}

void ParseArgs(int argc, char* argv[],
               double& mass_re, double& mass_im, double& tol, int& nrhs, int& reps){
  static struct option long_opts[] = {
    {"mass-re", required_argument, nullptr, 'r'},
    {"mass-im", required_argument, nullptr, 'i'},
    {"tol",     required_argument, nullptr, 't'},
    {"nrhs",    required_argument, nullptr, 'n'},
    {"reps",    required_argument, nullptr, 'p'},
    {"help",    no_argument,       nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
  };
  int opt, idx;
  while((opt = getopt_long(argc, argv, "r:i:t:n:p:h", long_opts, &idx)) != -1){
    switch(opt){
    case 'r': mass_re = std::stod(optarg); break;
    case 'i': mass_im = std::stod(optarg); break;
    case 't': tol     = std::stod(optarg); break;
    case 'n': nrhs    = std::stoi(optarg); break;
    case 'p': reps    = std::stoi(optarg); break;
    case 'h':
    case '?':
    default:  PrintHelp(); break;
    }
  }
}
//------------------------------------------

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);

  double mass_re=0.0, mass_im=0.0, tol=Comp::TOL_OUTER;
  int nrhs=1, reps=10;
  ParseArgs(argc, argv, mass_re, mass_im, tol, nrhs, reps);
  const Complex valence_mass(mass_re, mass_im);

  std::cout << "# test_mrhs: mass=("<<mass_re<<","<<mass_im<<")  tol="<<tol
            << "  nrhs="<<nrhs<<"  reps="<<reps << std::endl;

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();

  int device;
  CUDA_CHECK(cudaGetDeviceCount(&device));
  cudaDeviceProp device_prop[device];
  cudaGetDeviceProperties(&device_prop[0], 0);
  std::cout << "# dev = " << device_prop[0].name << std::endl;
  CUDA_CHECK(cudaSetDevice(0));

  constexpr Idx N  = Comp::N;

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;
  using Rng=ParallelRngExt<Base,Comp::Nt>;
  using Fermion=OverlapWMass<WilsonDirac>;

  // ---- operators (lifted from jj_conn_tpproj_claude.cu, free field) ----------------
  Base base(Comp::N_REFINE);
  const double M5 = -1.0;
  const double at = 0.2;
  const double nu1 = 1.0;
  if(Comp::Nt!=1) assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu1);

  Gauge U(base);                       // free field (U = 1)
  Rng rng(base, 1234);

  Fermion Dm(DW, valence_mass, 21);
  Dm.update(U);
  std::cout << "# overlap D_m set: lambda_min/max="<<Dm.lambda_min<<"/"<<Dm.lambda_max<<std::endl;

  // D_m = D_ov + m:  op_DmH forms the RHS (D_m^dag eta), op_Dmsq = D_m^dag D_m is the CG operator
  auto f_Dm   = std::bind(&Fermion::mult_deviceAsyncLaunch, &Dm, std::placeholders::_1, std::placeholders::_2);
  auto f_DmH  = std::bind(&Fermion::adj_deviceAsyncLaunch,  &Dm, std::placeholders::_1, std::placeholders::_2);
  auto f_Dmsq = std::bind(&Fermion::DDH_deviceAsyncLaunch,  &Dm, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_Dm(f_Dm), M_DmH(f_DmH), M_Dmsq(f_Dmsq);
  MatPoly op_DmH;  op_DmH.push_back(cplx(1.0), {&M_DmH});
  MatPoly op_Dmsq; op_Dmsq.push_back(cplx(1.0), {&M_Dmsq});

  // multi-shift CG operator: same op_Dmsq = D_m^dag D_m but the inner Zolotarev pole
  // solves go through solve_multishift (DDH_ms instead of DDH). A/B against op_Dmsq.
  auto f_Dmsq_ms = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms, &Dm, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_Dmsq_ms(f_Dmsq_ms);
  MatPoly op_Dmsq_ms; op_Dmsq_ms.push_back(cplx(1.0), {&M_Dmsq_ms});

  // ---- buffers (host-pinned FermionVector, same access pattern as jj) --------------
  // tmp[c]    = D_m^dag eta_c  (the CG RHS, exactly as in jj's forward solve)
  // x_ref[c]  = op_Dmsq^{-1} tmp[c]  via host-sync solve()   (reference baseline)
  std::vector<FermionVector> eta(nrhs), tmp(nrhs), x_ref(nrhs), x_ms(nrhs), res(nrhs);
  for(int c=0; c<nrhs; c++){
    eta[c].fill_z2_source(rng);
    op_DmH.from_cpu<N>(tmp[c].field, eta[c].field);
  }

  // ===== PREVIOUS TEST (baseline, reference-only) -- kept commented per workflow; =====
  // ===== superseded by the DDH vs DDH_ms A/B beneath. =================================
  // std::cout << "\n# === correctness (reference solve residual) ===" << std::endl;
  // double worst_res=0.0;
  // for(int c=0; c<nrhs; c++){
  //   op_Dmsq.solve<N>(x_ref[c].field, tmp[c].field, tol);
  //   op_Dmsq.from_cpu<N>(res[c].field, x_ref[c].field);
  //   double rn=0.0, bn=0.0;
  //   for(Idx i=0;i<N;i++){ rn += std::norm(res[c].field[i]-tmp[c].field[i]); bn += std::norm(tmp[c].field[i]); }
  //   const double rr = std::sqrt(rn) / (bn>0.0 ? std::sqrt(bn) : 1.0);
  //   worst_res=std::max(worst_res,rr);
  //   std::cout << "#  rhs "<<c<<": resid="<<rr << std::endl;
  // }
  // std::cout << "# WORST resid="<<worst_res<<"  ("<<(worst_res<10.0*tol?"PASS":"CHECK")<<")" << std::endl;
  // double sum=0.0;
  // for(int rp=0; rp<reps; rp++){
  //   CUDA_CHECK(cudaDeviceSynchronize());
  //   const auto t0 = std::chrono::steady_clock::now();
  //   for(int c=0; c<nrhs; c++) op_Dmsq.solve<N>(x_ref[c].field, tmp[c].field, tol);
  //   CUDA_CHECK(cudaDeviceSynchronize());
  //   sum += 1.0e3*std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
  // }
  // std::cout << "# reference solve() : "<<(sum/reps)<<" ms / batch"<<std::endl;
  // ===================================================================================

  // ---- correctness: reference DDH vs multi-shift DDH_ms, per RHS -------------------
  // op_Dmsq drives the inner Zolotarev pole solves via the original per-pole loop;
  // op_Dmsq_ms drives them via multi-shift CG (solve_multishift). Both solve the SAME
  // system to the same tol => x_ms == x_ref to ~tol and both residuals sit at ~tol.
  std::cout << "\n# === correctness (reference DDH vs multishift DDH_ms) ===" << std::endl;
  double worst_res=0.0, worst_resms=0.0, worst_diff=0.0;
  for(int c=0; c<nrhs; c++){
    op_Dmsq.solve   <N>(x_ref[c].field, tmp[c].field, tol);
    op_Dmsq_ms.solve<N>(x_ms[c].field,  tmp[c].field, tol);
    // residual ||op_Dmsq x - tmp|| / ||tmp|| (apply via op_Dmsq; same operator)
    auto resid = [&](const FermionVector& x)->double{
      op_Dmsq.from_cpu<N>(res[c].field, x.field);
      double rn=0.0, bn=0.0;
      for(Idx i=0;i<N;i++){ rn += std::norm(res[c].field[i]-tmp[c].field[i]); bn += std::norm(tmp[c].field[i]); }
      return std::sqrt(rn) / (bn>0.0 ? std::sqrt(bn) : 1.0);
    };
    const double rr = resid(x_ref[c]);
    const double rm = resid(x_ms[c]);
    double dabs=0.0, nref=0.0;
    for(Idx i=0;i<N;i++){ dabs=std::max(dabs, std::abs(x_ref[c].field[i]-x_ms[c].field[i])); nref+=std::norm(x_ref[c].field[i]); }
    const double drel = (nref>0.0 ? dabs/std::sqrt(nref) : dabs);
    worst_res=std::max(worst_res,rr); worst_resms=std::max(worst_resms,rm); worst_diff=std::max(worst_diff,drel);
    std::cout << "#  rhs "<<c<<": resid(ref)="<<rr<<"  resid(ms)="<<rm<<"  rel|x_ms-x_ref|="<<drel<<std::endl;
  }
  std::cout << "# WORST resid(ref)="<<worst_res<<"  resid(ms)="<<worst_resms<<"  rel(ms-ref)="<<worst_diff
            << "  ("<<((worst_resms<10.0*tol && worst_diff<1.0e-6)?"PASS":"CHECK")<<")" << std::endl;

  // ===== PREVIOUS TEST (timing-only benchmark) -- kept commented per workflow; =====
  // ===== superseded by the benchmark + per-rep diff beneath. =======================
  // std::cout << "\n# === benchmark (avg over "<<reps<<" reps of the "<<nrhs<<"-solve batch) ===" << std::endl;
  // auto bench = [&](bool ms)->double{
  //   double sum=0.0;
  //   for(int rp=0; rp<reps; rp++){
  //     CUDA_CHECK(cudaDeviceSynchronize());
  //     const auto t0 = std::chrono::steady_clock::now();
  //     for(int c=0; c<nrhs; c++){
  //       if(ms) op_Dmsq_ms.solve<N>(x_ms[c].field,  tmp[c].field, tol);
  //       else   op_Dmsq.solve   <N>(x_ref[c].field, tmp[c].field, tol);
  //     }
  //     CUDA_CHECK(cudaDeviceSynchronize());
  //     sum += 1.0e3*std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
  //   }
  //   return sum/reps;
  // };
  // const double t_ref = bench(false);
  // const double t_ms  = bench(true);
  // std::cout << "# reference DDH     : "<<t_ref<<" ms / batch  ("<<t_ref/nrhs<<" ms/solve)"<<std::endl;
  // std::cout << "# multishift DDH_ms : "<<t_ms <<" ms / batch  ("<<t_ms /nrhs<<" ms/solve)"<<std::endl;
  // std::cout << "# speedup (ref/ms)  : "<<(t_ms>0.0?t_ref/t_ms:0.0)<<"x"<<std::endl;
  // =================================================================================

  // ---- benchmark + PER-REP diff: reference DDH vs multishift DDH_ms ----------------
  // Re-check x_ms vs x_ref EVERY rep (not just once): a memory issue (stale/reused
  // d_Zblock/d_Yblock, a race, uninitialized scratch) would surface as a later rep whose
  // diff exceeds tol even though rep 0 passed. ref and ms are timed in separate synced
  // windows; the diff is computed on the host afterwards so it does not pollute timing.
  std::cout << "\n# === benchmark + per-rep diff (ref DDH vs multishift DDH_ms, "<<reps<<" reps) ===" << std::endl;
  double tsum_ref=0.0, tsum_ms=0.0, worst_rep_diff=0.0;
  bool all_reps_pass=true;
  for(int rp=0; rp<reps; rp++){
    CUDA_CHECK(cudaDeviceSynchronize());
    auto t0 = std::chrono::steady_clock::now();
    for(int c=0; c<nrhs; c++) op_Dmsq.solve<N>(x_ref[c].field, tmp[c].field, tol);
    CUDA_CHECK(cudaDeviceSynchronize());
    tsum_ref += 1.0e3*std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();

    CUDA_CHECK(cudaDeviceSynchronize());
    auto t1 = std::chrono::steady_clock::now();
    for(int c=0; c<nrhs; c++) op_Dmsq_ms.solve<N>(x_ms[c].field, tmp[c].field, tol);
    CUDA_CHECK(cudaDeviceSynchronize());
    tsum_ms += 1.0e3*std::chrono::duration<double>(std::chrono::steady_clock::now()-t1).count();

    double rep_diff=0.0;
    for(int c=0; c<nrhs; c++){
      double dabs=0.0, nref=0.0;
      for(Idx i=0;i<N;i++){ dabs=std::max(dabs, std::abs(x_ref[c].field[i]-x_ms[c].field[i])); nref+=std::norm(x_ref[c].field[i]); }
      rep_diff = std::max(rep_diff, nref>0.0 ? dabs/std::sqrt(nref) : dabs);
    }
    const bool ok = rep_diff < 1.0e-6;
    all_reps_pass = all_reps_pass && ok;
    worst_rep_diff = std::max(worst_rep_diff, rep_diff);
    std::cout << "#  rep "<<std::setw(2)<<rp<<": rel|x_ms-x_ref|="<<rep_diff<<"  ("<<(ok?"ok":"FAIL")<<")"<<std::endl;
  }
  const double t_ref = tsum_ref/reps, t_ms = tsum_ms/reps;
  std::cout << "# reference DDH     : "<<t_ref<<" ms / batch  ("<<t_ref/nrhs<<" ms/solve)"<<std::endl;
  std::cout << "# multishift DDH_ms : "<<t_ms <<" ms / batch  ("<<t_ms /nrhs<<" ms/solve)"<<std::endl;
  std::cout << "# speedup (ref/ms)  : "<<(t_ms>0.0?t_ref/t_ms:0.0)<<"x"<<std::endl;
  std::cout << "# all "<<reps<<" reps within tol: "<<(all_reps_pass?"PASS":"FAIL")
            << "  (worst rel|x_ms-x_ref|="<<worst_rep_diff<<")"<<std::endl;

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();
  return 0;
}
