// eig_lanczos_claude.cu
// _claude: entry point for the IRL + Chebyshev low-mode eigensolver (includes/lanczos_claude.h).
// Computes the Nk SMALLEST eigenvalues (and eigenvectors) of A = (D_ov + m)^\dagger (D_ov + m) via the
// Implicitly Restarted Lanczos algorithm with Chebyshev acceleration.  This is the eigensolver for the
// L4 deflation / LMA program (deflation_lma_qed3_impl_plan_claude.md).  REVISED chunk 1 = ADOPT + VALIDATE
// the pre-existing IRL (lanczos_claude.h, 2026-05-22) rather than rewrite; the dense cross-check is
// eig_wmass_val_claude.cu (geev of the SAME A).
//
// Algorithm: Sorensen implicitly restarted Lanczos (Bai et al., "Templates for the solution of algebraic
// eigenvalue problems", SIAM 2000; Grid ImplicitlyRestartedLanczos.h consulted).  Chebyshev filter maps the
// A-eigenvalue window [0, alpha^2] to the amplified band so IRL converges the SMALLEST modes first.
//
// Build: -DLREF=<1|2|4> selects the lattice (N_REFINE); -DNT_CLI=<Nt>.  See tmp_eig_validate_claude.sh.
// Usage: ./bin [mass_re] [mass_im] [alpha] [beta] [cheb_degree] [Nk] [Nm]
//   Chebyshev window per lanczos-2.pdf Eq (3.3), alpha > beta:  alpha = spectral UPPER bound (auto 2+|m|;
//   D_ov singular values in [0,2] on the GW circle),  beta = upper edge of the wanted LOW band (spectrum-
//   dependent; pass explicitly, default 0.6*alpha).  cheb_degree EVEN (T_n filters |q|>1 only for even n).
//   Free-field (cold gauge); config-read added chunk 2.
// Output: eig_lanczos_L<L>_nt<Nt>_claude.dat  (lowest Nk eigenvalues of A, ascending; cols: i  eval).

#include <typeinfo>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <filesystem>
#include <cstdint>
#include <complex>
#include <array>
#include <vector>
#include <map>
#include <memory>
#include <string>
#include <cmath>
#include <random>
#include <utility>
#include <functional>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

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

#ifndef LREF
#define LREF 1
#endif
#ifndef NT_CLI
#define NT_CLI 32
#endif
  constexpr int N_REFINE=LREF;
  constexpr int NS=2;
  constexpr int Nt=NT_CLI;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;
  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-8;
}

// _fermilab: ABSOLUTE geometry path (FNAL) -- run CWD is the /lustre2 output dir, relative breaks.
const std::string dir = "/project/affine/nmatsum/qed3/geometry/data/";

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

#include "overlap_wmass_claude.h"
#include "frozen_window_claude.h"       // FIXED Zolotarev window (lmin,lmax) per L -- same source of truth as HMC
#include "includes/lanczos_claude.h"    // ChebyshevFilter, IRLState, calc(...) -- the pre-existing IRL

// _claude: Hermitian H = D_W^dag D_W applied device-side, where the BARE Wilson D_W = M_DW + M (M5=-1 is
// baked into the overlap's M_DW, so M_DW = D_W - M).  Used by the DW_LOWMODES path to reach the low Wilson
// modes at L2/L4 (where dense geev is infeasible) via Hermitian IRL -- CSR matvecs only, no Zolotarev solves.
struct LinOpWilsonSq : public LinOp {
  using T = CuC;
  const CSR<Comp::N>& mdw;
  const CSR<Comp::N>& mdwh;
  const double Mshift;
  CuC* d_tmp;
  LinOpWilsonSq(const CSR<Comp::N>& mdw_, const CSR<Comp::N>& mdwh_, const double M_)
    : mdw(mdw_), mdwh(mdwh_), Mshift(M_)
  {
    CUDA_CHECK(cudaMalloc(&d_tmp, Comp::N * CD));
  }
  ~LinOpWilsonSq()
  {
    CUDA_CHECK(cudaFree(d_tmp));
  }
  void operator()(T* d_res, const T* d_v) const {
    constexpr Idx N = Comp::N;
    mdw(d_tmp, d_v);                                                                    // d_tmp = M_DW v
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_tmp, Mshift, d_v, d_tmp);   // d_tmp += M v = D_W v
    CUDA_CHECK(cudaDeviceSynchronize());
    mdwh(d_res, d_tmp);                                                                 // d_res = M_DWH (D_W v)
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_res, Mshift, d_tmp, d_res); // d_res += M(D_W v) = D_W^dag D_W v
    CUDA_CHECK(cudaDeviceSynchronize());
  }
  void Async(T* d_res, const T* d_v, const cudaStream_t stream) const {
    (void)stream;
    (*this)(d_res, d_v);
  }
};

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(12);
  std::clog << std::scientific << std::setprecision(12);

  int device;
  CUDA_CHECK(cudaGetDeviceCount(&device));
  cudaDeviceProp dp;
  cudaGetDeviceProperties(&dp, 0);
  std::cout << "# dev = " << dp.name << std::endl;
  CUDA_CHECK(cudaSetDevice(0));

  constexpr Idx N = Comp::N;

  // Chebyshev convention = lanczos-2.pdf Eq (3.3): q(l;alpha,beta) = (2 l^2 - (alpha^2+beta^2))/(alpha^2-beta^2),
  // ALPHA > BETA, with ALPHA = spectral upper bound (>= largest singular value) and BETA = upper edge of the
  // WANTED low band (~ largest of the m smallest singular values).  Wanted modes (l < beta) map to q < -1;
  // T_n filters |q|>1, so the degree n MUST be EVEN (T_n(q)>0 for |q|>1 only when n even -- see the PDF).
  // NOTE: the ChebyshevFilter comment in lanczos_claude.h has alpha/beta NAMED backwards; the FORMULA it
  // evaluates is Eq (3.3), so pass alpha=upper bound, beta=lower(wanted) edge as here.
  double mass_re = 0.0;
  double mass_im = 0.0;
  double alpha_cheb = -1.0;   // <=0 -> auto = 1 + |m|  (massless overlap spectral upper limit = 1)
  double beta_cheb  = -1.0;   // <=0 -> auto = 2*Dov.lambda_min  (upper edge of the wanted low band)
  int    cheb_degree = 12;    // EVEN
  int    Nk = 12;             // few genuine low modes at L1/L2
  int    Nm = 48;
  if(argc>1) mass_re    = atof(argv[1]);
  if(argc>2) mass_im    = atof(argv[2]);
  if(argc>3) alpha_cheb = atof(argv[3]);
  if(argc>4) beta_cheb  = atof(argv[4]);
  if(argc>5) cheb_degree= atoi(argv[5]);
  if(argc>6) Nk         = atoi(argv[6]);
  if(argc>7) Nm         = atoi(argv[7]);
  const Complex mass = Complex(mass_re, mass_im);

  if(cheb_degree % 2 != 0){
    std::cout << "# NOTE: cheb_degree "<<cheb_degree<<" is ODD; bumping to "<<(cheb_degree+1)
              << " (EVEN required so T_n filters |q|>1, per lanczos-2.pdf)" << std::endl;
    cheb_degree += 1;
  }

  std::cout << "# L(N_REFINE)="<<Comp::N_REFINE<<" Nt="<<Comp::Nt<<" N="<<N
            << "  mass="<<mass<<"  Nk="<<Nk<<" Nm="<<Nm<<" deg="<<cheb_degree<<std::endl;
  assert(Nm > Nk && Nm <= N && "need Nk < Nm <= N");

  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;
  using Fermion=OverlapWMass<WilsonDirac>;

  Base base(Comp::N_REFINE);
  const double r = 1.0;
  const double M5 = -1.0;
  const double at = 0.2;         // MUST match the ensemble (dir "at0.200000") and production (jj/hmc use
                                 // at=0.2): the frozen window edges below were chosen from the at=0.2 Wilson
                                 // scan, so at=0.03125 would push D_W's spectrum outside [lmin,lmax] -> wrong
                                 // sign.  (Earlier "lambda_max blows up" note was about the per-config window,
                                 // which we no longer use -- the window is FROZEN.)
  const double nu0 = 1.0;
  WilsonDirac DW(base, 0.0, r, M5, at, nu0);

  Gauge U(base);                 // cold (free field) by default; CONFIG_LAT env overrides with a real config
  const char* cfg = std::getenv("CONFIG_LAT");
  if(cfg != nullptr && std::string(cfg).size() > 0){
    assert(std::filesystem::exists(cfg) && "CONFIG_LAT path not found");
    U.read(cfg);
    std::cout << "# read gauge config: " << cfg << std::endl;
  }
  else{
    std::cout << "# free field (cold gauge)" << std::endl;
  }
  Fermion Dov(DW, mass, 21);
  double lmin_frozen, lmax_frozen;
  frozen_window(Comp::N_REFINE, lmin_frozen, lmax_frozen);   // FIXED Zolotarev window per L (matches HMC)
  Dov.set_lambda(lmin_frozen, lmax_frozen);                 // rebuild Zolotarev on frozen k + FREEZE (no per-config recompute)
  Dov.update(U);
  std::cout << "# Dov (FROZEN window): lambda_min/max = " << Dov.lambda_min << " / " << Dov.lambda_max
            << "  k = " << Dov.lambda_min/Dov.lambda_max << "  Delta = " << Dov.Delta()
            << "  [Wilson sign-function window, FROZEN per L]" << std::endl;

  // ---- optional INVERSE ITERATION check (env INV_ITER) -- independent of the Chebyshev/IRL machinery ----
  // v_{k+1} = A^{-1} v_k / ||.|| converges to the SMALLEST eigenvector of A = D_m^dag D_m; the Rayleigh
  // quotient <w,A w>/<w,w> = <w,v>/<w,w> (since A w = v) -> smallest eigenvalue.  Uses the SAME frozen-window
  // Dov the IRL uses, via the multishift solve.  Ground truth: does A actually have small eigenvalues?
  if(std::getenv("INV_ITER") != nullptr){
    auto f_Dmsq = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms, &Dov, std::placeholders::_1, std::placeholders::_2);
    LinOpWrapper M_Dmsq(f_Dmsq);
    MatPoly op_Dmsq;
    op_Dmsq.push_back(cplx(1.0), {&M_Dmsq});

    std::vector<Complex> v(N), w(N);
    std::mt19937 g(4321u);
    std::uniform_real_distribution<double> uu(-0.5, 0.5);
    for(Idx i=0;i<N;i++) v[i] = Complex(uu(g), uu(g));

    const int n_inv = 25;
    const double inv_tol = 1.0e-8;
    std::cout << "# INVERSE ITERATION on A=D_m^dag D_m (frozen window); lambda -> smallest eigenvalue" << std::endl;
    double lam = 0.0;
    for(int it=0; it<n_inv; it++){
      double nn = 0.0;
      for(Idx i=0;i<N;i++) nn += std::norm(v[i]);
      nn = std::sqrt(nn);
      for(Idx i=0;i<N;i++) v[i] /= nn;
      op_Dmsq.solve<N>(w.data(), v.data(), inv_tol);          // w = A^{-1} v
      Complex wv(0.0, 0.0);
      double ww = 0.0;
      for(Idx i=0;i<N;i++){ wv += std::conj(w[i]) * v[i]; ww += std::norm(w[i]); }
      lam = wv.real() / ww;                                    // <w,Aw>/<w,w> = <w,v>/<w,w>
      std::cout << "#   iter " << std::setw(2) << it << "  lambda_min ~ " << lam
                << "  singular ~ " << std::sqrt(std::max(0.0, lam)) << std::endl;
      v = w;
    }
    std::cout << "# INVERSE ITERATION result: smallest eig(A) ~ " << lam
              << "  (singular value ~ " << std::sqrt(std::max(0.0, lam)) << ")" << std::endl;
    for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
    return 0;
  }

  cublasHandle_t cublas_handle;
  CUBLAS_CHECK(cublasCreate(&cublas_handle));

  // ---- optional WILSON LOW-MODE spectrum (env DW_LOWMODES) -- the low complex eigenvalues of the BARE
  // Wilson D_W = M_DW + M at L2/L4 (dense geev infeasible).  IRL on the Hermitian H = D_W^dag D_W finds the
  // smallest-|D_W| (near-origin) modes; Rayleigh-Ritz B = V^dag D_W V (small dense geev) recovers the complex
  // D_W eigenvalues near the origin -> plot with the domain wall M=1 to see if they cross it at finer L.
  // "lambda<2" = |D_W| < 2 -> beta = 2.  H is CSR-only (cheap; no Zolotarev), so L4 is feasible.
  if(std::getenv("DW_LOWMODES") != nullptr){
    const double M_wall = 1.0;                                // D_W = M_DW + M (M5=-1 baked into M_DW)
    LinOpWilsonSq HermOpW(Dov.M_DW, Dov.M_DWH, M_wall);

    CuC* d_srcW;
    CUDA_CHECK(cudaMalloc(&d_srcW, N * CD));
    {
      std::vector<Complex> h(N);
      std::mt19937 g(777u);
      std::uniform_real_distribution<double> uu(-0.5, 0.5);
      for(Idx i=0;i<N;i++) h[i] = Complex(uu(g), uu(g));
      CUDA_CHECK(cudaMemcpy(d_srcW, reinterpret_cast<const CuC*>(h.data()), N * CD, H2D));
    }

    // alpha = spectral bound on |D_W| by 40-step power iteration on H (safe upper bound); beta = 2.
    double a_dw = alpha_cheb;
    if(a_dw <= 0.0){
      CuC* du;
      CuC* dw2;
      CUDA_CHECK(cudaMalloc(&du,  N * CD));
      CUDA_CHECK(cudaMalloc(&dw2, N * CD));
      CUDA_CHECK(cudaMemcpy(du, d_srcW, N * CD, D2D));
      double maxH = 0.0;
      for(int it=0; it<40; it++){
        CuC nn;
        CUBLAS_CHECK(cublasZdotc(cublas_handle, N, du, 1, du, 1, &nn));
        const double inv = 1.0/std::sqrt(cuCreal(nn));
        CUBLAS_CHECK(cublasZdscal(cublas_handle, N, &inv, du, 1));
        HermOpW(dw2, du);
        CuC rq;
        CUBLAS_CHECK(cublasZdotc(cublas_handle, N, du, 1, dw2, 1, &rq));
        maxH = cuCreal(rq);
        CUDA_CHECK(cudaMemcpy(du, dw2, N * CD, D2D));
      }
      a_dw = std::sqrt(std::max(0.0, maxH)) * 1.05;
      CUDA_CHECK(cudaFree(du));
      CUDA_CHECK(cudaFree(dw2));
    }
    const double b_dw = (beta_cheb > 0.0) ? beta_cheb : 2.0;
    std::cout << "# DW_LOWMODES: IRL on D_W^dag D_W  alpha="<<a_dw<<" beta="<<b_dw<<" deg="<<cheb_degree
              << " Nk="<<Nk<<" Nm="<<Nm<<"  M="<<M_wall << std::endl;

    ChebyshevFilter ChebOpW(HermOpW, cheb_degree, a_dw, b_dw);
    IRLState sw(Nk, Nm);
    Timer tW;
    // Spectrum-mapping study (NOT deflation) -- positions only.  The eigenVALUES stabilize in a few restarts
    // but the residual test (eigenVECTOR quality) never hits 1e-6 for the clustered low modes, so all_conv
    // never fires and a long run breaks down (beta->0 -> NaN).  So we STOP after a few restarts (MaxIter
    // small, well before the breakdown) and use the finite, eigenvalue-converged Ritz vectors.  eresid loose.
    calc(ChebOpW, HermOpW, sw, d_srcW, 1.0e-6, 8, cublas_handle);
    const int nk = (int)sw.eval.size();
    std::cout << "# IRL(D_W) done in "<<tW.currentSeconds()<<" s; subspace dim="<<nk << std::endl;

    // Rayleigh-Ritz for the (non-Hermitian) D_W on the converged subspace: B_ij = <v_i, D_W v_j>,
    // D_W v = M_DW v + M v.  Eigenvalues of the nk x nk B are the complex D_W Ritz values.
    std::vector<std::vector<Complex>> Vh(nk, std::vector<Complex>(N));
    std::vector<std::vector<Complex>> Wh(nk, std::vector<Complex>(N));
    CuC* dwv;
    CUDA_CHECK(cudaMalloc(&dwv, N * CD));
    for(int j=0;j<nk;j++){
      CUDA_CHECK(cudaMemcpy(Vh[j].data(), sw.d_evec[j], N * CD, D2H));
      Dov.M_DW(dwv, sw.d_evec[j]);                                     // dwv = M_DW v_j
      CUDA_CHECK(cudaDeviceSynchronize());
      CUDA_CHECK(cudaMemcpy(Wh[j].data(), dwv, N * CD, D2H));
      for(Idx i=0;i<N;i++) Wh[j][i] += M_wall * Vh[j][i];             // + M v_j -> D_W v_j
    }
    CUDA_CHECK(cudaFree(dwv));

    Eigen::MatrixXcd B(nk, nk);
    for(int i=0;i<nk;i++){
      for(int j=0;j<nk;j++){
        Complex acc(0.0, 0.0);
        for(Idx kk=0;kk<N;kk++) acc += std::conj(Vh[i][kk]) * Wh[j][kk];
        B(i,j) = acc;
      }
    }
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces(B, false);
    const Eigen::VectorXcd evs = ces.eigenvalues();

    const std::string outname = "eig_dw_irl_L"+std::to_string(Comp::N_REFINE)
                              + "_nt"+std::to_string(Comp::Nt)+"_claude.dat";
    std::ofstream out(outname);
    out << std::scientific << std::setprecision(14);
    out << "# IRL+Rayleigh-Ritz low spectrum of BARE Wilson D_W (= M_DW + M, M="<<M_wall<<")  L="
        << Comp::N_REFINE<<" Nt="<<Comp::Nt<<"  subspace="<<nk << std::endl;
    out << "# i   Re(D_W)   Im(D_W)   |D_W|" << std::endl;
    for(int i=0;i<nk;i++){
      const double re = evs(i).real();
      const double im = evs(i).imag();
      out << i << "   " << re << "   " << im << "   " << std::sqrt(re*re + im*im) << std::endl;
    }
    out.close();
    std::cout << "# wrote "<<outname<<" ("<<nk<<" D_W Ritz values near the origin)" << std::endl;

    CUBLAS_CHECK(cublasDestroy(cublas_handle));
    CUDA_CHECK(cudaFree(d_srcW));
    for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
    return 0;
  }

  LinOpDHDWrapper<Fermion> HermOp(Dov);                       // A = (D+m)^\dagger (D+m)

  // Chebyshev window (SINGULAR-value edges of D_m; A-eigenvalue is the square).
  //   alpha = spectral UPPER bound.  D_ov eigenvalues sit on the GW circle |z-1|=1, so its singular values
  //           are in [0, 2]; with mass D_m=D_ov+m they are in [|m|, 2+|m|].  EXACT bound: alpha = 2 + |m|.
  //           (Confirmed from the dense geev spectrum: max singular = 1.9993 for massless.  NOT 1.)
  //   beta  = upper edge of the wanted LOW band (a singular value); SPECTRUM-DEPENDENT.  Free-field is GAPPED
  //           with NO near-zero modes (smallest singular ~0.98), so beta ~ 1.1-1.4 there; the interacting L4
  //           near-zero tail wants a SMALL beta.  Dov.lambda_min is the WILSON D_W bound (~1.36 here), a
  //           DIFFERENT scale -- do not use it as an overlap edge.  Pass beta explicitly; default 0.6*alpha.
  if(alpha_cheb <= 0.0) alpha_cheb = 2.0 + std::abs(mass);
  if(beta_cheb  <= 0.0) beta_cheb  = 0.6 * alpha_cheb;
  assert(alpha_cheb > beta_cheb && "Cheby: need alpha (spectral bound) > beta (wanted-band edge)");
  std::cout << "# Chebyshev: alpha(spectral bound)="<<alpha_cheb
            << "  beta(wanted-band edge)="<<beta_cheb<<"  degree="<<cheb_degree<<" (even)"<<std::endl;

  ChebyshevFilter ChebOp(HermOp, cheb_degree, alpha_cheb, beta_cheb);
  IRLState s(Nk, Nm);

  // Random initial vector (fixed seed -> reproducible).
  CuC* d_src;
  CUDA_CHECK(cudaMalloc(&d_src, N * CD));
  {
    std::vector<Complex> h_src(N);
    std::mt19937 gen(20260715u);
    std::uniform_real_distribution<double> uni(-0.5, 0.5);
    for(Idx i=0;i<N;i++) h_src[i] = Complex(uni(gen), uni(gen));
    CUDA_CHECK(cudaMemcpy(d_src, reinterpret_cast<const CuC*>(h_src.data()), N * CD, H2D));
  }

  Timer timer;
  calc(ChebOp, HermOp, s, d_src, 1.0e-8, 200, cublas_handle);
  std::cout << "# IRL done in " << timer.currentSeconds() << " s;  converged evals = "
            << s.eval.size() << std::endl;

  // SELF-VALIDATION (replaces the dense geev cross-check): for each returned Ritz pair recompute the
  // Rayleigh quotient lam = <u,Au>/<u,u> and the RELATIVE residual ||A u - lam u|| / (lam ||u||).  This is
  // the mathematical DEFINITION of an eigenpair, so a small residual PROVES (lam,u) solves A u = lam u to
  // that tolerance -- no O(N^2) dense reference needed.  (The cheap free-field dense anchors that these are
  // the SMALLEST eigenvalues; the Chebyshev window + true bound alpha=2 guarantees the low end otherwise.)
  const int nconv = (int)s.eval.size();
  std::vector<std::pair<double,double>> lam_res(nconv);   // (lambda, relative residual)
  {
    CuC* d_w;
    CUDA_CHECK(cudaMalloc(&d_w, N * CD));
    for(int i=0;i<nconv;i++){
      HermOp(d_w, s.d_evec[i]);                            // w = A u_i
      CuC uu_c, uw_c;
      CUBLAS_CHECK(cublasZdotc(cublas_handle, N, s.d_evec[i], 1, s.d_evec[i], 1, &uu_c));
      CUBLAS_CHECK(cublasZdotc(cublas_handle, N, s.d_evec[i], 1, d_w, 1, &uw_c));
      const double uu  = cuCreal(uu_c);
      const double lam = cuCreal(uw_c) / uu;
      Taxpy_gen<CuC, double, N><<<NBlocks, NThreadsPerBlock>>>(d_w, -lam, s.d_evec[i], d_w);   // w -= lam u
      CUDA_CHECK(cudaDeviceSynchronize());
      CuC ww_c;
      CUBLAS_CHECK(cublasZdotc(cublas_handle, N, d_w, 1, d_w, 1, &ww_c));
      const double rres = std::sqrt(cuCreal(ww_c) / uu) / std::max(std::abs(lam), 1.0e-30);
      lam_res[i] = std::make_pair(lam, rres);
    }
    CUDA_CHECK(cudaFree(d_w));
  }
  std::sort(lam_res.begin(), lam_res.end());              // ascending in lambda (smallest = near-zero)

  const double RES_TOL = 1.0e-7;
  double max_res = 0.0;
  int n_ok = 0;
  for(int i=0;i<nconv;i++){
    max_res = std::max(max_res, lam_res[i].second);
    if(lam_res[i].second < RES_TOL) n_ok++;
  }
  std::cout << "# SELF-VALIDATION: " << n_ok << "/" << nconv << " modes with rel-residual < " << RES_TOL
            << "  (max rel-residual = " << max_res << ")" << std::endl;
  std::cout << (n_ok==nconv ? "# RESULT: PASS -- all returned pairs are eigenpairs of A"
                            : "# RESULT: FAIL -- some residuals too large (see .dat)") << std::endl;

  const std::string outname = "eig_lanczos_L"+std::to_string(Comp::N_REFINE)
                            + "_nt"+std::to_string(Comp::Nt)+"_claude.dat";
  std::ofstream out(outname);
  out << std::scientific << std::setprecision(14);
  out << "# IRL low spectrum of A=(D_ov+m)^dag(D_ov+m)   L="<<Comp::N_REFINE<<" Nt="<<Comp::Nt
      << " mass="<<mass<<" Nk="<<Nk<<" Nm="<<Nm<<" deg="<<cheb_degree<<std::endl;
  out << "# lambda_min/max(D_m singular) = "<<Dov.lambda_min<<" / "<<Dov.lambda_max<<std::endl;
  out << "# i   eval(A)   sqrt(eval)=singular_value   rel_residual" << std::endl;
  for(int i=0;i<nconv;i++){
    const double ev = lam_res[i].first;
    out << i << "   " << ev << "   " << std::sqrt(std::max(0.0, ev)) << "   " << lam_res[i].second << std::endl;
  }
  out.close();
  std::cout << "# wrote " << outname << " (" << nconv << " evals + residuals)" << std::endl;

  CUBLAS_CHECK(cublasDestroy(cublas_handle));
  CUDA_CHECK(cudaFree(d_src));
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
  return 0;
}
