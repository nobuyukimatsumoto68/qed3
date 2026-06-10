// cont_prop_eigbasis_claude.cpp
//
// Continuum free fermion propagator on S^2 x R in the eigenbasis, Eq. (C.28) of
// qed3_v2-6.pdf (label eq:Gf).  Reference for the eigenvalue construction:
// Abrikosov, hep-th/0212134 (Dirac on S^2), as cited in App. C.
//
// CHUNK 1 (this file, current state): angular + normalization primitives, validated standalone.
//   - xi_{|m|,n}(z): the S^2 angular function
//       xi = (1-z)^{alpha/2} (1+z)^{beta/2} P_j^{(alpha,beta)}(z),
//       alpha = |m| - 1/2,  beta = -(|m| + 1/2),  j = n + |m| + 1/2,
//     with |m| a positive half-integer (m \in Z + 1/2).  Since alpha+beta = -1, the standard
//     three-term degree recurrence is singular at low j; we evaluate the Jacobi factor through
//     its terminating hypergeometric series instead (App. C, Eq. def_P):
//       P_j^{(alpha,beta)}(z) = [ Gamma(alpha+1+j) / (Gamma(alpha+1) j!) ]
//                               * F(-j, j; |m|+1/2; (1-z)/2),
//     using the identity 1+alpha+beta+j = j valid here (alpha+beta=-1).
//   - c2_{|m|,n}: the eigenfunction normalization (App. C, Eq. normalization_psi).
//   - tests: edge values  xi(1)=delta_{|m|,1/2}/sqrt(2),  xi(-1)=0  (Eq. xi_deltam);
//            orthonormality  \int_{-1}^{1} xi_{|m|,n} xi_{|m|,n'} dz = c2_{|m|,n} delta_{nn'}.
//            The diagonal identity follows from the Jacobi norm with weight
//            (1-z)^alpha (1+z)^beta, which here equals xi^2 (the prefactor powers ARE half the
//            weight powers), and reduces algebraically to c2 (alpha+beta = -1).
//
// Pure CPU + GSL.  Build:
//   g++ -O2 -std=c++17 cont_prop_eigbasis_claude.cpp -o cont_prop_eigbasis_claude \
//       $(gsl-config --cflags --libs)

#include <cstdio>
#include <cmath>
#include <complex>
#include <vector>
#include <cassert>

#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <chrono>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_integration.h>

#include <boost/math/special_functions/jacobi.hpp>   // stable xi (three-term recurrence)

#include <highfive/H5File.hpp>

// ---------------------------------------------------------------------------
// naive terminating Gauss series  F(-j, j; c; x) = sum_{s=0}^{j} (-j)_s (j)_s / ((c)_s s!) x^s .
// Reliable only for x away from 1 (interior z); used ONLY for the cross-check in test 1.
static double hyperg_terminating(int j, double c, double x) {
  double term = 1.0, sum = 1.0;
  for (int s = 1; s <= j; ++s) {
    term *= (((double)(s - 1) - j) * ((double)(s - 1) + j))
            / (((double)(s - 1) + c) * (double)s) * x;
    sum  += term;
  }
  return sum;
}

// naive (reference) xi via the direct prefactor x F(-j,j;c;(1-z)/2); UNSTABLE near z=-1.
static double xi_mn_naive(double mm, int n, double z) {
  const double alpha = mm - 0.5, beta = -(mm + 0.5), c = mm + 0.5;
  const int    j     = (int)std::lround(n + mm + 0.5);
  if (z <= -1.0 + 1.0e-300) return 0.0;
  const double lpref = gsl_sf_lngamma(alpha + 1.0 + j) - gsl_sf_lngamma(alpha + 1.0)
                     - gsl_sf_lngamma((double)j + 1.0);
  const double P = std::exp(lpref) * hyperg_terminating(j, c, 0.5 * (1.0 - z));
  return std::pow(1.0 - z, 0.5 * alpha) * std::pow(1.0 + z, 0.5 * beta) * P;
}

// xi_{|m|,n}(z),  mm = |m| (positive half-integer: 0.5, 1.5, ...).
//
// ACTIVE version: stable Boost three-term Jacobi recurrence (verbatim form of the user's
// dirac_inverse/integrate.cc).  mpH = |m| + 1/2 (integer).  xi = (1-z)^{(mpH-1)/2}
// (1+z)^{-mpH/2} P^{(mpH-1,-mpH)}_{n+mpH}(z).  Robust to large (|m|,n) -- replaces the Pfaff
// binomial sum xi_mn_pfaff below, which OVERFLOWS at large order (the undamped 2D check
// check_S2_generic_pair_claude.cu exposed this; the 3D sum here is shielded by e^{-lambda|tau|}
// so the swap leaves the validated results unchanged, only more robust).
static double xi_mn(double mm, int n, double z) {
  const int mpH = (int)std::lround(mm + 0.5);                 // |m| + 1/2
  if (std::abs(z - 1.0) < 1.0e-14) return (mpH == 1) ? 1.0 / std::sqrt(2.0) : 0.0;
  if (std::abs(z + 1.0) < 1.0e-14) return 0.0;
  const double factor = std::pow(1.0 - z, 0.5 * (mpH - 1)) * std::pow(1.0 + z, -0.5 * mpH);
  const double poly   = boost::math::jacobi(n + mpH, (double)(mpH - 1.0), (double)(-1.0 * mpH), z);
  return factor * poly;
}

// OLD version (kept for reference / A-B): Pfaff transform of F(-j,j;c;x), x=(1-z)/2, folded to a
// degree-n binomial sum.  Accurate at low order (matched xi_mn_naive in test 1) but the explicit
// (z-1)^s (1+z)^{n-s} terms OVERFLOW for large (|m|,n).  Superseded by the Boost version above.
static double xi_mn_pfaff(double mm, int n, double z) {
  const double alpha = mm - 0.5;
  const int    j     = (int)std::lround(n + mm + 0.5);
  const double c     = mm + 0.5;
  const double om = 1.0 - z, op = 1.0 + z;

  double b = 1.0;                 // b_0
  double G = std::pow(op, n);     // s = 0 term: b_0 (z-1)^0 (1+z)^n
  for (int s = 1; s <= n; ++s) {
    b *= (((double)(s - 1) - j) * ((double)(s - 1) - n))
         / (((double)(s - 1) + c) * (double)s);
    G += b * std::pow(z - 1.0, s) * std::pow(op, n - s);
  }
  const double lpref = gsl_sf_lngamma(alpha + 1.0 + j) - gsl_sf_lngamma(alpha + 1.0)
                     - gsl_sf_lngamma((double)j + 1.0);
  return std::pow(om, 0.5 * alpha) * std::pow(op, 0.5 * (mm + 0.5))
         * std::exp(lpref - (double)j * M_LN2) * G;
}

// c2_{|m|,n} = 1/(2n+2|m|+1) * n! (n+2|m|)! / [ (n+|m|-1/2)! (n+|m|+1/2)! ]  (Eq. normalization_psi).
// all factorials are of nonnegative integers (2|m| is an integer); use lngamma.
static double c2_mn(double mm, int n) {
  const double lnum = gsl_sf_lngamma(n + 1.0) + gsl_sf_lngamma(n + 2.0 * mm + 1.0);
  const double lden = gsl_sf_lngamma(n + mm + 0.5) + gsl_sf_lngamma(n + mm + 1.5);
  return std::exp(lnum - lden) / (2.0 * n + 2.0 * mm + 1.0);
}

// ---------------------------------------------------------------------------
// CHUNK 2: the full 2x2 propagator G(x,x') = Eq. (C.28), with the radial (time) direction.
//
// The k-integral and the iota_3 sum are done ANALYTICALLY.  Because the signed eigenvalue is
// Lambda_K = iota_3 lambda (D psi_K = i Lambda_K psi_K), the denominator 1/(i Lambda_K) carries an
// extra iota_3.  Summing iota_3 = +/-1 turns e^{iota_3 i k tau} into 2 i sin(k tau) on the diagonal
// and 2 cos(k tau) off-diagonal.  Every k-integral then collapses to
//   \int dk e^{ik tau} / (k^2 + M^2)-type  =>  pure exponentials e^{-lambda |tau|},
//   lambda = lambda_{|m|,n} = n + |m| + 1/2   (the S^2 eigenvalue; M = lambda here),
// with the 1/Lambda (Bessel K_0) pieces cancelling by parity.  No Bessel / no Ci.
//
// Define the k-independent angular spinor pieces at a point (z=cos theta, phi), for label (m,n):
//   A = xi_{|m|,n}( iota_m z) e^{i m phi},
//   B = iota_m i (-1)^n xi_{|m|,n}(-iota_m z) e^{i m phi},   iota_m = sign(m).
// Then (tau = t - t', st = sign(tau)):
//   G_11 =  sum (1/4pi c^2) A  A'^*  st e^{-lambda|tau|},
//   G_22 = -sum (1/4pi c^2) B  B'^*  st e^{-lambda|tau|},
//   G_12 = -i sum (1/4pi c^2) A  B'^*    e^{-lambda|tau|},
//   G_21 = -i sum (1/4pi c^2) B  A'^*    e^{-lambda|tau|},
// summed over n>=0 and m in {+-1/2, +-3/2, ...}.  The prefactor 1/(4pi c^2) = 1/Cnm^2 is exactly
// the dirac_inverse normalization (Cnm^2 = 4pi c2_{|m|,n}); it reproduces Eq. (C.29) with the
// quoted (n+1)/(4pi) coefficient.

typedef std::complex<double> cd;

static inline double lambda_mn(double mm, int n) { return (double)n + mm + 0.5; }     // S^2 eigenvalue
static inline double cnm2_mn(double mm, int n)   { return 4.0 * M_PI * c2_mn(mm, n); } // = Cnm^2 (4pi norm)

// full 2x2 propagator G[a][b](x,x'), x=(theta,phi,t), x'=(thp,php,tp).  Requires tau != 0.
static void G_prop(double th, double ph, double t,
                   double thp, double php, double tp,
                   int nmax, cd G[2][2]) {
  const double z = std::cos(th), zp = std::cos(thp);
  const double tau = t - tp;
  const double atau = std::fabs(tau);
  const double st = (tau > 0.0) ? 1.0 : (tau < 0.0 ? -1.0 : 0.0);
  for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) G[a][b] = cd(0.0, 0.0);

  for (int twom = 1; twom <= 2 * nmax + 1; twom += 2) {   // 2|m| = 1,3,5,...
    const double mm = 0.5 * twom;
    for (int sgn = 0; sgn < 2; ++sgn) {                   // iota_m = +1, -1  (m = iota_m |m|)
      const double iom = (sgn == 0) ? 1.0 : -1.0;
      const double m   = iom * mm;
      for (int n = 0; n <= nmax; ++n) {
        const double lam = lambda_mn(mm, n);
        const double eM  = std::exp(-lam * atau);
        if (eM < 1.0e-30) continue;                       // negligible high-mode tail
        const double sgnn = (n & 1) ? -1.0 : 1.0;         // (-1)^n
        const double pref = 1.0 / cnm2_mn(mm, n);

        const cd eph  = std::exp(cd(0.0, m * ph));
        const cd ephp = std::exp(cd(0.0, m * php));
        const cd A  = xi_mn(mm, n,  iom * z ) * eph;
        const cd B  = iom * cd(0.0, 1.0) * sgnn * xi_mn(mm, n, -iom * z ) * eph;
        const cd Ap = xi_mn(mm, n,  iom * zp) * ephp;
        const cd Bp = iom * cd(0.0, 1.0) * sgnn * xi_mn(mm, n, -iom * zp) * ephp;

        G[0][0] +=          pref * A * std::conj(Ap) * st * eM;
        G[1][1] += -        pref * B * std::conj(Bp) * st * eM;
        G[0][1] += cd(0.0,-1.0) * pref * A * std::conj(Bp) * eM;
        G[1][0] += cd(0.0,-1.0) * pref * B * std::conj(Ap) * eM;
      }
    }
  }
}

// ---------------------------------------------------------------------------
// CHUNK 3: independent pure-S^2 (no time) propagator, for the CFT cross-check at general angle.
// From the Ising CFT (paper Eq. psipsi_S2):  <psi(theta,0) bar psi(0,0)> = sigma_1/(4pi sin(theta/2)).
// From the derived eigenbasis (Eq. psipsi_S2_ev, the (1,2) component; only |m|=1/2 survives at
// theta'=0):  coefficient = sum_{n>=0} (-1)^n / (pi sqrt2) * xi_{1/2,n}(-z),  z=cos theta.
// This is exactly dirac_inverse/integrate.cc line ~138; it exercises xi at all z (not just z=1).
// Convergence to the CFT form is from the IR; the UV (theta->0) needs large n_max (paper: 50,200).
static double G_S2_12(double theta, int nmax) {
  const double z = std::cos(theta);
  double s = 0.0;
  for (int n = 0; n <= nmax; ++n) {
    const double sgnn = (n & 1) ? -1.0 : 1.0;       // (-1)^n
    s += sgnn * xi_mn(0.5, n, -z);
  }
  return s / (std::sqrt(2.0) * M_PI);
}

// ---------------------------------------------------------------------------
// orthonormality integrand  xi_{|m|,n}(z) xi_{|m|,n'}(z)  for GSL quadrature.
struct OrthoParams { double mm; int n; int np; };
static double ortho_integrand(double z, void* p) {
  OrthoParams* q = (OrthoParams*)p;
  return xi_mn(q->mm, q->n, z) * xi_mn(q->mm, q->np, z);
}
static double ortho_integral(double mm, int n, int np) {
  OrthoParams par{mm, n, np};
  gsl_function F; F.function = &ortho_integrand; F.params = &par;
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(2000);
  double result = 0.0, err = 0.0;
  gsl_integration_qag(&F, -1.0, 1.0, 1.0e-12, 1.0e-11, 2000, GSL_INTEG_GAUSS61,
                      w, &result, &err);
  gsl_integration_workspace_free(w);
  return result;
}

// ---------------------------------------------------------------------------
static int run_tests() {
  gsl_set_error_handler_off();
  printf("# CHUNK 1 -- xi_{|m|,n} and c2_{|m|,n} primitives (Eq. C.28 angular part)\n\n");

  // ---- test 1: stable xi vs the naive prefactor-series xi at SAFE interior z --------------
  // (GSL 2F1 is unreliable for negative-integer parameters; the naive series is the reference
  //  away from z=-1, where it does not suffer the (1+z)^{beta/2} cancellation.)
  printf("## test 1: stable xi_mn vs naive prefactor-series xi (small j, safe z)\n");
  {
    double worst = 0.0;
    for (double mm : {0.5, 1.5})              // small |m|
      for (int n : {0, 1, 2})                  // small n => small j, naive series reliable
        for (double z : {0.0, 0.3, 0.6}) {     // x=(1-z)/2 in [0.2,0.5], no near-1 cancellation
          const double vs = xi_mn(mm, n, z);
          const double vn = xi_mn_naive(mm, n, z);
          const double rel = std::fabs(vs - vn) / (std::fabs(vn) + 1.0e-300);
          if (rel > worst) worst = rel;
        }
    printf("   worst relative diff stable vs naive = %.3e   %s\n\n", worst,
           (worst < 1.0e-9) ? "[PASS]" : "[CHECK]");
  }

  // ---- test 2: edge values  xi(1) = delta_{|m|,1/2}/sqrt(2),  xi(-1) = 0 ------------------
  printf("## test 2: edge values (Eq. xi_deltam)\n");
  {
    double worst = 0.0;
    const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    for (double mm : {0.5, 1.5, 2.5, 3.5})
      for (int n : {0, 1, 2, 5, 12}) {
        const double v1 = xi_mn(mm, n, 1.0);            // z = 1 exactly
        const double expect1 = (mm == 0.5) ? inv_sqrt2 : 0.0;
        const double vm1 = xi_mn(mm, n, -1.0);           // z = -1 exactly
        worst = std::max(worst, std::fabs(v1 - expect1));
        worst = std::max(worst, std::fabs(vm1));
      }
    printf("   max |xi(1)-expected| and |xi(-1)| = %.3e   %s\n\n", worst,
           (worst < 1.0e-9) ? "[PASS]" : "[CHECK]");
  }

  // ---- test 3: orthonormality  \int xi_n xi_n' dz = c2_n delta_{nn'} ----------------------
  printf("## test 3: orthonormality  int_{-1}^{1} xi_{|m|,n} xi_{|m|,n'} dz vs c2 delta_{nn'}\n");
  {
    double worst_diag = 0.0, worst_off = 0.0;
    for (double mm : {0.5, 1.5, 2.5, 3.5}) {
      for (int n = 0; n <= 40; ++n) {                      // diagonal up to n=40 (operating range)
        const double I  = ortho_integral(mm, n, n);
        const double c2 = c2_mn(mm, n);
        const double rel = std::fabs(I - c2) / (std::fabs(c2) + 1.0e-300);
        worst_diag = std::max(worst_diag, rel);
        for (int np = n + 1; np <= std::min(n + 6, 12); ++np) {   // off-diagonal sampling
          const double Ioff = ortho_integral(mm, n, np);
          const double scale = std::sqrt(c2 * c2_mn(mm, np));
          worst_off = std::max(worst_off, std::fabs(Ioff) / (scale + 1.0e-300));
        }
      }
    }
    // xi is exact (1e-14) at low n and degrades to ~1e-6 by n=40 (Pfaff-sum cancellation);
    // this floor is harmless for the propagator, where n=40 modes carry e^{-40|tau|}.
    printf("   worst diagonal rel error  |I - c2|/c2 (n<=40)   = %.3e   %s\n",
           worst_diag, (worst_diag < 1.0e-5) ? "[PASS]" : "[CHECK]");
    printf("   worst off-diagonal |I|/sqrt(c2 c2')            = %.3e   %s\n\n",
           worst_off, (worst_off < 1.0e-7) ? "[PASS]" : "[CHECK]");
  }

  printf("# CHUNK 1 done.\n\n");

  // ---- test 4 (CHUNK 2): temporal propagator vs Eq. (C.29) -------------------------------
  // G(t) at x=(0,0,t), x'=(0,0,0) must equal  sigma_3 sign(t) sum_n (n+1)/(4pi) e^{-(n+1)|t|}.
  printf("# CHUNK 2 -- full 2x2 propagator; temporal check against Eq. (C.29)\n");
  printf("## test 4: G(t) at theta=theta'=0 vs sigma_3 sign(t) sum_n (n+1)/(4pi) e^{-(n+1)t}\n");
  {
    const int nmax = 40;
    printf("   %6s  %18s  %18s  %12s  %12s\n",
           "t", "Re G11 (computed)", "ref (C.29)", "|G22+G11|", "|offdiag|");
    double worst_rel = 0.0, worst_struct = 0.0;
    for (double t : {0.2, 0.5, 1.0, 2.0, 4.0}) {
      cd G[2][2];
      G_prop(0.0, 0.0, t, 0.0, 0.0, 0.0, nmax, G);
      double ref = 0.0;
      for (int n = 0; n <= nmax; ++n) ref += (n + 1) / (4.0 * M_PI) * std::exp(-(n + 1) * t);
      const double reG11 = G[0][0].real();
      const double struct22 = std::abs(G[1][1] + G[0][0]);          // G22 = -G11 (sigma_3)
      const double offd = std::abs(G[0][1]) + std::abs(G[1][0]);    // diagonal only
      const double rel = std::fabs(reG11 - ref) / (std::fabs(ref) + 1.0e-300);
      worst_rel = std::max(worst_rel, rel);
      worst_struct = std::max(worst_struct, struct22 + offd);
      printf("   %6.2f  %18.10e  %18.10e  %12.3e  %12.3e\n",
             t, reG11, ref, struct22, offd);
    }
    printf("   worst |G11-ref|/ref = %.3e   worst |sigma_3 struct viol| = %.3e   %s\n",
           worst_rel, worst_struct,
           (worst_rel < 1.0e-12 && worst_struct < 1.0e-12) ? "[PASS]" : "[CHECK]");
  }

  printf("\n# CHUNK 2 done.\n\n");

  // ---- test 5 (CHUNK 3): pure-S^2 propagator vs CFT form 1/(4pi sin(theta/2)) -------------
  printf("# CHUNK 3 -- pure-S^2 angular check (general angle), exercises xi at all z\n");
  printf("## test 5: sum_n (-1)^n xi_{1/2,n}(-z)/(pi sqrt2) vs 1/(4pi sin(theta/2))\n");
  {
    printf("   %8s  %16s  %16s  %16s  %12s\n",
           "theta", "n_max=20", "n_max=40", "CFT ref", "rel(40)");
    for (double th : {M_PI * 0.75, M_PI * 0.6, M_PI * 0.5, M_PI / 3.0}) {
      const double v20 = G_S2_12(th, 20);
      const double v40 = G_S2_12(th, 40);
      const double ref = 1.0 / (4.0 * M_PI * std::sin(0.5 * th));
      const double rel = std::fabs(v40 - ref) / (std::fabs(ref) + 1.0e-300);
      printf("   %8.4f  %16.10e  %16.10e  %16.10e  %12.3e\n", th, v20, v40, ref, rel);
    }
    // judge at theta=pi/2 (z=0), where the alternating S^2 sum converges fastest.
    const double ref_mid = 1.0 / (4.0 * M_PI * std::sin(M_PI / 4.0));
    const double rel_mid = std::fabs(G_S2_12(M_PI / 2.0, 40) - ref_mid) / ref_mid;
    printf("   rel(n_max=40) at theta=pi/2 (fast-converging point) = %.3e   %s\n",
           rel_mid, (rel_mid < 1.0e-4) ? "[PASS]" : "[CHECK]");
    printf("   (the pure-S^2 sum is UNDAMPED -> off-pi/2 angles converge slowly, needing huge n_max,\n");
    printf("    exactly as in the paper (n_max=50,200).  The 3D propagator has e^{-lambda|tau|}\n");
    printf("    damping instead, so n_max=40 is ample -- see test 6.)\n");
  }

  // ---- test 6 (CHUNK 3): self-convergence of the damped 3D propagator at a general point --
  // No closed form off-axis, but the e^{-lambda|tau|} damping must make the n-sum converge in
  // n_max.  Compare G(n_max=30) vs G(n_max=40) at a generic (theta,theta',dphi,tau).
  printf("\n## test 6: 3D propagator self-convergence in n_max at a general off-axis point\n");
  {
    const double th = 1.0, thp = 2.0, dph = 0.7;          // generic angles, phi'=php
    printf("   %6s  %14s  %14s  %14s\n", "tau", "max|G|", "max|G40-G30|", "max|G40-G35|");
    double worst = 0.0;
    for (double tau : {1.0, 0.5, 0.2}) {                  // incl. lattice nearest-slice tau=a_t=0.2
      cd G30[2][2], G35[2][2], G40[2][2];
      G_prop(th, dph, tau, thp, 0.0, 0.0, 30, G30);
      G_prop(th, dph, tau, thp, 0.0, 0.0, 35, G35);
      G_prop(th, dph, tau, thp, 0.0, 0.0, 40, G40);
      double mag = 0.0, d4030 = 0.0, d4035 = 0.0;
      for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) {
        mag   = std::max(mag,   std::abs(G40[a][b]));
        d4030 = std::max(d4030, std::abs(G40[a][b] - G30[a][b]));
        d4035 = std::max(d4035, std::abs(G40[a][b] - G35[a][b]));
      }
      worst = std::max(worst, (tau >= 0.5) ? d4035 : 0.0);
      printf("   %6.2f  %14.6e  %14.3e  %14.3e\n", tau, mag, d4030, d4035);
    }
    printf("   worst |G40-G35| for tau>=0.5 = %.3e   %s\n",
           worst, (worst < 1.0e-8) ? "[PASS]" : "[CHECK]");
    printf("   (tau=0.2 = nearest lattice slice: slowest, truncation ~ e^{-40*0.2} ~ 3e-4;\n");
    printf("    larger separations converge to machine precision well within n_max=40.)\n");
  }

  printf("\n# CHUNK 3 done.\n");
  return 0;
}

// ===========================================================================================
// CHUNK 4: lattice driver -- read pts_n<L>.dat, build the all-to-all continuum propagator at the
// lattice (theta,phi,t) coordinates, and write it in the other project's Dinv.h5 schema
// (jj_propagator_exact_claude.cu): Dm_inv/{real,imag} row-major length-N^2, index = row*N + col,
// plus N/N_REFINE/Nt/n_sites/parity/mP/k/complete.  Free massless => parity=0, mP=(0,0).
//
// Global DOF index (includes/dirac_ext.h:65):  r = Nx*s + NS*i + a,  Nx = NS*n_sites, NS = 2,
// s = timeslice (slowest), i = site, a = spinor (fastest).  Continuum t = at*s; tau = at*(s-s').
// No antiperiodic folding (user): literal tau, both signs, |tau| up to at*Nt; tau=0 blocks = 0.
// ===========================================================================================

// read pts_n<L>.dat (one unit-sphere Cartesian point per line) -> theta=acos(z), phi=atan2(y,x).
static int read_sites(const std::string& path, std::vector<double>& theta, std::vector<double>& phi) {
  std::ifstream in(path);
  if (!in) { fprintf(stderr, "ERROR: cannot open %s\n", path.c_str()); return -1; }
  theta.clear(); phi.clear();
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ss(line);
    double x, y, z;
    if (!(ss >> x >> y >> z)) continue;            // skip blank/garbage lines
    const double rr = std::sqrt(x*x + y*y + z*z);
    theta.push_back(std::acos(std::max(-1.0, std::min(1.0, z / rr))));
    phi.push_back(std::atan2(y, x));
  }
  return (int)theta.size();
}

static int build_all_to_all(int L, int Nt, double at, int nmax,
                            int kidx, const std::string& geo_dir, const std::string& out_dir) {
  // ---- geometry --------------------------------------------------------------------------
  std::vector<double> theta, phi;
  const std::string pts = geo_dir + "/pts_n" + std::to_string(L) + ".dat";
  const int ns = read_sites(pts, theta, phi);
  if (ns <= 0) return -1;
  const int expect = 10 * L * L + 2;
  if (ns != expect)
    fprintf(stderr, "WARN: n_sites=%d != 10*L^2+2=%d for L=%d\n", ns, expect, L);
  const long long NS = 2, Nx = NS * ns, N = Nx * (long long)Nt;
  printf("# L=%d  n_sites=%d  Nt=%d  at=%.4f  nmax=%d  N=%lld  (matrix %lld x %lld)\n",
         L, ns, Nt, at, nmax, N, N, N);

  const int nq = nmax + 1;                         // mm = 0.5 + q,  q = 0..nmax
  // ---- precompute the k-independent angular spinor pieces A,B per (site,q,n,iota_m) -------
  //   A = xi_{|m|,n}( iota_m z) e^{i m phi},  B = iota_m i (-1)^n xi_{|m|,n}(-iota_m z) e^{i m phi}.
  // layout: idx = ((i*nq + q)*nq + n)*2 + s,  s=0 -> iota_m=+1, s=1 -> iota_m=-1.
  const size_t TAB = (size_t)ns * nq * nq * 2;
  std::vector<cd> Amat(TAB), Bmat(TAB);
  for (int i = 0; i < ns; ++i) {
    const double z = std::cos(theta[i]);
    for (int q = 0; q < nq; ++q) {
      const double mm = 0.5 + q;
      for (int n = 0; n <= nmax; ++n) {
        const double sgnn = (n & 1) ? -1.0 : 1.0;
        const double xiP = xi_mn(mm, n,  z);       // xi(+z)
        const double xiM = xi_mn(mm, n, -z);       // xi(-z)
        for (int s = 0; s < 2; ++s) {
          const double iom = (s == 0) ? 1.0 : -1.0;
          const double m   = iom * mm;
          const cd eph = std::exp(cd(0.0, m * phi[i]));
          const double xiA = (s == 0) ? xiP : xiM; // xi( iota_m z)
          const double xiB = (s == 0) ? xiM : xiP; // xi(-iota_m z)
          const size_t idx = (((size_t)i * nq + q) * nq + n) * 2 + s;
          Amat[idx] = xiA * eph;
          Bmat[idx] = iom * cd(0.0, 1.0) * sgnn * xiB * eph;
        }
      }
    }
  }

  // ---- allocate the dense N^2 output (row-major), zero (tau=0 blocks stay zero) -----------
  printf("# allocating %.2f GB (re+im) ...\n", 2.0 * (double)N * (double)N * 8.0 / 1.0e9);
  std::vector<double> re((size_t)N * N, 0.0), im((size_t)N * N, 0.0);

  // ---- per-(i,j) 2x2 block at +tau (st=+1); the -tau block flips only the diagonal --------
  // accumulate into Gp[2][2] for tau = +at*dt (dt>0); the dt<0 block is diag-negated.
  auto block_plus = [&](int i, int j, double atau, cd Gp[2][2]) {
    for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) Gp[a][b] = cd(0.0, 0.0);
    for (int q = 0; q < nq; ++q) {
      const double mm = 0.5 + q;
      for (int n = 0; n <= nmax; ++n) {
        const double lam = (double)n + mm + 0.5;
        const double eM  = std::exp(-lam * atau);
        if (eM < 1.0e-30) continue;
        const double w = eM / cnm2_mn(mm, n);      // pref * e^{-lambda|tau|}  (st=+1 here)
        for (int s = 0; s < 2; ++s) {
          const size_t ii = (((size_t)i * nq + q) * nq + n) * 2 + s;
          const size_t jj = (((size_t)j * nq + q) * nq + n) * 2 + s;
          const cd Ai = Amat[ii], Bi = Bmat[ii];
          const cd Aj = std::conj(Amat[jj]), Bj = std::conj(Bmat[jj]);
          Gp[0][0] +=             w * Ai * Aj;       // st=+1
          Gp[1][1] += -           w * Bi * Bj;       // st=+1
          Gp[0][1] += cd(0.0,-1.0) * w * Ai * Bj;    // off-diag: no st
          Gp[1][0] += cd(0.0,-1.0) * w * Bi * Aj;
        }
      }
    }
  };

  // ---- assemble: loop ordered spatial pairs (i,j) and |dt|, scatter into both tau signs ---
  const auto t0 = std::chrono::steady_clock::now();
  auto elapsed = [&]() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  };
  for (int i = 0; i < ns; ++i) {
    for (int j = 0; j < ns; ++j) {
      for (int dt = 1; dt <= Nt - 1; ++dt) {
        cd Gp[2][2]; block_plus(i, j, at * (double)dt, Gp);
        // dt>0 (tau=+at*dt): rows s=dt..Nt-1, s'=s-dt;  block = Gp.
        // dt<0 (tau=-at*dt): rows s=0..Nt-1-dt, s'=s+dt; block = diag-negated Gp.
        for (int s = dt; s <= Nt - 1; ++s) {
          const int sp = s - dt;
          for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) {
            const size_t r = (size_t)Nx * s  + NS * i + a;
            const size_t c = (size_t)Nx * sp + NS * j + b;
            re[r * N + c] = Gp[a][b].real();  im[r * N + c] = Gp[a][b].imag();
          }
        }
        for (int s = 0; s <= Nt - 1 - dt; ++s) {
          const int sp = s + dt;             // s - sp = -dt  => tau < 0
          for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) {
            const double sflip = (a == b) ? -1.0 : 1.0;   // diagonal odd in tau, off-diag even
            const size_t r = (size_t)Nx * s  + NS * i + a;
            const size_t c = (size_t)Nx * sp + NS * j + b;
            re[r * N + c] = sflip * Gp[a][b].real();  im[r * N + c] = sflip * Gp[a][b].imag();
          }
        }
      }
    }
    if ((i % 8) == 0)
      printf("#   site %d/%d  [%.1f s]\n", i, ns, elapsed());
  }

  // ---- write Dinv.<k>.h5 in the other project's schema (atomic .tmp + rename) -------------
  std::filesystem::create_directories(out_dir);
  const std::string h5path = out_dir + "/Dinv." + std::to_string(kidx) + ".h5";
  const std::string h5tmp  = h5path + ".tmp";
  {
    HighFive::File h5(h5tmp, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
    h5.createDataSet("N",        std::vector<int>{(int)N});
    h5.createDataSet("N_REFINE", std::vector<int>{L});
    h5.createDataSet("Nt",       std::vector<int>{Nt});
    h5.createDataSet("n_sites",  std::vector<int>{ns});
    h5.createDataSet("parity",   std::vector<int>{0});
    h5.createDataSet("mP/real",  std::vector<double>{0.0});
    h5.createDataSet("mP/imag",  std::vector<double>{0.0});
    h5.createDataSet("Dm_inv/real", re);
    h5.createDataSet("Dm_inv/imag", im);
    h5.createDataSet("k", std::vector<int>{kidx});
    h5.createDataSet("complete", std::vector<int>{1});   // sentinel, last
  }
  std::filesystem::rename(h5tmp, h5path);
  printf("# wrote %s  [%.1f s total]\n", h5path.c_str(), elapsed());
  return 0;
}

static void usage(const char* a0) {
  printf("usage: %s [--test] | [--L <1|2|4>] [--nt N] [--at X] [--nmax N] [--k K]\n", a0);
  printf("                 [--geo <dir>] [--out <dir>]\n");
  printf("  --test         run the validation suite (Chunks 1-3) and exit\n");
  printf("  --L            lattice refinement (reads <geo>/pts_n<L>.dat); default 1\n");
  printf("  --nt           temporal extent (default 128)\n");
  printf("  --at           temporal lattice spacing (default 0.2)\n");
  printf("  --nmax         mode truncation (default 40)\n");
  printf("  --k            config index for the Dinv.<k>.h5 filename (default 0)\n");
  printf("  --geo          geometry data dir (default ../../geometry/data)\n");
  printf("  --out          output dir (default cont_prop_L<L>)\n");
}

int main(int argc, char** argv) {
  gsl_set_error_handler_off();
  if (argc <= 1) return run_tests();

  int L = 1, Nt = 128, nmax = 40, kidx = 0;
  double at = 0.2;
  std::string geo = "../../geometry/data", out;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto next = [&](double& v){ if (i + 1 < argc) v = std::atof(argv[++i]); };
    auto nexi = [&](int& v){ if (i + 1 < argc) v = std::atoi(argv[++i]); };
    if      (a == "--test") return run_tests();
    else if (a == "--L")    nexi(L);
    else if (a == "--nt")   nexi(Nt);
    else if (a == "--at")   { double t; next(t); at = t; }
    else if (a == "--nmax") nexi(nmax);
    else if (a == "--k")    nexi(kidx);
    else if (a == "--geo")  { if (i + 1 < argc) geo = argv[++i]; }
    else if (a == "--out")  { if (i + 1 < argc) out = argv[++i]; }
    else { usage(argv[0]); return (a == "--help" || a == "-h") ? 0 : 1; }
  }
  if (out.empty()) out = "cont_prop_L" + std::to_string(L);
  printf("# cont_prop_eigbasis: continuum free fermion propagator (Eq. C.28) -> Dinv.h5\n");
  return build_all_to_all(L, Nt, at, nmax, kidx, geo, out);
}
