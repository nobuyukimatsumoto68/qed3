// check_C54_s2_prop_claude.cu
//
// Standalone analytic check of Appendix C.3 (App. C of qed3_v2-6.pdf), the pure-S^2
// ("N_t = 1", single time-slice) free fermion propagator.  We verify the derived
// eigenfunction sum, Eq. (C.55), against the Ising-CFT closed form, Eq. (C.54):
//
//   (C.54)   <psi(theta,0) bar psi(0,0)> = sigma_1 / (4 pi sin(theta/2)),
//   (C.55)   <psi(theta,0) bar psi(0,0)> = sigma_1 * sum_{n>=0} (-1)^n/(pi sqrt2) xi_{1/2,n}(-z),
//
// with z = cos theta and the source at the north pole (theta' = 0), where only the
// |m| = 1/2 sector survives.  Both sides are proportional to sigma_1, so we compare
// the scalar (1,2) coefficient (the paper takes the (1,2) component, see Fig. 13).
//
// Reference for xi (Jacobi/hypergeometric angular function): G. 't Hooft style spinor
// harmonics on S^2; here xi follows Abrikosov, hep-th/0212134, as cited in App. C.  The
// xi_mn / c2_mn routines below are copied verbatim (same Pfaff-transform form) from
// cont_prop_eigbasis_claude.cpp, where they are validated (orthonormality to 1e-14,
// edge values xi(1) = delta_{|m|,1/2}/sqrt2, xi(-1) = 0).
//
// The truncated sum (C.55) converges to (C.54) from the infrared; the ultraviolet
// (theta -> 0, where the CFT form diverges) is restored only as n_max -> infinity
// (paper uses n_max = 50, 200).  This file reproduces that Fig. 13 behavior.
//
// Pure host math (GSL for lngamma); nvcc compiles the host code.  Single-thread.
// Build/run: see run_check_C54_claude.sh in this directory.

#include <cstdio>
#include <cmath>
#include <vector>

#include <gsl/gsl_sf_gamma.h>
#include <boost/math/special_functions/jacobi.hpp>   // stable xi (three-term recurrence)

// xi_{|m|,n}(z),  mm = |m| (positive half-integer: 0.5, 1.5, ...).
// ACTIVE version: stable Boost three-term Jacobi recurrence (verbatim form of the user's
// dirac_inverse/integrate.cc).  mpH = |m| + 1/2 (integer).  Robust to large (|m|,n); replaces
// the Pfaff binomial sum (xi_mn_pfaff below) which overflows at large order.
static double xi_mn(double mm, int n, double z) {
  const int mpH = (int)std::lround(mm + 0.5);                 // |m| + 1/2
  if (std::abs(z - 1.0) < 1.0e-14) return (mpH == 1) ? 1.0 / std::sqrt(2.0) : 0.0;
  if (std::abs(z + 1.0) < 1.0e-14) return 0.0;
  const double factor = std::pow(1.0 - z, 0.5 * (mpH - 1)) * std::pow(1.0 + z, -0.5 * mpH);
  const double poly   = boost::math::jacobi(n + mpH, (double)(mpH - 1.0), (double)(-1.0 * mpH), z);
  return factor * poly;
}

// OLD version (kept for reference): Pfaff transform folded to a degree-n binomial sum.  Accurate
// at low order but the explicit (z-1)^s (1+z)^{n-s} terms OVERFLOW for large (|m|,n).
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

// Derived eigenfunction sum, Eq. (C.55): the scalar (1,2) coefficient of <psi(theta,0) bar psi(0,0)>,
//   sum_{n=0}^{nmax} (-1)^n / (pi sqrt2) xi_{1/2,n}(-z),   z = cos theta.
static double G_S2_12_sum(double theta, int nmax) {
  const double z = std::cos(theta);
  double s = 0.0;
  for (int n = 0; n <= nmax; ++n) {
    const double sgnn = (n & 1) ? -1.0 : 1.0;        // (-1)^n
    s += sgnn * xi_mn(0.5, n, -z);
  }
  return s / (std::sqrt(2.0) * M_PI);
}

// Ising-CFT closed form, Eq. (C.54): the (1,2) coefficient = 1 / (4 pi sin(theta/2)).
static double G_S2_12_cft(double theta) {
  return 1.0 / (4.0 * M_PI * std::sin(0.5 * theta));
}

int main(int argc, char** argv) {
  // theta-sweep, log-spaced from the deep UV to theta = pi (antipode), like Fig. 13.
  const double th_lo = 0.02;
  const double th_hi = M_PI;
  const int    npts  = 41;

  const int nmax_list[] = {50, 200, 1000};
  const int n_nmax = (int)(sizeof(nmax_list) / sizeof(nmax_list[0]));

  printf("# Appendix C.3 check: S^2 free fermion propagator (C.54) vs derived sum (C.55)\n");
  printf("# (1,2) coefficient; CFT = 1/(4 pi sin(theta/2)).  rel = |sum - cft| / cft\n");
  printf("#\n");
  printf("# %10s %14s", "theta", "cft");
  for (int k = 0; k < n_nmax; ++k) printf("   sum(nmax=%-5d)   rel", nmax_list[k]);
  printf("\n");

  std::vector<double> max_rel(n_nmax, 0.0);

  for (int i = 0; i < npts; ++i) {
    const double frac  = (double)i / (double)(npts - 1);
    const double theta = th_lo * std::pow(th_hi / th_lo, frac);   // geometric spacing
    const double cft   = G_S2_12_cft(theta);

    printf("  %10.6f %14.6e", theta, cft);
    for (int k = 0; k < n_nmax; ++k) {
      const double s   = G_S2_12_sum(theta, nmax_list[k]);
      const double rel = std::fabs(s - cft) / std::fabs(cft);
      printf("  %14.6e %9.2e", s, rel);
      if (rel > max_rel[k]) max_rel[k] = rel;
    }
    printf("\n");
  }

  printf("#\n");
  printf("# max rel error over the full sweep (incl. deep UV theta = %.3f):\n", th_lo);
  for (int k = 0; k < n_nmax; ++k)
    printf("#   nmax = %-5d : %.3e\n", nmax_list[k], max_rel[k]);

  // IR-window verdict: away from the UV divergence (theta >= pi/4) the sum must match the CFT.
  // This is the clean pass/fail; the UV (small theta) only restores as nmax -> infinity.
  const int    nmax_v = 1000;
  double max_rel_ir = 0.0;
  for (int i = 0; i < npts; ++i) {
    const double frac  = (double)i / (double)(npts - 1);
    const double theta = th_lo * std::pow(th_hi / th_lo, frac);
    if (theta < 0.25 * M_PI) continue;
    const double cft = G_S2_12_cft(theta);
    const double s   = G_S2_12_sum(theta, nmax_v);
    const double rel = std::fabs(s - cft) / std::fabs(cft);
    if (rel > max_rel_ir) max_rel_ir = rel;
  }
  printf("#\n");
  printf("# IR window (theta >= pi/4), nmax = %d : max rel = %.3e  -> %s\n",
         nmax_v, max_rel_ir, (max_rel_ir < 1.0e-3) ? "PASS" : "FAIL");

  return 0;
}
