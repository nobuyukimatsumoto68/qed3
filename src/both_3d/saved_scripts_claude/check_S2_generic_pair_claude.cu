// check_S2_generic_pair_claude.cu
//
// Off-diagonal / off-pole check of the continuum free fermion propagator, complementing
// the source-at-pole check (C.54)/(C.55) in check_C54_s2_prop_claude.cu.
//
// MOTIVATION.  (C.29) (temporal-at-pole) and (C.54) (spatial, source-at-pole) both restrict
// to |m| = 1/2.  Same-site homogeneity validates only the DIAGONAL addition theorem.  The
// OFF-DIAGONAL reassembly between two DISTINCT off-pole points -- the cross terms
// A(theta) conj(B(theta')) carrying the iota_m, (-1)^n, -i, (+/- iota_m z) sign machinery --
// is what feeds the JJ correlators and is otherwise untested (it cancels at same-site, which
// is WHY the same-site block is pure sigma_3).
//
// OBJECT.  The pure-S^2 (2D, "N_t = 1") free fermion propagator from the eigenbasis.  The S^2
// Dirac operator anticommutes with sigma_3, so the +/- lambda eigenspinors are
// Psi_- = sigma_3 Psi_+ with Psi_+ = (A, B)^T,
//   A = xi_{|m|,n}( iota_m z ) e^{i m phi},   B = iota_m i (-1)^n xi_{|m|,n}( -iota_m z ) e^{i m phi}.
// G = sum_K Psi_K Psi_K^dagger / (i Lambda_K Cnm^2) collapses to a PURELY OFF-DIAGONAL matrix,
//   G_12 = sum_{m,n} ( 2/(i lambda Cnm^2) ) A conj(B'),   G_21 = sum ( 2/(i lambda Cnm^2) ) B conj(A'),
// Cnm^2 = 4 pi c^2_{|m|,n}.  Source-at-pole reduces this to (C.55), fixing the (C.54) normalization.
//
// COMPARISON.  The 2D Ising/free-fermion CFT two-point on the round S^2 is, by conformal+S^2
// symmetry, a gamma-structure of magnitude 1/(4 pi sin(gamma/2)), gamma = geodesic angle
// ((C.54) carried off the pole by isometry).  Independent SU(2) frame on each point
// (G -> U_x G U_y^dagger) => only SINGULAR VALUES are comparable; here G is off-diagonal so
// they are {|G_12|, |G_21|}.  NOTE: |G_12| = |G_21| is an EXACT identity here
// (G_21 = -conj(G_12), from the +/-m pairing), so it is NOT a discriminating test -- only the
// MAGNITUDE |G_12| -> 1/(4 pi sin(gamma/2)) discriminates the off-diagonal machinery.
//
// NUMERICS.  This is an UNDAMPED 2D sum (no e^{-lambda|tau|}), so it stresses xi at large order.
// The project Pfaff xi (xi_mn below, from cont_prop_eigbasis_claude.cpp) OVERFLOWS for large
// (|m|,n).  We therefore evaluate the propagator with the user's Boost three-term-recurrence
// xi (xi_boost, copied verbatim from dirac_inverse/integrate.cc), which is stable, and keep
// xi_mn only for the divergence diagnostic (xi_pfaff_vs_boost) that locates the j where Pfaff
// breaks and shows it sits under the e^{-lambda|tau|} damping that protects the 3D propagator.
//
// Pure host math (GSL lngamma + Boost Jacobi); nvcc not needed -- build with g++ -x c++.
// Single-thread.

#include <cstdio>
#include <cmath>
#include <complex>
#include <algorithm>

#include <gsl/gsl_sf_gamma.h>
#include <boost/math/special_functions/jacobi.hpp>

typedef std::complex<double> cd;
typedef double Real;

// ---- xi via Boost Jacobi recurrence (stable) -- verbatim from dirac_inverse/integrate.cc ----
// mpH = |m| + 1/2 = 1, 2, 3, ...  (integer);  n = 0, 1, 2, ...
static Real xi_boost(const int mpH, const int n, const double z) {
  if (std::abs(z - 1.0) < 1.0e-14) return mpH == 1 ? 1.0 / std::sqrt(2.0) : 0.0;
  if (std::abs(z + 1.0) < 1.0e-14) return 0.0;
  const Real factor = std::pow(1.0 - z, 0.5 * (mpH - 1)) * std::pow(1.0 + z, -0.5 * mpH);
  const Real poly   = boost::math::jacobi(n + mpH, (Real)(mpH - 1.0), (Real)(-1.0 * mpH), z);
  return factor * poly;
}

// ---- xi via Pfaff binomial sum (project version; OVERFLOWS at large order) ------------------
// mm = |m| (half-integer).  Kept ONLY for the divergence diagnostic.
static double xi_mn(double mm, int n, double z) {
  const double alpha = mm - 0.5;
  const int    j     = (int)std::lround(n + mm + 0.5);
  const double c     = mm + 0.5;
  const double om = 1.0 - z, op = 1.0 + z;

  double b = 1.0;
  double G = std::pow(op, n);
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

// c2_{|m|,n}  (Eq. normalization_psi); identical to cont_prop_eigbasis_claude.cpp.
static double c2_mn(double mm, int n) {
  const double lnum = gsl_sf_lngamma(n + 1.0) + gsl_sf_lngamma(n + 2.0 * mm + 1.0);
  const double lden = gsl_sf_lngamma(n + mm + 0.5) + gsl_sf_lngamma(n + mm + 1.5);
  return std::exp(lnum - lden) / (2.0 * n + 2.0 * mm + 1.0);
}
static inline double cnm2_mn(double mm, int n) { return 4.0 * M_PI * c2_mn(mm, n); }

// Pure-S^2 propagator off-diagonal entries, summed shell by shell in lambda = n + |m| + 1/2.
// Uses xi_boost (stable).  Each shell: |m| = 1/2 .. lam-1/2, n = lam - |m| - 1/2, both iota_m.
static void G_S2_offdiag(double theta, double phi, double thp, double php,
                         int lammax, cd* G12, cd* G21) {
  const double z  = std::cos(theta), zp = std::cos(thp);
  cd s12(0.0, 0.0), s21(0.0, 0.0);

  for (int lam = 1; lam <= lammax; ++lam) {
    for (int twom = 1; twom <= 2 * lam - 1; twom += 2) {
      const double mm  = 0.5 * twom;
      const int    mpH = (twom + 1) / 2;               // = |m| + 1/2
      const int    n   = lam - mpH;                     // >= 0
      const double sgnn = (n & 1) ? -1.0 : 1.0;        // (-1)^n
      const double w  = 2.0 / (cnm2_mn(mm, n) * (double)lam);
      const cd     iw = cd(0.0, -1.0) * w;             // 2/(i lam Cnm^2)

      for (int sgn = 0; sgn < 2; ++sgn) {              // iota_m = +1, -1
        const double iom = (sgn == 0) ? 1.0 : -1.0;
        const double m   = iom * mm;

        const cd eph  = std::exp(cd(0.0, m * phi));
        const cd ephp = std::exp(cd(0.0, m * php));
        // const cd A = xi_mn(mm,n, iom*z)*eph;  // (Pfaff -- overflows; replaced by xi_boost)
        const cd A   = xi_boost(mpH, n,  iom * z ) * eph;
        const cd B   = iom * cd(0.0, 1.0) * sgnn * xi_boost(mpH, n, -iom * z ) * eph;
        const cd Ap  = xi_boost(mpH, n,  iom * zp) * ephp;
        const cd Bp  = iom * cd(0.0, 1.0) * sgnn * xi_boost(mpH, n, -iom * zp) * ephp;

        s12 += iw * A * std::conj(Bp);
        s21 += iw * B * std::conj(Ap);
      }
    }
  }
  *G12 = s12;
  *G21 = s21;
}

static double geodesic_angle(double theta, double phi, double thp, double php) {
  const double c = std::cos(theta) * std::cos(thp)
                 + std::sin(theta) * std::sin(thp) * std::cos(phi - php);
  return std::acos(std::max(-1.0, std::min(1.0, c)));
}

// ---- diagnostic: where does Pfaff xi diverge from Boost, vs the 3D damping weight? ----------
// For each shell lambda, scan all splits (|m|,n) with |m|+n+1/2 = lambda at an interior z,
// take the worst relative Pfaff-vs-Boost error, and print it beside e^{-lambda a_t} (a_t=0.2,
// the smallest 3D time separation).  This shows the 3D propagator only ever uses these shells
// multiplied by that (tiny) weight, so the Pfaff degradation is harmless there.
static void xi_pfaff_vs_boost(double z, double at) {
  printf("# Pfaff xi_mn vs Boost xi at z = %.3f : worst rel error over each lambda-shell,\n", z);
  printf("# beside the 3D damping weight e^{-lambda a_t} (a_t = %.2f, nearest time slice).\n", at);
  printf("#  %6s  %14s  %14s\n", "lambda", "worst_rel_err", "e^{-lam a_t}");
  for (int lam = 2; lam <= 130; ++lam) {
    double worst = 0.0;
    for (int twom = 1; twom <= 2 * lam - 1; twom += 2) {
      const double mm  = 0.5 * twom;
      const int    mpH = (twom + 1) / 2;
      const int    n   = lam - mpH;
      const double vb = xi_boost(mpH, n, z);
      const double vp = xi_mn(mm, n, z);
      const double rel = std::fabs(vp - vb) / (std::fabs(vb) + 1.0e-300);
      if (rel > worst) worst = rel;
    }
    if (lam % 8 == 0 || worst > 1.0e-6)
      printf("   %6d  %14.3e  %14.3e\n", lam, worst, std::exp(-lam * at));
  }
  printf("#\n");
}

int main(int argc, char** argv) {
  // --- diagnostic first: locate the Pfaff overflow vs the 3D damping ---
  printf("# ===== diagnostic: Pfaff xi vs Boost xi, and the 3D damping that hides it =====\n");
  xi_pfaff_vs_boost(0.30, 0.2);

  // --- main check: generic off-pole pairs, propagator via stable Boost xi ---
  const double thx = 1.0, phx = 0.3;
  struct Pair { double thp, php; };
  const Pair ys[] = {
    {0.7,  1.9}, {1.3,  2.6}, {2.0,  0.9}, {2.4,  3.5},
    {1.7, -1.1}, {0.9,  4.2}, {2.2,  1.6}, {1.1,  2.2}
  };
  const int npair = (int)(sizeof(ys) / sizeof(ys[0]));

  const int lam_list[] = {50, 100, 200, 400};
  const int n_lam = (int)(sizeof(lam_list) / sizeof(lam_list[0]));

  printf("# ===== generic off-pole S^2 two-point check (Boost xi) =====\n");
  printf("# G off-diagonal; |G12|=|G21| is an exact identity.  CFT: |G| = 1/(4 pi sin(gamma/2)).\n");
  printf("# discriminating test: |G12|/CFT -> 1 as lambda_max grows.\n#\n");

  double worst_magerr = 0.0, worst_asym = 0.0;

  for (int p = 0; p < npair; ++p) {
    const double thp = ys[p].thp, php = ys[p].php;
    const double gamma = geodesic_angle(thx, phx, thp, php);
    const double ref   = 1.0 / (4.0 * M_PI * std::sin(0.5 * gamma));

    printf("## pair %d:  x=(%.2f,%.2f)  y=(%.2f,%.2f)   gamma=%.4f  CFT|G|=%.6e\n",
           p, thx, phx, thp, php, gamma, ref);
    printf("#   %8s  %14s %14s  %10s  %12s\n",
           "lam_max", "|G12|", "|G21|", "asym", "|G12|/CFT");

    for (int k = 0; k < n_lam; ++k) {
      cd G12, G21;
      G_S2_offdiag(thx, phx, thp, php, lam_list[k], &G12, &G21);
      const double a12 = std::abs(G12), a21 = std::abs(G21);
      const double asym = std::fabs(a12 - a21) / (0.5 * (a12 + a21));
      printf("    %8d  %14.6e %14.6e  %10.3e  %12.6f\n",
             lam_list[k], a12, a21, asym, a12 / ref);
      if (k == n_lam - 1) {
        if (asym > worst_asym) worst_asym = asym;
        const double me = std::fabs(a12 / ref - 1.0);
        if (me > worst_magerr) worst_magerr = me;
      }
    }
    printf("#\n");
  }

  printf("# ===== verdict (lambda_max = %d) =====\n", lam_list[n_lam - 1]);
  printf("#  worst off-diagonal asymmetry (exact-identity check) : %.3e\n", worst_asym);
  printf("#  worst magnitude error |G12|/CFT - 1                 : %.3e  -> %s\n",
         worst_magerr, (worst_magerr < 1.0e-2) ? "PASS" : "CHECK");

  return 0;
}
