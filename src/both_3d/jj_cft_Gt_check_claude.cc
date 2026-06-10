// jj_cft_Gt_check_claude.cc
//
// Check qed3int_v2-11 Eq. (4.31): the temporally projected conserved-current two-point function
//   G_t(t; n,n) = -C_J e^{-Delta t} (1 - e^{-t})^{-2 Delta},   Delta = 2,
// (defined in Eq. (4.30)) by contracting the validated continuum eigenbasis propagator G = D^-1
// into the conserved VECTOR current  j_V^a = eta^dag sigma^a xi - xi^dag sigma^a eta  (Eqs. 3.6-3.7),
// with the action S_2 = eta^dag D xi - xi^dag D eta (Eq. 1.24).  See jj_cft_check_equations_claude.md.
//
// Wick contraction (props: <xi eta^dag>=G, <eta xi^dag>=-G, <xi xi^dag>=<eta eta^dag>=0; the ++ and
// -- pieces are equal, +- cross terms vanish):
//   f_cyl^{ab}(x1,x2) = - tr[ sigma^a G(x1,x2) sigma^b G(x2,x1) ]   (overall const absorbed into C_J).
// Temporal projection e^a_3 e^b_3 -> the a=b=3 (sigma^3) component (the global time axis in this
// frame; the same-site block is pure sigma_3):  G_t = f^{33},  G_s = f^{11}+f^{22}.
//
// Checks: (i) shape  G_t / [e^{-2t}(1-e^{-t})^{-4}] = const;  (ii) n-independence;
//   (iii) G_s/G_t = -(D-1) = -2;  (iv) agreement with the closed form
//   G_t = 2 c_3(t)^2 = (1/8pi^2) e^{-2t}(1-e^{-t})^{-4},  c_3(t) = (1/4pi) e^{-t}/(1-e^{-t})^2.
//
// xi_mn / c2_mn / G_prop are copied from cont_prop_eigbasis_claude.cpp (Boost-recurrence xi).
// Pure host (GSL + Boost Jacobi).  NO CUDA.  Single-thread.

#include <cstdio>
#include <cmath>
#include <complex>

#include <gsl/gsl_sf_gamma.h>
#include <boost/math/special_functions/jacobi.hpp>

typedef std::complex<double> cd;

// ---- xi via Boost Jacobi recurrence (stable); mm = |m| (half-integer) -----------------------
static double xi_mn(double mm, int n, double z) {
  const int mpH = (int)std::lround(mm + 0.5);                 // |m| + 1/2
  if (std::abs(z - 1.0) < 1.0e-14) return (mpH == 1) ? 1.0 / std::sqrt(2.0) : 0.0;
  if (std::abs(z + 1.0) < 1.0e-14) return 0.0;
  const double factor = std::pow(1.0 - z, 0.5 * (mpH - 1)) * std::pow(1.0 + z, -0.5 * mpH);
  const double poly   = boost::math::jacobi(n + mpH, (double)(mpH - 1.0), (double)(-1.0 * mpH), z);
  return factor * poly;
}

static double c2_mn(double mm, int n) {
  const double lnum = gsl_sf_lngamma(n + 1.0) + gsl_sf_lngamma(n + 2.0 * mm + 1.0);
  const double lden = gsl_sf_lngamma(n + mm + 0.5) + gsl_sf_lngamma(n + mm + 1.5);
  return std::exp(lnum - lden) / (2.0 * n + 2.0 * mm + 1.0);
}

static inline double lambda_mn(double mm, int n) { return (double)n + mm + 0.5; }
static inline double cnm2_mn(double mm, int n)   { return 4.0 * M_PI * c2_mn(mm, n); }

// ---- full 2x2 propagator G[a][b](x,x'), Eq. (C.28); copied from cont_prop_eigbasis_claude.cpp --
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
    for (int sgn = 0; sgn < 2; ++sgn) {                   // iota_m = +1, -1
      const double iom = (sgn == 0) ? 1.0 : -1.0;
      const double m   = iom * mm;
      for (int n = 0; n <= nmax; ++n) {
        const double lam = lambda_mn(mm, n);
        const double eM  = std::exp(-lam * atau);
        if (eM < 1.0e-30) continue;
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

// ---- Pauli matrices and 2x2 helpers ---------------------------------------------------------
static const cd SIG[3][2][2] = {
  { { cd(0,0), cd(1,0) }, { cd(1,0), cd(0,0) } },          // sigma^1
  { { cd(0,0), cd(0,-1)}, { cd(0,1), cd(0,0) } },          // sigma^2
  { { cd(1,0), cd(0,0) }, { cd(0,0), cd(-1,0)} }           // sigma^3
};

static void mat2_mul(const cd A[2][2], const cd B[2][2], cd C[2][2]) {
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      C[i][j] = A[i][0] * B[0][j] + A[i][1] * B[1][j];
}

static cd mat2_trace(const cd A[2][2]) { return A[0][0] + A[1][1]; }

// f_cyl^{ab}(t; n,n) = - tr[ sigma^a G(x1,x2) sigma^b G(x2,x1) ],  x1=(th,ph,t), x2=(th,ph,0).
static void fcyl_diag(double th, double ph, double t, int nmax, cd F[3][3]) {
  cd G12[2][2], G21[2][2];
  G_prop(th, ph, t, th, ph, 0.0, nmax, G12);     // G(x1, x2)
  G_prop(th, ph, 0.0, th, ph, t, nmax, G21);     // G(x2, x1)
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      cd t1[2][2], t2[2][2], t3[2][2];
      mat2_mul(SIG[a], G12, t1);                 // sigma^a G(x1,x2)
      mat2_mul(t1, SIG[b], t2);                  // sigma^a G sigma^b
      mat2_mul(t2, G21, t3);                      // sigma^a G sigma^b G(x2,x1)
      F[a][b] = -mat2_trace(t3);
    }
  }
}

// closed-form same-site temporal coefficient c_3(t) = (1/4pi) e^{-t}/(1-e^{-t})^2  (Eq. C.29 sum).
static double c3_analytic(double t) {
  const double e = std::exp(-t);
  return e / (4.0 * M_PI * (1.0 - e) * (1.0 - e));
}

int main(int argc, char** argv) {
  // nmax=40 (the validated cont_prop default).  Do NOT push much higher: at a GENERIC site the
  // same-site sum runs the full |m|-tower and |xi(-z)|^2 at mm ~ 100 overflows double (>1e308)
  // before the e^{-lambda t} damping tames it (worst at small t and z near the poles).  nmax=40
  // converges c_3(t) for t >~ 0.8 without overflow; smaller t needs the closed form, not the sum.
  const int nmax = (argc > 1) ? std::atoi(argv[1]) : 40;

  const double thref = 1.0, phref = 0.3;          // reference n-hat (z = 0.54, away from poles)
  const double ts[] = {0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0};
  const int nt = (int)(sizeof(ts) / sizeof(ts[0]));

  printf("# Check Eq. (4.31): G_t(t;n,n) = -C_J e^{-2t}(1-e^{-t})^{-4}  (Delta=2), nmax=%d\n", nmax);
  printf("# G_t = f^33, G_s = f^11+f^22, f^{ab} = -tr[sig^a G(x1,x2) sig^b G(x2,x1)].\n");
  printf("# ref(t) = e^{-2t}(1-e^{-t})^{-4}.  G_t/ref should be CONST = 1/(8pi^2) = %.6e.\n",
         1.0 / (8.0 * M_PI * M_PI));
  printf("#\n# %6s %14s %14s %10s %14s %12s %14s\n",
         "t", "G_t", "G_s", "G_s/G_t", "ref", "G_t/ref", "2*c3^2(anal)");

  double rmin = 1e300, rmax = -1e300, ratio_sum = 0.0;
  double gsgt_min = 1e300, gsgt_max = -1e300;
  for (int i = 0; i < nt; ++i) {
    const double t = ts[i];
    cd F[3][3];
    fcyl_diag(thref, phref, t, nmax, F);
    const double Gt = F[2][2].real();
    const double Gs = F[0][0].real() + F[1][1].real();
    const double ref = std::exp(-2.0 * t) * std::pow(1.0 - std::exp(-t), -4.0);
    const double ratio = Gt / ref;
    const double c3 = c3_analytic(t);
    const double Gt_anal = 2.0 * c3 * c3;
    const double gsgt = Gs / Gt;
    printf("  %6.2f %14.6e %14.6e %10.5f %14.6e %12.6e %14.6e\n",
           t, Gt, Gs, gsgt, ref, ratio, Gt_anal);
    if (ratio < rmin) rmin = ratio;
    if (ratio > rmax) rmax = ratio;
    ratio_sum += ratio;
    if (gsgt < gsgt_min) gsgt_min = gsgt;
    if (gsgt > gsgt_max) gsgt_max = gsgt;
  }
  const double rmean = ratio_sum / nt;
  printf("#\n# G_t/ref: mean=%.6e  spread=%.3e (rel %.2e)  [const => shape matches (4.31)]\n",
         rmean, rmax - rmin, (rmax - rmin) / std::fabs(rmean));
  printf("# G_s/G_t in [%.6f, %.6f]  [-> -2 = -(D-1), Eqs. (4.28)/(4.31)]\n", gsgt_min, gsgt_max);

  // n-independence: G_t at several n-hat, fixed t.
  const double tcheck = 1.5;
  // n-hat points kept near the equator (|z| <~ 0.6) to avoid the near-pole |xi(-z)| blow-up.
  const double nn[][2] = { {1.0, 0.3}, {2.0, 1.5}, {M_PI / 2, 0.0}, {1.3, 2.2} };
  const int nn_n = (int)(sizeof(nn) / sizeof(nn[0]));
  printf("#\n# n-independence of G_t at t=%.2f:\n", tcheck);
  double gmin = 1e300, gmax = -1e300;
  for (int k = 0; k < nn_n; ++k) {
    cd F[3][3];
    fcyl_diag(nn[k][0], nn[k][1], tcheck, nmax, F);
    const double Gt = F[2][2].real();
    printf("#   (theta,phi)=(%.3f,%.3f)  G_t=%.10e\n", nn[k][0], nn[k][1], Gt);
    if (Gt < gmin) gmin = Gt;
    if (Gt > gmax) gmax = Gt;
  }
  printf("#   G_t spread across n-hat: %.3e (rel %.2e)  [-> 0 = isotropic]\n",
         gmax - gmin, (gmax - gmin) / std::fabs(0.5 * (gmax + gmin)));

  // structure of f^{ab} at one (t, n-hat): should be diag(-2,-2,+2) c3^2.
  cd F[3][3];
  fcyl_diag(thref, phref, 1.0, nmax, F);
  printf("#\n# f^{ab} at t=1.0 (expect diagonal ~ (-2,-2,+2) c3^2, c3^2=%.6e):\n",
         c3_analytic(1.0) * c3_analytic(1.0));
  for (int a = 0; a < 3; ++a)
    printf("#   [ %12.5e %12.5e %12.5e ]\n",
           F[a][0].real(), F[a][1].real(), F[a][2].real());

  return 0;
}
