// jj_cft_ylm_check_claude.cc
//
// PART A: check the Appendix C real-spherical-harmonic identities (qed3int_v2-11):
//   (C.5) orthonormality, and the selection rules (C.15)-(C.17)
//     (1/4pi) int Y_{l1m1}(n1) Y_{l2m2}(n2) {1, n12, n12^2}
//       = { d_{l1,0}d_{l2,0}d_{m1,0}d_{m2,0},  (1/3)d_{l1,1}d_{l2,1}d_{m1m2},
//           (1/3)d_{l1,0}d_{l2,0} + (2/15)d_{l1,2}d_{l2,2}d_{m1m2} },
//   verified via the factorized single-sphere Cartesian moments int Y, int Y n^a, int Y n^a n^b.
//
// PART B: check the Eq. (4.34)/(4.35) descendant tower.  G_t(t;n1,n2)=f^33 is isotropic (depends
//   only on x=n12), so Funk-Hecke gives  G^t_{ll}(t) = (1/2) int_{-1}^1 G_t(t;x) P_l(x) dx, diagonal
//   in l,m.  Prediction (Delta=2, our Cj=-1/(8pi^2)):  G^t_11 -> e^{-2t}/(24pi^2),
//   G^t_22 -> e^{-3t}/(10pi^2), G^t_00 -> 0;  decay rates (2,3); ratio [G22 e^{3t}]/[G11 e^{2t}] -> 12/5.
//
// G_t from the propagator: f^33 = -tr[sigma3 G(x1,x2) sigma3 G(x2,x1)].  Equatorial pairs (z=0,
// overflow-safe), nmax=40, IR.  G_prop/xi_mn(Boost)/c2_mn copied from cont_prop_eigbasis_claude.cpp.
// Pure host (GSL + Boost).  NO CUDA.

#include <cstdio>
#include <cmath>
#include <complex>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>
#include <boost/math/special_functions/jacobi.hpp>

typedef std::complex<double> cd;

// ---- propagator machinery (copied) ----------------------------------------------------------
static double xi_mn(double mm, int n, double z) {
  const int mpH = (int)std::lround(mm + 0.5);
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

static void G_prop(double th, double ph, double t, double thp, double php, double tp,
                   int nmax, cd G[2][2]) {
  const double z = std::cos(th), zp = std::cos(thp);
  const double tau = t - tp, atau = std::fabs(tau);
  const double st = (tau > 0.0) ? 1.0 : (tau < 0.0 ? -1.0 : 0.0);
  for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) G[a][b] = cd(0.0, 0.0);
  for (int twom = 1; twom <= 2 * nmax + 1; twom += 2) {
    const double mm = 0.5 * twom;
    for (int sgn = 0; sgn < 2; ++sgn) {
      const double iom = (sgn == 0) ? 1.0 : -1.0;
      const double m = iom * mm;
      for (int n = 0; n <= nmax; ++n) {
        const double lam = lambda_mn(mm, n);
        const double eM = std::exp(-lam * atau);
        if (eM < 1.0e-30) continue;
        const double sgnn = (n & 1) ? -1.0 : 1.0;
        const double pref = 1.0 / cnm2_mn(mm, n);
        const cd eph = std::exp(cd(0.0, m * ph)), ephp = std::exp(cd(0.0, m * php));
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

// G_t = f^33 = -tr[sigma3 G(x1,x2) sigma3 G(x2,x1)],  sigma3=diag(1,-1).
static double Gt_of(double th1, double ph1, double th2, double ph2, double t, int nmax) {
  cd G12[2][2], G21[2][2];
  G_prop(th1, ph1, t, th2, ph2, 0.0, nmax, G12);
  G_prop(th2, ph2, 0.0, th1, ph1, t, nmax, G21);
  // sigma3 G12 sigma3 = [[G12_00,-G12_01],[-G12_10,G12_11]];  then * G21, trace.
  cd M[2][2];
  M[0][0] =  G12[0][0]; M[0][1] = -G12[0][1];
  M[1][0] = -G12[1][0]; M[1][1] =  G12[1][1];
  cd tr = (M[0][0]*G21[0][0] + M[0][1]*G21[1][0]) + (M[1][0]*G21[0][1] + M[1][1]*G21[1][1]);
  return -tr.real();
}

// ---- real spherical harmonics Y_{l,m} (C.1) via GSL orthonormalized sphPlm ------------------
static double Ylm(int l, int m, double th, double ph) {
  const int am = std::abs(m);
  const double P = gsl_sf_legendre_sphPlm(l, am, std::cos(th));   // sqrt((2l+1)/4pi (l-m)!/(l+m)!) P_l^m
  if (m > 0) return std::sqrt(2.0) * P * std::cos(am * ph);
  if (m < 0) return std::sqrt(2.0) * P * std::sin(am * ph);
  return P;
}

// ---- Part A: Appendix C identities ----------------------------------------------------------
struct LM { int l, m; };
static const LM HARM[] = { {0,0}, {1,-1},{1,0},{1,1}, {2,-2},{2,-1},{2,0},{2,1},{2,2} };
static const int NH = 9;

static void run_partA() {
  const int Nu = 32, Np = 32;
  gsl_integration_glfixed_table* gt = gsl_integration_glfixed_table_alloc(Nu);

  // single-sphere moments: Iv0[h], Iv1[h][a], Iv2[h][a][b]  (a,b = x,y,z).
  double Iv0[NH], Iv1[NH][3], Iv2[NH][3][3], ortho[NH][NH];
  for (int h = 0; h < NH; ++h) { Iv0[h] = 0.0; for (int a=0;a<3;++a){Iv1[h][a]=0.0; for(int b=0;b<3;++b)Iv2[h][a][b]=0.0;} }
  for (int h1 = 0; h1 < NH; ++h1) for (int h2 = 0; h2 < NH; ++h2) ortho[h1][h2] = 0.0;

  for (int iu = 0; iu < Nu; ++iu) {
    double u, wu;
    gsl_integration_glfixed_point(-1.0, 1.0, iu, &u, &wu, gt);
    const double th = std::acos(u), s = std::sqrt(1.0 - u*u);
    for (int ip = 0; ip < Np; ++ip) {
      const double ph = 2.0 * M_PI * (ip + 0.5) / Np;
      const double wphi = 2.0 * M_PI / Np;
      const double w = wu * wphi;
      const double na[3] = { s*std::cos(ph), s*std::sin(ph), u };
      double Y[NH];
      for (int h = 0; h < NH; ++h) Y[h] = Ylm(HARM[h].l, HARM[h].m, th, ph);
      for (int h = 0; h < NH; ++h) {
        Iv0[h] += w * Y[h];
        for (int a = 0; a < 3; ++a) { Iv1[h][a] += w * Y[h] * na[a];
          for (int b = 0; b < 3; ++b) Iv2[h][a][b] += w * Y[h] * na[a] * na[b]; }
      }
      for (int h1 = 0; h1 < NH; ++h1) for (int h2 = 0; h2 < NH; ++h2) ortho[h1][h2] += w * Y[h1] * Y[h2];
    }
  }
  gsl_integration_glfixed_table_free(gt);

  // (C.5) orthonormality: ortho == identity.
  double e5 = 0.0;
  for (int h1 = 0; h1 < NH; ++h1) for (int h2 = 0; h2 < NH; ++h2)
    e5 = std::max(e5, std::fabs(ortho[h1][h2] - (h1 == h2 ? 1.0 : 0.0)));

  // (C.15)-(C.17): build (1/4pi) sum over the contracted Cartesian moments, compare to RHS.
  const double inv4pi = 1.0 / (4.0 * M_PI);
  double e15 = 0.0, e16 = 0.0, e17 = 0.0;
  for (int h1 = 0; h1 < NH; ++h1) for (int h2 = 0; h2 < NH; ++h2) {
    const int l1=HARM[h1].l, m1=HARM[h1].m, l2=HARM[h2].l, m2=HARM[h2].m;
    const double r15 = inv4pi * Iv0[h1] * Iv0[h2];
    double s16 = 0.0; for (int a=0;a<3;++a) s16 += Iv1[h1][a]*Iv1[h2][a];
    const double r16 = inv4pi * s16;
    double s17 = 0.0; for (int a=0;a<3;++a) for (int b=0;b<3;++b) s17 += Iv2[h1][a][b]*Iv2[h2][a][b];
    const double r17 = inv4pi * s17;
    const double p15 = (l1==0&&l2==0&&m1==0&&m2==0) ? 1.0 : 0.0;
    const double p16 = (l1==1&&l2==1&&m1==m2) ? 1.0/3.0 : 0.0;
    const double p17 = ((l1==0&&l2==0)?1.0/3.0:0.0) + ((l1==2&&l2==2&&m1==m2)?2.0/15.0:0.0);
    e15 = std::max(e15, std::fabs(r15 - p15));
    e16 = std::max(e16, std::fabs(r16 - p16));
    e17 = std::max(e17, std::fabs(r17 - p17));
  }
  printf("# ===== PART A: Appendix C identities (sphere quadrature %dx%d) =====\n", Nu, Np);
  printf("#  (C.5) orthonormality  max|<Yl1m1,Yl2m2> - delta| : %.3e\n", e5);
  printf("#  (C.15) int YY*1        max dev from RHS           : %.3e\n", e15);
  printf("#  (C.16) int YY*n12      max dev from RHS (1/3 l=1) : %.3e\n", e16);
  printf("#  (C.17) int YY*n12^2    max dev from RHS (1/3,2/15): %.3e\n", e17);
}

// ---- Part B: Eq. (4.35) tower via Legendre projection of G_t -----------------------------
static void run_partB(int nmax) {
  const int N = 48;
  gsl_integration_glfixed_table* gt = gsl_integration_glfixed_table_alloc(N);
  const double ts[] = {2.0, 3.0, 4.0, 5.0};
  const int nt = 4;
  const double th1 = M_PI / 2, ph1 = 0.0;

  printf("# ===== PART B: Eq.(4.35) tower  G^t_ll(t) = (1/2) int G_t(t;x) P_l(x) dx =====\n");
  printf("# predicted: G11 -> e^{-2t}/(24pi^2)=%.4e*e^{-2t}; G22 -> e^{-3t}/(10pi^2)=%.4e*e^{-3t}; G00->0\n",
         1.0/(24.0*M_PI*M_PI), 1.0/(10.0*M_PI*M_PI));
  printf("# %5s   %13s %13s %13s %13s   %10s %10s\n",
         "t", "G^t_00", "G^t_11", "G^t_22", "G^t_33", "G11*e2t", "G22*e3t");

  double Gprev[4] = {0,0,0,0}, G11e2[nt], G22e3[nt];
  for (int it = 0; it < nt; ++it) {
    const double t = ts[it];
    double Gl[4] = {0,0,0,0};
    for (int k = 0; k < N; ++k) {
      double x, w;
      gsl_integration_glfixed_point(-1.0, 1.0, k, &x, &w, gt);
      const double gamma = std::acos(x);
      const double g = Gt_of(th1, ph1, M_PI / 2, gamma, t, nmax);   // equatorial pair, n12=cos gamma
      for (int l = 0; l < 4; ++l) Gl[l] += 0.5 * w * g * gsl_sf_legendre_Pl(l, x);
    }
    G11e2[it] = Gl[1] * std::exp(2.0 * t);
    G22e3[it] = Gl[2] * std::exp(3.0 * t);
    printf("  %5.2f   %13.5e %13.5e %13.5e %13.5e   %10.4e %10.4e\n",
           t, Gl[0], Gl[1], Gl[2], Gl[3], G11e2[it], G22e3[it]);
    for (int l = 0; l < 4; ++l) {
      if (it > 0 && Gprev[l] != 0.0 && Gl[l] != 0.0)
        ; // effective rates printed below
      Gprev[l] = Gl[l];
    }
  }
  gsl_integration_glfixed_table_free(gt);

  // effective decay rates between the two largest t (IR), and the convention-free ratio.
  // recompute the last two t for rates:
  printf("#\n# IR readout (t=%.0f -> %.0f):\n", ts[nt-2], ts[nt-1]);
  printf("#   G11*e^{2t} : %.5e -> %.5e   (-> 1/(24pi^2)=%.5e)\n",
         G11e2[nt-2], G11e2[nt-1], 1.0/(24.0*M_PI*M_PI));
  printf("#   G22*e^{3t} : %.5e -> %.5e   (-> 1/(10pi^2)=%.5e)\n",
         G22e3[nt-2], G22e3[nt-1], 1.0/(10.0*M_PI*M_PI));
  printf("#   ratio [G22 e3t]/[G11 e2t] = %.5f  (-> 12/5 = 2.4)\n",
         G22e3[nt-1] / G11e2[nt-1]);
}

int main(int argc, char** argv) {
  const int nmax = (argc > 1) ? std::atoi(argv[1]) : 40;
  run_partA();
  printf("#\n");
  run_partB(nmax);
  return 0;
}
