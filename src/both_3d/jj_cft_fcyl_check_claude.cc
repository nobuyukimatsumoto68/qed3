// jj_cft_fcyl_check_claude.cc
//
// Check qed3int_v2-11 Eq. (4.26): the general (n1 != n2) cylinder conserved-current correlator
//   f_cyl^{ab}(t;n1,n2) = C_j Lambda^a_alpha(n1) Lambda^b_beta(n2) e^{-Delta t}
//                          [ delta^{ab} - 2 v^a v^b / D ] D^{-Delta},
//   v = n1 - e^{-t} n2,   D = 1 - 2 n12 e^{-t} + e^{-2t} = |v|^2,   Delta = 2,
// by contracting the validated eigenbasis propagator into the vector current:
//   f^{ab}(x1,x2) = - tr[ sigma^a G(x1,x2) sigma^b G(x2,x1) ],  x1=(n1,t), x2=(n2,0).
//
// FRAME-FREE checks (the bracket is the Householder reflection  delta - 2 vhat vhat , vhat=v/sqrt(D),
// so T = C_j e^{-Delta t} D^{-Delta} * reflection, and f = Lambda1 T Lambda2^T with Lambda in SO(3)):
//   (1) all THREE singular values of f equal  s0 = |C_j| e^{-2t} D^{-2}   [encodes Delta=2, D^{-Delta},
//       AND the reflection tensor];  (2) det f < 0  (reflection sign);  (3) the frame-invariant
//       temporal component  G_t = f^{33} = C_j e^{-2t} D^{-2} [ n12 - 2 (n1.vhat)(n2.vhat) ].
// With f = -tr[...] (N=1):  |C_j| = 1/(8 pi^2),  C_j = -1/(8 pi^2)  (fixed by the n1=n2 case, Eq.4.31).
// Singular values via the invariants of M=f^dag f (I1=tr, I2, I3=det); all-equal iff I2=3(I1/3)^2,
// I3=(I1/3)^3, then s0^2=I1/3.  No eigensolver needed.
//
// G_prop / xi_mn(Boost) / c2_mn copied from cont_prop_eigbasis_claude.cpp.  Pure host.  NO CUDA.

#include <cstdio>
#include <cmath>
#include <complex>

#include <gsl/gsl_sf_gamma.h>
#include <boost/math/special_functions/jacobi.hpp>

typedef std::complex<double> cd;

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

// full 2x2 propagator G[a][b](x,x'); copied from cont_prop_eigbasis_claude.cpp.
static void G_prop(double th, double ph, double t,
                   double thp, double php, double tp,
                   int nmax, cd G[2][2]) {
  const double z = std::cos(th), zp = std::cos(thp);
  const double tau = t - tp;
  const double atau = std::fabs(tau);
  const double st = (tau > 0.0) ? 1.0 : (tau < 0.0 ? -1.0 : 0.0);
  for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) G[a][b] = cd(0.0, 0.0);

  for (int twom = 1; twom <= 2 * nmax + 1; twom += 2) {
    const double mm = 0.5 * twom;
    for (int sgn = 0; sgn < 2; ++sgn) {
      const double iom = (sgn == 0) ? 1.0 : -1.0;
      const double m   = iom * mm;
      for (int n = 0; n <= nmax; ++n) {
        const double lam = lambda_mn(mm, n);
        const double eM  = std::exp(-lam * atau);
        if (eM < 1.0e-30) continue;
        const double sgnn = (n & 1) ? -1.0 : 1.0;
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

static const cd SIG[3][2][2] = {
  { { cd(0,0), cd(1,0) }, { cd(1,0), cd(0,0) } },
  { { cd(0,0), cd(0,-1)}, { cd(0,1), cd(0,0) } },
  { { cd(1,0), cd(0,0) }, { cd(0,0), cd(-1,0)} }
};
static void mat2_mul(const cd A[2][2], const cd B[2][2], cd C[2][2]) {
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      C[i][j] = A[i][0] * B[0][j] + A[i][1] * B[1][j];
}
static cd mat2_trace(const cd A[2][2]) { return A[0][0] + A[1][1]; }

// f^{ab}(t; n1,n2) = - tr[ sigma^a G(x1,x2) sigma^b G(x2,x1) ], x1=(th1,ph1,t), x2=(th2,ph2,0).
static void fcyl_full(double th1, double ph1, double th2, double ph2, double t, int nmax,
                      cd F[3][3]) {
  cd G12[2][2], G21[2][2];
  G_prop(th1, ph1, t, th2, ph2, 0.0, nmax, G12);
  G_prop(th2, ph2, 0.0, th1, ph1, t, nmax, G21);
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b) {
      cd t1[2][2], t2[2][2], t3[2][2];
      mat2_mul(SIG[a], G12, t1);
      mat2_mul(t1, SIG[b], t2);
      mat2_mul(t2, G21, t3);
      F[a][b] = -mat2_trace(t3);
    }
}

// 3x3 complex helpers ------------------------------------------------------------------------
static double max_imag(const cd F[3][3]) {
  double m = 0.0;
  for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) m = std::max(m, std::fabs(F[a][b].imag()));
  return m;
}
static cd det3(const cd A[3][3]) {
  return A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
       - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
       + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
}
// M = F^dag F (Hermitian).  Returns I1=tr M, I2=sum of principal 2x2 minors, I3=det M (all real).
static void ff_invariants(const cd F[3][3], double* I1, double* I2, double* I3) {
  cd M[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      cd s(0.0, 0.0);
      for (int k = 0; k < 3; ++k) s += std::conj(F[k][i]) * F[k][j];
      M[i][j] = s;
    }
  const double i1 = (M[0][0] + M[1][1] + M[2][2]).real();
  double trM2 = 0.0;                              // tr(M^2) = sum_ij |M_ij|^2 (M Hermitian)
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) trM2 += std::norm(M[i][j]);
  *I1 = i1;
  *I2 = 0.5 * (i1 * i1 - trM2);
  *I3 = det3(M).real();
}

int main(int argc, char** argv) {
  const int nmax = (argc > 1) ? std::atoi(argv[1]) : 40;
  const double Cj = 1.0 / (8.0 * M_PI * M_PI);     // |C_j| for the N=1 contraction

  printf("# Check Eq. (4.26) general n1!=n2 via f^{ab} = -tr[sig^a G(x1,x2) sig^b G(x2,x1)], nmax=%d\n", nmax);
  printf("# D = 1 - 2 n12 e^{-t} + e^{-2t}.  Prediction: all 3 singular values = |C_j| e^{-2t} D^{-2},\n");
  printf("#   |C_j| = 1/(8pi^2) = %.6e;  det f < 0;  G_t = f^33 = closed form.\n#\n", Cj);

  // n1 fixed on the equator; n2 swept on the equator so n12 = cos(gamma), z1=z2=0 (overflow-safe).
  const double th1 = M_PI / 2, ph1 = 0.0;
  const double gammas[] = {0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8};
  const int ng = (int)(sizeof(gammas) / sizeof(gammas[0]));
  const double ts[] = {0.8, 1.2, 2.0, 3.0};
  const int nt = (int)(sizeof(ts) / sizeof(ts[0]));

  printf("# %5s %7s %9s   %13s %10s %10s   %12s   %13s %13s\n",
         "t", "n12", "D", "s0(=sqrt I1/3)", "SV_eqresid", "ratio/Cj", "sign(det f)",
         "G_t", "G_t_closed");
  double ratio_min = 1e300, ratio_max = -1e300, imag_max = 0.0, eqresid_max = 0.0;
  for (int it = 0; it < nt; ++it) {
    const double t = ts[it];
    const double et = std::exp(-t);
    for (int ig = 0; ig < ng; ++ig) {
      const double g = gammas[ig];
      const double n12 = std::cos(g);
      const double th2 = M_PI / 2, ph2 = g;        // both equatorial -> n1.n2 = cos(gamma)
      cd F[3][3];
      fcyl_full(th1, ph1, th2, ph2, t, nmax, F);

      const double D = 1.0 - 2.0 * n12 * et + et * et;
      double I1, I2, I3;
      ff_invariants(F, &I1, &I2, &I3);
      const double s0sq = I1 / 3.0;
      // all-equal residual: how well I2=3 s0sq^2 and I3=s0sq^3 hold (rel).
      const double r2 = std::fabs(I2 - 3.0 * s0sq * s0sq) / (std::fabs(I2) + 1e-300);
      const double r3 = std::fabs(I3 - s0sq * s0sq * s0sq) / (std::fabs(I3) + 1e-300);
      const double eqresid = std::max(r2, r3);
      const double s0 = std::sqrt(s0sq);
      const double ref_sv = Cj * et * et / (D * D);  // |C_j| e^{-2t} D^{-2}
      const double ratio = s0 / ref_sv;              // -> 1 (i.e. |C_j| const)
      const double detf = det3(F).real();

      // G_t closed form (frame-invariant), C_j = -1/(8pi^2):
      const double sqrtD = std::sqrt(D);
      const double n1v = (1.0 - et * n12) / sqrtD;
      const double n2v = (n12 - et) / sqrtD;
      const double Gt_closed = -Cj * et * et / (D * D) * (n12 - 2.0 * n1v * n2v);
      const double Gt = F[2][2].real();

      printf("  %5.2f %7.4f %9.4f   %13.5e %10.2e %10.5f   %12.1f   %13.5e %13.5e\n",
             t, n12, D, s0, eqresid, ratio, (detf < 0 ? -1.0 : 1.0), Gt, Gt_closed);
      if (ratio < ratio_min) ratio_min = ratio;
      if (ratio > ratio_max) ratio_max = ratio;
      const double im = max_imag(F) / (s0 + 1e-300);
      if (im > imag_max) imag_max = im;
      if (eqresid > eqresid_max) eqresid_max = eqresid;
    }
    printf("#\n");
  }
  printf("# s0/ref(=|C_j|/Cj) over all (t,n12): [%.6f, %.6f]  spread %.2e  [const=1 => (4.26) shape OK]\n",
         ratio_min, ratio_max, ratio_max - ratio_min);
  printf("# worst equal-singular-value residual: %.2e  [-> 0 = reflection structure of (4.26)]\n",
         eqresid_max);
  printf("# worst |Im f|/s0: %.2e  [-> 0 = real Euclidean correlator]\n", imag_max);

  // isotropy spot-check: off-equator pair with the SAME n12 must give the same s0 and G_t.
  printf("#\n# isotropy spot-check (same n12, off-equator), t=1.2:\n");
  const double tt = 1.2, ett = std::exp(-tt);
  // pair A: equatorial gamma=1.2 (n12=cos1.2).  pair B: tilt both by building a rotated pair.
  const double n12tgt = std::cos(1.2);
  // pair B: n1=(0.6,0.0), n2 chosen so that n1.n2 = n12tgt.
  // n1=(sin.6,0,cos.6); pick n2=(sin t2 cos p2, sin t2 sin p2, cos t2) with dot = n12tgt; use t2=1.7,
  // solve p2: cos.6 cos1.7 + sin.6 sin1.7 cos p2 = n12tgt.
  const double a1 = std::cos(0.6) * std::cos(1.7);
  const double b1 = std::sin(0.6) * std::sin(1.7);
  double cosp = (n12tgt - a1) / b1;
  if (cosp > 1.0) cosp = 1.0; if (cosp < -1.0) cosp = -1.0;
  const double p2 = std::acos(cosp);
  cd FA[3][3], FB[3][3];
  fcyl_full(M_PI / 2, 0.0, M_PI / 2, 1.2, tt, nmax, FA);
  fcyl_full(0.6, 0.0, 1.7, p2, tt, nmax, FB);
  double IA[3], IB[3];
  ff_invariants(FA, &IA[0], &IA[1], &IA[2]);
  ff_invariants(FB, &IB[0], &IB[1], &IB[2]);
  const double DD = 1.0 - 2.0 * n12tgt * ett + ett * ett;
  printf("#   n12=%.4f D=%.4f : s0(eq)=%.6e  s0(tilt)=%.6e  G_t(eq)=%.6e  G_t(tilt)=%.6e\n",
         n12tgt, DD, std::sqrt(IA[0] / 3.0), std::sqrt(IB[0] / 3.0), FA[2][2].real(), FB[2][2].real());

  return 0;
}
