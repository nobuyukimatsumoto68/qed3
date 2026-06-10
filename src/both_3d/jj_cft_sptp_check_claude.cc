// jj_cft_sptp_check_claude.cc
//
// Spatial (sp) and temporal (tp) projections of the general n1 != n2 cylinder current correlator
// f^{ab}(t;n1,n2) = -tr[sigma^a G(x1,x2) sigma^b G(x2,x1)], for random IR t, vs the canonical-frame
// (4.26) prediction.  Companion to jj_cft_fcyl_check_claude.cc (which did the frame-free singular
// values + G_t).
//
// tp:  G_t = f^33  (our spin frame has the radial/time leg = sigma_3 globally, so f^33 = the
//      temporal-temporal projection; canonical pred = nhat_1 . T . nhat_2).
// sp:  G_s = f^11 + f^22  (our spin frame's sigma_1,sigma_2 spatial legs).  Canonical pred uses the
//      (4.23)-(4.24) frame:  G_s = thetahat_1 . T . thetahat_2 + phihat_1 . T . phihat_2.
//
// Embedding (4.26) tensor:  T^{ab} = Cj e^{-Delta t} D^{-Delta} (delta^{ab} - 2 vhat^a vhat^b),
//   v = nhat1 - e^{-t} nhat2,  D = |v|^2 = 1 - 2 n12 e^{-t} + e^{-2t},  Delta = 2,  Cj = -1/(8 pi^2)
//   (signed; fixes the n1=n2 case to Eq.4.31).  a.T.b = Cj e^{-2t} D^{-2} ( a.b - 2 (a.vhat)(b.vhat) ).
//
// TEST: does our f^11+f^22 equal G_s^pred (=> spin frame spatial legs align with canonical theta,phi,
//   and G_s is the clean vierbein-independent sp projection)?  Plus an isotropy check (matched n12,
//   different absolute orientation -> same G_s, G_t).
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

// ---- 3-vector helpers and the canonical (4.23)-(4.24) frame ---------------------------------
struct V3 { double x, y, z; };
static double dot(const V3& a, const V3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static V3 nhat_of(double th, double ph)   { return { std::sin(th)*std::cos(ph), std::sin(th)*std::sin(ph), std::cos(th) }; }
static V3 phihat_of(double th, double ph) { return { -std::sin(ph), std::cos(ph), 0.0 }; }
static V3 thetahat_of(double th, double ph){ return { std::cos(th)*std::cos(ph), std::cos(th)*std::sin(ph), -std::sin(th) }; }

// a . T . b  with  T = Cj e^{-2t} D^{-2} (I - 2 vhat vhat^T).
static double aTb(const V3& a, const V3& b, const V3& vh, double K) {
  return K * (dot(a, b) - 2.0 * dot(a, vh) * dot(b, vh));
}

// ---- 3x3 real helpers, the canonical Lambda rows, and the embedding f / projector trace ------
static void mat3_mul(const double A[3][3], const double B[3][3], double C[3][3]) {
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      double s = 0.0;
      for (int k = 0; k < 3; ++k) s += A[i][k] * B[k][j];
      C[i][j] = s;
    }
}
static double mat3_tr(const double A[3][3]) { return A[0][0] + A[1][1] + A[2][2]; }
static void proj_tan(const V3& n, double P[3][3]) {        // P = I - n n^T (tangent projector)
  const double v[3] = { n.x, n.y, n.z };
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      P[i][j] = (i == j ? 1.0 : 0.0) - v[i] * v[j];
}
static void comp3(const V3& a, double out[3]) { out[0] = a.x; out[1] = a.y; out[2] = a.z; }

// embedding f^{alpha beta} = Lambda1^T f_canon Lambda2 (rows of Lambda = thetahat, phihat, nhat).
static void femb_from_canon(const cd F[3][3], double th1, double ph1, double th2, double ph2,
                            double fe[3][3]) {
  double L1[3][3], L2[3][3];                                 // L[a][alpha] = a-th frame vector comp
  comp3(thetahat_of(th1, ph1), L1[0]); comp3(phihat_of(th1, ph1), L1[1]); comp3(nhat_of(th1, ph1), L1[2]);
  comp3(thetahat_of(th2, ph2), L2[0]); comp3(phihat_of(th2, ph2), L2[1]); comp3(nhat_of(th2, ph2), L2[2]);
  for (int al = 0; al < 3; ++al)
    for (int be = 0; be < 3; ++be) {
      double s = 0.0;
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          s += L1[a][al] * L2[b][be] * F[a][b].real();
      fe[al][be] = s;
    }
}
// frame-free transverse-transverse projection  tr[ P1 M P2 ].
static double tr_P1_M_P2(const double M[3][3], const V3& n1, const V3& n2) {
  double P1[3][3], P2[3][3], tmp[3][3], out[3][3];
  proj_tan(n1, P1);
  proj_tan(n2, P2);
  mat3_mul(P1, M, tmp);
  mat3_mul(tmp, P2, out);
  return mat3_tr(out);
}

static void projections(double th1, double ph1, double th2, double ph2, double t, int nmax,
                        double* Gt_our, double* Gs_our, double* Gt_pred, double* Gs_pred,
                        double* Gs_inv_our, double* Gs_inv_pred, double* n12_out) {
  cd F[3][3];
  fcyl_full(th1, ph1, th2, ph2, t, nmax, F);
  *Gt_our = F[2][2].real();
  *Gs_our = F[0][0].real() + F[1][1].real();

  const V3 n1 = nhat_of(th1, ph1), n2 = nhat_of(th2, ph2);
  const double n12 = dot(n1, n2);
  *n12_out = n12;
  const double et = std::exp(-t);
  const double D  = 1.0 - 2.0 * n12 * et + et * et;
  const double Cj = -1.0 / (8.0 * M_PI * M_PI);          // signed
  const double K  = Cj * et * et / (D * D);
  const double sD = std::sqrt(D);
  const V3 vh = { (n1.x - et*n2.x)/sD, (n1.y - et*n2.y)/sD, (n1.z - et*n2.z)/sD };

  const V3 th1v = thetahat_of(th1, ph1), ph1v = phihat_of(th1, ph1);
  const V3 th2v = thetahat_of(th2, ph2), ph2v = phihat_of(th2, ph2);
  *Gt_pred = aTb(n1, n2, vh, K);
  *Gs_pred = aTb(th1v, th2v, vh, K) + aTb(ph1v, ph2v, vh, K);

  // frame-free spatial projection  G_s^inv = tr[P1 f_emb P2]  (isotropic; vierbein-independent).
  double fe[3][3];
  femb_from_canon(F, th1, ph1, th2, ph2, fe);
  *Gs_inv_our = tr_P1_M_P2(fe, n1, n2);
  double vha[3];
  comp3(vh, vha);
  double Tm[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      Tm[i][j] = K * ((i == j ? 1.0 : 0.0) - 2.0 * vha[i] * vha[j]);
  *Gs_inv_pred = tr_P1_M_P2(Tm, n1, n2);
}

int main(int argc, char** argv) {
  const int nmax = (argc > 1) ? std::atoi(argv[1]) : 40;
  printf("# sp/tp projections of f^{ab}(t;n1,n2) for n1!=n2, vs canonical-frame (4.26), nmax=%d\n", nmax);
  printf("# tp: G_t = f^33 vs nhat1.T.nhat2 ;  sp: G_s = f^11+f^22 vs theta1.T.theta2 + phi1.T.phi2.\n#\n");

  // n1 fixed equatorial; n2 swept equatorial (n12 = cos gamma).  Two IR times.
  const double th1 = M_PI / 2, ph1 = 0.0;
  const double gammas[] = {0.5, 1.0, 1.5, 2.0, 2.5};
  const int ng = (int)(sizeof(gammas) / sizeof(gammas[0]));
  const double ts[] = {1.7, 2.3};
  const int nt = (int)(sizeof(ts) / sizeof(ts[0]));

  printf("# %5s %8s   %13s %8s   %13s %8s   %13s %13s %8s\n",
         "t", "n12", "G_t(our)", "relerr", "G_s(our)", "relerr", "Gs_inv(our)", "Gs_inv(pred)", "relerr");
  for (int it = 0; it < nt; ++it) {
    for (int ig = 0; ig < ng; ++ig) {
      const double th2 = M_PI / 2, ph2 = gammas[ig];
      double Gt_our, Gs_our, Gt_pred, Gs_pred, Gsi_our, Gsi_pred, n12;
      projections(th1, ph1, th2, ph2, ts[it], nmax,
                  &Gt_our, &Gs_our, &Gt_pred, &Gs_pred, &Gsi_our, &Gsi_pred, &n12);
      const double et = std::fabs(Gt_pred) > 0 ? std::fabs(Gt_our - Gt_pred)/std::fabs(Gt_pred) : 0.0;
      const double es = std::fabs(Gs_pred) > 0 ? std::fabs(Gs_our - Gs_pred)/std::fabs(Gs_pred) : 0.0;
      const double ei = std::fabs(Gsi_pred) > 0 ? std::fabs(Gsi_our - Gsi_pred)/std::fabs(Gsi_pred) : 0.0;
      printf("  %5.2f %8.4f   %13.5e %8.1e   %13.5e %8.1e   %13.5e %13.5e %8.1e\n",
             ts[it], n12, Gt_our, et, Gs_our, es, Gsi_our, Gsi_pred, ei);
    }
    printf("#\n");
  }

  // Isotropy: matched n12 (=cos 1.2) at two absolute orientations; G_t and G_s must agree.
  printf("# isotropy (matched n12=cos(1.2)=%.4f, t=1.7): equatorial vs tilted pair\n", std::cos(1.2));
  const double tt = 1.7;
  double a_Gt, a_Gs, a_Gtp, a_Gsp, a_Gsi, a_Gsip, a_n12;
  projections(M_PI / 2, 0.0, M_PI / 2, 1.2, tt, nmax,
              &a_Gt, &a_Gs, &a_Gtp, &a_Gsp, &a_Gsi, &a_Gsip, &a_n12);
  // tilted pair with same n12: n1=(0.6,0), n2=(1.7,p2) solving dot = cos(1.2).
  const double n12t = std::cos(1.2);
  const double A1 = std::cos(0.6) * std::cos(1.7);
  const double B1 = std::sin(0.6) * std::sin(1.7);
  double cp = (n12t - A1) / B1;
  if (cp > 1.0) cp = 1.0;
  if (cp < -1.0) cp = -1.0;
  const double p2 = std::acos(cp);
  double b_Gt, b_Gs, b_Gtp, b_Gsp, b_Gsi, b_Gsip, b_n12;
  projections(0.6, 0.0, 1.7, p2, tt, nmax,
              &b_Gt, &b_Gs, &b_Gtp, &b_Gsp, &b_Gsi, &b_Gsip, &b_n12);
  printf("#   equatorial: n12=%.4f  G_t=%.6e  G_s=%.6e  Gs_inv=%.6e\n", a_n12, a_Gt, a_Gs, a_Gsi);
  printf("#   tilted    : n12=%.4f  G_t=%.6e  G_s=%.6e  Gs_inv=%.6e\n", b_n12, b_Gt, b_Gs, b_Gsi);
  printf("#   |dG_t|/G_t=%.2e   |dG_s|/G_s=%.2e   |dGs_inv|/Gs_inv=%.2e\n",
         std::fabs(a_Gt - b_Gt)/std::fabs(a_Gt), std::fabs(a_Gs - b_Gs)/std::fabs(a_Gs),
         std::fabs(a_Gsi - b_Gsi)/std::fabs(a_Gsi));
  printf("#   [G_t and Gs_inv -> 0 (isotropic, frame-free); G_s (canonical theta-phi) stays O(1e-2)]\n");

  return 0;
}
