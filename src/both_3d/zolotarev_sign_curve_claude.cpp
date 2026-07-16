// zolotarev_sign_curve_claude.cpp
// _claude: dump the Zolotarev sign-function approximation curve R(x) ~ sign(x) for the ACTION operator D
// (n=21) and the FORCE operator Df (n_f in {5,7,9,11}) that share the SAME frozen window per L. The two
// operators differ ONLY in the pole count n; R(x) is evaluated on the scaled variable x = lambda/lambda_max
// in [k, 1], k = lambda_min/lambda_max. Fewer poles -> larger equioscillation error Delta and a gentler
// roll-off just below lambda_min (the "tamer force / avoid infinite slope" of the split-pole scheme).
// Zolotarev struct copied VERBATIM from includes/overlap_wmass_claude.h (matches D.Delta()); rational
// approx: A.D. Kennedy, hep-lat/0402038. Windows: frozen_window_claude.h / frozen_window_claude.md.
//
// Build: g++ -O2 -std=c++17 zolotarev_sign_curve_claude.cpp -o zolotarev_sign_curve_claude.x
// Run:   ./zolotarev_sign_curve_claude.x       (writes zolotarev_sign_L{1,2,4}_claude.dat + a Delta table)
// Plot:  gnuplot zolotarev_sign_curve_claude.gp (-> zolotarev_sign_L{1,2,4}_claude.png)

#include <cstdio>
#include <cmath>
#include <vector>

struct Zolotarev{
  using FLOAT = double;
  static constexpr FLOAT ONE = 1.0;
  static constexpr FLOAT TWO = 2.0;
  static constexpr FLOAT HALF = ONE/TWO;

  const int n;
  double k;
  double kp;
  const int size;

  std::vector<double> c;
  std::vector<double> cp;
  double M;
  double lambda_inv;

  double E;
  double C;
  std::vector<double> A;

  Zolotarev(const double k_=0.01,
	    const int n_=21 )
    : n(n_)
    , k(k_)
    , kp(std::sqrt(1.0-k*k))
    , size( int(n/2)+1 )
    , c(size, 0.0)
    , cp(size, 0.0)
    , A(size, 0.0)
  {
    get_coeffs();
    partial_fraction();
    C = E * 2.0 / (1.0+lambda_inv) / (k*M);
  }

  // A.D.Kennedy 2004; https://arxiv.org/abs/hep-lat/0402038
  static void sncndnK(const FLOAT z, const FLOAT k,
		      FLOAT& sn, FLOAT& cn, FLOAT& dn,
		      FLOAT& K) {
    FLOAT agm;
    int sgn;
    sn = arithgeom(z, ONE, std::sqrt(ONE - k*k), &agm);
    K = M_PI / (TWO * agm);
    sgn = ((int) (std::abs(z) / K)) % 4;
    sgn ^= sgn >> 1;
    sgn = 1 - ((sgn & 1) << 1);
    cn = ((FLOAT) sgn) * std::sqrt(ONE - sn * sn);
    dn = std::sqrt(ONE - k*k* sn * sn);
  }

  static FLOAT arithgeom(const FLOAT z, FLOAT a, FLOAT b, FLOAT* agm) {
    static FLOAT pb = -ONE;
    FLOAT xi;
    if (b <= pb) {
      pb = -ONE;
      *agm = a;
      return std::sin(z * a);
    }
    pb = b;
    xi = arithgeom(z, HALF*(a+b), std::sqrt(a*b), agm);
    return 2*a*xi / ((a+b) + (a-b)*xi*xi);
  }

  void get_coeffs() {
    double Kp = 1.0, xibar;
    for(int m=0; m<size; m++){
      double sn, cn, dn;
      double z = 2.0 * Kp * m / n;
      sncndnK( z, kp, sn, cn, dn, Kp );
      if(m==0) continue;
      c[m] = - std::pow( cn / sn, 2 );
      z = 2.0 * Kp * (m-0.5) / n;
      sncndnK( z, kp, sn, cn, dn, Kp );
      cp[m] = - std::pow( cn / sn, 2 );
      if(m==1) xibar = 1.0/dn;
    }

    M = 1.0;
    for(int m=1; m<size; m++) M *= (1.0-c[m]) / (1.0-cp[m]);

    lambda_inv = xibar / M;
    for(int m=1; m<size; m++) lambda_inv *= (1.0-c[m]*xibar*xibar) / (1.0-cp[m]*xibar*xibar);
  }

  void partial_fraction() {
    E = 1.0;
    for(int m=1; m<size; m++) {
      double numer = 1.0, denom = 1.0;
      for(int ell=1; ell<size; ell++) {
	numer *= k*k/cp[m] - k*k/c[ell];
	if(m!=ell) denom *= k*k/cp[m] - k*k/cp[ell];
      }
      A[m] = numer/denom;
      E*=c[m]/cp[m];
    }
  }

  // partial-fraction sign approximation R(x) ~ sign(x) on x in [k,1] (copied from OverlapWMass::operator()).
  double operator()( const double x ) const {
    double res = 1.0;
    for(int m=1; m<size; m++) res += A[m] / (x*x - k*k/cp[m]);
    res *= C * x;
    return res;
  }

  // analytic first derivative w.r.t. the SCALED variable x: R(x)=C x (1 + sum A_m/(x^2-b_m)), b_m=k^2/cp[m],
  // so R'(x) = C[ 1 + sum A_m/(x^2-b_m) - 2 x^2 sum A_m/(x^2-b_m)^2 ]. (dR/dlambda = R'(x)/lambda_max.)
  double deriv( const double x ) const {
    double s1 = 0.0, s2 = 0.0;
    for(int m=1; m<size; m++){
      const double b = k*k/cp[m];
      const double d = x*x - b;
      s1 += A[m]/d;
      s2 += A[m]/(d*d);
    }
    return C * ( 1.0 + s1 - 2.0*x*x*s2 );
  }

  inline double Delta() const {
    return (lambda_inv - 1.0) / (lambda_inv + 1.0);
  }
};

int main(){
  // FINALIZED split (NM 2026-07-14): ACTION op D = n=31 on the full frozen window [lambda_min, lambda_max];
  // FORCE op Df = n_f=11 on the NARROWED window [2*lambda_min, lambda_max] (lambda_min*=2). Same lambda_max
  // and same D_W, so both are evaluated at the same scaled x = lambda/lambda_max. Frozen windows per L from
  // frozen_window_claude.h.
  struct Pick {
    int L;
    double lmin;
    double lmax;
    int nact;    // per-L action pole count (hasenbusch_naction)
  };
  std::vector<Pick> picks = {
    { 1, 0.1,   13.0, 25 },
    { 2, 0.06,  8.0,  25 },
    { 4, 0.008, 5.0,  31 },
  };
  const int    N_FRC   = 11;   // force pole count
  const double WIN_FRC = 2.0;  // force lambda_min = WIN_FRC * lambda_min (window narrowing)

  std::printf("# Zolotarev sign approx: ACTION n=n_act(L) (window [lmin,lmax]) vs FORCE n_f=%d (window [%g*lmin,lmax])\n",
              N_FRC, WIN_FRC);
  std::printf("# %-4s %-6s %-10s %-10s | %-13s %-13s | %-13s %-13s %-13s\n",
              "L", "n_act", "lambda_min", "lambda_max", "k_action", "Delta_action", "lmin_force", "k_force", "Delta_force");

  for(const Pick& p : picks){
    const double k_act = p.lmin / p.lmax;
    const double lmin_f = WIN_FRC * p.lmin;
    const double k_frc = lmin_f / p.lmax;

    Zolotarev z_act( k_act, p.nact );
    Zolotarev z_frc( k_frc, N_FRC );   // force n_f, NARROWED window (w=2): k = 2*lmin/lmax
    Zolotarev z_frc1( k_act, N_FRC );  // force n_f, FULL window (w=1): k = lmin/lmax (same window as the action)

    std::printf("  %-4d %-6d %-10.5f %-10.3f | %-13.6e %-13.6e | %-13.6f %-13.6e %-13.6e | Delta_force(w=1,n=%d)=%-13.6e\n",
                p.L, p.nact, p.lmin, p.lmax, k_act, z_act.Delta(), lmin_f, k_frc, z_frc.Delta(), N_FRC, z_frc1.Delta());

    // curve file: physical eigenvalue lambda (log-spaced, from ~0 to lmax), scaled x=lambda/lmax.
    // Columns: lambda, R_action, R_force(w2), R_force(w1), dRdlam_action, dRdlam_force(w2), dRdlam_force(w1).
    char fname[128];
    std::snprintf( fname, sizeof(fname), "zolotarev_sign_L%d_claude.dat", p.L );
    FILE* fp = std::fopen( fname, "w" );
    std::fprintf( fp, "# L=%d  lambda_min=%.6f  lambda_max=%.6f  lmin_force(w=2)=%.6f\n", p.L, p.lmin, p.lmax, lmin_f );
    std::fprintf( fp, "# k_action=%.6e Delta_action(n=%d)=%.6e  k_force(w2)=%.6e Delta_force(w2,n=%d)=%.6e  Delta_force(w1,n=%d)=%.6e\n",
                  k_act, p.nact, z_act.Delta(), k_frc, N_FRC, z_frc.Delta(), N_FRC, z_frc1.Delta() );
    std::fprintf( fp, "# lambda   R_action   R_force_w2   R_force_w1   dRdlam_action   dRdlam_force_w2   dRdlam_force_w1\n" );

    const int NP = 4000;
    const double lam_lo = 1.0e-3 * p.lmin;   // start ~at 0 (log-spaced; the tiny gap is invisible on a linear axis)
    const double lam_hi = p.lmax;
    const double dlog = std::log( lam_hi / lam_lo ) / (NP-1);
    for(int i=0; i<NP; i++){
      const double lam = lam_lo * std::exp( dlog * i );
      const double x = lam / p.lmax;                 // scaled variable
      std::fprintf( fp, "%.8e  %.8e  %.8e  %.8e  %.8e  %.8e  %.8e\n",
                    lam, z_act(x), z_frc(x), z_frc1(x),
                    z_act.deriv(x)/p.lmax, z_frc.deriv(x)/p.lmax, z_frc1.deriv(x)/p.lmax );
    }
    std::fclose( fp );
    std::printf("#   wrote %s (%d points, lambda in [%.4f, %.3f])\n", fname, NP, lam_lo, lam_hi);
  }
  return 0;
}
