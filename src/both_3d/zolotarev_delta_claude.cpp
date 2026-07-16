// zolotarev_delta_claude.cpp
// _claude: Zolotarev sign-function approximation error Delta(n,k) for the chosen overlap windows.
// Delta = (lambda_inv-1)/(lambda_inv+1) depends ONLY on the pole count n and the window ratio
// k = lambda_min/lambda_max (NOT on the gauge config -- the config only decides whether its own ratio
// falls inside [k,1]). Zolotarev struct copied VERBATIM from includes/overlap_wmass_claude.h so the
// number matches the operator's D.Delta(). Rational-approx algorithm: A.D. Kennedy, hep-lat/0402038.
//
// Build: g++ -O2 -std=c++17 zolotarev_delta_claude.cpp -o zolotarev_delta_claude.x
// Run:   ./zolotarev_delta_claude.x

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

  inline double Delta() const {
    return (lambda_inv - 1.0) / (lambda_inv + 1.0);
  }
};

int main(){
  const int n = 21;

  struct Pick {
    const char* L;
    double lmin;
    double lmax;
  };
  std::vector<Pick> picks = {
    { "L=1", 0.1,   13.0 },
    { "L=2", 0.06,  8.0  },
    { "L=4", 0.008, 5.0  },
  };

  std::printf("# Zolotarev sign-approx error Delta(n=%d, k=lambda_min/lambda_max)\n", n);
  std::printf("# %-5s %-10s %-10s %-14s %-14s\n", "L", "lambda_min", "lambda_max", "k", "Delta");
  for(const Pick& p : picks){
    const double k = p.lmin / p.lmax;
    Zolotarev z( k, n );
    std::printf("  %-5s %-10.5f %-10.3f %-14.6e %-14.6e\n", p.L, p.lmin, p.lmax, k, z.Delta());
  }

  // reference: the current fixed production window
  {
    const double k = 0.001;
    Zolotarev z( k, n );
    std::printf("  %-5s %-10s %-10s %-14.6e %-14.6e   (current fixed window)\n", "ref", "-", "-", k, z.Delta());
  }
  return 0;
}
