// test_multishift_claude.cpp
// C1 (host only, no GPU): validate the multi-shift CG recurrence (Jegerlehner,
// hep-lat/9612014) against independent per-shift CG solves, on a small random SPD
// matrix with positive shifts mimicking the Zolotarev pole structure
//   (A + \sigma_m) Z_m = b,   A SPD,   \sigma_m = -k^2/c'_m > 0.
// This de-risks the \zeta_m recurrence and the seed (= smallest shift) choice before
// any CUDA. The recurrence is shift-value-agnostic, so a representative positive
// shift set suffices to prove algorithmic equivalence.
//
// Compile:
//   g++ -O2 -std=c++17 -I/opt/eigen-3.4.0 test_multishift_claude.cpp -o test_multishift_claude.out
// Run:
//   ./test_multishift_claude.out

#include <Eigen/Dense>
#include <vector>
#include <random>
#include <iostream>
#include <iomanip>
#include <cmath>

using Vec = Eigen::VectorXd;
using Mat = Eigen::MatrixXd;

// Plain CG solving (A + sigma I) x = b. Same scalar convention as matpoly_claude.h
// solve(): al = mu/gam, bet = mu_new/mu. Returns x; reports iteration count (= matvecs).
static Vec cg_single(const Mat& A, double sigma, const Vec& b, double tol, int& iters, int maxit=100000){
  const int n = b.size();
  Vec x = Vec::Zero(n);
  Vec r = b;                 // r = b - (A+sigma) x, x = 0
  Vec p = r;
  double mu = r.squaredNorm();
  const double mu_crit = tol*tol*b.squaredNorm();
  int k=0;
  for(; k<maxit; ++k){
    if(mu < mu_crit) break;
    Vec q = A*p + sigma*p;   // (A + sigma) p
    const double gam = p.dot(q);
    const double al  = mu/gam;
    x += al*p;
    r -= al*q;
    const double mu_new = r.squaredNorm();
    const double bet = mu_new/mu;
    p = r + bet*p;
    mu = mu_new;
  }
  iters = k;
  return x;
}

// Multi-shift CG: seed = smallest shift sigma_s; ONE Krylov pass (one matvec/iter)
// solves (A + sigma_m) x_m = b for all m. Relative shifts \hat\sigma_m = sigma_m - sigma_s >= 0.
//   zeta_m^{j+1} = (zeta_m^j zeta_m^{j-1} al_{j-1}) /
//                  (al_j bet_{j-1}(zeta_m^{j-1}-zeta_m^j) + zeta_m^{j-1} al_{j-1}(1 + \hat\sigma_m al_j))
//   al_m = al_j  zeta_m^{j+1}/zeta_m^j ;  x_m += al_m p_m
//   bet_m = bet_j (zeta_m^{j+1}/zeta_m^j)^2 ;  p_m = zeta_m^{j+1} r_{j+1} + bet_m p_m
// (the seed m=s has \hat\sigma=0, so the formula yields zeta=1 identically -- no special case).
static std::vector<Vec> cg_multishift(const Mat& A, const std::vector<double>& sigma, const Vec& b,
                                      double tol, int& iters, int maxit=100000){
  const int n = b.size();
  const int npole = (int)sigma.size();
  int seed = 0; for(int m=1;m<npole;m++) if(sigma[m] < sigma[seed]) seed = m;
  const double sig0 = sigma[seed];
  std::vector<double> hat(npole);
  for(int m=0;m<npole;m++) hat[m] = sigma[m] - sig0;   // >= 0

  // seed CG state on M0 = A + sig0
  Vec r  = b;          // seed residual r_j
  Vec p0 = r;          // seed direction
  double mu = r.squaredNorm();
  const double mu_crit = tol*tol*b.squaredNorm();

  // per-shift state
  std::vector<Vec> x(npole, Vec::Zero(n));
  std::vector<Vec> p(npole, b);                 // p_m^0 = r_0 = b
  std::vector<double> zeta(npole, 1.0), zeta_old(npole, 1.0);
  double al_old = 1.0, bet_old = 0.0;           // al_{-1}=1, bet_{-1}=0

  int k=0;
  for(; k<maxit; ++k){
    if(mu < mu_crit) break;
    // seed alpha
    Vec q = A*p0 + sig0*p0;
    const double gam = p0.dot(q);
    const double al  = mu/gam;

    // x_m update (uses zeta^{j+1}, computed here)
    std::vector<double> zeta_new(npole);
    for(int m=0;m<npole;m++){
      const double denom = al*bet_old*(zeta_old[m]-zeta[m])
                         + zeta_old[m]*al_old*(1.0 + hat[m]*al);
      zeta_new[m] = (zeta[m]*zeta_old[m]*al_old)/denom;
      const double al_m = al * zeta_new[m]/zeta[m];
      x[m] += al_m * p[m];
    }

    // seed residual + beta
    r -= al*q;
    const double mu_new = r.squaredNorm();
    const double bet = mu_new/mu;

    // p_m update
    for(int m=0;m<npole;m++){
      const double ratio = zeta_new[m]/zeta[m];
      const double bet_m = bet * ratio * ratio;
      p[m] = zeta_new[m]*r + bet_m*p[m];
    }
    // seed direction + history roll
    p0 = r + bet*p0;
    zeta_old = zeta; zeta = zeta_new;
    al_old = al; bet_old = bet;
    mu = mu_new;
  }
  iters = k;
  return x;
}

int main(){
  std::cout << std::scientific << std::setprecision(6);
  std::mt19937 gen(12345);
  std::normal_distribution<double> nd(0.0,1.0);
  std::uniform_real_distribution<double> ud(0.0,1.0);

  const int n = 64;
  const int npole = 20;
  const double tol = 1.0e-9;   // mimic Comp::TOL_INNER

  // SPD A = Q diag(d) Q^T with eigenvalues in [dmin,dmax] (mimics normalized D_W^dag D_W in (0,1]).
  const double dmin = 0.02, dmax = 1.0;
  Mat B(n,n); for(int i=0;i<n;i++) for(int j=0;j<n;j++) B(i,j)=nd(gen);
  Eigen::HouseholderQR<Mat> qr(B);
  Mat Q = qr.householderQ();
  Vec d(n); for(int i=0;i<n;i++) d(i) = dmin + (dmax-dmin)*ud(gen);
  Mat A = Q * d.asDiagonal() * Q.transpose();
  A = 0.5*(A + A.transpose());   // symmetrize against round-off

  // positive shifts, geometric over [smin,smax] (representative of -k^2/c'_m)
  const double smin = 1.0e-3, smax = 1.0e2;
  std::vector<double> sigma(npole);
  for(int m=0;m<npole;m++) sigma[m] = smin * std::pow(smax/smin, (double)m/(npole-1));

  Vec b(n); for(int i=0;i<n;i++) b(i) = nd(gen);

  // reference: independent CG per shift
  std::vector<Vec> x_ref(npole);
  int single_total_iters = 0, single_max_iters = 0;
  for(int m=0;m<npole;m++){ int it=0; x_ref[m] = cg_single(A, sigma[m], b, tol, it);
    single_total_iters += it; single_max_iters = std::max(single_max_iters, it); }

  // multi-shift: one pass
  int ms_iters = 0;
  std::vector<Vec> x_ms = cg_multishift(A, sigma, b, tol, ms_iters);

  // compare + true residual ||(A+sigma_m) x_m - b||/||b|| for the multi-shift solutions
  std::cout << "# n="<<n<<"  npole="<<npole<<"  tol="<<tol << std::endl;
  std::cout << "# m   sigma_m        max|x_ms-x_ref|   rel            resid(ms)" << std::endl;
  double worst_diff=0.0, worst_rel=0.0, worst_res=0.0;
  const double bn = b.norm();
  for(int m=0;m<npole;m++){
    const double diff = (x_ms[m]-x_ref[m]).cwiseAbs().maxCoeff();
    const double rel  = diff / (x_ref[m].norm()>0 ? x_ref[m].norm() : 1.0);
    const double res  = ((A*x_ms[m] + sigma[m]*x_ms[m]) - b).norm() / bn;
    worst_diff=std::max(worst_diff,diff); worst_rel=std::max(worst_rel,rel); worst_res=std::max(worst_res,res);
    std::cout << "  "<<std::setw(2)<<m<<"  "<<sigma[m]<<"   "<<diff<<"   "<<rel<<"   "<<res<<std::endl;
  }
  std::cout << "# WORST  diff="<<worst_diff<<"  rel="<<worst_rel<<"  resid(ms)="<<worst_res << std::endl;

  // matvec accounting: independent does sum(iters_m) matvecs; multi-shift does iters(seed) matvecs.
  std::cout << "# matvecs: independent(sum)="<<single_total_iters
            << "  independent(max/pole)="<<single_max_iters
            << "  multishift="<<ms_iters
            << "  => matvec speedup ~"<<(double)single_total_iters/ms_iters<<"x" << std::endl;

  const bool pass = (worst_rel < 1.0e-6) && (worst_res < 1.0e-6);
  std::cout << "# C1 multi-shift recurrence: " << (pass?"PASS":"FAIL") << std::endl;
  return pass ? 0 : 1;
}
