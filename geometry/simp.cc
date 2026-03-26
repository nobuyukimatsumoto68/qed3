#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <map>
#include <limits>
#include <stdfloat>

#include <Eigen/Dense>

#include <gsl/gsl_integration.h>

using Double = std::float64_t;
//#include "geodesic.h"
// #include "integral.h"

#include "s2.h"


// using namespace Geodesic;
Double Mod(Double a, Double b=2.0*M_PI){
  int p = int(std::floor(a / b));
  Double r = a - p*b;
  return r;
}


using Idx = std::int32_t;

using Link = std::array<Idx,2>; // <Idx,Idx>;
using Face = std::vector<Idx>;
using VD=Eigen::Matrix<Double, 2, 1>;
using VE=Eigen::Matrix<Double, 3, 1>;


using I2=Eigen::Vector2i;

constexpr int DIM = 2;
using V2=Eigen::Matrix<Double, 1, 2>;

constexpr int EDIM = 3;
using V3=Eigen::Matrix<Double, 1, 3>;



Double TOL=1.0e-13;

// constexpr Double TOLLOOSE=1.0e-6;
// constexpr Double TOLLOOSE=1.0e-12;

// const int limit = 4000; // 1000;
// Double epsabs = 1.0e-13; // 0.;
// Double epsrel = 1.0e-13; // TOLLOOSE;
// int key = 5;

const int limit = 10000; // 1000;
double epsabs = 1.0e-15; // 0.;
double epsrel = 1.0e-13; // TOLLOOSE;
int key = 6;


using SpinMatrix = Eigen::Matrix2cd;


using Func=std::function<double(const double)>;
double unwrap(double x, void *p) {
  auto fp = static_cast<Func*>(p);
  return (*fp)(x);
}



namespace Sphere{
  //  using namespace Sphere;
V3 ptilde(const V3& r1, const V3& r2, const double stilde) {
  assert(0.<=stilde && stilde<=1.);
  assert( std::abs(r1.norm()-1.0)<1.0e-14 );
  assert( std::abs(r2.norm()-1.0)<1.0e-14 );
  return r1 + stilde*(r2-r1);
}

V3 p(const V3& r1, const V3& r2, const double stilde) {
  const V3 tmp = ptilde(r1,r2,stilde);
  return tmp/tmp.norm();
}

V3 etilde(const V3& r1, const V3& r2, const double stilde) {
  const V3 dp = r2-r1;
  const V3 p_ = p(r1,r2,stilde);
  const V3 ptilde_ = ptilde(r1,r2,stilde);
  const double tmp = p_.dot(dp);
  return (dp - tmp*p_)/ptilde_.norm();
}

V3 dtheta_dp(const V3& r) {
  const double sqrt_xsq_ysq = std::sqrt( r(0)*r(0)+r(1)*r(1) );
  const double normsq = r.squaredNorm();

  const double dtheta_dx = r(0)*r(2)/ normsq / sqrt_xsq_ysq;
  const double dtheta_dy = r(1)*r(2)/ normsq / sqrt_xsq_ysq;
  const double dtheta_dz = -sqrt_xsq_ysq / normsq;
  V3 res;
  res << dtheta_dx, dtheta_dy, dtheta_dz;
  return res;
}

V3 dphi_dp(const V3& r) {
  const double dphi_dx = -r(1)/(r(0)*r(0)+r(1)*r(1));
  const double dphi_dy = r(0)/(r(0)*r(0)+r(1)*r(1));
  V3 res;
  res << dphi_dx, dphi_dy, 0.0;
  return res;
}

V3 e(const V3& r1, const V3& r2, const double stilde) {
  const V3 tmp = etilde(r1,r2,stilde);
  return tmp/tmp.norm();
}


inline double theta(const V3& r) {
  assert( std::abs(r.norm()-1.0)<1.0e-14 );
  const double res = std::acos(r[2]);
  assert( 0<=res && res<=M_PI);
  return res;
}

inline double phi(const V3& r) {
  assert( std::abs(r.norm()-1.0)<1.0e-14 );
  const double res = std::atan2(r[1], r[0]);
  assert( -M_PI<=res && res<=M_PI);
  return res;
}


double dtheta_dstilde(const V3& r1, const V3& r2, const double stilde) {
  const V3 p_ = p(r1,r2,stilde);
  const V3 dp_dstilde_ = etilde(r1,r2,stilde);
  const V3 dtheta_dp_ = dtheta_dp( p_ );
  return dp_dstilde_.dot(dtheta_dp_);
}

double dphi_dstilde(const V3& r1, const V3& r2, const double stilde) {
  const V3 p_ = p(r1,r2,stilde);
  const V3 dp_dstilde_ = etilde(r1,r2,stilde);
  const V3 dphi_dp_ = dphi_dp( p_ );
  return dp_dstilde_.dot(dphi_dp_);
}

  V2 e2(const V3& r1, const V3& r2, const double stilde) {
  const double xx = dtheta_dstilde(r1,r2,stilde);
  const double yy = std::sin(theta(p(r1,r2,stilde))) * dphi_dstilde(r1,r2,stilde);
  V2 res;
  res << xx, yy;
  return res;
}


inline double alpha(const V3& r1, const V3& r2, const double stilde) {
  const V2 v = e2(r1,r2,stilde);
  return std::atan2( v[1], v[0] );
}
};






std::string dir = "./data/";

int main(int argc, char* argv[]){
  int n_refine=1;
  if(argc==2) n_refine = atoi(argv[1]);

  QfeLatticeS2 lattice(5, n_refine);


  Double tht = 1.0/n_refine;
  Eigen::Matrix<Double, 3, 3> rot;
  rot << std::cos(tht), 0.0, std::sin(tht),
    0.0, 1.0, 0.0,
    -std::sin(tht), 0.0, std::cos(tht);

  std::vector<VE> sites;
  for(Idx ix=0; ix<lattice.n_sites; ix++) {
    VE x = rot*lattice.r[ix];
    sites.push_back(x);
  }




  std::vector<Double> omegas;
  std::vector<Double> alpha0s, alpha1s;
  // std::vector<Sol> sols;
  int counter=0;
  for(const auto& link : lattice.links){
    counter++;
    Idx ix1 = link.sites[0];
    Idx ix2 = link.sites[1];

    V3 x1( sites[ix1] );
    V3 x2( sites[ix2] );

    // Sol sol = SolveGeodesics( x1, x2 );
    // gsl_integration_workspace * w = gsl_integration_workspace_alloc (limit);
    // double result, error;

    double omega12;
    {
      gsl_integration_workspace * w = gsl_integration_workspace_alloc (limit);
      double result, error;

      Func f1 = [&](const double stilde){
        const V3 p_ = Sphere::p(x1, x2, stilde);
        const double dphi_dstilde_ = Sphere::dphi_dstilde(x1, x2, stilde);
        const double theta_ = Sphere::theta( p_ );

        return -dphi_dstilde_ * std::cos( theta_ );
      };
      gsl_function F;
      F.function = &unwrap;
      F.params = &f1;
      gsl_integration_qag(&F, 0., 1., epsabs, epsrel, limit,
                          key, w, &result, &error);

      omega12 = -result;
      gsl_integration_workspace_free (w);
    }

    //{
    const double alpha_1 = Sphere::alpha( x1,x2,0.0);
    const double alpha_2 = Sphere::alpha( x1,x2,1.0) + M_PI;

    // double diff = alpha_2 - omega12 - alpha_1;
    double diff = alpha_2 + M_PI + omega12 - alpha_1;
    while(diff>M_PI) diff -= 2.0*M_PI;
    while(diff<=-M_PI) diff += 2.0*M_PI;
    if(std::abs(diff)>TOL){
      std::cout << "debug. diff = " << std::abs(diff) << std::endl;
      std::cout << "debug. x1 = " << x1.transpose() << std::endl;
      std::cout << "debug. x2 = " << x2.transpose() << std::endl;
      assert( std::abs(diff)<TOL );
    }
    // }

    const double phi1 = Sphere::phi(x1);
    const double phi2 = Sphere::phi(x2);
    //
    const double phi0 = M_PI/48.0;

    double dphi1 = phi1-phi0;
    double dphi2 = phi2-phi0;
    double dphi = phi1-phi2;

    if(dphi1 >= M_PI) dphi1 -= 2.0*M_PI;
    else if(dphi1 < -M_PI) dphi1 += 2.0*M_PI;
    if(dphi2 >= M_PI) dphi2 -= 2.0*M_PI;
    else if(dphi2 < -M_PI) dphi2 += 2.0*M_PI;
    if(dphi >= M_PI) dphi -= 2.0*M_PI;
    else if(dphi < -M_PI) dphi += 2.0*M_PI;
    if( std::abs(dphi1) + std::abs(dphi2) - std::abs(dphi) < 1.0e-14 ) omega12 -= 2.0*M_PI;

    // omegas[dual.directedlinkidx(if1, df)] = omega12;
    omegas.push_back( omega12 );

    alpha0s.push_back( alpha_1 );
    // alpha1s.push_back( Mod(alpha_2 + M_PI) );
    alpha1s.push_back( alpha_2 );


    // if(!sol.is_split){
    //   F f1 = [&](Double s){ return sol.Dphi(s) * std::cos( sol.theta(s) ); };
    //   gsl_function F;
    //   F.function = &unwrap;
    //   F.params = &f1;
    //   gsl_integration_qag(&F, 0., sol.ell, epsabs, epsrel, limit,
    //                       key, w, &result, &error);
    // }
    // else{
    //   double result1, result2;
    //   double error1, error2;
    //   {
    //     F f1 = [&](Double s){ return sol.Dphi(s, 0) * std::cos( sol.theta(s, 0) ); };
    //     gsl_function F;
    //     F.function = &unwrap;
    //     F.params = &f1;
    //     gsl_integration_qag(&F, 0.0, sol.sE, epsabs, epsrel, limit,
    //                         key, w, &result1, &error1);
    //   }
    //   {
    //     F f1 = [&](Double s){ return sol.Dphi(s, 1) * std::cos( sol.theta(s, 1) ); };
    //     gsl_function F;
    //     F.function = &unwrap;
    //     F.params = &f1;
    //     gsl_integration_qag(&F, sol.sE, sol.ell, epsabs, epsrel, limit,
    //                         key, w, &result2, &error2);
    //   }
    //   result = result1+result2;
    //   error=error1+error2;
    // }
    // omegas.push_back( result );
    // sols.push_back( sol );
  }

  // {
  //   const Double phi0 = M_PI/48.0;
  //   for(Idx ix=0; ix<lattice.n_sites; ix++) assert( std::abs(Mod(Pt(sites[ix])i[1]) - phi0)>TOLLOOSE );
  //   std::vector<Idx> section;
  //   for(Idx il=0; il<lattice.n_links; il++){
  //     const auto& link = lattice.links[il];
  //     Idx ix1 = link.sites[0];
  //     Idx ix2 = link.sites[1];

  //     const Double phi1 = Pt(sites[ix1])i[1];
  //     const Double phi2 = Pt(sites[ix2])i[1];
  //     Double dphi1 = phi1-phi0;
  //     Double dphi2 = phi2-phi0;
  //     Double dphi = phi1-phi2;

  //     if(dphi1 >= M_PI) dphi1 -= 2.0*M_PI;
  //     else if(dphi1 < -M_PI) dphi1 += 2.0*M_PI;
  //     if(dphi2 >= M_PI) dphi2 -= 2.0*M_PI;
  //     else if(dphi2 < -M_PI) dphi2 += 2.0*M_PI;
  //     if(dphi >= M_PI) dphi -= 2.0*M_PI;
  //     else if(dphi < -M_PI) dphi += 2.0*M_PI;

  //     if( std::abs(dphi1) + std::abs(dphi2) - std::abs(dphi) < TOLLOOSE ){
  //       section.push_back( il );
  //     }
  //   }
  //   for(Idx il : section) omegas[il] -= 2.0*M_PI;
  // }

  std::map<const Link, Double> omega;
  for(Idx il=0; il<lattice.n_links; il++){
    const auto& link = lattice.links[il];
    Idx ix1 = link.sites[0];
    Idx ix2 = link.sites[1];

    omega.insert( { Link{ix1, ix2},  omegas[il] } );
    omega.insert( { Link{ix2, ix1}, -omegas[il] } );
  }


  // std::vector<Double> alpha0s, alpha1s;
  // for( auto& sol : sols){
  //   Double alpha0, alpha1;
  //   if(!sol.is_split){
  //     const Double theta0 = sol.theta(0);
  //     const Double dtheta_ds0 = sol.Dtheta(0);
  //     const Double dphi_ds0 = sol.Dphi(0);
  //     assert( std::abs(dtheta_ds0*dtheta_ds0 + std::sin(theta0)*std::sin(theta0)*dphi_ds0*dphi_ds0-1.0)<TOLLOOSE );

  //     VD e0( dtheta_ds0, std::sin(theta0)*dphi_ds0 );
  //     alpha0 = _acos( e0[0] );
  //     if(e0[1]<0) alpha0 *= -1.0;
  //     assert( std::abs(std::cos(alpha0)-e0[0]) + std::abs(std::sin(alpha0)-e0[1]) < TOLLOOSE );

  //     const Double theta1 = sol.theta(sol.ell);
  //     const Double dtheta_ds1 = sol.Dtheta(sol.ell);
  //     const Double dphi_ds1 = sol.Dphi(sol.ell);
  //     assert( std::abs(dtheta_ds1*dtheta_ds1 + std::sin(theta1)*std::sin(theta1)*dphi_ds1*dphi_ds1-1.0)<TOLLOOSE );

  //     VD e1( dtheta_ds1, std::sin(theta1)*dphi_ds1 );
  //     alpha1 = _acos( e1[0] );
  //     if(e1[1]<0) alpha1 *= -1.0;
  //     assert( std::abs(std::cos(alpha1)-e1[0]) + std::abs(std::sin(alpha1)-e1[1]) < TOLLOOSE );
  //   }
  //   else{
  //     Double theta0, dtheta_ds0, dphi_ds0;
  //     if( !isModdable(sol.sE) ){
  //       theta0 = sol.theta(0, 0);
  //       dtheta_ds0 = sol.Dtheta(0, 0);
  //       dphi_ds0 = sol.Dphi(0, 0);
  //       // std::cout << "debug. " <<
  //       //   dtheta_ds0 << " " << theta0 << " " << dphi_ds0 << std::endl <<
  //       //   std::abs(dtheta_ds0*dtheta_ds0 + std::sin(theta0)*std::sin(theta0)*dphi_ds0*dphi_ds0-1.0)
  //       //           << std::endl;
  //       assert( std::abs(dtheta_ds0*dtheta_ds0 + std::sin(theta0)*std::sin(theta0)*dphi_ds0*dphi_ds0-1.0)<TOLLOOSE );
  //     }
  //     else{
  //       theta0 = sol.theta(sol.sE, 1);
  //       dtheta_ds0 = sol.Dtheta(sol.sE, 1);
  //       dphi_ds0 = sol.Dphi(sol.sE, 1);
  //       assert( std::abs(dtheta_ds0*dtheta_ds0 + std::sin(theta0)*std::sin(theta0)*dphi_ds0*dphi_ds0-1.0)<TOLLOOSE );
  //     }

  //     VD e0( dtheta_ds0, std::sin(theta0)*dphi_ds0 );
  //     alpha0 = _acos( e0[0] );
  //     if(e0[1]<0) alpha0 *= -1.0;
  //     assert( std::abs(std::cos(alpha0)-e0[0]) + std::abs(std::sin(alpha0)-e0[1]) < TOLLOOSE );

  //     const Double theta1 = sol.theta(sol.ell, 1);
  //     const Double dtheta_ds1 = sol.Dtheta(sol.ell, 1);
  //     const Double dphi_ds1 = sol.Dphi(sol.ell, 1);
  //     assert( std::abs(dtheta_ds1*dtheta_ds1 + std::sin(theta1)*std::sin(theta1)*dphi_ds1*dphi_ds1-1.0)<TOLLOOSE );

  //     VD e1( dtheta_ds1, std::sin(theta1)*dphi_ds1 );
  //     alpha1 = _acos( e1[0] );
  //     if(std::abs(e1[1]) > TOLLOOSE && e1[1]<0) alpha1 *= -1.0;
  //     // std::cout << "debug. " << std::abs(std::cos(alpha1)-e1[0]) + std::abs(std::sin(alpha1)-e1[1])
  //     //           << std::endl;
  //     assert( std::abs(std::cos(alpha1)-e1[0]) + std::abs(std::sin(alpha1)-e1[1]) < TOLLOOSE );
  //   }
  //   alpha0s.push_back( alpha0 );
  //   alpha1s.push_back( Mod(alpha1 + M_PI) );
  // }

  std::map<const Link, Double> alpha;
  for(Idx il=0; il<lattice.n_links; il++){
    const auto& link = lattice.links[il];
    Idx ix1 = link.sites[0];
    Idx ix2 = link.sites[1];

    alpha.insert( { Link{ix1, ix2}, alpha0s[il] } );
    alpha.insert( { Link{ix2, ix1}, alpha1s[il] } );
  }

  {
    std::clog << "# checking spin structure" << std::endl;
    for(Idx ix=0; ix<sites.size(); ix++){
      const auto x = lattice.sites[ix];
      for(int iw=0; iw<x.nn; iw++){
        Idx iy = x.neighbors[iw];
        const Double alpha1 = alpha.at(Link{ix,iy});
        Double alpha2 = alpha.at(Link{iy,ix});
        Double omega12 = omega.at(Link{ix,iy});

        Double diff = (alpha2 + M_PI + omega12) - alpha1;
        // assert( isModdable(diff, TOL) );
        if(diff>M_PI) diff -= 2.0*M_PI;
        else if(diff<=-M_PI) diff += 2.0*M_PI;
        if(std::abs(diff)>TOL){
          std::cout << "debug. diff = " << std::abs(diff) << std::endl;
          assert( std::abs(diff)<TOL );
        }


        // Double om = alpha1 - (alpha2 + M_PI);
        // const int br = decide_branch( om-omega12 );
        // om -= M_PI*br;
        // omega[Link({ix, iy})] = om;
      }}
  }


  {
    std::clog << "# checking deficits" << std::endl;
    int counter=0;
    for(int ia=0; ia<lattice.n_faces; ia++){
      int sign = 1;
      {
        VE x0 = sites[ lattice.faces[ia].sites[0] ];
        VE x1 = sites[ lattice.faces[ia].sites[1] ];
        VE x2 = sites[ lattice.faces[ia].sites[2] ];
        VE sum = x0+x1+x2;
        if((x1-x0).cross(x2-x0).dot(sum) < 0) sign = -1;
      }

      Double sum = 0.0;
      for(int i=0; i<3; i++){
        Idx ix = lattice.faces[ia].sites[i];
        Idx jx = lattice.faces[ia].sites[(i+1)%3];

        sum += omega.at( Link{ix,jx} );
      }
      sum *= sign;
      // std::cout << "sum = " << sum << std::endl;
      Double mod = Mod(sum, 4.0*M_PI);
      // std::cout << "sum (mod4pi) = " << mod << std::endl;
      if(mod>2.0*M_PI) mod -= 4.0*M_PI;
      // std::clog << "# sum (mod4pi, repr) = " << mod << std::endl;
      assert( (-1.5 * 4.0*M_PI/lattice.n_faces < mod && mod < 0.0) );
      counter++;
    }
  }



  {
    std::ofstream ofs(dir+"omega_n"+std::to_string(n_refine)+".dat");
    ofs << std::scientific << std::setprecision(15);
    for(Idx il=0; il<lattice.n_links; il++){
      const auto& link = lattice.links[il];
      Idx ix1 = link.sites[0];
      Idx ix2 = link.sites[1];

      ofs << std::setw(25) << ix1 << " ";
      ofs << std::setw(25) << ix2 << " ";
      ofs << std::setw(25) << omegas[il] << std::endl;
    }
  }

  {
    std::ofstream ofs(dir+"alpha_n"+std::to_string(n_refine)+".dat");
    ofs << std::scientific << std::setprecision(15);
    for(Idx il=0; il<lattice.n_links; il++){
      const auto& link = lattice.links[il];
      Idx ix1 = link.sites[0];
      Idx ix2 = link.sites[1];

      ofs << std::setw(25) << ix1 << " ";
      ofs << std::setw(25) << ix2 << " ";
      ofs << std::setw(25) << alpha0s[il] << std::endl;

      ofs << std::setw(25) << ix2 << " ";
      ofs << std::setw(25) << ix1 << " ";
      ofs << std::setw(25) << alpha1s[il] << std::endl;
    }
  }

  return 0;
}

