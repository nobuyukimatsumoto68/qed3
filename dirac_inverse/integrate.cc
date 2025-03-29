#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <map>
#include <limits>
#include <functional>
#include <stdfloat>

#include <vector>
#include <numbers>



#include <Eigen/Dense>
#include "boost/math/special_functions/jacobi.hpp"

using Real = double; // std::float64_t;
using Complex = std::complex<Real>;
using F=std::function<Real(const Real)>;
using VC=Eigen::Vector2cd;

constexpr Complex I = Complex(0.0, 1.0);

#include "integral.h"


const int limit = 10000; // 1000;

// Real epsabs = 1.0e-15; // 0.;
// Real epsrel = 1.0e-13; // TOLLOOSE;
Real epsabs = 1.0e-2; // 0.;
Real epsrel = 1.0e-2; // TOLLOOSE;

int key = 5;

// Real jacobi(unsigned n, Real alpha, Real beta, Real x);
// mpH = 1, 2, 3, ...
// n = 0, 1, 2, ...
Real xi( const int mpH, const int n, const int sm, const int s3, const double z){
  const Real factor = std::pow( 1-z, -0.5*(mpH-0.5-0.5*s3) ) * std::pow( 1+z, -0.5*(mpH-0.5+0.5*s3) ) * std::pow( 0.5*(1.0-s3*z), mpH-1 );
  const Real poly = boost::math::jacobi(n+mpH, mpH-1.0, -1.0*mpH, s3*z);

  return factor * poly;
}



Real xiPP( const int mpH, const int n, const double z){
  Real res = 1.0;
  // if(std::abs(z+1.0)<8.0e-3){
  //   res *= std::pow( 2.0, -mpH-1 ) * std::pow( 1.0-z, 0.5*(mpH-1) );
  //   res *= std::pow( 1.0+z, 0.5*mpH );
  //   res *= std::pow( -1.0, n ) * std::pow( 2.0, -mpH ) * tgamma(n+2*mpH-1)/tgamma(mpH)/tgamma(n+mpH-1);
  //   return res;
  // }
  // else{
  const Real factor = std::pow( 1-z, -0.5*(mpH-1) ) * std::pow( 1+z, -0.5*mpH ) * std::pow( 0.5*(1.0-z), mpH-1 );
  const Real poly = boost::math::jacobi(n+mpH, mpH-1.0, -1.0*mpH, z);
  return factor * poly;
  // }
}

Real Cnm( const int mpH, const int n ){
  Real tmp = 4.0*M_PI;
  tmp /= 2.0*n + 2.0*mpH;
  tmp *= std::tgamma( n+2.0*mpH ) * std::tgamma( n+1.0 );
  tmp /= std::tgamma( n+mpH ) * std::tgamma( n+mpH+1.0 );
  return std::sqrt(tmp);
}



Complex eta( const int mpH, const int n, const int sm, const int s3, const double z, const double phi){
  const Complex phase = phi * sm * (mpH-0.5);
  const Complex pf = std::exp( I * phase );
  const Complex Xi = xi(mpH, n, sm, s3, z);
  return pf * Xi;
}


VC psi( const int mpH, const int n, const int s, const double z, const double phi){
  const Complex first = eta(mpH,n, 1, 1,z,phi) + std::pow(-1,n)*s*I * eta(mpH,n,-1,-1,z,phi);
  const Complex secon = eta(mpH,n,-1, 1,z,phi) + std::pow(-1,n)*s*I * eta(mpH,n, 1,-1,z,phi);
  VC res;
  res << first, secon;
  return res / Cnm(mpH, n);
}

VC psi2( const int mpH, const int n, const int s, const double z, const double phi){
  const Complex first = eta(mpH,n, 1, 1,z,phi) - std::pow(-1,n)*s*I * eta(mpH,n,-1,-1,z,phi);
  const Complex secon = -eta(mpH,n,-1, 1,z,phi) + std::pow(-1,n)*s*I * eta(mpH,n, 1,-1,z,phi);
  VC res;
  res << first, secon;
  return res / Cnm(mpH, n);
}



const int mpH_max = 40;
const int n_max = mpH_max;


int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  // const double z = std::cos(0.2), phi = 0.3;
  // std::vector<int> mpHs;
  // std::vector<int> ns;
  // std::vector<Real> Zs;

  // for(int mpH=1; mpH<=mpH_max; mpH++){
  //   for(int n=0; n<=n_max-mpH; n++){
  //     // if(mpH>=8 && n!=0 ) continue;

  //     // const int mpH=2, n=0;
  //     // const int sm=1, s3=1;
  //     // std::cout << "eta = " << eta(mpH, n, sm, s3, z, phi) << std::endl;

  //     // std::cout << "psi+ = " << std::endl
  //     //           << psi(mpH, n, +1, z, phi) << std::endl
  //     //           << "psi- = " << std::endl
  //     //           << psi(mpH, n, -1, z, phi) << std::endl;
  //     // gsl_integration_workspace * w = gsl_integration_workspace_alloc (limit);
  //     // double result, error;

  //     // F fp = [&](Real z){
  //     //   const Real psiPP= xiPP(mpH, n, z);
  //     //   // const Real psiMM= xi(mpH, n,-1,-1, z);
  //     //   // const Real psiPM= xi(mpH, n, 1,-1, z);
  //     //   // const Real psiMP= xi(mpH, n,-1, 1, z);
  //     //   return psiPP*psiPP;
  //     // };

  //     // gsl_function F;
  //     // F.function = &unwrap;
  //     // F.params = &fp;

  //     // double pts[2];
  //     // pts[0] = -1.0;
  //     // pts[1] =  1.0;
  //     // // pts[0] = -0.99;
  //     // // pts[1] =  0.99;
  //     // // gsl_integration_qagp(&F, pts, 2, epsabs, epsrel, limit,
  //     // //                      w, &result, &error);
  //     // // gsl_integration_qags(&F, pts[0], pts[1], epsabs, epsrel, limit,
  //     // //                      w, &result, &error);
  //     // gsl_integration_qag(&F, pts[0], pts[1], epsabs, epsrel, limit,
  //     //                     key, w, &result, &error);
  //     // result *= 8.0*M_PI;
  //     // error *= 8.0*M_PI;

  //     // gsl_integration_workspace_free (w);

  //     mpHs.push_back(mpH);
  //     ns.push_back(n);
  //     Zs.push_back(result);

  //     // std::cout << mpH-0.5 << " "
  //     //           << n << " "
  //     //           << result << " "
  //     //           << error << " "
  //     //           << std::endl;
  //   }}

  const int nint=100;
  const Real delta=M_PI/nint;

  const int i0 = 0, i1=1;

  for(int i=1; i<nint; i++){
    const Real theta = delta*i;
    const Real phi = 0.0;

    const Real z = std::cos(theta);
    Complex wf=0.0;

    // Complex pfn = std::exp( I*3.0*M_PI/2.0 );

    for(int mpH=1; mpH<=mpH_max; mpH++){
      for(int n=0; n<=n_max-mpH; n++){
    // for(int j=0; j<mpHs.size(); j++){
    //   const int mpH = mpHs[j];
    //   const int n = ns[j];
        if(n>10 && mpH != 1) continue;

      // VC spinor = psi(mpH, n, 1, z, phi) - psi(mpH, n, -1, z, phi);
      // spinor += pfn * (psi2(mpH, n, 1, z, phi) - psi2(mpH, n, -1, z, phi));
      // VC spinor0 = psi(mpH, n, 1, 0.99, phi) - psi(mpH, n, -1, 0.99, phi);
      // spinor0 += pfn * (psi2(mpH, n, 1, 0.99, phi) - psi2(mpH, n, -1, 0.99, phi));
      // wf += -I/(1.0*n+mpH) / Zs[j] * spinor[0] * std::conj(spinor0[i1]);

      VC spinor = psi(mpH, n, 1, z, phi);
      VC spinor0 = psi(mpH, n, 1, 0.99, phi);
      Complex test = spinor[i0] * std::conj(spinor0[i1]);
      // std::cout << "debug. test = " << test << std::endl;
      wf += -I/(1.0*n+mpH) * spinor[i0] * std::conj(spinor0[i1]);

      spinor = psi(mpH, n, -1, z, phi);
      spinor0 = psi(mpH, n, -1, 0.99, phi);
      test = spinor[i0] * std::conj(spinor0[i1]);
      // std::cout << "debug. test = " << test << std::endl;
      wf += I/(1.0*n+mpH) * spinor[i0] * std::conj(spinor0[i1]);

      spinor = psi2(mpH, n, 1, z, phi);
      spinor0 = psi2(mpH, n, 1, 0.99, phi);
      test = spinor[i0] * std::conj(spinor0[i1]);
      // std::cout << "debug. test = " << test << std::endl;
      wf += -I/(1.0*n+mpH) * spinor[i0] * std::conj(spinor0[i1]);

      spinor = psi2(mpH, n, -1, z, phi);
      spinor0 = psi2(mpH, n, -1, 0.99, phi);
      test = spinor[i0] * std::conj(spinor0[i1]);
      // std::cout << "debug. test = " << test << std::endl;
      wf += I/(1.0*n+mpH) * spinor[i0] * std::conj(spinor0[i1]);

    }
    std::cout << theta << " " << std::abs(wf) << " " << std::sin(theta) << std::endl;
    }}

  return 0;
}

