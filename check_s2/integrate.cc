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
// Real epsabs = 1.0e-2; // 0.;
// Real epsrel = 1.0e-2; // TOLLOOSE;
Real epsabs = 1.0e-5; // 0.;
Real epsrel = 1.0e-5; // TOLLOOSE;

int key = 5;



// Real jacobi( const int n, const double alpha, const double beta, const double z, const int nthreshold=97, const int M=10 ){
Real jacobi( const int n, const double alpha, const double beta, const double z, const int nthreshold=200, const int M=10 ){
  return boost::math::jacobi(n, 1.0*alpha, 1.0*beta, z);
  // else return jacobi_asymp(n, alpha, beta, z, M);
  // else return jacobi_asymp_darboux(n, alpha, beta, z);
}



// Real jacobi(unsigned n, Real alpha, Real beta, Real x);
// mpH = 1, 2, 3, ...
// n = 0, 1, 2, ...
Real xi( const int mpH, const int n, const double z){
  assert( mpH >= 1);
  if(std::abs(z-1.0)<1.0e-14) return mpH==1 ? 1.0/std::sqrt(2.0) : 0.0;
  else if( std::abs(z+1.0)<1.0e-14 ) return 0.0;

  // const Real factor = std::pow( 1-z, 0.5*(mpH-1) ) * std::pow( 1+z, -0.5*mpH ) * std::pow(0.5, mpH-1);
  const Real factor = std::pow( 1-z, 0.5*(mpH-1) ) * std::pow( 1+z, -0.5*mpH ); // * std::pow(0.5, mpH-1);
  // const Real poly = boost::math::jacobi(n+mpH, mpH-1.0, -1.0*mpH, z);
  const Real poly = jacobi(n+mpH, mpH-1.0, -1.0*mpH, z);

  return factor * poly;
}



Real Cnm( const int mpH, const int n ){
  assert( mpH >= 1);
  Real tmp = 4.0*M_PI;
  tmp /= 2.0*n + 2.0*mpH;

  if(mpH!=1){
    tmp *= std::tgamma( n+2.0*mpH ) * std::tgamma( n+1.0 );
    tmp /= std::tgamma( n+mpH ) * std::tgamma( n+mpH+1.0 );
  }
  return std::sqrt(tmp);
}



// VC psiP( const int sign_m, const int mpH, const int n, const int s, const double z, const double phi){
//   VC res;
//   const Complex first = xi(mpH, n, z);
//   const Complex secon = I * std::pow(-1.0, n) * xi( mpH, n, -z);

//   if(sign_m>0) res << first, secon;
//   else res << secon, first;

//   res *= std::exp( I*(1.0*sign_m) * (mpH-0.5) * phi );

//   return res / Cnm(mpH, n);
// }


// VC psiM( const int sign_m, const int mpH, const int n, const int s, const double z, const double phi){
//   VC res;
//   const Complex first = xi(mpH, n, z);
//   const Complex secon = I * std::pow(-1.0, n) * xi( mpH, n, -z);

//   if(sign_m>0) res << first, secon;
//   else res << secon, first;

//   res[1] *= -1.0;
//   res *= std::exp( I*sign_m * (mpH-0.5) * phi );

//   return res / Cnm(mpH, n);
// }




const int mpH_max = 200;
const int n_max = mpH_max;


int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  const int nint=10000;
  const Real delta=M_PI/nint;

  // const int i0 = 0, i1=1;

  for(int i=1; i<nint; i++){
    const Real theta = delta*i;
    const Real phi = 0.0;

    const Real z = std::cos(theta);
    Complex wf=0.0;

    {
      const int mpH=1;
      for(int n=0; n<=n_max-mpH; n++){
        // Real z0 = 1.0;
        wf += std::pow(-1.0, n) * xi(1,n,-z) / (std::sqrt(2.0)*M_PI);
      }
      std::cout << theta << " " << std::abs(wf) << std::endl;
    }
  }

  return 0;
}

