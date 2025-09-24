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
}



Real f0( const int n, const double z){
  assert( n >= 1);
  return ( 1.0-z ) * jacobi(n, 1.0, -1.0, z);
}

Real df0( const int n, const double z){
  assert( n >= 1);
  return - jacobi(n, 1.0, -1.0, z) + 0.5 * (1.0+n) * ( 1.0-z ) *jacobi(n-1.0, 2.0, 0.0, z);
}



// Real Cnm( const int mpH, const int n ){
//   assert( mpH >= 1);
//   Real tmp = 4.0*M_PI;
//   tmp /= 2.0*n + 2.0*mpH;

//   if(mpH!=1){
//     tmp *= std::tgamma( n+2.0*mpH ) * std::tgamma( n+1.0 );
//     tmp /= std::tgamma( n+mpH ) * std::tgamma( n+mpH+1.0 );
//   }
//   return std::sqrt(tmp);
// }




const int mpH_max = 20;
const int n_max = mpH_max;


int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  const int nint=1000;
  const Real delta=M_PI/nint;

  // const int i0 = 0, i1=1;

  for(int i=1; i<nint; i++){
    const Real theta = delta*i;
    const Real phi = 0.0;

    const Real z = std::cos(theta);
    Real wf=0.0;

    {
      for(int n=1; n<=n_max; n++){
        wf -= (2.0*n+1.0)/(n+1.0) * df0(n,z) / (4.0*M_PI);
      }
      std::cout << z << " " << wf << std::endl;
    }
  }

  return 0;
}

