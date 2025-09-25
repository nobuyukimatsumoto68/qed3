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
constexpr Complex I = Complex(0.0, 1.0);


Real jacobi( const int n, const double alpha, const double beta, const double z){
  return boost::math::jacobi(n, 1.0*alpha, 1.0*beta, z);
}



// Real jacobi(unsigned n, Real alpha, Real beta, Real x);
// mpH = 1, 2, 3, ...
// n = 0, 1, 2, ...
Real xi( const int mpH, const int n, const double z){
  assert( mpH >= 1);
  if(std::abs(z-1.0)<1.0e-14) return mpH==1 ? 1.0/std::sqrt(2.0) : 0.0;
  else if( std::abs(z+1.0)<1.0e-14 ) return 0.0;

  const Real factor = std::pow( 1-z, 0.5*(mpH-1) ) * std::pow( 1+z, -0.5*mpH );
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


int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  int n_max = 200;
  if(argc==2) n_max=atoi(argv[1]);
  std::cout << "# n_max = " << n_max << std::endl;

  const int nint=10000;
  const Real delta=M_PI/nint;

  for(int i=1; i<nint; i++){
    const Real theta = delta*i;
    const Real z = std::cos(theta);

    {
      Complex wf=0.0;
      for(int n=0; n<=n_max; n++){
        wf += std::pow(-1.0, n) * xi(1,n,-z) / (std::sqrt(2.0)*M_PI);
      }
      std::cout << theta << " " << real(wf) << " " << imag(wf) << std::endl;
    }
  }

  return 0;
}

