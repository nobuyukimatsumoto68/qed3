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


Real jacobi( const int n, const double alpha, const double beta, const double z ){
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




int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  int n_max = 20;
  if(argc==2) n_max=atoi(argv[1]);
  std::cout << "# n_max = " << n_max << std::endl;

  const int nint=1000;
  const Real delta=M_PI/nint;

  for(int i=1; i<nint; i++){
    const Real theta = delta*i;

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

