#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <map>
#include <limits>
#include <stdfloat>

#include <Eigen/Dense>

// #include "geodesic.h"
// #include "integral.h"


#include "integral.h"


using Real = double;
using Complex = std::complex<Real>;

// using Idx = std::int32_t;
// using F=std::function<Real(const Real)>;
// using Complex = std::complex<double>;
// const Complex I = Complex(0.0, 1.0);
constexpr Complex I = Complex(0.0, 1.0);





Real TOL=1.0e-14;

constexpr Real TOLLOOSE=1.0e-6;

const int limit = 4000; // 1000;
// Real epsabs = 1.0e-15; // 0.;
// Real epsrel = 1.0e-13; // TOLLOOSE;
Real epsabs = 1.0e-7; // 0.;
Real epsrel = 1.0e-7; // TOLLOOSE;
int key = 5;


// template<typename T>
// T w( T x, double alpha, double beta ){
//   return std::pow( 1.0-x, alpha ) * std::pow( 1.0+x, beta );
// }


double integrand (double theta, void * params) {
  std::vector<double> par = *(std::vector<double> *) params;
  double n = par[0];
  double alpha = par[1];
  double beta = par[2];
  double x = par[3];

  double eps = 0.5*(1.0-std::abs(x));
  Complex z = x + eps*std::exp( I*theta );

  Complex f = std::pow( 1.0-z, n+alpha ) * std::pow( 1.0+z, n+beta ) * std::exp( -n*I*theta );

  return f.real();
}


double integrand_cos (double theta, void * params) {
  std::vector<double> par = *(std::vector<double> *) params;
  double n = par[0];
  double alpha = par[1];
  double beta = par[2];
  double x = par[3];

  double eps = 0.5*(1.0-std::abs(x));
  Complex z = x + eps*std::exp( I*theta );

  Complex f = std::pow( 1.0-z, n+alpha ) * std::pow( 1.0+z, n+beta ) * std::exp( -n*I*theta );

  return f.real();
}

double integrand_sin (double theta, void * params) {
  std::vector<double> par = *(std::vector<double> *) params;
  double n = par[0];
  double alpha = par[1];
  double beta = par[2];
  double x = par[3];

  double eps = 0.5*(1.0-std::abs(x));
  Complex z = x + eps*std::exp( I*theta );

  Complex f = std::pow( 1.0-z, n+alpha ) * std::pow( 1.0+z, n+beta ) * std::exp( -n*I*theta );

  return f.real();
}




double jacobi_integral( gsl_integration_workspace* w,
                        const double n, const double alpha, const double beta, const double x ){
  double result, error;
  std::vector<double> par(4);
  par[0] = n; // n
  par[1] = alpha; // alpha
  par[2] = beta; // beta
  par[3] = x; // x

  double eps = 0.5*(1.0-std::abs(x));

  gsl_function F;
  F.function = &integrand;
  F.params = &par;

  gsl_integration_qag (&F, 0, 2.0*M_PI, epsabs, epsrel, 1000, 6,
                       w, &result, &error);

  double factor = std::pow(-0.5/eps, n);
  factor /= std::pow( 1.0-x, alpha ) * std::pow( 1.0+x, beta );
  factor /= 2.0*M_PI;

  return factor * result;
}



int main(int argc, char* argv[]){

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  const int n=21;
  const int alpha=0;
  const int beta=-2;
  // const int n=1;
  // const int alpha=2;
  // const int beta=-2;
  const double x=-0.5;
  std::cout << jacobi_integral(w, n, alpha, beta, x) << std::endl;

  gsl_integration_workspace_free (w);

  // double result, error;
  // double expected = 1.25;

  // double alpha = 1.0;
  // std::vector<double> par(4);
  // par[0] = 1.0; // n
  // par[1] = 0.0; // alpha
  // par[2] = -10.0; // beta
  // par[3] = 0.4; // x

  // gsl_function F;
  // F.function = &f;
  // F.params = &par;

  // gsl_integration_qag (&F, 0, 2.0*M_PI, 0, 1e-10, 1000, 6,
  //                      w, &result, &error);

  // printf ("result          = % .18f\n", result);
  // printf ("exact result    = % .18f\n", expected);
  // printf ("estimated error = % .18f\n", error);
  // printf ("actual error    = % .18f\n", result - expected);
  // printf ("intervals =  %d\n", w->size);

  return 0;
}

