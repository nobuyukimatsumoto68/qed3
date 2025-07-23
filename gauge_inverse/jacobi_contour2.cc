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
Real epsabs = 1.0e-15; // 0.;
Real epsrel = 1.0e-15; // TOLLOOSE;
int key = 5;


// template<typename T>
// T w( T x, double alpha, double beta ){
//   return std::pow( 1.0-x, alpha ) * std::pow( 1.0+x, beta );
// }

double rho = 0.5;

double integrand (double theta, void * params) {
  std::vector<double> par = *(std::vector<double> *) params;
  double n = par[0];
  double alpha = par[1];
  double beta = par[2];
  double x = par[3];

  double eps = rho*(1.0-std::abs(x));
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

  double eps = rho*(1.0-std::abs(x));
  Complex z = x + eps*std::exp( I*theta );

  Complex f = std::pow( 1.0-z, n+alpha ) * std::pow( 1.0+z, n+beta );

  return f.real();
}

double integrand_sin (double theta, void * params) {
  std::vector<double> par = *(std::vector<double> *) params;
  double n = par[0];
  double alpha = par[1];
  double beta = par[2];
  double x = par[3];

  double eps = rho*(1.0-std::abs(x));
  Complex z = x + eps*std::exp( I*theta );

  Complex f = std::pow( 1.0-z, n+alpha ) * std::pow( 1.0+z, n+beta );

  return f.imag();
}




double jacobi_integral( gsl_integration_workspace* w,
                        const double n, const double alpha, const double beta, const double x ){
  double result, error;
  std::vector<double> par(4);
  par[0] = n; // n
  par[1] = alpha; // alpha
  par[2] = beta; // beta
  par[3] = x; // x

  double eps = rho*(1.0-std::abs(x));

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

  const int n=41;
  const int alpha=0;
  const int beta=-2;
  // const int n=1;
  // const int alpha=2;
  // const int beta=-2;
  const double x=-0.9;

  // std::cout << jacobi_integral(w, n, alpha, beta, x) << std::endl;


  double omega = n;
  double L = 2.0*M_PI;

  size_t nnn = 5000;

  gsl_integration_qawo_table *tc = gsl_integration_qawo_table_alloc(omega, L, GSL_INTEG_COSINE, nnn);
  gsl_integration_qawo_table *ts = gsl_integration_qawo_table_alloc(omega, L, GSL_INTEG_SINE,   nnn);

  gsl_integration_qawo_table_set(tc, omega, L, GSL_INTEG_COSINE);
  gsl_integration_qawo_table_set(ts, omega, L, GSL_INTEG_SINE);


  std::vector<double> par(4);
  par[0] = n; // n
  par[1] = alpha; // alpha
  par[2] = beta; // beta
  par[3] = x; // x


  gsl_function Fc;
  Fc.function = &integrand_cos;
  Fc.params = &par;

  double result_c, error_c;
  gsl_integration_qawo(&Fc, 0.0, epsabs, epsrel, 1000, w, tc, &result_c, &error_c);
  std::cout << result_c << " " << error_c << std::endl;

  gsl_function Fs;
  Fs.function = &integrand_sin;
  Fs.params = &par;

  double result_s, error_s;
  gsl_integration_qawo(&Fs, 0.0, epsabs, epsrel, 1000, w, ts, &result_s, &error_s);
  std::cout << result_s << " " << error_s << std::endl;

  //
  double eps = rho*(1.0-std::abs(x));

  double factor = std::pow(-0.5/eps, n);
  factor /= std::pow( 1.0-x, alpha ) * std::pow( 1.0+x, beta );
  factor /= 2.0*M_PI;

  std::cout << factor * (result_c+result_s) << " " << factor * error_c << std::endl;
  // std::cout << factor * result_s << " " << factor * error_s << std::endl;

  gsl_integration_qawo_table_free(tc);
  gsl_integration_qawo_table_free(ts);


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

