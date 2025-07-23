#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <map>
#include <limits>
#include <functional>
#include <stdfloat>

#include <vector>
#include <numbers>

#include <gsl/gsl_sf_hyperg.h>

#include <Eigen/Dense>
#include "boost/math/special_functions/jacobi.hpp"

using Real = double; // std::float64_t;
using Complex = std::complex<Real>;
using F=std::function<Real(const Real)>;
using VC=Eigen::Vector2cd;


using Int = int64_t;

constexpr Complex I = Complex(0.0, 1.0);

#include "integral.h"


const Int limit = 10000; // 1000;

// Real epsabs = 1.0e-15; // 0.;
// Real epsrel = 1.0e-13; // TOLLOOSE;
Real epsabs = 1.0e-2; // 0.;
Real epsrel = 1.0e-2; // TOLLOOSE;

int key = 5;

double hypergeometric( double a, double b, double c, double x )
{
  assert( c>0 );
  assert( std::abs(x)<1 );
  const double TOLERANCE = 1.0e-12;
  double term = a * b * x / c;
  double value = 1.0 + term;
  Int n = 1;

  while ( abs( term ) > TOLERANCE )
    {
      a++, b++, c++, n++;
      term *= a * b * x / c / n;
      value += term;
    }

  return value;
}

std::vector<std::vector<Int>> binom;

Int binomialCoeff(Int n, Int k) {
  if (k > n || k<0) return 0;
  if (k == 0 || k == n) return 1;
  return binom[n][k];
}

// Int binomialCoeff_gen(Int n, Int k) {
//   if (k > n || k<0) return 0;
//   if (k == 0 || k == n) return 1;
//   return binom[n-1][k-1] + binom[n-1][k];
// }


Real jacobi_( const Int n, const Int alpha, const Int beta, const Real z ){
  if(n==0) return 1.0;

  Real res = 0.0;
  Real pow1 = 1;
  // Real pow2 = std::pow( 0.5*(z+1.0), n );
  for(Int s=0; s<=n; s++){
    // std::cout << n+alpha << " " << n-s << " " << binomialCoeff( n+alpha, n-s ) << std::endl;
    // std::cout << n+beta << " " << s << " " << binomialCoeff( n+beta, s ) << std::endl;
    res += binomialCoeff( n+alpha, n-s ) * binomialCoeff( n+beta, s ) * pow1 * std::pow( 0.5*(z+1.0), n-s );
    pow1 *= 0.5*(z-1.0);
    // pow2 /= 0.5*(z+1.0);
  }
  return res;
}

Real jacobi_asymp_darboux(const Int n, const Int alpha, const Int beta, const double z){
  Real pre_inv = std::pow( std::sin(0.5*std::acos(z)), alpha+0.5 ) * std::pow( std::cos(0.5*std::acos(z)), beta+0.5 );
  pre_inv *= std::sqrt( M_PI * n );
  const Real angle = 0.5*( 2.0*n+alpha+beta+1.0 )*std::acos(z) - 0.25*(2.0*alpha+1.0)*M_PI;
  return std::cos(angle) / pre_inv;
}



// Real poch( const double x, const int n, const int M=100 ){
//   if(x>M) return std::exp( (x+n-1.0)*std::log(x+n-1.0) - (x-1.0)*std::log(x-1.0) - n );
//   return std::tgamma(x+n)/std::tgamma(x);
// }


// Real Theta(const double theta, const int n, const int alpha, const int beta, const int m, const int l){
//   Real res = 0.5*(2.0*n+alpha+beta+m+1.0) * theta;
//   res -= 0.5*(alpha+l+0.5)*M_PI;
//   return res;
// }


// inline Real tc( const int m, const int l, const int alpha, const int beta){
//   return poch( 0.5+alpha, l) * poch( 0.5-alpha, l) * poch( 0.5+beta, m-l) * poch( 0.5-beta, m-l);
// }

// Real f(const int m, const int n, const int alpha, const int beta, const double theta ){
//   Real sum=0.0;
//   for(int l=0; l<=m; l++){
//     Real factor = tc(m,l,alpha,beta);
//     factor /= std::tgamma( l+1 ) * std::tgamma( m-l+1 );

//     Real numer = std::cos( Theta(theta,n,alpha,beta,m,l) );
//     Real denom = std::pow( std::sin(0.5*theta), l ) * std::pow( std::cos(0.5*theta), m-l );

//     sum += factor * numer / denom;
//   }
//   return sum;
// }



// Real jacobi_asymp(const Int n, const Int alpha, const Int beta, const Int z, const Int M=20){
//   Real pre_inv = std::pow( std::sin(0.5*std::acos(z)), alpha+0.5 ) * std::pow( std::cos(0.5*std::acos(z)), beta+0.5 );
//   pre_inv *= M_PI;

//   Real pre2 = std::pow( 2, 2*n+alpha+beta+1 );
//   pre2 *= std::beta( n+alpha+1, n+beta+1 );

//   Real sum=0.0;
//   for(Int m=0; m<=M; m++){
//     Real tmp = f(m,n,alpha,beta,std::acos(z));
//     tmp /= std::pow(2, m);
//     tmp /= poch(2*n+alpha+beta+2, m);
//     if(!std::isnormal(tmp)){
//       std::cout << "f = " << f(m,n,alpha,beta,std::acos(z)) << std::endl;
//       std::cout << "poch = " << poch(2*n+alpha+beta+2, m) << std::endl;
//       std::cout << "n = " << n << std::endl;
//       std::cout << "m = " << m << std::endl;
//       assert(false);
//     }
//     sum += tmp;
//   }
//   return pre2 * sum / pre_inv;
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
  double factor = std::pow(-0.5/eps, n);
  factor /= std::pow( 1.0-x, alpha ) * std::pow( 1.0+x, beta );
  factor /= 2.0*M_PI;

  return factor * f.real();
}


double jacobi( gsl_integration_workspace* w,
               const double n, const double alpha, const double beta, const double x ){
  double result, error;
  std::vector<double> par(4);
  par[0] = n; // n
  par[1] = alpha; // alpha
  par[2] = beta; // beta
  par[3] = x; // x

  gsl_function F;
  F.function = &integrand;
  F.params = &par;

  gsl_integration_qag (&F, 0, 2.0*M_PI, 0, 1e-10, 1000, 6,
                       w, &result, &error);

  return result;
}






Real jacobi( const Int n, const Int alpha, const Int beta, const double z, const int nthreshold=50 ){
  if(n<nthreshold) return jacobi_(n, alpha, beta, z);
  // else return jacobi_asymp(n, alpha, beta, z);
  else return jacobi_asymp_darboux(n, alpha, beta, z);
}





// Real jacobi( const Int n, const double alpha, const double beta, const double z ){
//   // return boost::math::jacobi(n, 1.0*alpha, 1.0*beta, z);
//   // return boost::math::jacobi(n, 1.0*alpha, 1.0*beta, z);
//   // gsl_sf_result result;
//   // gsl_sf_hyperg_2F1_e(-n, 1.0+alpha+beta+n, alpha+1.0, 0.5*(1.0-z), &result);
//   // return result.val;
//   // return hypergeometric(-n, 1.0+alpha+beta+n, alpha+1.0, 0.5*(1.0-z));
//   return jacobi_( n, 1.0+alpha+beta+n, alpha+1.0, 0.5*(1.0-z));
// }


// n = 1, 2, ...
Real fn0_( const Int n, const Real z ){
  assert(n>0);
  return (1.0-z) * jacobi( n, 1, -1, z );
}

// n = 1, 2, ...
Real Dtilfn0_( const Int n, const Real z){
  assert(n>0);
  // std::cout << "jacobi1: " << jacobi( n, 1, -1, z ) << std::endl;
  // std::cout << "jacobi2: " << jacobi( n-1, 2, 0, z ) << std::endl;
  return (1.0-z*z) * ( -jacobi( n, 1, -1, z ) + (1.0-z)*0.5*(n+1.0) * jacobi( n-1, 2, 0, z ) );
}


// n = 0, 1, 2, ...
Real fnm_( const Int n, const Int m, const Real z){
  assert(m!=0);
  assert(n>=0);
  const Int sign = m>0 ? 1 : -1;
  const Int absm = sign*m;
  return jacobi( n+absm+1, absm-1, -absm-1, -sign*z );
}


// n = 0, 1, 2, ...
Real Dtilfnm_( const Int n, const Int m, const Real z){
  assert(m!=0);
  assert(n>=0);

  const Int sign = m>0 ? 1 : -1;
  const Int absm = sign*m;

  // const Real first = -z * jacobi( n+absm+1, absm-1, -absm-1, -sign*z );
  const Real secon = -sign * 0.5 * (n+absm) * jacobi( n+absm, absm, -absm, -sign*z );
  return (1.0-z*z) * secon;
}


Real fnm( const Int n, const Int m, const Real z){
  if(m==0) return fn0_( n, z );
  else return fnm_( n, m, z );
}

Real Dtilfnm( const Int n, const Int m, const Real z){
  if(m==0) return Dtilfn0_( n, z );
  else return Dtilfnm_( n, m, z );
}




Real Cn0_( const Int n ){
  assert(n>0);
  Real tmp = 8.0*M_PI;
  tmp /= 2.0*n+1.0;
  tmp *= n+1.0;
  tmp /= 1.0*n;
  return std::sqrt(tmp);
}


Real Cnm_( const Int n, const Int m ){
  assert(m!=0);
  assert(n>=0);

  const Int absm = std::abs(m);

  Real tmp = 2.0*M_PI;
  tmp /= 2.0*n+2.0*absm+1.0;
  tmp *= std::tgamma(n+2*absm+1);
  tmp *= std::tgamma(n+1);
  tmp /= std::tgamma(n+absm);
  tmp /= std::tgamma(n+absm+2);
  return std::sqrt(tmp);
}

Real Cnm( const Int n, const Int m ){
  if(m==0) return Cn0_( n );
  else return Cnm_( n, m );
  // return 1.0;
}


Real xinm( const Int n, const Int m, const Real z ){
  if(std::abs(z-1.0)<1.0e-15){
    if(m==-1) return 0.5 / Cnm( n, m );
    else return 0.0;
  }

  const Real factor1 = std::pow( 1.0-z, -0.5*( m+1) );
  const Real factor2 = std::pow( 1.0+z, -0.5*(-m+1) );
  return factor1 * factor2 * fnm( n, m, z ) / Cnm( n, m );
}

Real Dxinm( const Int n, const Int m, const Real z ){
  if(std::abs(z-1.0)<1.0e-15){
    if(m==-1) return -0.5 / Cnm( n, m );
    else return 0.0;
  }

  const Real factor1 = std::pow( 1.0-z, -0.5*( m+1) );
  const Real factor2 = std::pow( 1.0+z, -0.5*(-m+1) );

  const Real derivA = m * factor1 * factor2 * fnm( n,m,z ) / Cnm( n, m );
  const Real derivB = factor1 * factor2 * Dtilfnm( n,m,z ) / Cnm( n, m );

  return derivA + derivB;
}


Int lambdanm( const Int n, const Int m ){
  const Int absm = std::abs(m);
  return (n+absm)*(n+absm+1);
}



// const Int m_max = 100;
// const Int n_max = 80;
const Int n_max = 20;
// const Int n_max = 15;

// const Int m_max = 6;
// const Int n_max = 8;


int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  const int nint=20;
  const Real delta=M_PI/nint;
  // Real z0 = 0.999;
  Real z0 = 1.0;

  binom.resize(10*n_max);
  for(Int n=0; n<10*n_max; n++){
    binom[n].resize(n+1);
    binom[n][0] = 1;
    binom[n][n] = 1;
    for(Int k=1; k<n; k++) binom[n][k] = binom[n-1][k-1] + binom[n-1][k];
  }

  // gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  // const int n=1;
  // std::cout << "fn0_( 1, z0 ) = " << fn0_( n, z0 ) << std::endl;
  // std::cout << "fnm_( 1, 1, 0.5 ) = " << fnm_( n, 1, 0.5 ) << std::endl;

  // // std::cout << "jacobi( n, 0, -2, 0.5 ) = " << jacobi( w, n, 0, -2, 0.5 ) << std::endl;
  // // std::cout << "jacobi( n, 1, -1, 0.5 ) = " << jacobi( w, n, 1, -1, 0.5 ) << std::endl;
  // // std::cout << "jacobi( n, 0, -2, -0.5 ) = " << jacobi( w, n, 0, -2, -0.5 ) << std::endl;
  // // std::cout << "jacobi( n, 1, -1, -0.5 ) = " << jacobi( w, n, 1, -1, -0.5 ) << std::endl;

  // gsl_integration_workspace_free (w);


  // std::cout << "fn0( 1, 0, 0.5 ) = " << fnm( n, 0, 0.5 ) << std::endl;
  // std::cout << "fnm( 1, 1, 0.5 ) = " << fnm( n, 1, 0.5 ) << std::endl;

  // std::cout << "Dtilfnm( 1, 0, 0.5 ) = " << Dtilfnm( n, 0, 0.5 ) << std::endl;
  // std::cout << "Dtilfn0_( 1, 0, 0.5 ) = " << Dtilfn0_( n, 0.5 ) << std::endl;
  // std::cout << "Dtilfnm( 1, 1, 0.5 ) = " << Dtilfnm( n, 1, 0.5 ) << std::endl;

  // std::cout << "xinm( 1, 0, 0.5 ) = " << xinm( n, 0, 0.5 ) << std::endl;
  // std::cout << "xinm( 1, 1, 0.5 ) = " << xinm( n, 1, 0.5 ) << std::endl;

  // std::cout << "Dxinm( 1, 0, 0.5 ) = " << Dxinm( n, 0, 0.5 ) << std::endl;
  // std::cout << "Dxinm( 1, 1, 0.5 ) = " << Dxinm( n, 1, 0.5 ) << std::endl;

  // {
  //   Real alpha = 0.0;
  //   Real beta = -2.0;
  //   Real z = -0.5;
  //   // gsl_sf_result result;
  //   // gsl_sf_hyperg_2F1_e(-n, 1.0+alpha+beta+n, alpha+1.0, 0.5*(1.0-z), &result);
  //   // std::cout << "Dxinm( 1, 1, 0.5 ) = " << result.val << std::endl;
  //   std::cout << "2F1( 1, 1, 0.5 ) = " << hypergeometric(-n, 1.0+alpha+beta+n, alpha+1.0, 0.5*(1.0-z)) << std::endl;
  //   // return ;
  //   //
  //   //gsl_sf_hyperg_2F1_e(-n, 1.0+alpha+beta+n, alpha+1.0, 0.5*(1.0-z), &result)
  // }

  // return 1;


  for(Int i=1; i<nint; i++){
    const Real theta = delta*i;

    const Real z = std::cos(theta);
    Complex wf=0.0;

    Real One=1.0;
    for(Int n=0; n<=n_max; n++){
      Int m = -1;
      // for(Int m=-m_max; m<=m_max; m++){
      // if(n==0 && m==0) continue;
      // {
      //   Real source = xinm(n, m, z0);
      //   Real sink = xinm(n, m, z);
      //   wf += source * sink / lambdanm(n,m);
      // }
      // std::cout << i << " " << n << std::endl;
      Complex tmp=0.0;
      {
        Real source = Dxinm(n, m, z0);
        Real sink = Dxinm(n, m, z);
        tmp += source * sink;// / lambdanm(n,m);
        // std::cout << source * sink << std::endl;
      }
      {
        Real source = Dxinm(n, m, z0);
        Real sink = xinm(n, m, z);
        tmp -= I*(One*m) * source * sink;// / lambdanm(n,m);
        // std::cout << -I*(One*m) *source * sink << std::endl;
      }
      {
        Real source = xinm(n, m, z0);
        Real sink = Dxinm(n, m, z);
        tmp += I*(One*m) * source * sink;// / lambdanm(n,m);
        // std::cout << I*(One*m) *source * sink << std::endl;
      }
      {
        Real source = xinm(n, m, z0);
        Real sink = xinm(n, m, z);
        tmp += (One*m)*(One*m) * source * sink;// / lambdanm(n,m);
        // std::cout << (One*m*m) *source * sink << std::endl;
      }
      tmp /= One*lambdanm(n,m);
      // std::cout << "tmp = " << tmp << std::endl;
      wf += tmp;
    }
    // }
  // std::cout << theta << " " << std::abs(wf) << std::endl;
  std::cout << theta << " " << std::real(wf) << " " << std::imag(wf) << std::endl;
  }

  return 0;
}

