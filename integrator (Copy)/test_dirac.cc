#include <iostream>
#include <fstream>
#include <cstdlib>

#include "dirac_s2.h"

constexpr double TOLLOOSE=1.0e-10;

double Mod(double a, double b=2.0*M_PI){
  int p = int(std::round(a / b));
  double r = a - p*b;
  return r;
}

bool isModdable( const double v, const double q, const double TOL=TOLLOOSE ){
  const double tmp1 = Mod(std::abs(v),q);
  const double tmp2 = Mod(-std::abs(v),q);
  return tmp1<TOL || tmp2<TOL;
}


int main(int argc, char* argv[]){
  

  return 0;
}

