#pragma once

#include <array>
#include <cmath>
#include <map>
#include <Eigen/Dense>

#include "lattice.h"


using Complex = std::complex<double>;
const Complex I = Complex(0.0, 1.0);
using Idx = int;

constexpr int NS = 2;
using MS=Eigen::Matrix2cd;

constexpr int DIM = 2;
using VD=Eigen::Vector2d;


double Mod(double a, double b=2.0*M_PI){
  int p = int(std::round(a / b));
  double r = a - p*b;
  return r;
}


Idx idx(const Idx x, const Idx y) { return x+N*y; }

void get_xy(Idx &x, Idx &y, const Idx i) {
  x = (i+N)%N;
  y = (i-x)/N;
}

void cshift(Idx& xp, Idx& yp, const Idx x, const Idx y, const int mu) {
  int dx, dy;
  if(mu%3==0){
    dx = -1;
    dy = 0;
  }
  else if(mu%3==1){
    dx = +1;
    dy = -1;
  }
  else if(mu%3==2){
    dx = 0;
    dy = +1;
  }
  else assert(false);

  if(mu>=3){
    dx *= -1;
    dy *= -1;
  }
  xp = (x+dx+N)%N;
  yp = (y+dy+N)%N;
}


double alpha(const int mu) {
  double res;
  if(mu%3==0){
    res = M_PI;
  }
  else if(mu%3==1){
    res = -1.0/3.0 * M_PI;
  }
  else if(mu%3==2){
    res = +1.0/3.0 * M_PI;
  }
  else assert(false);

  if(mu>=3) res = Mod(res+M_PI);

  return res;
}


Idx cshift(const Idx i, const int mu) {
  Idx x, y, xp, yp;
  get_xy(x, y, i);
  cshift(xp, yp, x, y, mu);
  return idx(xp, yp);
}


struct Dirac1fonFlat {
  QfeLattice& lattice;

  const double m;
  const double r;

  std::array<MS, 4> sigma;

  double ell;
  double ellstar;
  double site_vol;

  Dirac1fonFlat()=delete;

  Dirac1fonFlat(QfeLattice& lattice_,
	      const double m_=0.0,
	      const double r_=1.0)
    : lattice(lattice_)
    , m(m_)
    , r(r_)
  {
    set_sigma();
    ellstar = 1.0/std::sqrt(3); // @@@
    site_vol = std::sqrt(3); // @@@
  }

  Dirac1fonFlat & operator=(const Dirac1fonFlat&) = delete;

  void set_sigma(){
    assert(NS==2);
    sigma[0] << 1,0,0,1;
    sigma[1] << 0,1,1,0;
    sigma[2] << 0,-I,I,0;
    sigma[3] << 1,0,0,-1;
  }

  MS gamma(const int mu) const {
    const double al = alpha(mu); // @@@
    return std::cos(al)*sigma[1] + std::sin(al)*sigma[2];
  }

  Eigen::MatrixXcd matrix_form() const {
    Eigen::MatrixXcd res = Eigen::MatrixXcd::Zero(NS*lattice.n_sites, NS*lattice.n_sites);

    for(int ix=0; ix<lattice.n_sites; ix++){
      for(int mu=0; mu<6; mu++){
	const int iy = cshift(ix, mu);
	res.block<NS,NS>(NS*ix,NS*iy) = - ellstar * ( r*sigma[0] - gamma(mu) );
      }
      res.block<NS,NS>(NS*ix,NS*ix) = site_vol * (m + DIM*r) * sigma[0];
    }

    return res;
  } // end matrix_form


};
