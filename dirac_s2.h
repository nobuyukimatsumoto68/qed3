#pragma once

#include <array>
#include <cmath>
#include <Eigen/Dense>
#include "s2.h"

// using Complex = std::complex<double>;
const Complex I = Complex(0.0, 1.0);

constexpr int NS = 2;
using MS=Eigen::Matrix2cd;

constexpr int DIM = 2;
using VD=Eigen::Vector2d;


VD projection(const Vec3& x){
  const double r = x.norm();
  VD res;

  res[0] = std::acos(x[2]/r);

  if( res[0]<1.0e-15 || M_PI-res[0]<1.0e-15){
    res[1]=0.0;
  }
  else{
    res[1]=std::asin(x[1] / (r*std::sin(res[0])));
  }

  return res;
}


struct Dirac1fonS2 {
  QfeLatticeS2& lattice;
  std::vector<std::vector<double>> alpha;
  std::vector<double> omega;
  const double kappa;
  const double m;
  const double r;
  std::array<MS, 4> sigma;

  Dirac1fonS2()=delete;

  Dirac1fonS2(QfeLatticeS2& lattice_,
	      const double kappa_=1.0,
	      const double m_=0.0,
	      const double r_=1.0)
    : lattice(lattice_)
    , kappa(kappa_)
    , m(m_)
    , r(r_)
  {
    set_sigma();
  }


  Dirac1fonS2( const Dirac1fonS2& other )
    : lattice(other.lattice)
    , kappa(other.kappa)
    , m(other.m)
    , r(other.r)
  {
    set_sigma();
  }


  Dirac1fonS2 & operator=(const Dirac1fonS2&) = delete;


  void set_sigma(){
    assert(NS==2);
    sigma[0] << 1,0,0,1;
    sigma[1] << 0,1,1,0;
    sigma[2] << 0,-I,I,0;
    sigma[3] << 1,0,0,-1;
  }


  bool is_nn(const int ix, const int iy) const {
    if(ix==iy) return false;

    bool res = false;
    const QfeSite x = lattice.sites[ix];

    for(int kk=0; kk<x.nn; kk++){
      const int iz = x.neighbors[kk];
      if(iy==iz) {
	res = true;
	break;
      }
    }
    return res;
  }

  // +1 if (ix, iy, iz) is RH
  int sign(const int ix, const int iy, const int iz) const {
    const Vec3 rx = lattice.r[ix];
    const Vec3 ry = lattice.r[iy];
    const Vec3 rz = lattice.r[iz];

    const Vec3 cross = (ry-rx).cross(rz-rx);
    const Vec3 sum = rx+ry+rz;
    const double inner = cross.dot(sum);

    int res = 0;
    if(inner>0) res = 1;
    else if(inner<0) res = -1;
    else assert(false);
    return res;
  }



  MS gamma(const int ix, const int iA) const {
    return std::cos(alpha[ix][iA])*sigma[1] + std::sin(alpha[ix][iA])*sigma[2];
  }


  MS Omega(const int il) const {
    return std::cos(omega[il])*sigma[0] + std::cos(omega[il])*sigma[3];
  }


  Eigen::MatrixXcd matrix_form() const {
    Eigen::MatrixXcd res(NS*lattice.n_sites, NS*lattice.n_sites);

    for(int ix=0; ix<lattice.n_sites; ix++){
      for(int jj=0; jj<lattice.sites[ix].nn; jj++){
	const int iy = lattice.sites[ix].neighbors[jj];

	// hopping
	if(ix<iy){
	  res.block<NS,NS>(NS*ix,NS*iy) = - 0.5*kappa * ( r*sigma[0] - gamma(ix,jj) ) * Omega(jj);
	  res.block<NS,NS>(NS*iy,NS*ix) = - 0.5*kappa * Omega(jj).adjoint() * ( r*sigma[0] + gamma(ix,jj) );
	}
	res.block<NS,NS>(NS*ix,NS*ix) = (m + DIM*r) * sigma[0];
      }}

    return res;
  } // end matrix_form


  // Eigen::MatrixXd matrix_form( const CompactU1onS2& U ) const {
  //   Eigen::MatrixXd res(NS*lattice.n_sites, NS*lattice.n_sites);

  //   for(int ix=0; ix<lattice.n_sites; ix++){
  //     for(int jj=0; jj<lattice.sites[ix].nn; jj++){
  // 	const int iy = lattice.sites[ix].neighbors[jj];
  // 	const int il = lattice.sites[ix].links[jj];

  // 	// hopping
  // 	if(ix<iy){
  // 	  res.block(NS*ix,NS*iy, NS,NS) = - 0.5*kappa*U[il] * ( r*sigma[0] - gamma(ix,jj) ) * Omega(jj);
  // 	  res.block(NS*iy,NS*ix, NS,NS) = - 0.5*kappa*std::conj(U[il]) * Omega(jj).adjoint() * ( r*sigma[0] + gamma(ix,jj) );
  // 	}
  // 	res.block(NS*ix,NS*ix, NS,NS) = (m + DIM*r) * sigma[0];
  //     }
  //   }


  // } // end matrix_form

  

};
