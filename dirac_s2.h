#pragma once

#include <array>
#include <Eigen/Dense>
#include "s2.h"

// using Complex = std::complex<double>;
const Complex I = Complex(0.0, 1.0);

constexpr int NS = 2;
using MS=Eigen::Matrix2cd;
constexpr int DIM = 2;

struct Dirac1fonS2 {
  QfeLatticeS2& lattice;
  const std::vector<std::vector<double>> alphas;
  const double kappa;
  const double m;
  const double r;
  std::array<MS, 4> sigma;

  Dirac1fonS2()=delete;

  Dirac1fonS2(QfeLatticeS2& lattice_,
	      const std::vector<std::vector<double>> alphas_,
	      const double kappa_,
	      const double m_=0.0,
	      const double r_=1.0)
    : lattice(lattice_)
    , alphas(alphas_)
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

  M2 gamma(const int ix, const int iA) const {
    const double alpha = alphas[ix][iA];
    return std::cos(alpha)*sigma[1] + std::sin(alpha)*sigma[2];
  }

  M2 Omega(const int il) const {
  }

  Eigen::MatrixXd matrix_form( const CompactU1onS2& U ) const {
    Eigen::MatrixXd res(NS*lattice.n_sites, NS*lattice.n_sites);

    for(int ix=0; ix<lattice.n_sites; ix++){
      for(int jj=0; jj<lattice.sites[ix].nn; jj++){
	const int iy = lattice.sites[ix].neighbors[jj];
	const int il = lattice.sites[ix].links[jj];

	// hopping
	if(ix<iy){
	  res.block(NS*ix,NS*iy, NS,NS) = - 0.5*kappa*U[il] * ( r*sigma[0] - gamma(ix,jj) ) * Omega(jj);
	  res.block(NS*iy,NS*ix, NS,NS) = - 0.5*kappa*std::conj(U[il]) * Omega(jj).adjoint() * ( r*sigma[0] + gamma(ix,jj) );
	}
	res.block(NS*ix,NS*ix, NS,NS) = (m + DIM*r) * sigma[0];
      }
    }


  } // end matrix_form

  

};
