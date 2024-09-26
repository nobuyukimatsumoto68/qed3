#pragma once

#include <array>
#include <cmath>
#include <map>
#include <Eigen/Dense>
#include "s2.h"

// using Complex = std::complex<double>;
const Complex I = Complex(0.0, 1.0);

constexpr int NS = 2;
using MS=Eigen::Matrix2cd;

constexpr int DIM = 2;
using VD=Eigen::Vector2d;


double Mod(double a, double b=2.0*M_PI){
  int p = int(std::round(a / b));
  double r = a - p*b;
  return r;
}


struct SpinStructure{

  using Link = std::array<int,2>; // <int,int>;

  std::map<const Link, const double> omega;
  std::map<const Link, const double> alpha;
  std::map<const int, const int> NM2EO;

  std::map<const Link, const double> omegaEO;
  std::map<const Link, const double> alphaEO;


  void set_omega() const {
  }

  SpinStructure()
  {
    // set_omega();
    {
      std::ifstream file("omega.dat");

      std::string str;
      std::string file_contents;
      while (std::getline(file, str)){
	std::istringstream iss(str);
	int i,j;
	double v;
	iss >> i;
	iss >> j;
	iss >> v;
	omega.insert( { Link{i,j}, v } );
	omega.insert( { Link{j,i}, -v } );
      }
    }

    {
      std::ifstream file("alpha.dat");

      std::string str;
      std::string file_contents;
      while (std::getline(file, str)){
	std::istringstream iss(str);
	int i,j;
	double v;
	iss >> i;
	iss >> j;
	iss >> v;
	alpha.insert( {Link{i,j}, v} );
      }
    }

    {
      NM2EO.insert( {10, 0} );

      NM2EO.insert( { 3, 3} );
      NM2EO.insert( { 9, 5} );
      NM2EO.insert( { 1, 8} );
      NM2EO.insert( { 7, 9} );
      NM2EO.insert( { 8,11} );

      NM2EO.insert( { 6,10} );
      NM2EO.insert( { 5, 1} );
      NM2EO.insert( { 2, 2} );
      NM2EO.insert( {12, 4} );
      NM2EO.insert( { 4, 7} );

      NM2EO.insert( {11, 6} );
    }

    {
      for(auto elem : alpha){
	int ix1 = elem.first[0];
	int iy1 = elem.first[1];
	int ix2 = NM2EO[ix1];
	int iy2 = NM2EO[iy1];
	alphaEO.insert( { Link{ix2,iy2}, alpha[Link{ix1,iy1}] } );
      }
      for(auto elem : omega){
	int ix1 = elem.first[0];
	int iy1 = elem.first[1];
	int ix2 = NM2EO[ix1];
	int iy2 = NM2EO[iy1];
	omegaEO.insert( { Link{ix2,iy2}, omega[Link{ix1,iy1}] } );
      }
    }
  }
};




struct Dirac1fonS2 {
  QfeLatticeS2& lattice;

  // sign for the ordering of Evan's face.sites; +1 for clockwise rotation from the origin
  std::vector<int> face_signs; // index: ia (Evan's label for faces)

  const double kappa;
  const double m;
  const double r;

  using Link = std::array<int,2>; // <int,int>;
  const std::map<const Link, const double> omega;
  const std::map<const Link, const double> alpha;

  std::array<MS, 4> sigma;

  Dirac1fonS2()=delete;

  Dirac1fonS2(QfeLatticeS2& lattice_,
	      const double kappa_=1.0,
	      const double m_=0.0,
	      const double r_=1.0,
	      const SpinStructure spin_=SpinStructure() )
    : lattice(lattice_)
    , kappa(kappa_)
    , face_signs(lattice.n_faces)
    , m(m_)
    , r(r_)
    , omega(spin_.omegaEO)
    , alpha(spin_.alphaEO)
  {
    set_sigma();
    set_face_signs();

    // check
    {
      for(int ix=0; ix<lattice.n_sites; ix++){
	for(int jj=0; jj<lattice.sites[ix].nn; jj++){
	  const int iy = lattice.sites[ix].neighbors[jj];

	  const double alpha1 = alpha.at(Link{ix,iy});
	  double alpha2 = alpha.at(Link{iy,ix});
	  double omega12 = omega.at(Link{ix,iy});

	  double diff = (alpha2 + M_PI + omega12) - alpha1;
	  assert( std::abs(Mod(diff))<1.0e-14 );
	}}
    }

    {
      for(int ia=0; ia<lattice.n_faces; ia++){
	double omega_sum = 0.0;

	for(int i=0; i<3; i++){
	  int ix = lattice.faces[ia].sites[i];
	  int iy = lattice.faces[ia].sites[(i+1)%3];
	  omega_sum += omega.at(Link{ix,iy});
	}

	double diff = Mod( face_signs[ia]*omega_sum ) + M_PI/5.0;
	assert( std::abs(diff)<1.0e-14 );
      }
    }
  }

  Dirac1fonS2 & operator=(const Dirac1fonS2&) = delete;

  void set_sigma(){
    assert(NS==2);
    sigma[0] << 1,0,0,1;
    sigma[1] << 0,1,1,0;
    sigma[2] << 0,-I,I,0;
    sigma[3] << 1,0,0,-1;
  }


  int face_sign(const int i_face) const {
    const QfeFace& face = lattice.faces[i_face];
    const Vec3 r0 = lattice.r[face.sites[0]];
    const Vec3 r1 = lattice.r[face.sites[1]];
    const Vec3 r2 = lattice.r[face.sites[2]];

    const Vec3 cross = (r1-r0).cross(r2-r0);
    const Vec3 sum = r0+r1+r2;

    const double inner = cross.dot(sum);

    int res = 1;
    if(inner<0) res = -1;
    return res;
  }


  void set_face_signs() { for(int i=0; i<lattice.n_faces; i++) face_signs[i] = face_sign(i); }


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


  MS gamma(const int ix, const int iy, const double shift = 0.0) const { // located at x
    const double al = alpha.at(Link{ix,iy}) + shift;
    return std::cos(al)*sigma[1] + std::sin(al)*sigma[2];
    // return std::cos(al)*sigma[1] - std::sin(al)*sigma[2];
  }

  MS Omega(const int ix, const int iy) const {
    const double om = omega.at(Link{ix,iy});
    return std::cos(0.5*om)*sigma[0] - I*std::sin(0.5*om)*sigma[3];
  }


  Eigen::MatrixXcd matrix_form() const {
    Eigen::MatrixXcd res(NS*lattice.n_sites, NS*lattice.n_sites);

    for(int ix=0; ix<lattice.n_sites; ix++){
      for(int jj=0; jj<lattice.sites[ix].nn; jj++){
	const int iy = lattice.sites[ix].neighbors[jj];
	res.block<NS,NS>(NS*ix,NS*iy) = - 0.5*kappa * ( r*sigma[0] - gamma(ix, iy) ) * Omega(ix, iy);
	// res.block<NS,NS>(NS*ix,NS*iy) = - 0.5*kappa * Omega(ix, iy) * ( r*sigma[0] - gamma(iy, ix, M_PI) );
      }
      res.block<NS,NS>(NS*ix,NS*ix) = (m + DIM*r) * sigma[0];
    }

    return res;
  } // end matrix_form

};
