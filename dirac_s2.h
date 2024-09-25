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


inline double double_mod(const double x, const double y=2.0*M_PI, const double z=-M_PI){
  double tmp = (x+z) + 100*y;
  return std::fmod(tmp, y) + z;
}


// VD projection(const Vec3& x){
//   const double r = x.norm();
//   VD res;

//   res[0] = std::acos(x[2]/r);

//   if( res[0]<1.0e-15 || M_PI-res[0]<1.0e-15){
//     res[1]=0.0;
//   }
//   else{
//     res[1]=std::asin(x[1] / (r*std::sin(res[0])));
//   }

//   return res;
// }






struct Dirac1fonS2 {
  QfeLatticeS2& lattice;

  // sign for the ordering of Evan's face.sites; +1 for clockwise rotation from the origin
  std::vector<int> face_signs; // index: ia (Evan's label for faces)

  // orientation: from a given direction (ilx=0), clockwise ilx=0,...,x.nn-1
  // The primary directions coincide between the nn_ and link_ lists
  // the function alpha( ilx ) gives the rotation angle for the ordered label ilx
  std::vector<std::vector<int> > nn_oriented; // first index: ix (Evan's label for sites)
  std::vector<std::vector<int> > link_oriented; // first index: ix (Evan's label for sites)

  // always outgoing from ix to ilx
  std::vector<double> alpha0; // index: ix (Evan's label for sites)
  // std::vector<std::vector<double> > alpha; // first index: ix (Evan's label for sites); second index: oriented ilx

  std::vector<double> omega; // index: il (Evan's label for links); direction x->y (ix<iy)

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
    , face_signs(lattice.n_faces)
    , alpha0(lattice.n_sites)
    // , alpha(lattice.n_sites)
    , omega(lattice.n_links)
    , m(m_)
    , r(r_)
  {
    set_sigma();
    set_face_signs();
    set_oriented_info();
    // set_omega();
  }


  // Dirac1fonS2( const Dirac1fonS2& other )
  //   : lattice(other.lattice)
  //   , kappa(other.kappa)
  //   , m(other.m)
  //   , r(other.r)
  // {
  //   set_sigma();
  // }
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


  void set_oriented_info(){
    for(int ix=0; ix<lattice.n_sites; ix++){
      auto x = lattice.sites[ix];
      std::vector<int> x_nn_oriented;
      std::vector<int> x_link_oriented;

      int iy = x.neighbors[0];
      x_nn_oriented.push_back(iy);
      x_link_oriented.push_back(x.links[0]);

      for(int ell=0; ell<x.nn-1; ell++){
	bool is_break = false;

	for(int kk=0; kk<x.nn; kk++){
	  const int iz = x.neighbors[kk];

	  if(this->is_nn(iy,iz)){
	    if(this->sign(ix, iy, iz)==1){
	      iy = iz;
	      x_nn_oriented.push_back(iy);
	      x_link_oriented.push_back(x.links[kk]);
	      is_break = true;
	    }
	  }

	  if(is_break) break;
	}
      }

      nn_oriented.push_back(x_nn_oriented);
      link_oriented.push_back(x_link_oriented);
    }
  }


  double alpha(const int ix, const int ixl) const {
    assert(0<=ixl && ixl<5);
    return alpha0[ix] + 2.0*M_PI/5.0 * ixl;
  }


  void set_omega(){
    for(int il=0; il<lattice.n_links; il++){
      const QfeLink link = lattice.links[il];
      const int ix = std::min(link.sites[0], link.sites[1]);
      const int iy = std::max(link.sites[0], link.sites[1]);

      const auto itx_ell = std::find(link_oriented[ix].begin(),
				     link_oriented[ix].end(),
				     il);
      const auto ity_ell = std::find(link_oriented[iy].begin(),
				     link_oriented[iy].end(),
				     il);
      const int ixl = std::distance(link_oriented[ix].begin(), itx_ell);
      const int iyl = std::distance(link_oriented[iy].begin(), ity_ell);

      // omega[il] = double_mod( this->alpha(ixl) - this->alpha(iyl) - M_PI, 2.0*M_PI );
      omega[il] = double_mod( this->alpha(ix, ixl) - this->alpha(iy, iyl) - M_PI, 2.0*M_PI );
    }
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
    // return std::cos(alpha[ix][iA])*sigma[1] + std::sin(alpha[ix][iA])*sigma[2];
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
