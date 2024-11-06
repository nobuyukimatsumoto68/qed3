#pragma once

#include <array>
#include <cmath>
#include <map>
#include <Eigen/Dense>

#include "s2n.h"




double Mod(double a, double b=2.0*M_PI){
  int p = int(std::round(a / b));
  double r = a - p*b;
  return r;
}


VE circumcenter(const VE& r0, const VE& r1, const VE& r2){
  const VE r10 = r1 - r0;
  const VE r20 = r2 - r0;

  const VE tmp1 = r10.squaredNorm() * r20 - r20.squaredNorm() * r10;
  const VE cross = r10.cross(r20);
  const VE numer = tmp1.cross(cross);
  const double denom = 2.0*cross.squaredNorm();

  return numer/denom + r0;
}




struct SpinStructure{
  std::map<const Link, const double> omega;
  std::map<const Link, const double> alpha;

  SpinStructure(const int n_refine)
  {
    {
      std::cout << "reading omega" << std::endl;
      std::ifstream file("./dats/omega_dual_n"+std::to_string(n_refine)+"_singlepatch.dat");

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
      std::cout << "reading alpha" << std::endl;
      std::ifstream file("./dats/alpha_dual_n"+std::to_string(n_refine)+"_singlepatch.dat");

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
  }
};




struct Dirac1fonS2 : public SpinStructure{
  const Lattice& lattice;

  const double m;
  const double r;
  double a = 1.0;

  // SpinStructure spin;
  std::array<MS, 4> sigma;

  std::vector<double> ell; // link label
  std::vector<double> link_volume; // link label
  // std::vector<std::array<VD, 3>> exM; // site, M=A,B,C

  Dirac1fonS2()=delete;

  Dirac1fonS2(const Lattice& lattice_,
	      const double m_=0.0,
	      const double r_=1.0)
    : SpinStructure(lattice_.n_refine)
    , lattice(lattice_)
    , m(m_)
    , r(r_)
    // , spin(lattice.n_refine)
    , ell(lattice.n_links)
    , link_volume(lattice.n_links)
    // , exM(lattice.n_sites)
  {
    set_sigma();
    set_ell_and_link_volumes();

    // check
    double TOL=1.0e-6;
    {
      std::cout << "checking spin structure" << std::endl;
      for(int ix=0; ix<lattice.n_sites; ix++){
	for(int iy : lattice.nns[ix]){
	  const double alpha1 = alpha.at(Link{ix,iy});
	  double alpha2 = alpha.at(Link{iy,ix});
	  double omega12 = omega.at(Link{ix,iy});

	  double diff = (alpha2 + M_PI + omega12) - alpha1;
	  assert( std::abs(Mod(diff))<TOL );
	}}
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
  

  // MS gamma(const int ix, const int iy, const double shift=0.0) const { // located at x
  //   const double al = alpha.at(Link{ix,iy}) + shift;
  //   // return std::cos(al)*sigma[1] + std::sin(al)*sigma[2];
  //   return std::cos(al)*sigma[1] + std::sin(al)*sigma[2];
  // }
  MS gamma(const int ix, const int jj) const { // located at x
    // const double al = alpha.at(Link{ix,iy}) + shift;
    // return std::cos(al)*sigma[1] + std::sin(al)*sigma[2];
    return lattice.v[ix][jj](0) * sigma[1] + lattice.v[ix][jj](1) * sigma[2];
  }

  MS Omega(const int ix, const int iy) const {
    const double om = omega.at(Link{ix,iy});
    return std::cos(0.5*om)*sigma[0] - I*std::sin(0.5*om)*sigma[3];
  }

  Eigen::MatrixXcd matrix_form() const {
    Eigen::MatrixXcd res = Eigen::MatrixXcd::Zero(NS*lattice.n_sites, NS*lattice.n_sites);

    for(int ix=0; ix<lattice.n_sites; ix++){
      // for(int iy : lattice.nns[ix]){
      for(int jj=0; jj<3; jj++){
	int iy = lattice.nns[ix][jj];
	// naive
	// res.block<NS,NS>(NS*ix,NS*iy) += 0.5*ellstar[il] * gamma(ix, iy) * Omega(ix, iy);
	// res.block<NS,NS>(NS*ix,NS*iy) += gamma(ix, iy) * Omega(ix, iy);

	// wilson // WITH GEODESIC
	// res.block<NS,NS>(NS*ix,NS*iy) -= 0.25 * (link_volume[il]/ell[il]) * (r*sigma[0] - gamma(ix, iy)) * Omega(ix, iy);
	// res.block<NS,NS>(NS*ix,NS*iy) -= (r*sigma[0] - gamma(ix, iy)) * Omega(ix, iy);
	// res.block<NS,NS>(NS*ix,NS*iy) -= (r*sigma[0] - gamma(ix, iy)) * Omega(ix, iy);
	res.block<NS,NS>(NS*ix,NS*iy) -= lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj]*sigma[0] - gamma(ix, jj)) * Omega(ix, iy);
	// res.block<NS,NS>(NS*ix,NS*iy) -= (r*sigma[0] - gamma(ix, jj)) * Omega(ix, iy);

	res.block<NS,NS>(NS*ix,NS*ix) += lattice.vol[ix]/lattice.mean_vol * r*lattice.u[ix][jj]*sigma[0];
      }

      // res.block<NS,NS>(NS*ix,NS*ix) = site_vol[ix] * (m + DIM*r) * sigma[0];
      // res.block<NS,NS>(NS*ix,NS*ix) = (m + 3*r) * sigma[0];
      // res.block<NS,NS>(NS*ix, NS*ix) += site_vol[ix] * m * sigma[0];
    }

    return res;
  } // end matrix_form


  void set_ell_and_link_volumes() {
    // for(int il=0; il<lattice.n_links; il++) {
    //   const auto link = lattice.links[il];
    //   const int iA = link.faces[0];
    //   const int iB = link.faces[1];

    //   double ellA=0.0, ellB=0.0;
    //   double areaA=0.0, areaB=0.0;
    //   {
    // 	const QfeFace& face = lattice.faces[iA];
    // 	VE r0, r1, r2; // r0,1: link
    // 	if(face.sites[0]==link.sites[0] && face.sites[1]==link.sites[1]){
    // 	  r0 = lattice.r[face.sites[0]];
    // 	  r1 = lattice.r[face.sites[1]];
    // 	  r2 = lattice.r[face.sites[2]];
    // 	}
    // 	else if(face.sites[1]==link.sites[0] && face.sites[2]==link.sites[1]){
    // 	  r0 = lattice.r[face.sites[1]];
    // 	  r1 = lattice.r[face.sites[2]];
    // 	  r2 = lattice.r[face.sites[0]];
    // 	}
    // 	else if(face.sites[2]==link.sites[0] && face.sites[0]==link.sites[1]){
    // 	  r0 = lattice.r[face.sites[2]];
    // 	  r1 = lattice.r[face.sites[0]];
    // 	  r2 = lattice.r[face.sites[1]];
    // 	} // reverse
    // 	else if(face.sites[1]==link.sites[0] && face.sites[0]==link.sites[1]){
    // 	  r0 = lattice.r[face.sites[1]];
    // 	  r1 = lattice.r[face.sites[0]];
    // 	  r2 = lattice.r[face.sites[2]];
    // 	}
    // 	else if(face.sites[2]==link.sites[0] && face.sites[1]==link.sites[1]){
    // 	  r0 = lattice.r[face.sites[2]];
    // 	  r1 = lattice.r[face.sites[1]];
    // 	  r2 = lattice.r[face.sites[0]];
    // 	}
    // 	else if(face.sites[0]==link.sites[0] && face.sites[2]==link.sites[1]){
    // 	  r0 = lattice.r[face.sites[0]];
    // 	  r1 = lattice.r[face.sites[2]];
    // 	  r2 = lattice.r[face.sites[1]];
    // 	}
    // 	else assert(false);

    // 	//
    // 	const VE p = circumcenter(r0, r1, r2).transpose();
    // 	assert( std::abs( (p-r0).norm() - (p-r1).norm() )<1.0e-14 );
    // 	assert( std::abs( (p-r0).norm() - (p-r2).norm() )<1.0e-14 );

    // 	double a_ = std::acos( r0.dot(p) /(r0.norm()* p.norm()) );
    // 	double b_ = std::acos( r1.dot(p) /(r1.norm()* p.norm()) );
    // 	double c_ = std::acos( r0.dot(r1)/(r0.norm()*r1.norm()) ); // ell

    // 	double s_ = 0.5*(a_+b_+c_);
    // 	double tmp = std::tan(0.5*s_) * std::tan(0.5*(s_-a_)) * std::tan(0.5*(s_-b_)) * std::tan(0.5*(s_-c_));
    // 	double area_ = 4.0*std::atan( std::sqrt( tmp ) );

    // 	ellA = c_;
    // 	areaA = area_;
    //   }
    //   {
    // 	const QfeFace& face = lattice.faces[iB];
    // 	VE r0, r1, r2; // r0,1: link
    // 	if(face.sites[0]==link.sites[0] && face.sites[1]==link.sites[1]){
    // 	  r0 = lattice.r[face.sites[0]];
    // 	  r1 = lattice.r[face.sites[1]];
    // 	  r2 = lattice.r[face.sites[2]];
    // 	}
    // 	else if(face.sites[1]==link.sites[0] && face.sites[2]==link.sites[1]){
    // 	  r0 = lattice.r[face.sites[1]];
    // 	  r1 = lattice.r[face.sites[2]];
    // 	  r2 = lattice.r[face.sites[0]];
    // 	}
    // 	else if(face.sites[2]==link.sites[0] && face.sites[0]==link.sites[1]){
    // 	  r0 = lattice.r[face.sites[2]];
    // 	  r1 = lattice.r[face.sites[0]];
    // 	  r2 = lattice.r[face.sites[1]];
    // 	} // reverse
    // 	else if(face.sites[1]==link.sites[0] && face.sites[0]==link.sites[1]){
    // 	  r0 = lattice.r[face.sites[1]];
    // 	  r1 = lattice.r[face.sites[0]];
    // 	  r2 = lattice.r[face.sites[2]];
    // 	}
    // 	else if(face.sites[2]==link.sites[0] && face.sites[1]==link.sites[1]){
    // 	  r0 = lattice.r[face.sites[2]];
    // 	  r1 = lattice.r[face.sites[1]];
    // 	  r2 = lattice.r[face.sites[0]];
    // 	}
    // 	else if(face.sites[0]==link.sites[0] && face.sites[2]==link.sites[1]){
    // 	  r0 = lattice.r[face.sites[0]];
    // 	  r1 = lattice.r[face.sites[2]];
    // 	  r2 = lattice.r[face.sites[1]];
    // 	}
    // 	else assert(false);

    // 	const VE p = circumcenter(r0, r1, r2).transpose();
    // 	assert( std::abs( (p-r0).norm() - (p-r1).norm() )<1.0e-14 );
    // 	assert( std::abs( (p-r0).norm() - (p-r2).norm() )<1.0e-14 );

    // 	double a_ = std::acos( r0.dot(p) /(r0.norm()* p.norm()) );
    // 	double b_ = std::acos( r1.dot(p) /(r1.norm()* p.norm()) );
    // 	double c_ = std::acos( r0.dot(r1)/(r0.norm()*r1.norm()) ); // ell

    // 	double s_ = 0.5*(a_+b_+c_);
    // 	double tmp = std::tan(0.5*s_) * std::tan(0.5*(s_-a_)) * std::tan(0.5*(s_-b_)) * std::tan(0.5*(s_-c_));
    // 	double area_ = 4.0*std::atan( std::sqrt( tmp ) );

    // 	ellB = c_;
    // 	areaB = area_;
    //   }

    //   assert( std::abs(ellA-ellB)<1.0e-14 );
    //   ell[il] = ellA;
    //   link_volume[il] = areaA + areaB;
    // }

    // int counter = 0;
    // a = 0.0;
    // for(int il=0; il<lattice.n_links; il++) {
    //   a += ell[il];
    //   counter++;
    //   std::cout << "ell[il] =" << ell[il] << std::endl;
    // }
    // a /= counter;
    // a *= 0.1;
    // std::cout << "a = " << a << std::endl;
  }

};
