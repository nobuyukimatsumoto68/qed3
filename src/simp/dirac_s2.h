#pragma once

#include <array>
#include <cmath>
#include <map>
#include <Eigen/Dense>

#include "s2n.h"

// using Complex = std::complex<double>;

// using VD=Eigen::Vector2d;


// double Mod(double a, double b=2.0*M_PI){
//   int p = int(std::round(a / b));
//   double r = a - p*b;
//   return r;
// }



std::string dir = "/mnt/hdd_barracuda/qed3/dats/";


struct SpinStructure{

  using Link = std::array<int,2>; // <int,int>;

  std::map<const Link, const double> omega;
  std::map<const Link, const double> alpha;
  // std::map<const int, const int> NM2EO;


  SpinStructure(const int n_refine)
  {
    {
      // std::ifstream file("./dats/omega_n"+std::to_string(n_refine)+".dat");
      std::ifstream file(dir+"/omega_n"+std::to_string(n_refine)+"_singlepatch.dat");

      std::string str;
      std::string file_contents;
      while (std::getline(file, str)){
	std::istringstream iss(str);
	Idx i,j;
	double v;
	iss >> i;
	iss >> j;
	iss >> v;
	omega.insert( { Link{i,j}, v } );
	omega.insert( { Link{j,i}, -v } );
      }
    }

    {
      // std::ifstream file("./dats/alpha_n"+std::to_string(n_refine)+".dat");
      std::ifstream file(dir+"/alpha_n"+std::to_string(n_refine)+"_singlepatch.dat");

      std::string str;
      std::string file_contents;
      while (std::getline(file, str)){
	std::istringstream iss(str);
	Idx i,j;
	double v;
	iss >> i;
	iss >> j;
	iss >> v;
	alpha.insert( {Link{i,j}, v} );
      }
    }
  }
};






template<typename Gauge>
struct DiracS2Simp {
  using Lattice=S2Simp;

  using MS=Eigen::Matrix2cd;
  using VD=Eigen::Vector2d;
  using VE=Eigen::Vector3d;
  using VC=Eigen::VectorXcd;

  static constexpr int NS = 2;
  // using MS=Eigen::Matrix2cd;
  static constexpr int DIM = 2;
  static constexpr Complex I = Complex(0.0, 1.0);



  Lattice& lattice;

  // sign for the ordering of Evan's face.sites; +1 for clockwise rotation from the origin
  std::vector<int> face_signs; // index: ia (Evan's label for faces)

  const double m;
  const double r;

  double a = 1.0;

  using Link = std::array<int,2>; // <int,int>;
  SpinStructure spin;
  const std::map<const Link, const double> omega;
  const std::map<const Link, const double> alpha;

  std::array<MS, 4> sigma;

  std::vector<double> ell; // evan's link label

  // std::vector<double> ellstar; // evan's link label
  std::vector<double> link_volume; // evan's link label
  // std::vector<double> site_vol; // evan's site label

  DiracS2Simp()=delete;

  DiracS2Simp(Lattice& lattice_,
	      const int n_refine,
	      const double m_=0.0,
	      const double r_=1.0)
    : lattice(lattice_)
    , face_signs(lattice.n_faces)
    , m(m_)
    , r(r_)
    , spin(n_refine)
    , omega(spin.omega)
    , alpha(spin.alpha)
    , ell(lattice.n_links)
      // , ellstar(lattice.n_links)
    , link_volume(lattice.n_links)
      // , site_vol(lattice.n_sites)
  {
    set_sigma();
    set_face_signs();
    set_ell_and_link_volumes();
    // set_ell_ellstar();
    // set_site_vol();

    // std::cout << "debug.2.1" << std::endl;
    // // check
    // double TOL=1.0e-6;
    // {
    //   for(int ix=0; ix<lattice.n_sites; ix++){
    //     for(int jj=0; jj<lattice.sites[ix].nn; jj++){
    //       const int iy = lattice.sites[ix].neighbors[jj];

    //       std::cout << "ix, iy = " << ix << " " << iy << std::endl;
    //       const double alpha1 = alpha.at(Link{ix,iy});
    //       std::cout << "debug, pt1" << std::endl;
    //       double alpha2 = alpha.at(Link{iy,ix});
    //       std::cout << "debug, pt2" << std::endl;
    //       double omega12 = omega.at(Link{ix,iy});
    //       std::cout << "debug, pt3" << std::endl;

    //       double diff = (alpha2 + M_PI + omega12) - alpha1;
    //       // assert( std::abs(Mod(diff))<TOL );
    //     }}
    // }
    // std::cout << "debug.2.2" << std::endl;
    // {
    //   for(int ia=0; ia<lattice.n_faces; ia++){
    //     double sum = 0.0;

    //     for(int i=0; i<3; i++){
    //       int ix = lattice.faces[ia].sites[i];
    //       int iy = lattice.faces[ia].sites[(i+1)%3];
    //       sum -= omega.at(Link{ix,iy});
    //       sum += alpha.at(Link{ix,iy});
    //       sum -= alpha.at(Link{iy,ix}) + M_PI;
    //     }

    //     // std::cout << sum << std::endl;
    //     // assert( std::abs(Mod(-std::abs(Mod(sum)))) < TOL );
    //   }
    // }
  }

  DiracS2Simp & operator=(const DiracS2Simp&) = delete;

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


  MS gamma(const int ix, const int iy, const double shift=0.0) const { // located at x
    const double al = alpha.at(Link{ix,iy}) + shift;
    return std::cos(al)*sigma[1] + std::sin(al)*sigma[2];
  }

  MS Omega(const int ix, const int iy) const {
    const double om = omega.at(Link{ix,iy});
    return std::cos(0.5*om)*sigma[0] - I*std::sin(0.5*om)*sigma[3];
  }


  Eigen::MatrixXcd matrix_form() const {
    Eigen::MatrixXcd res = Eigen::MatrixXcd::Zero(NS*lattice.n_sites, NS*lattice.n_sites);

    for(int ix=0; ix<lattice.n_sites; ix++){
      for(int jj=0; jj<lattice.nn(ix); jj++){
	const int iy = lattice.nns[ix][jj];
	const int il = lattice.sites[ix].links[jj];

        // wilson // BROWER ET AL.
	res.block<NS,NS>(NS*ix,NS*iy) += 0.5/a * (link_volume[il]/ell[il]) * gamma(ix, iy) * Omega(ix, iy);
	res.block<NS,NS>(NS*ix,NS*ix) += 0.5 * r * (link_volume[il]/(ell[il]*ell[il])) * sigma[0];
	res.block<NS,NS>(NS*ix,NS*iy) -= 0.5 * r * (link_volume[il]/(ell[il]*ell[il])) * Omega(ix, iy);
      }
    }

    return res;
  } // end matrix_form


  void coo_format( std::vector<Complex>& v,
		   const Gauge& U ) const {
    const Idx N = lattice.n_sites * NS;
    for(Idx i=0; i<N; i++) v[i] = 0.0;

    Idx counter=0;
// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL)
// #endif
    for(Idx ix=0; ix<lattice.n_sites; ix++){
      // Idx counter = 3*8*ix;
      for(int jj=0; jj<lattice.nn(ix); jj++){
	const Idx iy = lattice.nns[ix][jj];
        const int il = lattice.sites[ix].links[jj];

        const MS tmp = 0.5/a * (link_volume[il]/ell[il]) * gamma(ix, iy) * Omega(ix, iy) - 0.5 * r * (link_volume[il]/(ell[il]*ell[il])) * Omega(ix, iy);
	const MS tmp2 = 0.5 * r * ( link_volume[il]/(ell[il]*ell[il]) ) * sigma[0];

	// res[NS*ix] += -tmp(0,0)*v[NS*iy] - tmp(0,1)*v[NS*iy+1];
	v[counter] = tmp(0,0); counter++;
	v[counter] = tmp(0,1); counter++;

	// res[NS*ix] += tmp(0,0)*v[NS*ix] + tmp(0,1)*v[NS*ix+1];
	v[counter] = tmp2(0,0); counter++;
	v[counter] = tmp2(0,1); counter++;

	// res[NS*ix+1] += -tmp(1,0)*v[NS*iy] - tmp(1,1)*v[NS*iy+1];
	v[counter] = tmp(1,0); counter++;
	v[counter] = tmp(1,1); counter++;

	// res[NS*ix+1] += tmp(1,0)*v[NS*ix] + tmp(1,1)*v[NS*ix+1];
	v[counter] = tmp2(1,0); counter++;
	v[counter] = tmp2(1,1); counter++;
      }
    }
  }




  Vec3 circumcenter(const Vec3& r0, const Vec3& r1, const Vec3& r2){
    const Vec3 r10 = r1 - r0;
    const Vec3 r20 = r2 - r0;

    const Vec3 tmp1 = r10.squaredNorm() * r20 - r20.squaredNorm() * r10;
    const Vec3 cross = r10.cross(r20);
    const Vec3 numer = tmp1.cross(cross);
    const double denom = 2.0*cross.squaredNorm();

    return numer/denom + r0;
  }



  void set_ell_and_link_volumes() {
    for(int il=0; il<lattice.n_links; il++) {
      const auto link = lattice.links[il];
      const int iA = link.faces[0];
      const int iB = link.faces[1];

      double ellA=0.0, ellB=0.0;
      double areaA=0.0, areaB=0.0;
      {
	const QfeFace& face = lattice.faces[iA];
	Vec3 r0, r1, r2; // r0,1: link
	if(face.sites[0]==link.sites[0] && face.sites[1]==link.sites[1]){
	  r0 = lattice.r[face.sites[0]];
	  r1 = lattice.r[face.sites[1]];
	  r2 = lattice.r[face.sites[2]];
	}
	else if(face.sites[1]==link.sites[0] && face.sites[2]==link.sites[1]){
	  r0 = lattice.r[face.sites[1]];
	  r1 = lattice.r[face.sites[2]];
	  r2 = lattice.r[face.sites[0]];
	}
	else if(face.sites[2]==link.sites[0] && face.sites[0]==link.sites[1]){
	  r0 = lattice.r[face.sites[2]];
	  r1 = lattice.r[face.sites[0]];
	  r2 = lattice.r[face.sites[1]];
	} // reverse
	else if(face.sites[1]==link.sites[0] && face.sites[0]==link.sites[1]){
	  r0 = lattice.r[face.sites[1]];
	  r1 = lattice.r[face.sites[0]];
	  r2 = lattice.r[face.sites[2]];
	}
	else if(face.sites[2]==link.sites[0] && face.sites[1]==link.sites[1]){
	  r0 = lattice.r[face.sites[2]];
	  r1 = lattice.r[face.sites[1]];
	  r2 = lattice.r[face.sites[0]];
	}
	else if(face.sites[0]==link.sites[0] && face.sites[2]==link.sites[1]){
	  r0 = lattice.r[face.sites[0]];
	  r1 = lattice.r[face.sites[2]];
	  r2 = lattice.r[face.sites[1]];
	}
	else assert(false);

	//
	const Vec3 p = circumcenter(r0, r1, r2).transpose();
	assert( std::abs( (p-r0).norm() - (p-r1).norm() )<1.0e-14 );
	assert( std::abs( (p-r0).norm() - (p-r2).norm() )<1.0e-14 );

	// double a_ = (r0-p).norm();
	// double b_ = (r1-p).norm();
	// double c_ = (r0-r1).norm(); // ell

	// double s_ = 0.5*(a_+b_+c_);
	// double tmp = s_ * (s_-a_) * (s_-b_) * (s_-c_);
	// double area_ = std::sqrt( tmp );

        double a_ = std::acos( r0.dot(p) /(r0.norm()* p.norm()) );
	double b_ = std::acos( r1.dot(p) /(r1.norm()* p.norm()) );
	double c_ = std::acos( r0.dot(r1)/(r0.norm()*r1.norm()) ); // ell

	double s_ = 0.5*(a_+b_+c_);
	double tmp = std::tan(0.5*s_) * std::tan(0.5*(s_-a_)) * std::tan(0.5*(s_-b_)) * std::tan(0.5*(s_-c_));
	double area_ = 4.0*std::atan( std::sqrt( tmp ) );

	ellA = c_;
	areaA = area_;
      }
      {
	const QfeFace& face = lattice.faces[iB];
	Vec3 r0, r1, r2; // r0,1: link
	if(face.sites[0]==link.sites[0] && face.sites[1]==link.sites[1]){
	  r0 = lattice.r[face.sites[0]];
	  r1 = lattice.r[face.sites[1]];
	  r2 = lattice.r[face.sites[2]];
	}
	else if(face.sites[1]==link.sites[0] && face.sites[2]==link.sites[1]){
	  r0 = lattice.r[face.sites[1]];
	  r1 = lattice.r[face.sites[2]];
	  r2 = lattice.r[face.sites[0]];
	}
	else if(face.sites[2]==link.sites[0] && face.sites[0]==link.sites[1]){
	  r0 = lattice.r[face.sites[2]];
	  r1 = lattice.r[face.sites[0]];
	  r2 = lattice.r[face.sites[1]];
	} // reverse
	else if(face.sites[1]==link.sites[0] && face.sites[0]==link.sites[1]){
	  r0 = lattice.r[face.sites[1]];
	  r1 = lattice.r[face.sites[0]];
	  r2 = lattice.r[face.sites[2]];
	}
	else if(face.sites[2]==link.sites[0] && face.sites[1]==link.sites[1]){
	  r0 = lattice.r[face.sites[2]];
	  r1 = lattice.r[face.sites[1]];
	  r2 = lattice.r[face.sites[0]];
	}
	else if(face.sites[0]==link.sites[0] && face.sites[2]==link.sites[1]){
	  r0 = lattice.r[face.sites[0]];
	  r1 = lattice.r[face.sites[2]];
	  r2 = lattice.r[face.sites[1]];
	}
	else assert(false);

	const Vec3 p = circumcenter(r0, r1, r2).transpose();
	assert( std::abs( (p-r0).norm() - (p-r1).norm() )<1.0e-14 );
	assert( std::abs( (p-r0).norm() - (p-r2).norm() )<1.0e-14 );

	double a_ = std::acos( r0.dot(p) /(r0.norm()* p.norm()) );
	double b_ = std::acos( r1.dot(p) /(r1.norm()* p.norm()) );
	double c_ = std::acos( r0.dot(r1)/(r0.norm()*r1.norm()) ); // ell

	double s_ = 0.5*(a_+b_+c_);
	double tmp = std::tan(0.5*s_) * std::tan(0.5*(s_-a_)) * std::tan(0.5*(s_-b_)) * std::tan(0.5*(s_-c_));
	double area_ = 4.0*std::atan( std::sqrt( tmp ) );
	// double a_ = (r0-p).norm();
	// double b_ = (r1-p).norm();
	// double c_ = (r0-r1).norm(); // ell

	// double s_ = 0.5*(a_+b_+c_);
	// double tmp = s_ * (s_-a_) * (s_-b_) * (s_-c_);
	// double area_ = std::sqrt( tmp );

	ellB = c_;
	areaB = area_;
      }

      assert( std::abs(ellA-ellB)<1.0e-14 );
      ell[il] = ellA;
      link_volume[il] = areaA + areaB;
    }

    int counter = 0;
    a = 0.0;
    for(int il=0; il<lattice.n_links; il++) {
      a += ell[il];
      counter++;
      // std::cout << "ell[il] =" << ell[il] << std::endl;
    }
    a /= counter;
    a *= 1.0;
    // std::cout << "a = " << a << std::endl;
  }





  // void set_ell_and_link_volumes() { // geodesic
  //   for(int il=0; il<lattice.n_links; il++) {
  //     const auto link = lattice.links[il];
  //     const int iA = link.faces[0];
  //     const int iB = link.faces[1];

  //     double ellA=0.0, ellB=0.0;
  //     double areaA=0.0, areaB=0.0;
  //     {
  // 	const QfeFace& face = lattice.faces[iA];
  // 	Vec3 r0, r1, r2; // r0,1: link
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
  // 	const Vec3 p = circumcenter(r0, r1, r2).transpose();
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
  //     }
  //     {
  // 	const QfeFace& face = lattice.faces[iB];
  // 	Vec3 r0, r1, r2; // r0,1: link
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

  // 	const Vec3 p = circumcenter(r0, r1, r2).transpose();
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
  //     }

  //     assert( std::abs(ellA-ellB)<1.0e-14 );
  //     ell[il] = ellA;
  //     link_volume[il] = areaA + areaB;
  //   }

  //   int counter = 0;
  //   a = 0.0;
  //   for(int il=0; il<lattice.n_links; il++) {
  //     a += ell[il];
  //     counter++;
  //     std::cout << "ell[il] =" << ell[il] << std::endl;
  //   }
  //   a /= counter;
  //   a *= 1.0;
  //   std::cout << "a = " << a << std::endl;
  // }
























  // void set_site_vol(){
  //   for(int i=0; i<lattice.n_sites; i++){
  //     site_vol[i] = 0.0;
  //     const auto x = lattice.sites[i];
  //     //for(const int il : x.links){
  //     for(int jj=0; jj<x.nn; jj++){
  // 	const int il = x.links[jj];
  // 	site_vol[i] += 0.25*ell[il]*ellstar[il];
  //     }
  //   }
  // }
  


};
