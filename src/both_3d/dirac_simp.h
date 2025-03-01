#pragma once

#include <array>
#include <cmath>
#include <map>
#include <Eigen/Dense>

#include "s2n_simp.h"

// using Complex = std::complex<double>;

// using VD=Eigen::Vector2d;


// double Mod(double a, double b=2.0*M_PI){
//   int p = int(std::round(a / b));
//   double r = a - p*b;
//   return r;
// }




struct SpinStructureSimp{

  using Link = std::array<Idx,2>; // <int,int>;
  using Face = std::vector<Idx>;

  std::map<const Link, const double> omega;
  std::map<const Link, const double> alpha;
  // std::map<const int, const int> NM2EO;


  SpinStructureSimp(const int n_refine)
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
struct DiracS2Simp : public SpinStructureSimp{
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
  std::vector<Idx> face_signs; // index: ia (Evan's label for faces)

  const double m;
  const double r;
  const double M5;

  // double a = 1.0;

  using Link = std::array<Idx,2>; // <int,int>;
  // SpinStructureSimp spin;
  // const std::map<const Link, const double> omega;
  // const std::map<const Link, const double> alpha;

  std::array<MS, 4> sigma;

  // std::vector<double> ell; // evan's link label
  // std::vector<double> ellstar; // evan's link label
  // std::vector<double> link_volume; // evan's link label
  // std::vector<double> site_vol; // evan's site label
  std::vector<double> kappa; // evan's link label

  DiracS2Simp()=delete;

  DiracS2Simp(Lattice& lattice_,
	      const double m_=0.0,
	      const double r_=1.0,
              const double M5_=0.0 )
    : SpinStructureSimp(Comp::N_REFINE)
    , lattice(lattice_)
    , face_signs(lattice.n_faces)
    , m(m_)
    , r(r_)
    , M5(M5_)
      // , ell(lattice.n_links)
      // , link_volume(lattice.n_links)
    , kappa(lattice.n_links)
  {
    set_sigma();
    set_face_signs();
    // set_ell_and_link_volumes();
    set_kappa();
  }

  DiracS2Simp & operator=(const DiracS2Simp&) = delete;

  void set_sigma(){
    assert(NS==2);
    sigma[0] << 1,0,0,1;
    sigma[1] << 0,1,1,0;
    sigma[2] << 0,-I,I,0;
    sigma[3] << 1,0,0,-1;
  }


  int face_sign(const Idx i_face) const {
    const Face& face = lattice.faces[i_face];
    const VE r0 = lattice.sites[face[0]];
    const VE r1 = lattice.sites[face[1]];
    const VE r2 = lattice.sites[face[2]];

    const VE cross = (r1-r0).cross(r2-r0);
    const VE sum = r0+r1+r2;

    const double inner = cross.dot(sum);

    int res = 1;
    if(inner<0) res = -1;
    return res;
  }


  void set_face_signs() { for(Idx i=0; i<lattice.n_faces; i++) face_signs[i] = face_sign(i); }


  bool is_nn(const Idx ix, const Idx iy) const {
    if(ix==iy) return false;

    bool res = false;
    // const VE x = lattice.sites[ix];

    for(Idx iz : lattice.nns[ix]){
      if(iy==iz) {
	res = true;
	break;
      }
    }
    return res;
  }


  // +1 if (ix, iy, iz) is RH
  int sign(const Idx ix, const Idx iy, const Idx iz) const {
    const VE rx = lattice.sites[ix];
    const VE ry = lattice.sites[iy];
    const VE rz = lattice.sites[iz];

    const VE cross = (ry-rx).cross(rz-rx);
    const VE sum = rx+ry+rz;
    const double inner = cross.dot(sum);

    int res = 0;
    if(inner>0) res = 1;
    else if(inner<0) res = -1;
    else assert(false);
    return res;
  }


  MS gamma(const Idx ix, const Idx iy) const { // located at x
    const double al = alpha.at(Link{ix,iy});
    return std::cos(al)*sigma[1] + std::sin(al)*sigma[2];
  }

  MS Omega(const Idx ix, const Idx iy) const {
    const double om = omega.at(Link{ix,iy});
    return std::cos(0.5*om)*sigma[0] - I*std::sin(0.5*om)*sigma[3];
  }


  // Eigen::MatrixXcd matrix_form() const {
  //   Eigen::MatrixXcd res = Eigen::MatrixXcd::Zero(NS*lattice.n_sites, NS*lattice.n_sites);

  //   for(Idx ix=0; ix<lattice.n_sites; ix++){
  //     for(int jj=0; jj<lattice.nn(ix); jj++){
  //       const Idx iy = lattice.nns[ix][jj];
  //       const Idx il = lattice.sites[ix].links[jj];

  //       // wilson // BROWER ET AL.
  //       res.block<NS,NS>(NS*ix,NS*iy) += 0.5/a * (link_volume[il]/ell[il]) * gamma(ix, iy) * Omega(ix, iy);
  //       res.block<NS,NS>(NS*ix,NS*ix) += 0.5 * r * (link_volume[il]/(ell[il]*ell[il])) * sigma[0];
  //       res.block<NS,NS>(NS*ix,NS*iy) -= 0.5 * r * (link_volume[il]/(ell[il]*ell[il])) * Omega(ix, iy);
  //     }
  //   }

  //   return res;
  // } // end matrix_form


  void coo_format( std::vector<Complex>& v,
		   const Gauge& u ) const {
    const Idx N = lattice.n_sites * NS;
    for(Idx i=0; i<N; i++) v[i] = 0.0;

    // Idx counter=0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL)
#endif
    for(Idx ix=0; ix<lattice.n_sites; ix++){
      Idx counter = lattice.counter_accum[ix];
      for(const Idx iy : lattice.nns[ix]){
        const Idx il = lattice.map2il.at(Link{ix,iy});

        // const MS tmp = 0.5/a * (link_volume[il]/ell[il]) * gamma(ix, iy) * Omega(ix, iy) - 0.5 * r * (link_volume[il]/(ell[il]*ell[il])) * Omega(ix, iy);
	// const MS tmp2 = 0.5 * r * ( link_volume[il]/(ell[il]*ell[il]) ) * sigma[0] + M5/lattice.nn(ix) * sigma[0];

        // const MS tmp = kappa[il] * ( - lattice.alat/ell[il] * r *sigma[0] + gamma(ix, iy) ) * Omega(ix, iy);
	// const MS tmp2 = lattice.alat/ell[il] * r*kappa[il] * sigma[0] + M5/lattice.nn(ix) * sigma[0];

        // const MS tmp = kappa[il] * ( -r *sigma[0] + gamma(ix, iy) ) * Omega(ix, iy);
	// const MS tmp2 = r*kappa[il] * sigma[0] + M5/lattice.nn(ix) * sigma[0];

        const MS tmp = 0.5 * kappa[il] * ( -r *sigma[0] + gamma(ix, iy) ) * std::exp( I* u(Link{ix,iy})) * Omega(ix, iy);
	const MS tmp2 = 0.5 * r*kappa[il] * sigma[0] + M5/lattice.nns[ix].size() * sigma[0];

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



  void d_coo_format( std::vector<COOEntry>& elem,
        	     const Gauge& u,
        	     const Link& el ) const {
    const Idx ix0 = el[0];
    const Idx iy0 = el[1];

    elem.clear();
    {
      // pos
      const Idx ix = ix0;
      // std::cout << "debug. pos. ix = " << ix << std::endl;
      for(int jj=0; jj<lattice.nns[ix].size(); jj++){
        const Idx iy = lattice.nns[ix][jj];
        if(iy!=iy0) continue;
        // std::cout << "debug. pos. iy = " << iy << std::endl;
        const Idx il = lattice.map2il.at(Link{ix,iy});
        const MS tmp = 0.5 * kappa[il] * ( -r *sigma[0] + gamma(ix, iy) ) * I*std::exp( I* u(Link{ix,iy})) * Omega(ix, iy);
        // std::cout << "debug. tmp = " << tmp << std::endl;

        // res[NS*ix] += -tmp(0,0)*v[NS*iy] - tmp(0,1)*v[NS*iy+1];
        elem.push_back(COOEntry(tmp(0,0),NS*ix,NS*iy));
        elem.push_back(COOEntry(tmp(0,1),NS*ix,NS*iy+1));

        // res[NS*ix+1] += -tmp(1,0)*v[NS*iy] - tmp(1,1)*v[NS*iy+1];
        elem.push_back(COOEntry(tmp(1,0),NS*ix+1,NS*iy));
        elem.push_back(COOEntry(tmp(1,1),NS*ix+1,NS*iy+1));
      }
    }

    {
      // neg
      const Idx iy = iy0;
      // std::cout << "debug. neg. iy = " << iy << std::endl;
      for(int jj=0; jj<lattice.nns[iy].size(); jj++){
        const Idx ix = lattice.nns[iy0][jj];
        if(ix!=ix0) continue;
        // std::cout << "debug. neg. ix = " << ix << std::endl;
        const Idx il = lattice.map2il.at(Link{ix,iy});
        const MS tmp = -0.5 * kappa[il] * ( -r *sigma[0] + gamma(iy, ix) ) * I*std::exp( I* u(Link{iy,ix})) * Omega(iy, ix);
        // std::cout << "debug. tmp = " << tmp << std::endl;

        // res[NS*iy] += -tmp(0,0)*v[NS*ix] - tmp(0,1)*v[NS*ix+1];
        elem.push_back(COOEntry(tmp(0,0),NS*iy,NS*ix));
        elem.push_back(COOEntry(tmp(0,1),NS*iy,NS*ix+1));

        // res[NS*iy+1] += -tmp(1,0)*v[NS*ix] - tmp(1,1)*v[NS*ix+1];
        elem.push_back(COOEntry(tmp(1,0),NS*iy+1,NS*ix));
        elem.push_back(COOEntry(tmp(1,1),NS*iy+1,NS*ix+1));
      }
    }
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



  void set_kappa() {
    for(Idx il=0; il<lattice.n_links; il++) {
      kappa[il] = lattice.link_volume[il]/lattice.ell[il]/(lattice.alat/std::sqrt(3.0));
    }
    // std::cout << "a = " << a << std::endl;
  }

  // void set_kappa() {
  //   for(Idx il=0; il<lattice.n_links; il++) {
  //     const Link link = lattice.links[il];
  //     const Idx ix = lattice.links[il][0];
  //     const Idx iy = lattice.links[il][1];

  //     const Idx iA = lattice.dual_links[il][0];
  //     const Idx iB = lattice.dual_links[il][1];

  //     // std::cout << "debug. il = " << il << std::endl
  //     //           << "debug. link = " << link[0] << " " << link[1] << std::endl
  //     //           << "debug. iA = " << iA << std::endl
  //     //           << "debug. iB = " << iB << std::endl;

  //     double ellA=0.0, ellB=0.0;
  //     double areaA=0.0, areaB=0.0;

  //     const VE x = lattice.sites[ix];
  //     const VE y = lattice.sites[iy];
  //     {
  //       const VE p = lattice.dual_sites[iA];

  //       // double a_ = (x-p).norm();
  //       // double b_ = (y-p).norm();
  //       // double c_ = (x-y).norm(); // ell

  //       // double s_ = 0.5*(a_+b_+c_);
  //       // double tmp = s_ * (s_-a_) * (s_-b_) * (s_-c_);
  //       // double area_ = std::sqrt( tmp );

  //       double a_ = std::acos( x.dot(p) /(x.norm()* p.norm()) );
  //       double b_ = std::acos( y.dot(p) /(y.norm()* p.norm()) );
  //       double c_ = std::acos( x.dot(y)/(x.norm()*y.norm()) ); // ell

  //       double s_ = 0.5*(a_+b_+c_);
  //       double tmp = std::tan(0.5*s_) * std::tan(0.5*(s_-a_)) * std::tan(0.5*(s_-b_)) * std::tan(0.5*(s_-c_));
  //       double area_ = 4.0*std::atan( std::sqrt( tmp ) );

  //       ellA = c_;
  //       areaA = area_;
  //     }
  //     {
  //       const VE p = lattice.dual_sites[iB];

  //       // double a_ = (x-p).norm();
  //       // double b_ = (y-p).norm();
  //       // double c_ = (x-y).norm(); // ell

  //       // double s_ = 0.5*(a_+b_+c_);
  //       // double tmp = s_ * (s_-a_) * (s_-b_) * (s_-c_);
  //       // double area_ = std::sqrt( tmp );

  //       double a_ = std::acos( x.dot(p) /(x.norm()* p.norm()) );
  //       double b_ = std::acos( y.dot(p) /(y.norm()* p.norm()) );
  //       double c_ = std::acos( x.dot(y)/(x.norm()*y.norm()) ); // ell

  //       double s_ = 0.5*(a_+b_+c_);
  //       double tmp = std::tan(0.5*s_) * std::tan(0.5*(s_-a_)) * std::tan(0.5*(s_-b_)) * std::tan(0.5*(s_-c_));
  //       double area_ = 4.0*std::atan( std::sqrt( tmp ) );

  //       ellB = c_;
  //       areaB = area_;
  //     }
  //     // {
  //     //   const Face& face = lattice.faces[iA];
  //     //   std::cout << "face " << iA << std::endl;
  //     //   for(Idx elem : face) std::cout << elem << std::endl;
  //     //   VE r0, r1, r2; // r0,1: link
  //     //   if(face[0]==link[0] && face[1]==link[1]){
  //     //     r0 = lattice.sites[face[0]];
  //     //     r1 = lattice.sites[face[1]];
  //     //     r2 = lattice.sites[face[2]];
  //     //   }
  //     //   else if(face[1]==link[0] && face[2]==link[1]){
  //     //     r0 = lattice.sites[face[1]];
  //     //     r1 = lattice.sites[face[2]];
  //     //     r2 = lattice.sites[face[0]];
  //     //   }
  //     //   else if(face[2]==link[0] && face[0]==link[1]){
  //     //     r0 = lattice.sites[face[2]];
  //     //     r1 = lattice.sites[face[0]];
  //     //     r2 = lattice.sites[face[1]];
  //     //   } // reverse
  //     //   else if(face[1]==link[0] && face[0]==link[1]){
  //     //     r0 = lattice.sites[face[1]];
  //     //     r1 = lattice.sites[face[0]];
  //     //     r2 = lattice.sites[face[2]];
  //     //   }
  //     //   else if(face[2]==link[0] && face[1]==link[1]){
  //     //     r0 = lattice.sites[face[2]];
  //     //     r1 = lattice.sites[face[1]];
  //     //     r2 = lattice.sites[face[0]];
  //     //   }
  //     //   else if(face[0]==link[0] && face[2]==link[1]){
  //     //     r0 = lattice.sites[face[0]];
  //     //     r1 = lattice.sites[face[2]];
  //     //     r2 = lattice.sites[face[1]];
  //     //   }
  //     //   else assert(false);

  //     //   //
  //     //   const VE p = circumcenter(r0, r1, r2).transpose();
  //     //   assert( std::abs( (p-r0).norm() - (p-r1).norm() )<1.0e-14 );
  //     //   assert( std::abs( (p-r0).norm() - (p-r2).norm() )<1.0e-14 );

  //     //   // double a_ = (r0-p).norm();
  //     //   // double b_ = (r1-p).norm();
  //     //   // double c_ = (r0-r1).norm(); // ell

  //     //   // double s_ = 0.5*(a_+b_+c_);
  //     //   // double tmp = s_ * (s_-a_) * (s_-b_) * (s_-c_);
  //     //   // double area_ = std::sqrt( tmp );

  //     //   double a_ = std::acos( r0.dot(p) /(r0.norm()* p.norm()) );
  //     //   double b_ = std::acos( r1.dot(p) /(r1.norm()* p.norm()) );
  //     //   double c_ = std::acos( r0.dot(r1)/(r0.norm()*r1.norm()) ); // ell

  //     //   double s_ = 0.5*(a_+b_+c_);
  //     //   double tmp = std::tan(0.5*s_) * std::tan(0.5*(s_-a_)) * std::tan(0.5*(s_-b_)) * std::tan(0.5*(s_-c_));
  //     //   double area_ = 4.0*std::atan( std::sqrt( tmp ) );

  //     //   ellA = c_;
  //     //   areaA = area_;
  //     // }
  //     // {
  //     //   const Face& face = lattice.faces[iB];
  //     //   std::cout << "face " << iB << std::endl;
  //     //   for(Idx elem : face) std::cout << elem << std::endl;

  //     //   VE r0, r1, r2; // r0,1: link
  //     //   if(face[0]==link[0] && face[1]==link[1]){
  //     //     r0 = lattice.sites[face[0]];
  //     //     r1 = lattice.sites[face[1]];
  //     //     r2 = lattice.sites[face[2]];
  //     //   }
  //     //   else if(face[1]==link[0] && face[2]==link[1]){
  //     //     r0 = lattice.sites[face[1]];
  //     //     r1 = lattice.sites[face[2]];
  //     //     r2 = lattice.sites[face[0]];
  //     //   }
  //     //   else if(face[2]==link[0] && face[0]==link[1]){
  //     //     r0 = lattice.sites[face[2]];
  //     //     r1 = lattice.sites[face[0]];
  //     //     r2 = lattice.sites[face[1]];
  //     //   } // reverse
  //     //   else if(face[1]==link[0] && face[0]==link[1]){
  //     //     r0 = lattice.sites[face[1]];
  //     //     r1 = lattice.sites[face[0]];
  //     //     r2 = lattice.sites[face[2]];
  //     //   }
  //     //   else if(face[2]==link[0] && face[1]==link[1]){
  //     //     r0 = lattice.sites[face[2]];
  //     //     r1 = lattice.sites[face[1]];
  //     //     r2 = lattice.sites[face[0]];
  //     //   }
  //     //   else if(face[0]==link[0] && face[2]==link[1]){
  //     //     r0 = lattice.sites[face[0]];
  //     //     r1 = lattice.sites[face[2]];
  //     //     r2 = lattice.sites[face[1]];
  //     //   }
  //     //   else assert(false);

  //     //   const VE p = circumcenter(r0, r1, r2).transpose();
  //     //   assert( std::abs( (p-r0).norm() - (p-r1).norm() )<1.0e-14 );
  //     //   assert( std::abs( (p-r0).norm() - (p-r2).norm() )<1.0e-14 );

  //     //   double a_ = std::acos( r0.dot(p) /(r0.norm()* p.norm()) );
  //     //   double b_ = std::acos( r1.dot(p) /(r1.norm()* p.norm()) );
  //     //   double c_ = std::acos( r0.dot(r1)/(r0.norm()*r1.norm()) ); // ell

  //     //   double s_ = 0.5*(a_+b_+c_);
  //     //   double tmp = std::tan(0.5*s_) * std::tan(0.5*(s_-a_)) * std::tan(0.5*(s_-b_)) * std::tan(0.5*(s_-c_));
  //     //   double area_ = 4.0*std::atan( std::sqrt( tmp ) );

  //     //   // double a_ = (r0-p).norm();
  //     //   // double b_ = (r1-p).norm();
  //     //   // double c_ = (r0-r1).norm(); // ell

  //     //   // double s_ = 0.5*(a_+b_+c_);
  //     //   // double tmp = s_ * (s_-a_) * (s_-b_) * (s_-c_);
  //     //   // double area_ = std::sqrt( tmp );

  //     //   ellB = c_;
  //     //   areaB = area_;
  //     // }

  //     assert( std::abs(ellA-ellB)<1.0e-14 );
  //     ell[il] = ellA;
  //     link_volume[il] = areaA + areaB;
  //     // kappa[il] = link_volume[il]/ell[il]/a

  //     // VE rA = lattice.sites[iA];
  //     // VE rB = lattice.sites[iB];
  //     // kappa[il] = std::acos( rA.dot(rB) /(rA.norm() * rB.norm()) );
  //   }

  //   // Idx counter = 0;
  //   // a = 0.0;
  //   // for(Idx il=0; il<lattice.n_links; il++) {
  //   //   a += ell[il];
  //   //   counter++;
  //   //   // std::cout << "ell[il] =" << ell[il] << std::endl;
  //   // }
  //   // a /= counter;
  //   // a *= 1.0;

  //   for(Idx il=0; il<lattice.n_links; il++) {
  //     kappa[il] = link_volume[il]/ell[il]/(lattice.alat/std::sqrt(3.0));
  //     // kappa[il] /= a;
  //   }
  //   // std::cout << "a = " << a << std::endl;
  // }



};
