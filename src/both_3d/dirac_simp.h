#pragma once

#include <array>
#include <cmath>

#include "s2n_simp.h"
#include "dirac_base.h"

// using Complex = std::complex<double>;

// using VD=Eigen::Vector2d;


// double Mod(double a, double b=2.0*M_PI){
//   int p = int(std::round(a / b));
//   double r = a - p*b;
//   return r;
// }




struct SpinStructureSimp{
  using Link = std::array<Idx,2>; // <int,int>;


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




struct DiracS2Simp : public DiracBase, public SpinStructureSimp{
  using Lattice=S2Simp;
  using Link = std::array<Idx,2>; // <int,int>;
  using Face = std::vector<Idx>;

  // using MS=Eigen::Matrix2cd;
  // using VD=Eigen::Vector2d;
  // using VE=Eigen::Vector3d;
  // using VC=Eigen::VectorXcd;

  // static constexpr int NS = 2;
  // static constexpr int DIM = 2;
  // static constexpr Complex I = Complex(0.0, 1.0);

  Lattice& lattice;

  // sign for the ordering of Evan's face.sites; +1 for clockwise rotation from the origin
  // std::vector<Idx> face_signs; // index: ia (Evan's label for faces)

  const double m;
  const double r;
  const double M5;

  // double a = 1.0;

  // SpinStructureSimp spin;
  // const std::map<const Link, const double> omega;
  // const std::map<const Link, const double> alpha;
  std::vector<double> kappa; // evan's link label

  DiracS2Simp()=delete;

  DiracS2Simp(Lattice& lattice_,
	      const double m_=0.0,
	      const double r_=1.0,
              const double M5_=0.0 )
    : SpinStructureSimp(Comp::N_REFINE)
    , lattice(lattice_)
    , m(m_)
    , r(r_)
    , M5(M5_)
    , kappa(lattice.n_links)
  {
    // set_sigma();
    set_kappa();
  }

  DiracS2Simp & operator=(const DiracS2Simp&) = delete;

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


  // void set_face_signs() { for(Idx i=0; i<lattice.n_faces; i++) face_signs[i] = face_sign(i); }

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


  void coo_structure( std::vector<Idx>& is,
                      std::vector<Idx>& js ) const {
    for(Idx ix=0; ix<lattice.n_sites; ix++){
      for(const Idx iy : lattice.nns[ix]){
	is.push_back( NS*ix ); js.push_back( NS*iy );
	is.push_back( NS*ix ); js.push_back( NS*iy+1 );
	is.push_back( NS*ix ); js.push_back( NS*ix );
	is.push_back( NS*ix ); js.push_back( NS*ix+1 );
	is.push_back( NS*ix+1 ); js.push_back( NS*iy );
	is.push_back( NS*ix+1 ); js.push_back( NS*iy+1 );
	is.push_back( NS*ix+1 ); js.push_back( NS*ix );
	is.push_back( NS*ix+1 ); js.push_back( NS*ix+1 );
      }
    }
    assert( is.size()==js.size() );
  }

  template<typename Gauge>
  void coo_format( std::vector<Complex>& v,
		   const Gauge& u ) const {
    const Idx N = lattice.n_sites * NS;
    for(Idx i=0; i<N; i++) v[i] = 0.0;

    // Idx counter=0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL)
#endif
    for(Idx ix=0; ix<lattice.n_sites; ix++){
      Idx counter = 8*lattice.counter_accum[ix];
      for(const Idx iy : lattice.nns[ix]){
        const Idx il = lattice.map2il.at(Link{ix,iy});

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



  template<typename Gauge>
  void d_coo_format( std::vector<COOEntry>& elem,
        	     const Gauge& u,
        	     const Link& el ) const {
    const Idx ix0 = el[0];
    const Idx iy0 = el[1];

    elem.clear();
    {
      // pos
      const Idx ix = ix0;
      for(int jj=0; jj<lattice.nns[ix].size(); jj++){
        const Idx iy = lattice.nns[ix][jj];
        if(iy!=iy0) continue;
        const Idx il = lattice.map2il.at(Link{ix,iy});
        const MS tmp = 0.5 * kappa[il] * ( -r *sigma[0] + gamma(ix, iy) ) * I*std::exp( I* u(Link{ix,iy})) * Omega(ix, iy);

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
      for(int jj=0; jj<lattice.nns[iy].size(); jj++){
        const Idx ix = lattice.nns[iy0][jj];
        if(ix!=ix0) continue;
        const Idx il = lattice.map2il.at(Link{ix,iy});
        const MS tmp = -0.5 * kappa[il] * ( -r *sigma[0] + gamma(iy, ix) ) * I*std::exp( I* u(Link{iy,ix})) * Omega(iy, ix);

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
    kappa.clear();
    kappa.resize(lattice.n_links);
    for(Idx il=0; il<lattice.n_links; il++) {
      kappa[il] = lattice.link_volume[il] / lattice.mean_link_volume * lattice.mean_ell / lattice.ell[il];
      // kappa[il] = lattice.link_volume[il] / lattice.ell[il] / (lattice.alat/std::sqrt(3.0));
    }
  }


};
