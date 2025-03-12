#pragma once

#include <cmath>
#include <map>
#include <complex>


double Mod(double a, double b=2.0*M_PI){
  int p = int(std::round(a / b));
  double r = a - p*b;
  return r;
}


// std::string dir = "/mnt/hdd_barracuda/qed3/dats/";
// std::string dir = "/mnt/hdd_barracuda/qed3/dats_saved2/";


struct SpinStructureDual{
  // using Link = std::array<Idx,2>; // <Idx,Idx>;
  // using Face = std::vector<Idx>;

  // static constexpr int NS = 2;
  // static constexpr int DIM = 2;
  // static constexpr int EDIM = 3;

  std::map<const Link, const double> omega;
  std::map<const Link, const double> alpha;

  SpinStructureDual(const int n_refine)
  {
    {
      std::cout << "# reading omega" << std::endl;
      std::ifstream file(dir+"omega_dual_n"+std::to_string(n_refine)+"_singlepatch.dat");
      assert(file.is_open());

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
      std::cout << "# reading alpha" << std::endl;
      std::ifstream file(dir+"alpha_dual_n"+std::to_string(n_refine)+"_singlepatch.dat");
      assert(file.is_open());

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


struct DiracS2Dual : public DiracBase, public SpinStructureDual{
  using Lattice=S2Trivalent;

  // static constexpr int NS = 2;
  // static constexpr int DIM = 2;
  // static constexpr Complex I = Complex(0.0, 1.0);

  const Lattice& lattice;

  const double m;
  const double r;
  const double M5;

  // std::array<MS, 4> sigma;
  std::vector<double> kappa; // link label

  DiracS2Dual()=delete;

  DiracS2Dual(const Lattice& lattice_,
	      const double m_=0.0,
	      const double r_=1.0,
              const double M5_=0.0
              )
    : SpinStructureDual(Comp::N_REFINE)
    , lattice(lattice_)
    , m(m_)
    , r(r_)
    , M5(M5_)
    , kappa(lattice.n_links)
  {
    double TOL=1.0e-6;
    {
      std::cout << "# checking spin structure" << std::endl;
      for(Idx ix=0; ix<lattice.n_sites; ix++){
	for(Idx iy : lattice.nns[ix]){
	  const double alpha1 = alpha.at(Link{ix,iy});
	  double alpha2 = alpha.at(Link{iy,ix});
	  double omega12 = omega.at(Link{ix,iy});

	  double diff = (alpha2 + M_PI + omega12) - alpha1;
	  assert( std::abs(Mod(diff))<TOL );
	}}
    }
    set_kappa();
  }

  DiracS2Dual & operator=(const DiracS2Dual&) = delete;

  // void set_sigma(){
  //   assert(NS==2);
  //   sigma[0] << 1,0,0,1;
  //   sigma[1] << 0,1,1,0;
  //   sigma[2] << 0,-I,I,0;
  //   sigma[3] << 1,0,0,-1;
  // }

  // MS gamma(const Idx ix, const int jj) const { // located at x
  //   const Idx iy = lattice.nns[ix][jj];
  //   const double al = alpha.at(Link{ix,iy});
  //   return std::cos(al)*sigma[1] + std::sin(al)*sigma[2];
  //   // return lattice.v[ix][jj](0) * sigma[1] + lattice.v[ix][jj](1) * sigma[2];
  // }

  MS gamma(const Idx ix, const Idx iy) const { // located at x
    const double al = alpha.at(Link{ix,iy});
    return std::cos(al)*sigma[1] + std::sin(al)*sigma[2];
    // return lattice.v[ix][jj](0) * sigma[1] + lattice.v[ix][jj](1) * sigma[2];
  }

  MS Omega(const Idx ix, const Idx iy) const {
    const double om = omega.at(Link{ix,iy});
    return std::cos(0.5*om)*sigma[0] - I*std::sin(0.5*om)*sigma[3];
  }

  // Eigen::MatrixXcd matrix_form() const {
  //   Eigen::MatrixXcd res = Eigen::MatrixXcd::Zero(NS*lattice.n_sites, NS*lattice.n_sites);

  //   for(Idx ix=0; ix<lattice.n_sites; ix++){
  //     for(int jj=0; jj<3; jj++){
  //       Idx iy = lattice.nns[ix][jj];
  //       // res.block<NS,NS>(NS*ix,NS*iy) -= lattice.vol[ix]/lattice.mean_vol * (r*sigma[0] - gamma(ix, jj)) * Omega(ix, iy);
  //       // res.block<NS,NS>(NS*ix,NS*ix) += lattice.vol[ix]/lattice.mean_vol * (r + m/3.0)*sigma[0];
  //       // res.block<NS,NS>(NS*ix,NS*iy) -= lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj]*sigma[0] - gamma(ix, jj)) * Omega(ix, iy);
  //       // res.block<NS,NS>(NS*ix,NS*ix) += lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj] + m/3.0)*sigma[0];
  //       res.block<NS,NS>(NS*ix,NS*ix) += M5/3.0*sigma[0];
  //     }
  //   }

  //   return res;
  // } // end matrix_form


  // Eigen::MatrixXcd matrix_form( const Gauge& U ) const {
  //   Eigen::MatrixXcd res = Eigen::MatrixXcd::Zero(NS*lattice.n_sites, NS*lattice.n_sites);

  //   for(Idx ix=0; ix<lattice.n_sites; ix++){
  //     for(int jj=0; jj<3; jj++){
  //       Idx iy = lattice.nns[ix][jj];
  //       res.block<NS,NS>(NS*ix,NS*iy) -= lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj]*sigma[0] - gamma(ix, jj)) * std::exp( I* U(Link{ix,iy})) * Omega(ix, iy);
  //       res.block<NS,NS>(NS*ix,NS*ix) += lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj] + m/3.0)*sigma[0];
  //       res.block<NS,NS>(NS*ix,NS*ix) += M5/3.0*sigma[0];
  //     }
  //   }

  //   return res;
  // } // end matrix_form



  template<typename Gauge>
  void coo_format( std::vector<Complex>& v,
		   const Gauge& U ) const {
    const Idx N = lattice.n_sites * NS;
    for(Idx i=0; i<N; i++) v[i] = 0.0;

    // int counter=0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL)
#endif
    for(Idx ix=0; ix<lattice.n_sites; ix++){
      Idx counter = 3*8*ix;
      for(int jj=0; jj<3; jj++){
	const Idx iy = lattice.nns[ix][jj];
        const Idx il = lattice.map2il.at(Link{ix,iy});

        const MS tmp = 0.5 * kappa[il] * ( -r *sigma[0] + gamma(ix, iy) ) * Omega(ix, iy);
	const MS tmp2 = 0.5 * r*kappa[il] * sigma[0] + M5/lattice.nn(ix) * sigma[0];

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

//   void H_coo_format( Complex* v,
// 		     const Gauge& U,
// 		     const double lambda_max=1.0 ) const {
//     const Idx N = lattice.n_sites * NS;
//     // for(Idx i=0; i<N; i++) v[i] = 0.0;

//     // Idx counter=0;
// #ifdef _OPENMP
// #pragma omp parallel for num_threads(CompilationConst::NPARALLEL)
// #endif
//     for(Idx ix=0; ix<lattice.n_sites; ix++){
//       Idx counter = 3*8*ix;
//       for(int jj=0; jj<3; jj++){
// 	Idx iy = lattice.nns[ix][jj];
// 	const MS tmp = lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj]*sigma[0] - gamma(ix, jj)) * std::exp( I* U(Link{ix,iy})) * Omega(ix, iy) / lambda_max;
// 	const MS tmp2 = lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj] + m/3.0)*sigma[0] / lambda_max;

// 	// res[NS*ix] += -tmp(0,0)*v[NS*iy] - tmp(0,1)*v[NS*iy+1];
// 	v[counter] = -tmp(0,0); counter++;
// 	v[counter] = -tmp(0,1); counter++;

// 	// res[NS*ix] += tmp(0,0)*v[NS*ix] + tmp(0,1)*v[NS*ix+1];
// 	v[counter] = tmp2(0,0); counter++;
// 	v[counter] = tmp2(0,1); counter++;

// 	// res[NS*ix+1] -= -tmp(1,0)*v[NS*iy] - tmp(1,1)*v[NS*iy+1];
// 	v[counter] = tmp(1,0); counter++;
// 	v[counter] = tmp(1,1); counter++;

// 	// res[NS*ix+1] -= tmp(1,0)*v[NS*ix] + tmp(1,1)*v[NS*ix+1];
// 	v[counter] = -tmp2(1,0); counter++;
// 	v[counter] = -tmp2(1,1); counter++;
//       }
//     }
//   }


  // void d_coo_format( std::vector<COOEntry>& elem,
  //       	     const Gauge& U,
  //       	     const Link& el ) const {
  //   const Idx ix0 = el[0];
  //   const Idx iy0 = el[1];

  //   elem.clear();
  //   {
  //     // pos
  //     const Idx ix = ix0;
  //     for(int jj=0; jj<3; jj++){
  //       const Idx iy = lattice.nns[ix][jj];
  //       if(iy!=iy0) continue;
  //       // const MS tmp = lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj]*sigma[0] - gamma(ix, jj)) * I*std::exp( I* U(Link{ix,iy})) * Omega(ix, iy);
  //       const Idx il = lattice.map2il.at(Link{ix,iy});
  //       const MS tmp = 0.5/a * (link_volume[il]/ell[il]) * gamma(ix, iy) * Omega(ix, iy) - 0.5 * r * (link_volume[il]/(ell[il]*ell[il])) * I*std::exp( I* U(Link{ix,iy})) * Omega(ix, iy);

  //       // res[NS*ix] += -tmp(0,0)*v[NS*iy] - tmp(0,1)*v[NS*iy+1];
  //       elem.push_back(COOEntry(tmp(0,0),NS*ix,NS*iy));
  //       elem.push_back(COOEntry(tmp(0,1),NS*ix,NS*iy+1));

  //       // res[NS*ix+1] += -tmp(1,0)*v[NS*iy] - tmp(1,1)*v[NS*iy+1];
  //       elem.push_back(COOEntry(tmp(1,0),NS*ix+1,NS*iy));
  //       elem.push_back(COOEntry(tmp(1,1),NS*ix+1,NS*iy+1));
  //     }
  //   }

  //   {
  //     // neg
  //     const Idx iy = iy0;
  //     for(int jj=0; jj<3; jj++){
  //       const Idx ix = lattice.nns[iy0][jj];
  //       if(ix!=ix0) continue;
  //       // const MS tmp = -lattice.vol[iy]/lattice.mean_vol * (r*lattice.u[iy][jj]*sigma[0] - gamma(iy, jj)) * I*std::exp( I* U(Link{iy,ix})) * Omega(iy, ix);
  //       const Idx il = lattice.map2il.at(Link{ix,iy});
  //       // const MS tmp = -0.5/a * (link_volume[il]/ell[il]) * gamma(iy, ix) * I*std::exp( I* U(Link{iy,ix}))y* Omega(iy, ix) - 0.5 * r * (link_volume[il]/(ell[il]*ell[il])) * I*std::exp( I* U(Link{iy,ix})) * Omega(iy, ix);

  //       // res[NS*iy] += -tmp(0,0)*v[NS*ix] - tmp(0,1)*v[NS*ix+1];
  //       elem.push_back(COOEntry(tmp(0,0),NS*iy,NS*ix));
  //       elem.push_back(COOEntry(tmp(0,1),NS*iy,NS*ix+1));

  //       // res[NS*iy+1] += -tmp(1,0)*v[NS*ix] - tmp(1,1)*v[NS*ix+1];
  //       elem.push_back(COOEntry(tmp(1,0),NS*iy+1,NS*ix));
  //       elem.push_back(COOEntry(tmp(1,1),NS*iy+1,NS*ix+1));
  //     }
  //   }
  // }




  void set_kappa() {
    for(Idx il=0; il<lattice.n_links; il++) {
      kappa[il] = lattice.link_volume[il] / lattice.mean_link_volume * lattice.mean_ell / lattice.ell[il];
    }
  }

};


template <typename T> // eigen
void matmulgam5( T* res, T* v, const Idx Nx) {
  for(Idx ix=0; ix<Nx; ix++){
    res[2*ix] = v[2*ix];
    res[2*ix+1] = -v[2*ix+1];
  }
}

template <typename T> // eigen
void mult_a( T* res, const T a, const Idx N) {
  for(Idx i=0; i<N; i++) res[i] *= a;
}


template <typename T>
T matmultgam5(const T& v) {
  T res = v;
  for(Idx i=0; i<res.rows(); i++) res.row(i) *= -2*(i%2) + 1;
  return res;
}


