#pragma once

#include <array>
#include <cmath>
#include <map>
#include <Eigen/Dense>
#include <complex>


double Mod(double a, double b=2.0*M_PI){
  int p = int(std::round(a / b));
  double r = a - p*b;
  return r;
}

struct SpinStructure{
  using Link = std::array<int,2>; // <int,int>;
  using Face = std::vector<int>;

  static constexpr int NS = 2;
  static constexpr int DIM = 2;
  static constexpr int EDIM = 3;

  std::map<const Link, const double> omega;
  std::map<const Link, const double> alpha;

  SpinStructure(const int n_refine)
  {
    {
      std::cout << "# reading omega" << std::endl;
      std::ifstream file("../../dats/omega_dual_n"+std::to_string(n_refine)+"_singlepatch.dat");
      assert(file.is_open());

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
      std::cout << "# reading alpha" << std::endl;
      std::ifstream file("../../dats/alpha_dual_n"+std::to_string(n_refine)+"_singlepatch.dat");
      assert(file.is_open());

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
  using Complex = std::complex<double>;
  using Gauge=U1onS2;

  using MS=Eigen::Matrix2cd;
  using VD=Eigen::Vector2d;
  using VE=Eigen::Vector3d;
  using VC=Eigen::VectorXcd;

  static constexpr Complex I = Complex(0.0, 1.0);

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
    , ell(lattice.n_links)
    , link_volume(lattice.n_links)
  {
    set_sigma();
    // set_ell_and_link_volumes();

    // check
    double TOL=1.0e-6;
    {
      std::cout << "# checking spin structure" << std::endl;
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
  

  MS gamma(const int ix, const int jj) const { // located at x
    return lattice.v[ix][jj](0) * sigma[1] + lattice.v[ix][jj](1) * sigma[2];
  }

  MS Omega(const int ix, const int iy) const {
    const double om = omega.at(Link{ix,iy});
    return std::cos(0.5*om)*sigma[0] - I*std::sin(0.5*om)*sigma[3];
  }

  Eigen::MatrixXcd matrix_form() const {
    Eigen::MatrixXcd res = Eigen::MatrixXcd::Zero(NS*lattice.n_sites, NS*lattice.n_sites);

    for(int ix=0; ix<lattice.n_sites; ix++){
      for(int jj=0; jj<3; jj++){
	int iy = lattice.nns[ix][jj];
	res.block<NS,NS>(NS*ix,NS*iy) -= lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj]*sigma[0] - gamma(ix, jj)) * Omega(ix, iy);
	res.block<NS,NS>(NS*ix,NS*ix) += lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj] + m/3.0)*sigma[0];
      }
    }

    return res;
  } // end matrix_form


  Eigen::MatrixXcd matrix_form( const Gauge& U ) const {
    Eigen::MatrixXcd res = Eigen::MatrixXcd::Zero(NS*lattice.n_sites, NS*lattice.n_sites);

    for(int ix=0; ix<lattice.n_sites; ix++){
      for(int jj=0; jj<3; jj++){
	int iy = lattice.nns[ix][jj];
	res.block<NS,NS>(NS*ix,NS*iy) -= lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj]*sigma[0] - gamma(ix, jj)) * std::exp( I* U(Link{ix,iy})) * Omega(ix, iy);
	res.block<NS,NS>(NS*ix,NS*ix) += lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj] + m/3.0)*sigma[0];
      }
    }

    return res;
  } // end matrix_form



  void coo_format( Complex* v,
		   const Gauge& U ) const {
    const int N = lattice.n_sites * NS;
    // for(int i=0; i<N; i++) v[i] = 0.0;

    // int counter=0;
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(nparallel)
#endif
    for(int ix=0; ix<lattice.n_sites; ix++){
      int counter = 3*8*ix;
      for(int jj=0; jj<3; jj++){
	int iy = lattice.nns[ix][jj];
	const MS tmp = lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj]*sigma[0] - gamma(ix, jj)) * std::exp( I* U(Link{ix,iy})) * Omega(ix, iy);
	const MS tmp2 = lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj] + m/3.0)*sigma[0];

	// res[NS*ix] += -tmp(0,0)*v[NS*iy] - tmp(0,1)*v[NS*iy+1];
	v[counter] = -tmp(0,0); counter++;
	v[counter] = -tmp(0,1); counter++;

	// res[NS*ix] += tmp(0,0)*v[NS*ix] + tmp(0,1)*v[NS*ix+1];
	v[counter] = tmp2(0,0); counter++;
	v[counter] = tmp2(0,1); counter++;

	// res[NS*ix+1] += -tmp(1,0)*v[NS*iy] - tmp(1,1)*v[NS*iy+1];
	v[counter] = -tmp(1,0); counter++;
	v[counter] = -tmp(1,1); counter++;

	// res[NS*ix+1] += tmp(1,0)*v[NS*ix] + tmp(1,1)*v[NS*ix+1];
	v[counter] = tmp2(1,0); counter++;
	v[counter] = tmp2(1,1); counter++;
      }
    }
  }

  void H_coo_format( Complex* v,
		     const Gauge& U,
		     const double lambda_max=1.0 ) const {
    const int N = lattice.n_sites * NS;
    // for(int i=0; i<N; i++) v[i] = 0.0;

    // int counter=0;
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(nparallel)
#endif
    for(int ix=0; ix<lattice.n_sites; ix++){
      int counter = 3*8*ix;
      for(int jj=0; jj<3; jj++){
	int iy = lattice.nns[ix][jj];
	const MS tmp = lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj]*sigma[0] - gamma(ix, jj)) * std::exp( I* U(Link{ix,iy})) * Omega(ix, iy) / lambda_max;
	const MS tmp2 = lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj] + m/3.0)*sigma[0] / lambda_max;

	// res[NS*ix] += -tmp(0,0)*v[NS*iy] - tmp(0,1)*v[NS*iy+1];
	v[counter] = -tmp(0,0); counter++;
	v[counter] = -tmp(0,1); counter++;

	// res[NS*ix] += tmp(0,0)*v[NS*ix] + tmp(0,1)*v[NS*ix+1];
	v[counter] = tmp2(0,0); counter++;
	v[counter] = tmp2(0,1); counter++;

	// res[NS*ix+1] -= -tmp(1,0)*v[NS*iy] - tmp(1,1)*v[NS*iy+1];
	v[counter] = tmp(1,0); counter++;
	v[counter] = tmp(1,1); counter++;

	// res[NS*ix+1] -= tmp(1,0)*v[NS*ix] + tmp(1,1)*v[NS*ix+1];
	v[counter] = -tmp2(1,0); counter++;
	v[counter] = -tmp2(1,1); counter++;
      }
    }
  }

  
  void d_coo_format( std::vector<Complex>& v,
		     std::vector<int>& is,
		     std::vector<int>& js,
		     const Gauge& U,
		     const Link& ell ) const {
    const int ix0 = ell[0];
    const int iy0 = ell[1];

    {
      // pos
      const int ix = ix0;
      for(int jj=0; jj<3; jj++){
	const int iy = lattice.nns[ix][jj];
	if(iy!=iy0) continue;
	const MS tmp = lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj]*sigma[0] - gamma(ix, jj)) * I*std::exp( I* U(Link{ix,iy})) * Omega(ix, iy);

	// res[NS*ix] += -tmp(0,0)*v[NS*iy] - tmp(0,1)*v[NS*iy+1];
	v.push_back(-tmp(0,0)); is.push_back(NS*ix); js.push_back(NS*iy);
	v.push_back(-tmp(0,1)); is.push_back(NS*ix); js.push_back(NS*iy+1);

	// res[NS*ix+1] += -tmp(1,0)*v[NS*iy] - tmp(1,1)*v[NS*iy+1];
	v.push_back(-tmp(1,0)); is.push_back(NS*ix+1); js.push_back(NS*iy);
	v.push_back(-tmp(1,1)); is.push_back(NS*ix+1); js.push_back(NS*iy+1);
      }
    }

    {
      // neg
      const int iy = iy0;
      for(int jj=0; jj<3; jj++){
	const int ix = lattice.nns[iy0][jj];
	if(ix!=ix0) continue;
	const MS tmp = -lattice.vol[iy]/lattice.mean_vol * (r*lattice.u[iy][jj]*sigma[0] - gamma(iy, jj)) * I*std::exp( I* U(Link{iy,ix})) * Omega(iy, ix);

	// res[NS*iy] += -tmp(0,0)*v[NS*ix] - tmp(0,1)*v[NS*ix+1];
	v.push_back(-tmp(0,0)); is.push_back(NS*iy); js.push_back(NS*ix);
	v.push_back(-tmp(0,1)); is.push_back(NS*iy); js.push_back(NS*ix+1);

	// res[NS*iy+1] += -tmp(1,0)*v[NS*ix] - tmp(1,1)*v[NS*ix+1];
	v.push_back(-tmp(1,0)); is.push_back(NS*iy+1); js.push_back(NS*ix);
	v.push_back(-tmp(1,1)); is.push_back(NS*iy+1); js.push_back(NS*ix+1);
      }
    }
  }


};


template <typename T> // eigen
void matmulgam5( T* res, T* v, const int Nx) {
  for(int ix=0; ix<Nx; ix++){
    res[2*ix] = v[2*ix];
    res[2*ix+1] = -v[2*ix+1];
  }
}

template <typename T> // eigen
void mult_a( T* res, const T a, const int N) {
  for(int i=0; i<N; i++) res[i] *= a;
}


template <typename T>
T matmultgam5(const T& v) {
  T res = v;
  for(int i=0; i<res.rows(); i++) res.row(i) *= -2*(i%2) + 1;
  return res;
}


