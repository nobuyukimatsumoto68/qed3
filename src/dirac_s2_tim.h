#pragma once

#include <array>
#include <cmath>
#include <iomanip>
// #include <map>
#include <Eigen/Dense>


// using Complex = std::complex<double>;
// const Complex I = Complex(0.0, 1.0);

// constexpr int NS = 2;
// using MS=Eigen::Matrix2cd;

// constexpr int DIM = 2;
// using VD=Eigen::Vector2d;


// double Mod(double a, double b=2.0*M_PI){
//   int p = int(std::round(a / b));
//   double r = a - p*b;
//   return r;
// }

template<typename Derived>
std::ostream& operator<<(std::ostream& os, const Eigen::MatrixBase<Derived>& W) {
  os << std::scientific << std::setprecision(15);
  for(int i=0; i<W.rows(); i++){
    for(int j=0; j<W.rows(); j++){
      os << std::setw(22) << W(i,j) << " ";
    }
    os << std::endl;
  }
  return os;
}

// std::ostream& operator<<(std::ostream& os, const MC& W) {
//   os << std::scientific << std::setprecision(15);
//   for(int i=0; i<W.rows(); i++){
//     for(int j=0; j<W.rows(); j++){
//       os << std::setw(22) << W(i,j).real() << " "
// 	 << std::setw(22) << W(i,j).imag() << " ";
//     }
//     os << std::endl;
//   }
//   return os;
// }

// std::ostream& operator<<(std::ostream& os, const VR& W) {
//   os << std::scientific << std::setprecision(15);
//   for(int i=0; i<W.size(); i++){
//     os << std::setw(22) << W(i) << " ";
//   }
//   os << std::endl;
//   return os;
// }

// std::ostream& operator<<(std::ostream& os, const VC& W) {
//   os << std::scientific << std::setprecision(15);
//   for(int i=0; i<W.size(); i++){
//     os << std::setw(22) << W(i).real() << " "
//        << std::setw(22) << W(i).imag() << " ";
//   }
//   os << std::endl;
//   return os;
// }


struct Dirac1fonS2Tim {
  const std::array<MS, 4> sigma;

  const std::vector<std::vector<double>> alpha;
  const std::vector<std::vector<double>> omega;
  const std::vector<std::vector<int>> nn;
  const std::vector<std::vector<double>> As;
  const double r;
  const std::vector<std::vector<double>> kappa;
  const std::vector<double> Bs;
  const std::vector<int> Tim2Evan;

  // Dirac1fonS2Tim()=delete;

  int find_jj( const int ix, const int iy ) const {
    for(int jj=0; jj<nn[ix].size(); jj++){
      if(nn[ix][jj]==iy) return jj;
    }
    assert(false);
  }

  Dirac1fonS2Tim()
    : sigma( get_sigma() )
    , alpha( get_alpha() )
    , omega( get_omega() )
    , nn( get_nn() )
    , As( get_A() )
    , r( get_r() )
    , kappa( get_kappa() )
    , Bs( get_B() )
    , Tim2Evan( get_Tim2Evan() )
  {
    // check
    double TOL=1.0e-6;
    {
      const int n_sites = nn.size();
      for(int ix=0; ix<n_sites; ix++){
	for(int jj=0; jj<nn[ix].size(); jj++){
	  const int iy = nn[ix][jj];

	  const double alpha1 = alpha[ix][jj];

	  // int kk;
	  // for(kk=0; kk<nn[iy].size(); kk++){
	  //   if(nn[iy][kk]==ix) break;
	  // }
	  // if(kk==nn[iy].size()) assert(false);
	  int kk = find_jj( iy, ix );

	  const double alpha2 = alpha[iy][kk];
	  double omega12 = omega[ix][jj];
	  double omega21 = omega[iy][kk];

	  double diff1 = (alpha2 + M_PI + omega12) - alpha1;
	  double diff2 = omega12 + omega21;
	  // std::cout << diff1 << " " << diff2 << std::endl;
	  assert( std::abs(Mod(diff1))<TOL );
	  assert( std::abs(Mod(diff2))<TOL );
	}}
    }

    // {
    //   for(int ia=0; ia<lattice.n_faces; ia++){
    // 	double sum = 0.0;

    // 	for(int i=0; i<3; i++){
    // 	  int ix = lattice.faces[ia].sites[i];
    // 	  int iy = lattice.faces[ia].sites[(i+1)%3];
    // 	  sum -= omega.at(Link{ix,iy});
    // 	  sum += alpha.at(Link{ix,iy});
    // 	  sum -= alpha.at(Link{iy,ix}) + M_PI;
    // 	}

    // 	assert( std::abs(Mod(-std::abs(Mod(sum)))) < TOL );
    //   }
    // }
  }

  // @@ alpha[ix][jj]
  MS gamma(const int ix, const int jj) const { // located at x
    const double al = alpha[ix][jj];
    return std::cos(al)*sigma[1] + std::sin(al)*sigma[2];
  }

  // @@ omega[ix][jj]
  MS Omega(const int ix, const int jj) const {
    const double om = omega[ix][jj];
    return std::cos(0.5*om)*sigma[0] - I*std::sin(0.5*om)*sigma[3];
  }

  // @@ nn[ix][jj] // vector<vector>
  // @@ A[ix][jj] //vector<vector>
  // @@ r
  // @@ kappa[ix][jj]
  // @@ B[ix]
  Eigen::MatrixXcd matrix_form() const {
    const int n_sites = nn.size();
    Eigen::MatrixXcd res = Eigen::MatrixXcd::Zero(NS*n_sites, NS*n_sites);

    for(int ix=0; ix<n_sites; ix++){
      for(int jj=0; jj<nn[ix].size(); jj++){
	const int iy = nn[ix][jj];
	res.block<NS,NS>(NS*ix, NS*iy) += 0.5*As[ix][jj]*( -r*sigma[0] + kappa[ix][jj]*gamma(ix, jj) ) * Omega(ix, jj);
      }
      res.block<NS,NS>(NS*ix,NS*ix) += 0.5 * Bs[ix] * sigma[0];
    }

    return res;
  } // end matrix_form


  Eigen::MatrixXcd matrix_form_alternative( const Dirac1fonS2& D0 ) const {
    const int n_sites = nn.size();
    Eigen::MatrixXcd res = Eigen::MatrixXcd::Zero(NS*n_sites, NS*n_sites);

    for(int ix=0; ix<n_sites; ix++){
      const int ixEvan = Tim2Evan[ix];

      for(int jj=0; jj<nn[ix].size(); jj++){
	const int iy = nn[ix][jj];
	const int iyEvan = Tim2Evan[iy];

	res.block<NS,NS>(NS*ix, NS*iy) += 0.5*As[ix][jj]*( -r*sigma[0] + kappa[ix][jj]*D0.gamma(ixEvan, iyEvan) ) * D0.Omega(ixEvan, iyEvan);
      }
      res.block<NS,NS>(NS*ix,NS*ix) += 0.5 * Bs[ix] * sigma[0];
    }

    return res;
  } // end matrix_form



  std::array<MS, 4> get_sigma(){
    std::array<MS, 4> res;
    assert(NS==2);
    res[0] << 1,0,0,1;
    res[1] << 0,1,1,0;
    res[2] << 0,-I,I,0;
    res[3] << 1,0,0,-1;
    return res;
  }

  std::vector<std::vector<double>> get_alpha(){
    std::ifstream file("./prev/alphaijTim.dat");
    std::vector<std::vector<double>> res;

    std::string str;
    while (std::getline(file, str)){
      std::istringstream iss(str);
      std::vector<double> row;
      double v;
      while( iss >> v ) row.push_back(v);
      res.push_back( row );
    }
    return res;
  }

  std::vector<std::vector<double>> get_omega(){
    std::ifstream file("./prev/phiijTim.dat");
    std::vector<std::vector<double>> res;

    std::string str;
    while (std::getline(file, str)){
      std::istringstream iss(str);
      std::vector<double> row;
      double v;
      while( iss >> v ) row.push_back(v);
      res.push_back( row );
    }
    return res;
  }

  std::vector<std::vector<int>> get_nn(){
    std::ifstream file("./prev/nntableTim.dat");
    std::vector<std::vector<int>> res;

    std::string str;
    while (std::getline(file, str)){
      std::istringstream iss(str);
      std::vector<int> row;
      int v;
      while( iss >> v ) row.push_back(v-1);
      res.push_back( row );
    }
    return res;
  }

  std::vector<std::vector<double>> get_A(){
    std::ifstream file("./prev/AijTim.dat");
    std::vector<std::vector<double>> res;

    std::string str;
    while (std::getline(file, str)){
      std::istringstream iss(str);
      std::vector<double> row;
      double v;
      while( iss >> v ) row.push_back(v);
      res.push_back( row );
    }
    return res;
  }

  double get_r(){
    std::ifstream file("./prev/rTim.dat");
    double res;

    std::string str;
    while (std::getline(file, str)){
      std::istringstream iss(str);
      double v;
      while( iss >> v ) res = v;
    }
    return res;
  }

  std::vector<std::vector<double>> get_kappa(){
    std::ifstream file("./prev/kappaijTim.dat");
    std::vector<std::vector<double>> res;

    std::string str;
    while (std::getline(file, str)){
      std::istringstream iss(str);
      std::vector<double> row;
      double v;
      while( iss >> v ) row.push_back(v);
      res.push_back( row );
    }
    return res;
  }

  std::vector<double> get_B(){
    std::ifstream file("./prev/BiTim.dat");
    std::vector<double> res;

    std::string str;
    while (std::getline(file, str)){
      std::istringstream iss(str);
      double v;
      while( iss >> v ) res.push_back(v);
    }
    return res;
  }

  std::vector<int> get_Tim2Evan(){
    std::ifstream file("./prev/Tim2Evan.dat");
    std::vector<int> res;

    std::string str;
    while (std::getline(file, str)){
      std::istringstream iss(str);
      int v;
      while( iss >> v ) res.push_back(v-1);
    }
    return res;
  }


};
