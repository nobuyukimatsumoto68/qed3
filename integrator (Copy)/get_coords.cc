#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <map>
#include <limits>

#include <Eigen/Dense>

// #include "geodesic.h"
// #include "integral.h"

#include "geodesic.h"
#include "integral.h"


using namespace Geodesic;

using Idx = std::int32_t;
// using Complex = std::complex<double>;
// const Complex I = Complex(0.0, 1.0);



using Link = std::array<Idx,2>; // <Idx,Idx>;
using Face = std::vector<Idx>;
// using MS=Eigen::Matrix2cd;
// using VD=Eigen::Vector2d;
// using VE=Eigen::Vector3d;
// using VC=Eigen::VectorXcd;
// using MS=Eigen::Matrix<Complex, 2, 2>;
using VD=Eigen::Matrix<Double, 2, 1>;
using VE=Eigen::Matrix<Double, 3, 1>;
// using VC=Eigen::Matrix<Complex, Eigen::Dynamic, 1>;


Double TOL=TOLLOOSE;

std::string dir = "/mnt/hdd_barracuda/qed3/dats/";


int main(int argc, char* argv[]){
  int n_refine=1;
  if(argc==2) n_refine = atoi(argv[1]);

  std::vector<VE> sites;
  {
    std::ifstream file(dir+"pts_dual_n"+std::to_string(n_refine)+"_singlepatch.dat");

    std::string str;
    while (std::getline(file, str)){
      std::istringstream iss(str);
      Double v1, v2, v3;
      iss >> v1;
      iss >> v2;
      iss >> v3;
      sites.push_back( VE(v1, v2, v3) );
    }
  }

  std::vector<Double> thetas, phis;
  for(const VE& site : sites){
    Pt x(site);
    thetas.push_back(x.xi[0]);
    phis.push_back(x.xi[1]);
  }

  {
    std::ofstream ofs(dir+"dual_coords_n"+std::to_string(n_refine)+"_singlepatch.dat");
    ofs << std::scientific << std::setprecision(15);
    for(Idx ix=0; ix<sites.size(); ix++){
      ofs << std::setw(25) << thetas[ix] << " "
          << std::setw(25) << phis[ix] << std::endl;
    }
  }


  return 0;
}

