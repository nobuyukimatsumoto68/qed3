#include <iostream>
#include <cstdlib>
#include <ctime>
// #include <boost/math/special_functions/bessel.hpp>

// #include "statistics.h"
#include "dirac_s2.h"





int main(int argc, char* argv[]){

  // geometry
  const int q=5; // icosahedron
  const int n_refine=1; // no refinement

  // ----------------------------------

  Eigen::IOFormat fmt(15, 0,
		      "\t",
		      "\n",
		      "",
		      "",
		      "",
		      "");

  // ----------------------------------

  QfeLatticeS2 lattice(q, n_refine);
  Dirac1fonS2 D(lattice);

  std::vector<std::vector<int> > nn_oriented;

  for(int ix=0; ix<lattice.n_sites; ix++){
    auto x = lattice.sites[ix];
    std::vector<int> x_nn_oriented;

    int iy = x.neighbors[0];
    x_nn_oriented.push_back(iy);

    for(int ell=0; ell<x.nn-1; ell++){
      bool is_break = false;

      for(int kk=0; kk<x.nn; kk++){
	const int iz = x.neighbors[kk];

	if(D.is_nn(iy,iz)){
	  if(D.sign(ix, iy, iz)==1){
	    iy = iz;
	    x_nn_oriented.push_back(iy);
	    is_break = true;
	  }
	}

	if(is_break) break;
      }
    }

    nn_oriented.push_back(x_nn_oriented);
  }

  std::cout << "nn_oriented = " << std::endl;
  for(auto elem : nn_oriented[3]){
    std::cout << elem << " ";
  }
  std::cout << std::endl;


  // std::cout << "# x y z" << std::endl;
  // std::cout << std::scientific << std::setprecision(15);

  // {
  //   auto vec = lattice.r[ix];
  //   for(auto elem : vec) {
  //     std::cout << std::setw(25) << elem << " ";
  //   }
  //   std::cout << std::endl;
  //   std::cout << projection(vec).format(fmt).transpose() << std::endl;
  // }
  // {
  //   auto vec = lattice.r[iy];
  //   for(auto elem : vec) {
  //     std::cout << std::setw(25) << elem << " ";
  //   }
  //   std::cout << std::endl;
  //   std::cout << projection(vec).format(fmt).transpose() << std::endl;
  // }


  return 0;
}
