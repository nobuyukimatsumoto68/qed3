#include <iostream>
#include <iomanip>
#include <fstream>

#include "s2.h"

int main(){

  const int q=5; // icosahedron
  const int n_refine=5;

  QfeLatticeS2 lattice(q, n_refine);

  {
    std::ofstream ofs("vertex_coordinates_n"+std::to_string(n_refine)+".dat");
    Eigen::IOFormat fmt(15, 0,
			"\t",
			"\n",
			"",
			"",
			"",
			"");
    ofs << "# x y z" << std::endl;
    ofs << std::scientific << std::setprecision(15);
    for(auto vec : lattice.r) {
      for(auto elem : vec) {
	ofs << std::setw(25) << elem << " ";
      }
      ofs << std::endl;
    }
  }

  // {
  //   int counter=0;
  //   std::ofstream ofs("nearest_neighbor_table.dat");
  //   ofs << "# ix iy il" << std::endl;
  //   ofs << std::scientific << std::setprecision(15);
  //   for(int ix=0; ix<lattice.n_sites; ix++) {
  //     const auto x = lattice.sites[ix];
  //     for(int iw=0; iw<x.nn; iw++){
  // 	ofs << std::setw(10) << ix << " ";
  // 	ofs << std::setw(10) << x.neighbors[iw] << " ";
  // 	ofs << std::setw(10) << x.links[iw] << " ";
  // 	ofs << std::endl;
  //     }
  //   }
  // }

  return 0;
}
