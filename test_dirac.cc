#include <iostream>
#include <fstream>
#include <cstdlib>

#include "dirac_s2.h"

int main(int argc, char* argv[]){

  // geometry
  const int q=5; // icosahedron
  const int n_refine=1; // no refinement

  QfeLatticeS2 lattice(q, n_refine);
  Dirac1fonS2 D(lattice);

  // ----------------------------------

  auto mat = D.matrix_form();
  auto ev = mat.eigenvalues();
  for(int i=0; i<ev.size(); i++){
    std::cout << ev[i].real() << " " << ev[i].imag() << std::endl;
  }

  // ----------------------------------

  return 0;
}



  // for(int ix=0; ix<lattice.n_sites; ix++){
  //   for(int jj=0; jj<lattice.sites[ix].nn; jj++){
  //     const int iy = lattice.sites[ix].neighbors[jj];
  //     auto mat1 = ( D.sigma[0] - D.gamma(ix, iy) ) * D.Omega(ix, iy);
  //     auto mat2 = D.Omega(ix, iy) * ( D.sigma[0] - D.gamma(iy, ix, M_PI) );
  //     std::cout << mat1-mat2 << std::endl;
  //   }}

  // ----------------------------------
