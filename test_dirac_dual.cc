#include <iostream>
#include <fstream>
#include <cstdlib>

#include "dirac_s2_dual.h"

int main(int argc, char* argv[]){

  // geometry
  const int n_refine=3;

  Lattice lattice(n_refine);
  // SpinStructure spin(n_refine);
  Dirac1fonS2 D(lattice);

  // ----------------------------------

  // for(int ix=0; ix<lattice.n_sites; ix++){
  //   for(int jj=0; jj<lattice.sites[ix].nn; jj++){
  // int ix=0;
  // int jj=0;
  // const int iy = lattice.nns[ix][jj];
  // const int il = lattice.sites[ix].links[jj];

  // {
  //   auto mat_xy = 0.125*D.ellstar[il] * ( D.r*D.sigma[0] - D.gamma(ix, iy) ) * D.Omega(ix, iy);
  //   std::cout << "mat_xy = " << std::endl
  // 	      << mat_xy << std::endl;
  // }
  // {
  //   auto mat_yx = 0.125*D.ellstar[il] * ( D.r*D.sigma[0] - D.gamma(iy, ix) ) * D.Omega(iy, ix);
  //   std::cout << "mat_yx = " << std::endl
  // 	      << mat_yx << std::endl;
  // }


  auto mat = D.matrix_form();
  // std::cout << mat.real() << std::endl;
  auto ev = mat.eigenvalues();
  for(int i=0; i<ev.size(); i++){
    std::clog << ev[i].real() << " " << ev[i].imag() << std::endl;
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
