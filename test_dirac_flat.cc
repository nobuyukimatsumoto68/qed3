#include <iostream>
#include <fstream>
#include <cstdlib>

constexpr int N = 32;

#include "dirac_flat.h"

int main(int argc, char* argv[]){

  QfeLattice lattice;

  // omp_set_num_threads(8);
  // Eigen::setNbThreads(8);

  int K1, K2, K3;
  K1=K2=K3=1;
  lattice.InitTriangle(N, K1, K2, K3);

  Dirac1fonFlat D(lattice);

  // ----------------------------------

  Eigen::MatrixXcd mat = D.matrix_form();
  Eigen::MatrixXcd sq = mat.adjoint() * mat;

  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver( mat );
  // const Eigen::MatrixXcd evec = solver.eigenvectors();
  Eigen::VectorXcd ev = solver.eigenvalues();
  // for(int i=0; i<evec.rows(); i++){
  //   const Eigen::VectorXcd check1 = sq * evec.col(i);
  //   const Eigen::VectorXcd check2 = eval[i] * evec.col(i);
  //   assert( (check1-check2).norm() < 1.0e-8 );

  //   const Eigen::VectorXcd MV = mat * evec.col(i);
  //   std::cout << ( MV.array() / evec.col(i).array() - 1.0).abs().maxCoeff() << std::endl;
  // }

  // auto ev = mat.eigenvalues();
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
