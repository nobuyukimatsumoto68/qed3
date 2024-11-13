#include "rng.h"
#include "s2n.h"
#include "u1_s2_dual.h"
#include "dirac_s2_dual.h"
#include "sparse.h"
#include "cg_cuda.h"
#include "metropolis.h"



int main(int argc, char* argv[]){
  using MS=Eigen::Matrix2cd;
  using VD=Eigen::Vector2d;
  using VE=Eigen::Vector3d;
  using VC=Eigen::VectorXcd;
  using Complex=std::complex<double>;

  // geometry
  const int n_refine=1;
  Lattice lattice(n_refine);
  Dirac1fonS2 D(lattice);

  const double gR = 0.4;
  const double width = 5.0 * gR / std::sqrt( lattice.n_faces );

  const bool is_compact = false;
  U1onS2 U(lattice);
  U1Wilson SW(gR, is_compact);
  Metropolis<U1Wilson, U1onS2> met(SW, width);

  for(int i=0; i<100; i++) double r = met( U );

  // ---------------------------------------

  {
    auto mat = D.matrix_form( U );
    auto gam5_D = matmultgam5( mat );
    std::cout << gam5_D.adjoint() - gam5_D << std::endl;
    std::cout << "det D = "  << mat.determinant() << std::endl;
    std::cout << "det HW = " << gam5_D.determinant() << std::endl;

    // std::cout << "det = " << mat.determinant() << std::endl;
    // const VC r = VC::Random( mat.cols() );
    // auto tmp1 = (mat.adjoint()*mat).inverse() * r;

    // // ----------------

    // const CGCUDA solver( lattice, D );

    // Complex v[solver.sparse.N];
    // for(int i=0; i<solver.sparse.N; i++) v[i] = r[i];
    // Complex res[solver.sparse.N];

    // solver( res, v, U );

    // double norm = 0.0;
    // for(int i=0; i<solver.sparse.N; i++) {
    //   // std::cout << "i = " << i << ", " << std::abs(tmp1[i] - res[i]) << std::endl;
    //   norm += std::abs(tmp1[i] - res[i]);
    // }
    // std::cout << "norm = " << norm << std::endl;
  }



  return 0; // EXIT_SUCCESS;

}
