#include <iostream>
#include <fstream>
#include <cstdlib>

#include "s2n.h"
#include "rng.h"

#include "u1_s2_dual.h"
#include "dirac_s2_dual.h"

#include "sparse.h"

#include "metropolis.h"


int main(int argc, char* argv[]){

  // geometry
  const int n_refine=1;

  Lattice lattice(n_refine);
  // SpinStructure spin(n_refine);
  Dirac1fonS2 D(lattice);

  const double gR = 0.4;
  const double width = 5.0 * gR / std::sqrt( lattice.n_faces );

  const bool is_compact = false;
  U1onS2 U(lattice);
  U1Wilson SW(gR, is_compact);
  Metropolis<U1Wilson, U1onS2> met(SW, width);

  for(int i=0; i<100; i++) double r = met( U );

  {
    auto mat = D.matrix_form( U );
    const VC r = VC::Random( mat.cols() );

    // ----------------

    auto tmp1 = mat * r;

    // ----------------

    Sparse sparse(lattice);
    Complex v_coo[ sparse.len ], v_csr[ sparse.len ];
    D.coo_format( v_coo, sparse.N, U );
    sparse.coo2csr( v_csr, v_coo );

    Complex res[sparse.N];
    Complex v[sparse.N];
    for(int i=0; i<sparse.N; i++) v[i] = r[i];
    sparse.mult( res, v, v_csr );

    double norm = 0.0;
    for(int i=0; i<sparse.N; i++)  norm += std::abs(tmp1[i] - res[i]);
    std::cout << "norm = " << norm << std::endl;
  }

  {
    auto mat = D.matrix_form( U );
    const VC r = VC::Random( mat.cols() );

    // ----------------

    auto tmp1 = mat.transpose() * r;

    // ----------------

    Sparse sparse(lattice);
    Complex v_coo[ sparse.len ], v_csrT[ sparse.len ];
    D.coo_format( v_coo, sparse.N, U );
    sparse.coo2csrT( v_csrT, v_coo );

    Complex res[sparse.N];
    Complex v[sparse.N];
    for(int i=0; i<sparse.N; i++) v[i] = r[i];
    sparse.multT( res, v, v_csrT );

    double norm = 0.0;
    for(int i=0; i<sparse.N; i++)  norm += std::abs(tmp1[i] - res[i]);
    std::cout << "norm = " << norm << std::endl;
  }

  {
    auto mat = D.matrix_form( U );
    const VC r = VC::Random( mat.cols() );

    // ----------------

    auto tmp1 = mat.adjoint() * r;

    // ----------------

    Sparse sparse(lattice);
    Complex v_coo[ sparse.len ], v_csrH[ sparse.len ];
    D.coo_format( v_coo, sparse.N, U );
    sparse.coo2csrH( v_csrH, v_coo );

    Complex res[sparse.N];
    Complex v[sparse.N];
    for(int i=0; i<sparse.N; i++) v[i] = r[i];
    sparse.multT( res, v, v_csrH );

    double norm = 0.0;
    for(int i=0; i<sparse.N; i++)  norm += std::abs(tmp1[i] - res[i]);
    std::cout << "norm = " << norm << std::endl;
  }

  {
    auto mat = D.matrix_form( U );
    const VC r = VC::Random( mat.cols() );

    // ----------------

    auto tmp1 = mat.adjoint() * mat * r;

    // ----------------

    Sparse sparse(lattice);
    Complex v_coo[ sparse.len ], v_csr[ sparse.len ], v_csrH[ sparse.len ];
    D.coo_format( v_coo, sparse.N, U );
    sparse.coo2csr_csrH( v_csr, v_csrH, v_coo );

    Complex res[sparse.N];
    Complex v[sparse.N];
    for(int i=0; i<sparse.N; i++) v[i] = r[i];
    sparse.mult( res, v, v_csr );
    for(int i=0; i<sparse.N; i++) v[i] = res[i];
    sparse.multT( res, v, v_csrH );

    double norm = 0.0;
    for(int i=0; i<sparse.N; i++)  norm += std::abs(tmp1[i] - res[i]);
    std::cout << "norm = " << norm << std::endl;
  }

  return 0;
}
