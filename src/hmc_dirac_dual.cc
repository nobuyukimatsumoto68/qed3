#include <iostream>

// #include "statistics.h"
#include "s2n.h"

#include "rng.h"

#include "u1_s2_dual.h"
#include "dirac_s2_dual.h"

#include "metropolis.h"

int main(int argc, char* argv[]){

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


  for(int i=0; i<100; i++) {
    double r = met( U );
  }

  // {
  //   auto mat = D.matrix_form();
  //   const VC r = VC::Random( mat.cols() );
  //   std::cout << mat * r - D( r ) << std::endl;
  // }
  {
    auto mat = D.matrix_form();
    const VC r = VC::Random( mat.cols() );
    // std::cout << mat * r - D( r ) << std::endl;

    auto tmp = D(U,r);
    // std::cout << mat * r - tmp << std::endl;

    const int N = 2*lattice.n_sites;
    Complex v[N];
    for(int i=0; i<N; i++) v[i] = r[i];
    // {
    //   Complex res[N];

    //   D(res, v, U);
    //   for(int i=0; i<N; i++) std::cout << tmp[i] - res[i] << std::endl;
    // }

    // {
    //   std::vector<Complex> val;
    //   std::vector<int> cols;
    //   std::vector<int> rows;
    //   D.csr( val, cols, rows, U );
    //   // for(auto elem : rows) std::cout << elem << std::endl;

    //   Complex res[N];
    //   matmul<Complex> ( res, v, val, cols, rows );
    //   for(int i=0; i<N; i++) std::cout << tmp[i] - res[i] << std::endl;
    // }

    {
      std::vector<Complex> val;
      std::vector<int> cols;
      std::vector<int> rows;
      D.coo( val, cols, rows, U );
      // for(auto elem : rows) std::cout << elem << std::endl;

      Complex res[N];
      matmulcoo<Complex> ( res, v, val, cols, rows, N );

      double norm = 0.0;
      // for(int i=0; i<N; i++) std::cout << tmp[i] - res[i] << std::endl;
      for(int i=0; i<N; i++) norm += std::abs( tmp[i] - res[i]);
      std::cout << "norm = " << norm << std::endl;
    }
  }

  // {
  //   auto mat = D.matrix_form();
  //   const VC r = VC::Random( mat.cols() );

  //   auto tmp = mat.adjoint() * mat * r;

  //   const int N = 2*lattice.n_sites;
  //   Complex v[N];
  //   for(int i=0; i<N; i++) v[i] = r[i];

  //   {
  //     std::vector<Complex> val;
  //     std::vector<int> cols;
  //     std::vector<int> rows;
  //     D.csr( val, cols, rows, U );
  //     // for(auto elem : rows) std::cout << elem << std::endl;

  //     Complex res[N];
  //     matmul<Complex> ( res, v, val, cols, rows );

  //     memcpy( v, res, N*sizeof(Complex) );
  //     matmulgam5<Complex> ( res, v, lattice.n_sites );

  //     memcpy( v, res, N*sizeof(Complex) );
  //     matmul<Complex> ( res, v, val, cols, rows );

  //     memcpy( v, res, N*sizeof(Complex) );
  //     matmulgam5<Complex> ( res, v, lattice.n_sites );

  //     for(int i=0; i<N; i++) std::cout << tmp[i] - res[i] << std::endl;
  //   }
  // }

  // {
  //   auto mat = D.matrix_form( U );
  //   const VC r = VC::Random( mat.cols() );

  //   auto tmp1 = mat * r;
  //   auto tmp = mat.adjoint() * mat * r;

  //   const int N = 2*lattice.n_sites;
  //   Complex v[N];
  //   for(int i=0; i<N; i++) v[i] = r[i];
  //   {
  //     std::vector<Complex> val;
  //     std::vector<int> cols;
  //     std::vector<int> rows;
  //     D.coo( val, cols, rows, U );

  //     Complex res[N];
  //     matmulcoo<Complex> ( res, v, val, cols, rows, N );
  //     for(int i=0; i<N; i++) std::cout << tmp1[i] - res[i] << std::endl;

  //     memcpy( v, res, N*sizeof(Complex) );
  //     matmuladjointcoo<Complex> ( res, v, val, cols, rows, N );

  //     double norm = 0.0;
  //     for(int i=0; i<N; i++) norm += std::abs( tmp[i] - res[i]);
  //     std::cout << "norm = " << norm << std::endl;
  //     // for(int i=0; i<N; i++) std::cout << tmp[i] - res[i] << std::endl;
  //   }
  // }

  return 0;
}
