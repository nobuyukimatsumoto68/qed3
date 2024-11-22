#include <iostream>
#include <cstdlib>
#include <ctime>
#include <boost/math/special_functions/bessel.hpp>

#include "statistics.h"
#include "s2util.h"
#include "u1_s2.h"
#include "metropolis.h"

int main(int argc, char* argv[]){

  // geometry
  const int q=5; // icosahedron

  // wilson gauge action
  const double gR = 0.4;
  int n_refine=1; // no refinement
  if(argc==2) n_refine = std::atoi(argv[1]);
  QfeLatticeS2 lattice(q, n_refine);
  std::cout << "n_refine = " << n_refine << std::endl;

  // metropolis
  const double width = 5.0 * gR / std::sqrt( lattice.n_faces );
  const int met_iter = 2e4;
  const int n_therm = 2e3;
  const int interval = 20;

  // random
  srand(time(NULL));
  const unsigned int seed = rand();

  // ----------------------------------

  std::cout << "gR = " << gR << std::endl;
  lattice.SeedRng(seed);

  const bool is_compact = false;
  U1onS2 U(lattice);
  U1Wilson SW(gR, is_compact);
  Metropolis<U1Wilson, U1onS2> met(SW, width);

  QfeMeasReal plaq;
  QfeMeasReal acceptance;

  for(int i=0; i<met_iter; i++) {
    double r = met( U );

    if(i>=n_therm && i%interval==0){
      double mean, square;
      U.dist_plaqsq( mean, square );
      plaq.Measure( mean );
      acceptance.Measure( r );
    }
  }

  std::cout << "acceptance : " << std::endl
	    << acceptance.Mean() << std::endl;
  std::cout << std::endl;

  std::cout << "F2 (mean, error) : " << std::endl;
  const double exact = gR*gR;
  std::cout << 1.0/std::sqrt(lattice.n_faces) << ", ";
  std::cout << plaq.Mean()/exact << ", ";
  std::cout << plaq.Error()/exact * std::sqrt( std::max(1.0, plaq.AutocorrBack()) )
	    << std::endl;  

  return 0;
}
