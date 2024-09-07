#include <iostream>
#include <cstdlib>
#include <ctime>

#include "statistics.h"
#include "u1_s2.h"

int main(int argc, char* argv[]){

  // geometry
  const int q=5; // icosahedron
  const int n_refine=1; // no refinement

  // wilson gauge action
  double beta = 1.0;
  if(argc==2) beta = std::atof(argv[1]);

  // metropolis
  const double width = 1.0;
  const int met_iter = 1e4;
  const int n_therm = 5e3;

  // random
  srand(time(NULL));
  const unsigned int seed = rand();

  // ----------------------------------

  QfeLatticeS2 lattice(q, n_refine);
  lattice.SeedRng(seed);

  U1Wilson SW(beta);
  CompactU1onS2 U(lattice, beta);
  Metropolis<U1Wilson, CompactU1onS2> met(SW, width);

  QfeMeasReal plaq;
  QfeMeasReal acceptance;

  for(int i=0; i<met_iter; i++) {
    double r = met( U );

    if(i>=n_therm){
      plaq.Measure( U.average_plaquette() );
      acceptance.Measure( r );
    }
  }

  std::cout << "acceptance : " << std::endl
	    << acceptance.Mean() << std::endl;
  std::cout << std::endl;

  std::cout << "plaquette (mean, error) : " << std::endl
	    << plaq.Mean() << " "
	    << plaq.Error() * std::sqrt( plaq.AutocorrBack() )
	    // << " "
	    // << plaq.AutocorrFront() << " "
	    // << plaq.AutocorrBack()
	    << std::endl;

  
  double denom = 0.0;
  double numer = 0.0;
  for(int n=0; n<lattice.n_faces; n++){
    const double bessel_n = std::cyl_bessel_i(n, beta);
    const double bessel_np1 = std::cyl_bessel_i(n+1, beta);
    const double bessel_nm1 = std::cyl_bessel_i(std::abs(n-1), beta);
    const double d_bessel_n = 0.5*(bessel_np1 + bessel_nm1);

    denom += std::pow(bessel_n, lattice.n_faces);
    numer += d_bessel_n * std::pow(bessel_n, lattice.n_faces-1);
  }
  const double exact = numer / denom;
  std::cout << "exact = " << exact << std::endl;
  
  return 0;
}
