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
  QfeLatticeS2 lattice(q, n_refine);
  Dirac1fonS2 D(lattice);
  
  return 0;
}
