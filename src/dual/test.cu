#include <iostream>
#include <iomanip>

#include <omp.h>
constexpr int nparallel=10;
#define IsVerbose



#include "s2n.h"
#include "rng.h"
#include "gauge.h"
#include "action.h"
#include "dirac.h"
#include "sparse.h"
// #include "pseudofermion.h"
#include "overlap.h"
// #include "hmc.h"



int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(14);
  std::clog << std::scientific << std::setprecision(14);


  Zolotarev f;
  std::cout << "# Delta = " << f.Delta() << std::endl;

  for(double x = -1.0; x<=1.0; x+=0.001){
    std::cout << x << " " << 1.0-std::abs(f(x)) << std::endl;
  }

  

  
  return 0; // EXIT_SUCCESS;
}
