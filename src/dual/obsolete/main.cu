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
#include "pseudofermion.h"
#include "cg_cuda.h"
#include "hmc.h"



int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(14);
  std::clog << std::scientific << std::setprecision(14);

  using Gauge=U1onS2;
  using Force=U1onS2;
  using Action=U1Wilson;
  using Fermion=Dirac1fonS2;
  using HMC=HMC<Force,Gauge,Action,Fermion>;
  using Rng=ParallelRng;

  // geometry
  const int n_refine=2;
  Lattice lattice(n_refine);
  Fermion D(lattice);
  Rng rng(lattice);
  std::cout << "# n_sites = " << lattice.n_sites << std::endl;

  const double gR = 0.4;

  Gauge U(lattice);
  U.gaussian( rng );

  const bool is_compact = false;
  Action SW(gR, is_compact);

  // ---------------------------------------

  double stot = 1.0;
  int nsteps = 40;
  HMC hmc(rng, SW, D, stot, nsteps);

  const int ntherm=100;
  const int ntot=100;
  for(int n=0; n<ntherm; n++){
    double r, dH;
    bool is_accept;
    hmc.run( U, r, dH, is_accept );
    std::cout << n << "; " << is_accept << "; dH = " << dH << std::endl;
  }
  double acceptance=0.0;
  for(int n=0; n<ntot; n++){
    double r, dH;
    bool is_accept;
    hmc.run( U, r, dH, is_accept );
    acceptance+=r;
    std::cout << n << "; " << is_accept << "; dH = " << dH << std::endl;
  }
  acceptance /= ntot;
  std::cout << "acceptance = " << acceptance << std::endl;


  return 0; // EXIT_SUCCESS;
}
