#include <iomanip>

#include <omp.h>
constexpr int nparallel=10;
// #define IsVerbose


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
  using MS=Eigen::Matrix2cd;
  using VD=Eigen::Vector2d;
  using VE=Eigen::Vector3d;
  using VC=Eigen::VectorXcd;
  using Complex=std::complex<double>;

  std::cout << std::scientific << std::setprecision(14);
  std::clog << std::scientific << std::setprecision(14);

  using Gauge=U1onS2;
  using Force=U1onS2;
  using Action=U1Wilson;
  using Fermion=Dirac1fonS2;
  using HMC=HMC<Force,Gauge,Action,Fermion>;
  using Rng=ParallelRng;

  // geometry
  const int n_refine=1;
  Lattice lattice(n_refine);
  Fermion D(lattice);
  Rng rng(lattice);

  const double gR = 0.4;

  Gauge U(lattice);
  U.gaussian( rng );

  const bool is_compact = false;
  Action SW(gR, is_compact);

  // ---------------------------------------

  double stot = 1.0;

  Force pi( lattice );
  pi.gaussian( rng );

  for(int nsteps=10; nsteps<100; nsteps+=10){
    rng.reseed( 1 );
    HMC hmc(rng, SW, D, stot, nsteps);

    Force pi1(pi);
    Gauge U1( U );
    hmc.phi.gen( U, rng );

    const double h0 = hmc.H(pi1, U1);
    hmc.leapfrog_explicit( pi1, U1 );
    const double h1 = hmc.H(pi1, U1);

    std::cout << nsteps << " " << h1-h0 << std::endl;
  }  


  // double stot = 1.0;
  // int nsteps = 10;
  // HMC hmc(rng, SW, D, stot, nsteps);

  // for(int n=0; n<10; n++){
  //   double r, dH;
  //   bool is_accept;
  //   hmc.run( U, r, dH, is_accept );
  //   std::cout << n << "; " << is_accept << "; dH = " << dH << std::endl;
  // }  



  return 0; // EXIT_SUCCESS;
}
