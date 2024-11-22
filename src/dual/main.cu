#include <iomanip>

#include "s2n.h"
#include "rng.h"
#include "gauge.h"
#include "action.h"
#include "dirac.h"
#include "sparse.h"
#include "pseudofermion.h"
#include "cg_cuda.h"
// #include "metropolis.h"
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
  const double width = 5.0 * gR / std::sqrt( lattice.n_faces );

  const bool is_compact = false;
  Gauge U(lattice);
  U.gaussian( rng );

  Action SW(gR, is_compact);

  double stot = 1.0;
  int nsteps = 10;

  // ---------------------------------------

  for(int nsteps=10; nsteps<100; nsteps+=10){
    HMC hmc(rng, SW, D, stot, nsteps);
    // HMC<GaugeForce,GaugeField,GaugeAction> hmc(rng, SW, stot, nsteps);

    GaugeField U1( U );
    GaugeForce pi1(pi);

    std::cout << "pi1 = " << std::endl;
    for(auto elem : pi1 ) std::cout << elem << " ";
    std::cout << std::endl;
    std::cout << "U1 = " << std::endl;
    for(auto elem : U1 ) std::cout << elem << " ";
    std::cout << std::endl;

    const double h0 = hmc.H(pi1, U1);
    hmc.leapfrog_explicit( pi1, U1 );
    const double h1 = hmc.H(pi1, U1);

    std::clog << nsteps << " " << h1-h0 << std::endl;
  }  



  return 0; // EXIT_SUCCESS;
}
