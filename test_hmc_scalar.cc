#include "s2n.h"
#include "rng.h"
#include "u1_s2_dual.h"
#include "dirac_s2_dual.h"
#include "sparse.h"
// #include "cg_cuda.h"
#include "metropolis.h"
#include "hmc.h"

#include <iomanip>

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(14);
  std::clog << std::scientific << std::setprecision(14);

  using GaugeField=U1onS2;
  using GaugeForce=U1onS2;
  using GaugeAction=U1Wilson;
  using FermionAction=Dirac1fonS2;

  using MS=Eigen::Matrix2cd;
  using VD=Eigen::Vector2d;
  using VE=Eigen::Vector3d;
  using VC=Eigen::VectorXcd;
  using Complex=std::complex<double>;

  // geometry
  const int n_refine=1;

  Lattice lattice(n_refine);
  FermionAction D(lattice);
  ParallelRng rng(lattice);

  const double gR = 0.4;
  const double width = 5.0 * gR / std::sqrt( lattice.n_faces );

  const bool is_compact = false;
  GaugeField U(lattice);
  GaugeAction SW(gR, is_compact);
  Metropolis<GaugeAction, GaugeField> met(rng, SW, width);

  for(int i=0; i<100; i++) double r = met( U );

  // {
  //   double stot = 1.0;
  //   int nsteps = 10;
  //   HMC<GaugeForce,GaugeField,GaugeAction> hmc(rng, SW, stot, nsteps);

  //   GaugeForce pi( lattice );
  //   pi.gaussian( rng );

  //   GaugeForce pi1(pi);
  //   hmc.leapfrog_explicit( pi1, U );
  //   GaugeForce pir(pi1);
  //   pir *= -1.0;
  //   hmc.leapfrog_explicit( pir, U );
  //   pir += pi;
  //   for(auto elem : pir ) std::cout << elem << " ";
  // }

  for(int nsteps=10; nsteps<100; nsteps+=10){
    HMC<GaugeForce,GaugeField,GaugeAction> hmc(rng, SW, stot, nsteps);

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

  // for(int nsteps=10; nsteps<40; nsteps+=5){
  //   HMC<GaugeForce,GaugeField,GaugeAction> hmc(rng, SW, stot, nsteps);
  //   double r, dH;
  //   bool is_accept;
  //   hmc.run( U, r, dH, is_accept );
  //   std::cout << "dH = " << dH << std::endl;
  // }  

  return 0;
}
