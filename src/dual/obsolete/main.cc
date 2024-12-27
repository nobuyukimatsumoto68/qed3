#include <iomanip>

#include "s2n.h"
#include "rng.h"
#include "gauge.h"
#include "action.h"
#include "dirac.h"
#include "sparse.h"
#include "pseudofermion.h"
#include "cg_cuda.h"
#include "metropolis.h"



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
  using GaugeAction=U1Wilson;
  using Fermion=Dirac1fonS2;
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
  GaugeAction SW(gR, is_compact);

  Metropolis<GaugeAction, Gauge> met(rng, SW, width);
  for(int i=0; i<100; i++) double r = met( U );

  // ---------------------------------------

  {
    using Link = std::array<int,2>; // <int,int>;
    PseudoFermion phi( D, U, rng );
    // phi.gen( U, rng );

    Force exact = phi.get_force( U );

    for(int il=0; il<U.lattice.n_links; il++){
      // Link link = lattice.links[il];
      // const int ix=link[0], iy=link[1];

      const double eps = 1.0e-5;
      // const int ell=lattice.map2il.at( link );
      Gauge UP(U);
      Gauge UM(U);

      UP[il] += eps;
      UM[il] -= eps;

      std::vector<Complex> etaP = phi.get_eta( UP );
      std::vector<Complex> etaM = phi.get_eta( UM );

      Complex SP = phi.dot( phi.phi, etaP );
      Complex SM = phi.dot( phi.phi, etaM );
      Complex numeric = ( SP - SM ) / (2.0*eps);
      std::cout << exact[il] << " " << numeric << std::endl;
    }
  }



  return 0; // EXIT_SUCCESS;
}
