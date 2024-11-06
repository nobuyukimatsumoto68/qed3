#include <iostream>
#include <cstdlib>
#include <ctime>
#include <boost/math/special_functions/bessel.hpp>

#include "statistics.h"
#include "s2util.h"
#include "u1_s2.h"

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
  const double width = 2.0 * gR / std::sqrt( lattice.n_faces );
  const int met_iter = 2e3;
  const int n_therm = 5e2;
  const int interval = 5;

  // random
  srand(time(NULL));
  const unsigned int seed = rand();

  // ----------------------------------

  std::cout << "gR = " << gR << std::endl;
  lattice.SeedRng(seed);

  const bool is_compact = true;
  U1onS2 U(lattice);
  U1Wilson SW(gR, is_compact);
  Metropolis<U1Wilson, U1onS2> met(SW, width);

  // std::cout << SW.gR << std::endl;
  // for(auto elem : U.info_p.vps) std::cout << elem << std::endl;
  // // FaceInfo info_p( lattice );
  // // for(auto elem : info_p.vps) std::cout << elem << std::endl;
  // std::cout << U.plaquette_angle( 0 ) << std::endl;
  // std::cout << SW( U ) << std::endl;

  QfeMeasReal plaq;
  QfeMeasReal var;
  QfeMeasReal acceptance;

  // const int iface = 1;
  // const double factor = 0.5 / U.info_p.vp(iface);
  // std::cout << "vp = " << U.info_p.vp(iface) << std::endl;
  // const double factor = 0.5;

  for(int i=0; i<met_iter; i++) {
    double r = met( U );

    if(i>=n_therm && i%interval==0){
      // if(is_compact) plaq.Measure( U.average_plaquette() );
      // else
      // double val = U.plaquette_angle(iface);
      // while(val>M_PI) val -= 2.0*M_PI;
      // while(val<-M_PI) val += 2.0*M_PI;
      // std::cout << val << std::endl;

      double mean, square;
      U.dist_plaqsq( mean, square );
      plaq.Measure( mean );
      // var.Measure( square - mean*mean );
      acceptance.Measure( r );
    }
  }

  std::cout << "acceptance : " << std::endl
	    << acceptance.Mean() << std::endl;
  std::cout << std::endl;

  std::cout << "F2 (mean, error) : " << std::endl;
  // if(is_compact) {
  //   std::cout << lattice.n_faces * ( 1.0-plaq.Mean() ) << " ";
  //   std::cout << lattice.n_faces * plaq.Error() * std::sqrt( plaq.AutocorrBack() )
  // 	      << std::endl;
  // }
  // else{
  const double exact = gR*gR;
  std::cout << plaq.Mean()/exact << " ";
  std::cout << plaq.Error()/exact * std::sqrt( plaq.AutocorrBack() )
	    << std::endl;
  // std::cout << var.Mean() << " ";
  // std::cout << var.Error() * std::sqrt( var.AutocorrBack() )
  // 	    << std::endl;
  // }
  // << " "
  // << plaq.AutocorrFront() << " "
  // << plaq.AutocorrBack()

  
  // double denom = 0.0;
  // double numer = 0.0;

  // double tmp = 0.0;
  // // const int n=2;
  // // for(int n=0; n<lattice.n_faces; n++){
  // //   std::cout << std::exp(-beta) * boost::math::cyl_bessel_i(n, beta) << std::endl;
  // //   const double bessel_n = ine(n, beta);
  // //   const double bessel_np1 = ine(n+1, beta);
  // //   const double bessel_nm1 = ine(std::abs(n-1), beta);
  // //   const double d_bessel_n = 0.5*(bessel_np1 + bessel_nm1);
  // //   // std::cout << "bessel_n = " << bessel_n << std::endl;
  // //   // std::cout << "bessel_np1 (boost)" << std::exp(-beta) * boost::math::cyl_bessel_i(n+1, beta) << std::endl;
  // //   // std::cout << "bessel_np1 = " << bessel_np1 << std::endl;
  // //   // std::cout << "bessel_nm1 = " << bessel_nm1 << std::endl;
  // //   // std::cout << "d_bessel_n = " << d_bessel_n << std::endl;

  // //   denom += std::pow(bessel_n / bessel_0, lattice.n_faces);
  // //   numer += d_bessel_n / bessel_0 * std::pow(bessel_n / bessel_0, lattice.n_faces-1);

  // //   std::cout << "denom = " << denom << std::endl;
  // //   std::cout << "numer = " << numer << std::endl;
  // // }
  // double beta = 4.0;
  // const double bessel_0 = std::exp(-beta) * boost::math::cyl_bessel_i(0, beta);
  // for(int n=0; n<lattice.n_faces; n++){
  //   const double bessel_n = std::exp(-beta) * boost::math::cyl_bessel_i(n, beta);
  //   const double bessel_np1 = std::exp(-beta) * boost::math::cyl_bessel_i(n+1, beta);
  //   const double bessel_nm1 = std::exp(-beta) * boost::math::cyl_bessel_i(std::abs(n-1), beta);
  //   const double d_bessel_n = 0.5*(bessel_np1 + bessel_nm1);
  //   // const double bessel_n = ine(n, beta);
  //   // std::cout << "bessel = " << bessel_n << std::endl;
  //   // const double bessel_np1 = ine(n+1, beta);
  //   // std::cout << "bessel = " << bessel_np1 << std::endl;
  //   // const double bessel_nm1 = ine(std::abs(n-1), beta);
  //   // std::cout << "bessel = " << bessel_nm1 << std::endl;
  //   // const double d_bessel_n = 0.5*(bessel_np1 + bessel_nm1);

  //   denom += std::pow(bessel_n / bessel_0, lattice.n_faces);
  //   numer += d_bessel_n / bessel_0 * std::pow(bessel_n / bessel_0, lattice.n_faces-1);
  //   // std::cout << "denom = " << denom << std::endl;
  //   // std::cout << "numer = " << numer << std::endl;
  //   // if( std::abs(numer / denom - tmp ) < 1.0e-10  ) break;
  //   // tmp = numer / denom;
  // }
  // const double exact = numer / denom;
  // std::cout << "exact = " << lattice.n_faces * ( 1.0 - exact ) << std::endl;
  
  return 0;
}
