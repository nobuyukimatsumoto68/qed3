#include <iostream>
// #include <fstream>
// #include <cstdlib>
#include <cmath>

#include <Eigen/Dense>

#include "geodesic.h"
#include "integral.h"


using Complex = std::complex<double>;
const Complex I = Complex(0.0, 1.0);



int main(int argc, char* argv[]){

  // monotonic
  Pt x1( 0.11594391015821674, 0.6448028608904725, 0.7555039909124803 );
  Pt x2( -0.00967100461303834, 0.5470756243846032, 0.8370273190726639 );
  // // altering
  // Pt x1( -0.042367955582481, 0.13039515939997406, 0.990556438949753 );
  // Pt x2( 0.11092074774878105, 0.08058864047764432, 0.990556438949753 );

  Sol sol = SolveGeodesics( x1, x2 );

  std::cout << "theta = " << sol.theta( 0.5*sol.sE ) << std::endl;
  std::cout << "phi = " << sol.phi( 0.5*sol.sE ) << std::endl;

  std::cout << "Dtheta = "
	    << sol.Dtheta( 0.5*sol.sE ) << std::endl;
  std::cout << "Dtheta = "
	    << (sol.theta( 0.5*sol.sE+EPSNUMDER )-sol.theta( 0.5*sol.sE-EPSNUMDER ))/(2.0*EPSNUMDER) << std::endl;
  std::cout << "Dphi = "
	    << sol.Dphi( 0.5*sol.sE ) << std::endl;
  std::cout << "Dphi = "
	    << (sol.phi( 0.5*sol.sE+EPSNUMDER )-sol.phi( 0.5*sol.sE-EPSNUMDER ))/(2.0*EPSNUMDER) << std::endl;

  std::cout << "Dphi = "
	    << sol.Dphi( 0.1 ) << std::endl;
  std::cout << "Dphi = "
	    << (sol.phi( 0.1+EPSNUMDER )-sol.phi( 0.1-EPSNUMDER ))/(2.0*EPSNUMDER) << std::endl;
  std::cout << "phi = " << sol.phi( 0.1 ) << std::endl;
  
  if(sol.is_split){
    double ss = 0.5*(sol.sE+sol.ell);
    std::cout << "theta = " << sol.theta( ss, 1) << std::endl;
    std::cout << "phi = " << sol.phi( ss, 1) << std::endl;

    std::cout << "Dtheta = "
	      << sol.Dtheta( ss, 1 ) << std::endl;
    std::cout << "Dtheta = "
	      << (sol.theta( ss+EPSNUMDER, 1)-sol.theta( ss-EPSNUMDER, 1))/(2.0*EPSNUMDER) << std::endl;
    std::cout << "Dphi = "
	      << sol.Dphi( ss ) << std::endl;
    std::cout << "Dphi = "
	      << (sol.phi( ss+EPSNUMDER, 1)-sol.phi( ss-EPSNUMDER, 1))/(2.0*EPSNUMDER) << std::endl;
  }

  std::cout << "ell = " << sol.ell << std::endl;

  {
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double result, error;
    double expected = 0.156635;

    F f1 = [&](double s){ return sol.Dphi(s) * std::cos( sol.theta(s) ); };

    std::cout << "debug. f(0.1) = " << f1(0.1) << std::endl;
    std::cout << "debug. dphi(0.1) = " << sol.Dphi(0.1) << std::endl;
    std::cout << "debug. theta(0.1) = " << sol.theta(0.1) << std::endl;

    gsl_function F;
    F.function = &unwrap;
    F.params = &f1;
     
    gsl_integration_qag(&F, 0., sol.ell, 0, 1e-7, 1000,
			3, w, &result, &error); 
     
    printf ("result          = % .18f\n", result);
    printf ("exact result    = % .18f\n", expected);
    printf ("estimated error = % .18f\n", error);
    printf ("actual error    = % .18f\n", result - expected);
    printf ("intervals =  %d\n", w->size);
  }

  
  return 0;
}

