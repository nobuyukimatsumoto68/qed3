#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cassert>
#include <cstdint>
#include <cmath>
#include <complex>
#include <array>
#include <vector>

using Double = double;
using Idx = std::int32_t;
using Complex = std::complex<double>;

static constexpr int NS = 2;
static constexpr int DIM = 2;
static constexpr Complex I = Complex(0.0, 1.0);


namespace Comp{
  constexpr int N_REFINE=16;
  constexpr int NS=2;

  constexpr int Nt=64;
  // constexpr int Nt=16;

#ifdef IS_DUAL
  constexpr Idx N_SITES=20*N_REFINE*N_REFINE;
#else
  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
#endif

  constexpr Idx Nx=NS*N_SITES; // matrix size of DW
  constexpr Idx N=Nx*Nt; // matrix size of DW
}

const int nmax=100;
const double c = 1.0;
const int nsteps = 1000;

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  for(Idx t=0; t<10*nsteps; t++){
    Double tmp = 0.0;
    // for(int n=0; n<nmax; n++) tmp += (n+1) * std::cyl_bessel_k( 0, 1.0*(n+1)*t/nsteps );
    for(int n=0; n<nmax; n++) tmp += (n+1) * std::exp( -1.0*(n+1)*t/nsteps );
    // for(int n=1; n<nmax; n++) tmp += (n+1) * std::exp( -1.0*(n+1)*t/nsteps );
    tmp /= 2.0*M_PI;
    std::cout << 1.0*t/nsteps << " " << tmp << std::endl;
  }

  return 0;
}
