#include <typeinfo>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <filesystem>
#include <chrono>
#include <cstdint>
#include <complex>
#include <array>
#include <vector>
#include <map>
#include <Eigen/Dense>

using Double = double;
using Idx = std::int32_t;
using Complex = std::complex<double>;

using Link = std::array<Idx,2>; // <int,int>;
using Face = std::vector<Idx>;

using MS=Eigen::Matrix2cd;
using VD=Eigen::Vector2d;
using VE=Eigen::Vector3d;
using VC=Eigen::VectorXcd;

static constexpr int NS = 2;
static constexpr int DIM = 2;
static constexpr Complex I = Complex(0.0, 1.0);


// #define IS_DUAL
#define IS_OVERLAP
// #define IS_DAGGER
// #undef _OPENMP

// #define GAUGE_TRSF


namespace Comp{
  constexpr bool is_compact=false;

  constexpr int NPARALLEL_DUPDATE=4;
  constexpr int NPARALLEL=NPARALLEL_DUPDATE; // 12
  constexpr int NSTREAMS=NPARALLEL_DUPDATE; // 4
  constexpr int NPARALLEL_GAUGE=NPARALLEL_DUPDATE; // 12
  constexpr int NPARALLEL_SORT=NPARALLEL_DUPDATE; // 12

  constexpr int N_REFINE=8;
  constexpr int NS=2;
  constexpr int Nt=96;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE; // 30, 120, 480

  constexpr Idx Nx=NS*N_SITES; // matrix size of DW
  constexpr Idx N=Nx*Nt; // matrix size of DW

  const double TOL_INNER=1.0e-15;
  const double TOL_OUTER=1.0e-14;
}

const std::string dir = "../../dats/";

// // #define IsVerbose
// #define IsVerbose2
// // #define InfoForce
// #define InfoDelta

#include "timer.h"

#include "../../integrator/geodesic.h"

#include "s2n_simp.h"
#include "s2n_dual.h"
#include "rng.h"


#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cublas_api.h>
#include <cusolverDn.h>
using CuC = cuDoubleComplex;
#include "gpu_header.h"



#include "valence.h"
#include "gauge_ext.h"
#include "action_ext.h"

// ======================================

#include "sparse_matrix.h"
#include "dirac_base.h"
#include "dirac_simp.h"
#include "dirac_dual.h"
#include "dirac_ext.h"
// // #include "pseudofermion.h"
// #include "dirac.h"

#include "sparse_dirac.h"
#include "matpoly.h"

// #include "dirac_pf.h"
#include "overlap.h"





//------------------------------------------
// https://gist.github.com/ashwin/d88184923c7161d368a9
#include <getopt.h>

void PrintHelp()
{
  std::cout <<
    "--gsq <gsq>"
    "--Nf <Nf>"
    "--nu0 <nu0>"
    "--nu1 <nu1>"
    "--nhits <nhits>"
    "--dt <dt>"
    "--ellmax <ellmax>"
            << std::endl;
    // "--em <em>"
    // "--ab <ab>";
  exit(1);
}

void ParseArgs(int argc, char** argv,
               double& gsq,
               int& Nf,
               double& nu0,
               double& nu1,
               int& nhits,
               int& dt,
               int& ellmax
               // int& em,
               // int& ab
               ){
  const char* const short_opts = "gnpqsdl";
  const option long_opts[] = {
    {"gsq", required_argument, nullptr, 'g'},
    {"Nf", required_argument, nullptr, 'n'},
    {"nu0", required_argument, nullptr, 'p'},
    {"nu1", required_argument, nullptr, 'q'},
    {"nhits", required_argument, nullptr, 's'},
    {"dt", required_argument, nullptr, 'd'},
    {"ellmax", required_argument, nullptr, 'l'},
    //  {"em", required_argument, nullptr, 'm'},
    // {"ab", required_argument, nullptr, 'a'},
    // {"help", no_argument, nullptr, 'h'},
    {nullptr, no_argument, nullptr, 0}
  };

  int option_index;
  int opt;
  while ((opt = getopt_long(argc, argv, short_opts, long_opts, &option_index)) != -1){
    // const auto opt = getopt_long(argc, argv, short_opts, long_opts, &option_index);
    // if (-1 == opt) break;

    switch (opt) {
    case 'g':
      gsq = std::stof(optarg);
      // std::cout << "Num set to: " << num << std::endl;
      break;

    case 'n':
      Nf = std::stoi(optarg);
      // std::cout << "Beep is set to true\n";
      break;

    case 'p':
      nu0 = std::stof(optarg);
      // std::cout << "Num set to: " << num << std::endl;
      break;

    case 'q':
      nu1 = std::stof(optarg);
      // std::cout << "Beep is set to true\n";
      break;

    case 's':
      nhits = std::stoi(optarg);
      // std::cout << "Num set to: " << num << std::endl;
      break;

    case 'd':
      dt = std::stoi(optarg);
      // std::cout << "Beep is set to true\n";
      break;

    case 'l':
      ellmax = std::stoi(optarg);
      // std::cout << "Beep is set to true\n";
      break;

    // case 'm':
    //   em = std::stoi(optarg);
    //   // std::cout << "Beep is set to true\n";
    //   break;

    // case 'a':
    //   ab = std::stoi(optarg);
    //   // std::cout << "Beep is set to true\n";
    //   break;

    case 'h': // -h or --help
    case '?': // Unrecognized option
    default: PrintHelp();
      break;
    }
  }
}
//------------------------------------------









// TODO: Cusparse for SparseMatrix::act_gpu, probably defining handle in matpoly.h
// make 2 streams in V Vdag in square in Overlap
// all the operation on GPU in Overlap::operator()
// pseudofermion
// 3d
// __m256 to vectorize with AVX2


int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  double gsq = 4.0;
  // if(argc>1) gsq = atof(argv[1]);
  int Nf = 6;
  // if(argc>2) Nf = atoi(argv[2]);
  // std::cout << "# gsq = " << gsq << " Nf = " << Nf << std::endl;
  double nu0 = 1.0;
  // if(argc>3) nu0 = atof(argv[3]);
  double nu1 = 1.0;
  // if(argc>4) nu1 = atof(argv[4]);
  int nhits = 64;
  // if(argc>5) nhits = atoi(argv[5]);
  int dt = 24;
  // if(argc>6) dt = atoi(argv[6]);
  int ellmax = 3;
  // if(argc>7) ell = atoi(argv[7]);
  // int em = 0;
  // if(argc>8) em = atoi(argv[8]);
  // int ab = 3;

  ParseArgs(argc, argv,
            gsq, Nf,
            nu0, nu1,
            nhits, dt,
            ellmax);

  std::cout << "# gsq = " << gsq << " Nf = " << Nf << " nu0 = " << nu0 << " nu1 = " << nu1 << " nhits = " << nhits << " dt = " << dt << " ellmax = " << ellmax << std::endl;

  // ---------------------------------------

  constexpr Idx N = Comp::N;
  constexpr int Nt = Comp::Nt;

  using Base=S2Simp;

  std::cout << "Ylm(0, 0, 0.5, 0.5) = " << Ylm_real(0, 0, 0.5, 0.5) << std::endl;
  std::cout << "Ylm(1, 1, 0.2, 0.4) = " << Ylm_real(1, 1, 0.2, 0.4) << std::endl;
  std::cout << "Ylm(2, -1, 0.5, 1.2) = " << Ylm_real(2, -1, 0.5, 1.2) << std::endl;
  std::cout << "Ylm(3, -2, 2.4, 4.4) = " << Ylm_real(3, -2, 2.4, 4.4) << std::endl;
  std::cout << "Ylm(4, -3, 1.1, 0.4) = " << Ylm_real(4, -3, 1.1, 0.4) << std::endl;
  std::cout << "Ylm(5, -5, 0.1, 2.9) = " << Ylm_real(5, -5, 0.1, 2.9) << std::endl;
  std::cout << "Ylm(6, 5, 0.8, 0.1) = " << Ylm_real(6, 5, 0.8, 0.1) << std::endl;
  std::cout << "Ylm(7, 7, 1.5, 5.9) = " << Ylm_real(7, 7, 1.5, 5.9) << std::endl;

  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;

// #ifdef GAUGE_TRSF
//   FermionVector gauge; // (base, Nt, rng);
//   gauge.set_random_gauge(rng);
// #endif

  FermionVector src1; // (base, Nt, rng);
  FermionVector DH_src1; // (base, Nt, rng);
  FermionVector Dinv_src1; // (base, Nt, rng);
  FermionVector src2; // (base, Nt, rng);
  FermionVector DH_src2; // (base, Nt, rng);
  FermionVector Dinv_src2; // (base, Nt, rng);


  std::cout << "# calculating sink" << std::endl;



  for(int ell1=0; ell1<=ellmax; ell1++){
    for(int em1=-ell1; em1<=ell1; em1++){
      for(int ell2=0; ell2<=ellmax; ell2++){
        for(int em2=-ell2; em2<=ell2; em2++){
          src1.fill_one();
          src1.mult_Ylm_real(ell1, em1, base);

          src2.fill_one();
          src2.mult_Ylm_real_nomeasure(ell2, em2, base);

          // for(Idx s=0; s<Comp::Nt; s++) {
          Idx s=0; {
            const VC psi1 = src1.slice((s)%Comp::Nt);
            const VC psi2 = src2.slice((s)%Comp::Nt);
            const Complex tr = psi1.dot(psi2);
            std::cout << ell1 << " " << em1 << " "
                      << ell2 << " " << em2 << " "
                      << tr.real() << std::endl;
          }
        }}
    }} // for ell, em, ab

  

  return 0;

}

