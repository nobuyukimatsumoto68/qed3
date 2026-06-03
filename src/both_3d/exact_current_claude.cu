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


// #include <stdfloat>
// #include <boost/multiprecision/float128.hpp>

// #if __STDCPP_FLOAT64_T__ != 1
//     #error "64-bit float type required"
// #endif


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

  constexpr int NPARALLEL_DUPDATE=1;
  constexpr int NPARALLEL=NPARALLEL_DUPDATE; // 12
  constexpr int NSTREAMS=NPARALLEL_DUPDATE; // 4
  constexpr int NPARALLEL_GAUGE=NPARALLEL_DUPDATE; // 12
  constexpr int NPARALLEL_SORT=NPARALLEL_DUPDATE; // 12

  constexpr int N_REFINE=1;
  constexpr int NS=2;
  constexpr int Nt=128;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE; // 30, 120, 480

  constexpr Idx Nx=NS*N_SITES; // matrix size of DW
  constexpr Idx N=Nx*Nt; // matrix size of DW

  const double TOL_INNER=1.0e-15;
  const double TOL_OUTER=1.0e-14;
}

// const std::string dir = "../../dats/";
const std::string dir = "../../geometry/data/";

// // #define IsVerbose
// #define IsVerbose2
// // #define InfoForce
// #define InfoDelta

#include "timer.h"

// #include "../../integrator/geodesic.h"

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
// #include "matpoly.h"
#include "matpoly_claude.h"

// #include "dirac_pf.h"
#include "overlap.h"





//------------------------------------------
// https://gist.github.com/ashwin/d88184923c7161d368a9
#include <getopt.h>

void PrintHelp()
{
  printf("Usage: ./a.out [options]\n");
  printf("  --gsq <gsq>        Wilson coupling squared (default: 2.0)\n");
  printf("  --Nf <Nf>          number of fermion flavors (default: 2)\n");
  printf("  --nu0 <nu0>        valence mass parameter 0 (default: 1.0)\n");
  printf("  --nu1 <nu1>        valence mass parameter 1 (default: 1.0)\n");
  printf("  --nhits <nhits>    number of stochastic hits (default: 1)\n");
  printf("  --dt <dt>          time slice separation (default: Nt/2)\n");
  printf("  --ellmax <ellmax>  max angular momentum ell (default: 2)\n");
  printf("  -h, --help         show this help\n");
  exit(0);
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
  const char* const short_opts = "gnpqsdlh";
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
    {"help", no_argument, nullptr, 'h'},
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




// template<class Base>
// double Dot( const VC& psi1, const VC& psi2, const Base& base ){
//   // for(Idx ix=0; ix<Comp::N_SITES; ix++){
//   //   const VE r = lattice.sites[ix];
//   //   const double area = lattice.dual_areas[ix];
//   //   const double ylm = Ylm_real(ell, em, r);
//   //   for(int s=0; s<Nt; s++){
//   //     (*this)(s, ix, 0) *= area * ylm;
//   //     (*this)(s, ix, 1) *= area * ylm;
//   //   }
//   // }
//   boost::multiprecision::float128 sum = 0.0;
//   for(Idx i=0; i<psi1.size(); i++){
//     // std::cout << std::setw(15) << i << " " << std::setw(20) << (std::conj(psi1[i])*psi2[i]).real() << std::endl;
//     sum += (std::conj(psi1[i])*psi2[i]).real();
//   }
//   return (double)sum;
// }







// TODO: Cusparse for SparseMatrix::act_gpu, probably defining handle in matpoly.h
// make 2 streams in V Vdag in square in Overlap
// all the operation on GPU in Overlap::operator()
// pseudofermion
// 3d
// __m256 to vectorize with AVX2


int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  double gsq = 2.0;
  // if(argc>1) gsq = atof(argv[1]);
  int Nf = 2;
  // if(argc>2) Nf = atoi(argv[2]);
  // std::cout << "# gsq = " << gsq << " Nf = " << Nf << std::endl;
  double nu0 = 1.0;
  // if(argc>3) nu0 = atof(argv[3]);
  double nu1 = 1.0;
  // if(argc>4) nu1 = atof(argv[4]);
  int nhits = 1;
  // if(argc>5) nhits = atoi(argv[5]);
  int dt = Comp::Nt/2;
  // if(argc>6) dt = atoi(argv[6]);
  int ellmax = 2;
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

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();


  int device;
  CUDA_CHECK(cudaGetDeviceCount(&device));
  cudaDeviceProp device_prop[device];
  cudaGetDeviceProperties(&device_prop[0], 0);
  std::cout << "# dev = " << device_prop[0].name << std::endl;
  CUDA_CHECK(cudaSetDevice(0));// "TITAN V"
  std::cout << "# (GPU device is set.)" << std::endl;

  // ---------------------------------------

  constexpr Idx N = Comp::N;
  constexpr int Nt = Comp::Nt;

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;

  using Force=GaugeExt<Base,Nt,Comp::is_compact>;
  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  using Action=U1WilsonExt<Base>;

  using Rng=ParallelRngExt<Base,Nt>;
  // using FermionVector = FermionVector<Base>;


  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;

  // ----------------------
  const double at = 0.2; // T/Comp::Nt;
  if(Nt!=1) assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);

  std::string dir3, dir4;
  dir3="Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
  dir4="data_Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
  // dir4="data_at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"nu1"+std::to_string(nu1)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";

  std::filesystem::create_directory(dir4);

  Gauge U(base);
  srand( time(NULL) );
  Rng rng(base, rand());


  const double M5 = -1.0;
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu1);
  std::cout << "# DW set. " << std::endl;
  using Fermion=Overlap<WilsonDirac>;
  // Fermion D(DW, 31);
  Fermion D(DW, 21);
  std::cout << "# D set. " << std::endl;

  D.update( U );
  std::cout << "# D updated. " << std::endl;

  // auto f_DH = std::bind(&Fermion::adj_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);
  // auto f_DHD = std::bind(&Fermion::DHD_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);

  // LinOpWrapper M_DH( f_DH );
  // MatPoly op_DH; op_DH.push_back ( cplx(1.0), {&M_DH} );
  // LinOpWrapper M_DHD( f_DHD );
  // MatPoly op_DHD; op_DHD.push_back ( cplx(1.0), {&M_DHD} );

  // ---------------------

  


  std::pair<int, Idx> link(0,0);
  // const CuC* d_eta;

  COO<N> coo;
  DW.d_coo_format(coo.en, U, link);
  coo.do_it();

  FermionVector eta; // (base, Nt, rng);
  FermionVector Keta; // (base, Nt, rng);

  Complex coeff=I;
  MatPoly op_DH; op_DH.push_back ( cplx(I), {&coo} );
  op_DH.from_cpu<N>( Keta.data(), eta.data() );
  // coo.Async( d_dDeta[m], d_eta, stream[m] );

  // COO<N> coo;
  // DW.d_coo_format(coo.en, U, link);
  // coo.do_it();



















  // ------------------

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();


  return 0;

}

