#include <typeinfo>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <highfive/H5File.hpp>
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

using Link = std::array<Idx,2>;
using Face = std::vector<Idx>;

using MS=Eigen::Matrix2cd;
using VD=Eigen::Vector2d;
using VE=Eigen::Vector3d;
using VC=Eigen::VectorXcd;

static constexpr int NS = 2;
static constexpr int DIM = 2;
static constexpr Complex I = Complex(0.0, 1.0);


namespace Comp{
  constexpr bool is_compact=false;

  constexpr int NPARALLEL_DUPDATE=1;
  constexpr int NPARALLEL=NPARALLEL_DUPDATE;
  constexpr int NSTREAMS=NPARALLEL_DUPDATE;
  constexpr int NPARALLEL_GAUGE=NPARALLEL_DUPDATE;
  constexpr int NPARALLEL_SORT=NPARALLEL_DUPDATE;

  constexpr int N_REFINE=1;
  constexpr int NS=2;
  constexpr int Nt=128;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;

  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-15;
  const double TOL_OUTER=1.0e-14;
}

const std::string dir = "../../geometry/data/";

#include "timer.h"

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

#include "valence_claude.h"
#include "gauge_ext.h"
#include "action_ext.h"

#include "sparse_matrix.h"
#include "dirac_base.h"
#include "dirac_simp.h"
#include "dirac_dual.h"
#include "dirac_ext.h"

#include "sparse_dirac.h"
#include "matpoly.h"

#include "overlap.h"

//------------------------------------------
#include <getopt.h>

void PrintHelp(){
  printf("Usage: ./a.out [options]\n");
  printf("  --gsq <gsq>          Wilson coupling squared (default: 8.0)\n");
  printf("  --Nf <Nf>            number of fermion flavors (default: 2)\n");
  printf("  --nu0 <nu0>          sea quark asymmetry parameter (default: 1.0)\n");
  printf("  --nu1 <nu1>          valence asymmetry parameter (default: 1.0)\n");
  printf("  --nhits <nhits>      number of stochastic hits (default: 1)\n");
  printf("  --t_block <t_block>  timeslices per dilution block (default: 8)\n");
  printf("  -h, --help           show this help\n");
  exit(0);
}

void ParseArgs(int argc, char** argv,
               double& gsq,
               int&    Nf,
               double& nu0,
               double& nu1,
               int&    nhits,
               int&    t_block){
  const char* const short_opts = "gNnvHth";
  const option long_opts[] = {
    {"gsq",     required_argument, nullptr, 'g'},
    {"Nf",      required_argument, nullptr, 'N'},
    {"nu0",     required_argument, nullptr, 'n'},
    {"nu1",     required_argument, nullptr, 'v'},
    {"nhits",   required_argument, nullptr, 'H'},
    {"t_block", required_argument, nullptr, 't'},
    {"help",    no_argument,       nullptr, 'h'},
    {nullptr,   no_argument,       nullptr, 0}
  };

  int option_index;
  int opt;
  while((opt = getopt_long(argc, argv, short_opts, long_opts, &option_index)) != -1){
    switch(opt){
    case 'g': gsq     = std::stod(optarg); break;
    case 'N': Nf      = std::stoi(optarg); break;
    case 'n': nu0     = std::stod(optarg); break;
    case 'v': nu1     = std::stod(optarg); break;
    case 'H': nhits   = std::stoi(optarg); break;
    case 't': t_block = std::stoi(optarg); break;
    case 'h':
    case '?':
    default:  PrintHelp(); break;
    }
  }
}
//------------------------------------------

void compute_disc_corr(const std::vector<Complex>& L,
                       std::vector<double>& G_real,
                       std::vector<double>& G_imag){
  std::fill(G_real.begin(), G_real.end(), 0.0);
  std::fill(G_imag.begin(), G_imag.end(), 0.0);
  for(int t0=0; t0<Comp::Nt; t0++){
    for(int dt=0; dt<Comp::Nt; dt++){
      const Complex contrib = L[t0] * std::conj(L[(t0+dt) % Comp::Nt]);
      G_real[dt] += contrib.real() / Comp::Nt;
      G_imag[dt] += contrib.imag() / Comp::Nt;
    }
  }
}

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  double gsq    = 8.0;
  int    Nf     = 2;
  double nu0    = 1.0;
  double nu1    = 1.0;
  int    nhits  = 1;
  int    t_block = 8;

  ParseArgs(argc, argv, gsq, Nf, nu0, nu1, nhits, t_block);
  std::cout << "# gsq = " << gsq
            << " Nf = "      << Nf
            << " nu0 = "     << nu0
            << " nu1 = "     << nu1
            << " nhits = "   << nhits
            << " t_block = " << t_block << std::endl;

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].allocate();

  int device;
  CUDA_CHECK(cudaGetDeviceCount(&device));
  cudaDeviceProp device_prop[device];
  cudaGetDeviceProperties(&device_prop[0], 0);
  std::cout << "# dev = " << device_prop[0].name << std::endl;
  CUDA_CHECK(cudaSetDevice(0));
  std::cout << "# (GPU device is set.)" << std::endl;

  constexpr Idx N  = Comp::N;
  constexpr int Nt = Comp::Nt;

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;

  using Gauge=GaugeExt<Base,Nt,Comp::is_compact>;
  using Action=U1WilsonExt<Base>;
  using Rng=ParallelRngExt<Base,Nt>;

  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;

  const double at = 0.2;
  if(Nt!=1) assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);

  Action SW(gsq, at, base);
  std::cout << "# alat = " << base.mean_ell << std::endl;

  std::string dir3, dir4;
  if(Nf==0){
    dir3 = "gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
    dir4 = "disc_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu1"+std::to_string(nu1)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"tb"+std::to_string(t_block)+"/";
  } else {
    dir3 = "Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"/";
    dir4 = "disc_Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)+"nu0"+std::to_string(nu0)+"nu1"+std::to_string(nu1)+"nt"+std::to_string(Comp::Nt)+"L"+std::to_string(Comp::N_REFINE)+"tb"+std::to_string(t_block)+"/";
  }
  std::cout << "# dir3 = " << dir3 << std::endl;
  std::cout << "# dir4 = " << dir4 << std::endl;
  std::filesystem::create_directory(dir4);

  Gauge U(base);
  srand(time(NULL));
  Rng rng(base, rand());

  const double M5 = -1.0;
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu1);
  std::cout << "# DW set. " << std::endl;

  using Fermion=Overlap<WilsonDirac>;
  Fermion D(DW, 21);
  std::cout << "# D set. " << std::endl;

  D.update(U);
  std::cout << "# D updated. " << std::endl;

  auto f_DH  = std::bind(&Fermion::adj_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);
  auto f_DHD = std::bind(&Fermion::DHD_deviceAsyncLaunch, &D, std::placeholders::_1, std::placeholders::_2);

  LinOpWrapper M_DH(f_DH);
  MatPoly op_DH; op_DH.push_back(cplx(1.0), {&M_DH});
  LinOpWrapper M_DHD(f_DHD);
  MatPoly op_DHD; op_DHD.push_back(cplx(1.0), {&M_DHD});

  FermionVector eta;
  FermionVector DH_eta;
  FermionVector phi;
  FermionVector Gamma_phi;
  std::vector<std::vector<Complex>> summed_trace(4, std::vector<Complex>(Comp::Nt, Complex(0.0, 0.0)));
  std::vector<std::vector<double>> G_real(4, std::vector<double>(Comp::Nt, 0.0));
  std::vector<std::vector<double>> G_imag(4, std::vector<double>(Comp::Nt, 0.0));

  const int k_ckpoint = 10;
  const int kmax      = 1000;

  int k_tmp = 0;
  {
    for(k_tmp=k_ckpoint; k_tmp<=kmax; k_tmp+=k_ckpoint){
      const std::string str_lat = dir3+"ckpoint_lat."+std::to_string(k_tmp);
      if(!std::filesystem::exists(str_lat)) break;
    }
    k_tmp -= k_ckpoint;
  }

  int k_init = 0;
  {
    for(k_init=k_ckpoint; k_init<=kmax; k_init+=k_ckpoint){
      const std::string path = dir4+"disc_corr."+std::to_string(k_init)+".h5";
      if(!std::filesystem::exists(path)) break;
    }
    if(k_init > k_ckpoint) k_init -= k_ckpoint;
  }

  for(int k=k_init; k<=k_tmp; k+=k_ckpoint){
    const std::string h5path = dir4+"disc_corr."+std::to_string(k)+".h5";
    if(std::filesystem::exists(h5path)){
      try {
        HighFive::File h5check(h5path, HighFive::File::ReadOnly);
        const std::string last_ds = std::to_string(nhits-1)+"/3/real";
        if(h5check.exist(last_ds)){
          std::cout << "# skipping k = " << k << " (complete)" << std::endl;
          continue;
        }
      } catch(...) {}
    }

    std::cout << "# k = " << k << std::endl;
    U.read(dir3+"ckpoint_lat."+std::to_string(k));
    D.update(U);

    HighFive::File h5file(h5path,
        HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

    for(int h=0; h<nhits; h++){
      for(int ab=0; ab<4; ab++) std::fill(summed_trace[ab].begin(), summed_trace[ab].end(), Complex(0.0, 0.0));
      const int interval = Comp::Nt / t_block;
      for(int t_s=0; t_s<interval; t_s++){
        for(int spin=0; spin<NS; spin++){
          eta.time_spin_dilution(rng, t_s, t_block, spin);
          op_DH.from_cpu<N>(DH_eta.field, eta.field);
          op_DHD.solve<N>(phi.field, DH_eta.field);
          for(int ab=0; ab<4; ab++){
            Gamma_phi = phi;
            Gamma_phi.mult_sigma(ab);
            eta.accumulate_loop_gamma(summed_trace[ab], Gamma_phi, t_s, t_block, spin);
          }
        }
      }
      for(int ab=0; ab<4; ab++){
        compute_disc_corr(summed_trace[ab], G_real[ab], G_imag[ab]);
        const std::string ds_prefix = std::to_string(h) + "/" + std::to_string(ab) + "/";
        h5file.createDataSet(ds_prefix + "real", G_real[ab]);
        h5file.createDataSet(ds_prefix + "imag", G_imag[ab]);
      }
    }
  } // end for k

  for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();

  return 0;
}
