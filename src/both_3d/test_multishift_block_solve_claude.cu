// test_multishift_block_solve_claude.cu
// C6b (device): validate MatPoly::solve_multishift_block (mrhs) against the single-RHS
// MatPoly::solve_multishift, column by column, on the REAL Zolotarev inner-pole system
// (seed A = (1/lambda_max^2) D_W^dag D_W, shifts sigma_m = -k^2/c'_m).
//
//   block:  Aseed.solve_multishift_block(d_Xblk, d_B, M_DW, M_DWH, 1/lmax^2, sigma, npole, nstack)
//   ref:    for each RHS c, Aseed.solve_multishift(d_ref, d_B + c*N, sigma, npole)
// By construction (per-column independence + column freeze-on-convergence), the block column c
// must equal the single-RHS solve on b_c BIT-FOR-BIT. Two checks: nstack=1 (reproduces
// solve_multishift) and nstack=NSTACK. PASS = max|block - single| < 1e-10 over every column.
// Free field (U=1): the solver math is gauge-independent.
//
// Compile: handled by the both_3d Makefile (auto-picks every *.cu).
// Run: ./test_multishift_block_solve_claude.o [--mass-re x] [--mass-im y]

#include <typeinfo>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <complex>
#include <array>
#include <vector>
#include <map>
#include <string>
#include <random>
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
  constexpr int NSTREAMS=4;
  constexpr int NPARALLEL_GAUGE=NPARALLEL_DUPDATE;
  constexpr int NPARALLEL_SORT=NPARALLEL_DUPDATE;

  constexpr int N_REFINE=1;
  constexpr int NS=2;
  constexpr int Nt=128;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;

  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-5;
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
#include "matpoly_claude.h"
#include "overlap_wmass_claude.h"
#include "blocked_mat_claude.h"   // C7: BlockedMat<N,NSTACK,Op>

#include <getopt.h>

// max | block - ref | over nlen complex entries (device blocks; ref single column-major [m*N+i])
static double cmp_block(const CuC* d_blk, const CuC* d_ref, size_t nlen){
  std::vector<Complex> a(nlen), b(nlen);
  CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(a.data()), d_blk, nlen*CD, D2H));
  CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(b.data()), d_ref, nlen*CD, D2H));
  double d=0.0; for(size_t i=0;i<nlen;++i) d=std::max(d, std::abs(a[i]-b[i]));
  return d;
}

int main(int argc, char* argv[]){
  double mass_re=0.0, mass_im=0.0;
  static struct option long_opts[] = {
    {"mass-re", required_argument, nullptr, 'r'},
    {"mass-im", required_argument, nullptr, 'i'},
    {nullptr,0,nullptr,0}
  };
  int opt, idx;
  while((opt=getopt_long(argc,argv,"r:i:",long_opts,&idx))!=-1){
    if(opt=='r') mass_re=std::stod(optarg);
    else if(opt=='i') mass_im=std::stod(optarg);
  }
  const Complex valence_mass(mass_re, mass_im);

  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();

  cudaDeviceProp dp; cudaGetDeviceProperties(&dp,0);
  std::cout << "# dev = " << dp.name << std::endl;
  CUDA_CHECK(cudaSetDevice(0));

  constexpr Idx N = Comp::N;

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;
  using Rng=ParallelRngExt<Base,Comp::Nt>;
  using Fermion=OverlapWMass<WilsonDirac>;

  Base base(Comp::N_REFINE);
  const double M5=-1.0, at=0.2, nu1=1.0;
  if(Comp::Nt!=1) assert(std::sqrt(3.0)*base.mean_ell/at - 4.0/std::sqrt(3.0) > -1.0e-14);
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu1);

  Gauge U(base);            // free field U=1
  Rng rng(base, 1234);

  Fermion Dm(DW, valence_mass, 21);
  Dm.update(U);
  std::cout << "# D_m updated: lambda_min/max="<<Dm.lambda_min<<"/"<<Dm.lambda_max
            << "  size="<<Dm.size<<"  k="<<Dm.k << std::endl;

  const int npole = Dm.size - 1;
  std::vector<double> sigma(npole);
  for(int m=1;m<Dm.size;m++) sigma[m-1] = -Dm.k*Dm.k/Dm.cp[m];
  const double coeff = 1.0/(Dm.lambda_max*Dm.lambda_max);
  MatPoly Aseed; Aseed.push_back(cplx(coeff), {&Dm.M_DW, &Dm.M_DWH});  // handle owner; vec_mats unused by _block

  std::cout << std::scientific << std::setprecision(3);
  std::cout << "# npole="<<npole<<"\n";

  bool all_pass = true;

  // preallocate the thread-owned block scratch once, sized for the largest batch we test (N_SITES).
  // A solve with NSTACK <= cap uses the FRONT portion -> one allocation serves nstack=1/4/12.

  // run the block solve for a COMPILE-TIME NSTACK and compare each column vs single solve_multishift
  auto run = [&]<int NSTACK>(){
    const size_t Nblk = (size_t)N*NSTACK*npole;
    CuC *d_B,*d_Xblk,*d_ref;
    CUDA_CHECK(cudaMalloc(&d_B,   (size_t)N*NSTACK*CD));
    CUDA_CHECK(cudaMalloc(&d_Xblk,Nblk*CD));
    CUDA_CHECK(cudaMalloc(&d_ref, (size_t)N*npole*CD));

    // NSTACK distinct random Z2 source columns
    for(int c=0;c<NSTACK;c++){ FermionVector eta; eta.fill_z2_source(rng);
      CUDA_CHECK(cudaMemcpy(d_B + (size_t)c*N, reinterpret_cast<CuC*>(eta.field), N*CD, H2D)); }

    BlockedMat<N,NSTACK,Fermion> blk(Dm);   // owns block scratch (RAII), wraps Dm
    auto t0=std::chrono::steady_clock::now();
    blk.solve_multishift_block(d_Xblk, d_B, Comp::TOL_INNER);
    CUDA_CHECK(cudaDeviceSynchronize());
    const double t_blk=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();

    double worst=0.0; double t_ref=0.0;
    for(int c=0;c<NSTACK;c++){
      auto tr=std::chrono::steady_clock::now();
      Aseed.solve_multishift<N>(d_ref, d_B + (size_t)c*N, sigma.data(), npole, Comp::TOL_INNER);
      CUDA_CHECK(cudaDeviceSynchronize());
      t_ref += std::chrono::duration<double>(std::chrono::steady_clock::now()-tr).count();
      const double d = cmp_block(d_Xblk + (size_t)c*npole*N, d_ref, (size_t)N*npole);
      worst=std::max(worst,d);
    }
    const bool ok = worst<1e-10; all_pass &= ok;
    std::cout << "# nstack="<<std::setw(2)<<NSTACK<<"  max|block-single|="<<worst
              <<"  ("<<(ok?"PASS":"FAIL")<<")   block="<<t_blk<<"s  ref(serial "<<NSTACK<<" solves)="<<t_ref<<"s\n";
    CUDA_CHECK(cudaFree(d_B)); CUDA_CHECK(cudaFree(d_Xblk)); CUDA_CHECK(cudaFree(d_ref));
  };

  run.operator()<1>();             // must reproduce solve_multishift bit-for-bit
  run.operator()<4>();             // mrhs batch
  run.operator()<Comp::N_SITES>(); // nstack = n_sites (the jj sink-loop batch at this L)

  std::cout << "# C6b RESULT: " << (all_pass ? "ALL PASS" : "FAIL") << std::endl;
  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
  return all_pass ? 0 : 1;
}
