// condensate_deter_free_claude.cu
// -----------------------------------------------------------------------------
// DETERMINISTIC (exact) free-field cross-check for the condensate traces computed stochastically in
// jj_corr_dilute_claude.cu (Eqs. 1.23 PS, 1.55 FS of qed3int_v2-14).  The exact area-weighted traces are
// the stochastic estimator with a COMPLETE unit-source basis {e_x} (the N_hits -> infinity limit, since
// sum_x e_x e_x^dag = identity), using the IDENTICAL overlap operators -> apples-to-apples.
//
// Spacetime average uses the simplicial site-area weights A_n = dual_areas[n] (diag(A_n) before the trace):
//   etadag_xi        = <eta^dag xi>             = -tr[A D_m^{-1}]                          (ALL mass cases)
//   xidag_eta        = <xi^dag eta>             : massless/m_F = conj(etadag_xi);
//                                                 m_P = -(1+m_P)^{-1} tr[A tilde D_{m_P}^{-dag}]
//   xidag_1mDdag_eta = <xi^dag (1-D_ov^dag) eta>: massless/m_F = -conj(tr[A (1-D_ov) D_m^{-1}]);
//                                                 m_P = -(1+m_P)^{-1} tr[A (1-D_ov^dag) tilde D_{m_P}^{-dag}]
// Analysis: sigma_PS = etadag_xi + xidag_eta ; sigma_FS = etadag_xi - xidag_1mDdag_eta ; spacetime-average
//   by /Vst, Vst = Nt * sum_n dual_areas[n].   (Mass-case propagators: Eqs. 3.50/3.51, 3.57/3.58, 3.60/3.61.)
//
// Method (free U=1): for x=0..N-1, psi_x = D_m^{-1} e_x; accumulate A_{n(x)} psi_x[x] and
//   A_{n(x)} [(1-D_ov)psi_x]_x.  For m_P also phimm_x = tilde D_{m_P}^{-dag} e_x and the (1-D_ov^dag) leg.
// SPEEDUP (free field): the propagator is t-translation-invariant, so the diagonal D_m^{-1}(x,x) depends
//   only on the SPATIAL site, not the timeslice.  We thus sweep only ONE timeslice (x=0..Nx-1) and multiply
//   the traces by Nt -- O(Nx) solves instead of O(N).  A t=1 diagonal is checked against t=0 to confirm.
// Output: data_free_vmRe<..>vmIm<..>/condensate_deter_L<L>/cond.h5 (same keys as the dilute h0/condensate/*).
// Use the SAME mass as the sea quarks.  Plan: condensate_impl_plan_claude.md (OC3).

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
#include <memory>
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
  const double TOL_OUTER=1.0e-10;   // tight: exact reference (vs the dilute's looser correlator tol)
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

#include "sparse_dirac_claude.h"
#include "matpoly_claude.h"

#include "overlap_wmass_claude.h"
#include "blocked_mat_claude.h"

#include <getopt.h>

static void write_scalar(HighFive::File& h5, const std::string& key, const Complex& c){
  h5.createDataSet(key+"/real", std::vector<double>{ c.real() });
  h5.createDataSet(key+"/imag", std::vector<double>{ c.imag() });
}

int main(int argc, char** argv){
  // ---- args: --mass-re, --mass-im, --nu0, --out-tag (positional gsq/Nf retained for the op set) ----
  double gsq=1.0; int Nf=2; double nu0=1.0, nu1=-1.0;
  double mass_re=0.0, mass_im=0.0;
  std::string out_tag="";
  int npole=51;   // Zolotarev pole count (CLI --npole); default 51 for the high-precision contact test
  static struct option longopts[] = {
    {"mass-re",required_argument,0,'r'},{"mass-im",required_argument,0,'i'},
    {"nu0",required_argument,0,'u'},{"out-tag",required_argument,0,'O'},
    {"npole",required_argument,0,'p'},{0,0,0,0}
  };
  int c, idx;
  while((c=getopt_long(argc,argv,"r:i:u:O:p:",longopts,&idx))!=-1){
    switch(c){
      case 'r': mass_re=std::stod(optarg); break;
      case 'i': mass_im=std::stod(optarg); break;
      case 'u': nu0=std::stod(optarg); break;
      case 'O': out_tag=optarg; break;
      case 'p': npole=std::stoi(optarg); break;
    }
  }
  if(nu1 < 0.0) nu1 = nu0;

  const Complex valence_mass(mass_re, mass_im);
  const bool parity = (std::abs(mass_im) > 1.0e-15) && (std::abs(mass_re) <= 1.0e-15);
  std::cout << "# condensate dense free check: mass=("<<mass_re<<","<<mass_im<<")"
            << (parity?"  [parity m_P]":"  [massless/m_F]") << std::endl;

  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();
  int device; CUDA_CHECK(cudaGetDeviceCount(&device));
  cudaDeviceProp dp; cudaGetDeviceProperties(&dp,0);
  std::cout << "# dev = " << dp.name << std::endl;
  CUDA_CHECK(cudaSetDevice(0));

  constexpr Idx N  = Comp::N;
  constexpr Idx Nx = Comp::Nx;

  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;
  using Fermion=OverlapWMass<WilsonDirac>;

  Base base(Comp::N_REFINE);
  const double M5=-1.0, at=0.2;
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu1);
  Gauge U(base);                       // free field U=1
  std::cout << "# lattice + DW set (free U=1)." << std::endl;

  // overlap: D = D_ov (massless GW), D_m = D_ov+m, tilde D_{m_P} = D_ov + m/(1-m)
  std::cout << "# Zolotarev npole = " << npole << std::endl;
  Fermion D   (DW, Complex(0.0), npole);
  Fermion Dm  (DW, valence_mass, npole);
  Fermion Dtil(DW, valence_mass / (Complex(1.0) - valence_mass), npole);

  auto f_D    = std::bind(&Fermion::mult_deviceAsyncLaunch_ms, &D,    std::placeholders::_1, std::placeholders::_2);
  auto f_DH   = std::bind(&Fermion::adj_deviceAsyncLaunch_ms,  &D,    std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_D(f_D), M_DH(f_DH);
  MatPoly op_oneMinusD; op_oneMinusD.push_back(cplx(1.0),{}); op_oneMinusD.push_back(cplx(-1.0),{&M_D});  // 1 - D_ov
  MatPoly op_DH; op_DH.push_back(cplx(1.0),{&M_DH});          // D_ov^dag

  auto f_DmH  = std::bind(&Fermion::adj_deviceAsyncLaunch_ms,  &Dm,   std::placeholders::_1, std::placeholders::_2);
  auto f_Dmsq = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms,  &Dm,   std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_DmH(f_DmH), M_Dmsq(f_Dmsq);
  MatPoly op_DmH; op_DmH.push_back(cplx(1.0),{&M_DmH});
  MatPoly op_Dmsq; op_Dmsq.push_back(cplx(1.0),{&M_Dmsq});

  auto f_tilDm   = std::bind(&Fermion::mult_deviceAsyncLaunch_ms, &Dtil, std::placeholders::_1, std::placeholders::_2);
  auto f_tilDmsq = std::bind(&Fermion::DDH_deviceAsyncLaunch_ms,  &Dtil, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_tilDm(f_tilDm), M_tilDmsq(f_tilDmsq);
  MatPoly op_tilDm;   op_tilDm.push_back(cplx(1.0),{&M_tilDm});
  MatPoly op_tilDmsq; op_tilDmsq.push_back(cplx(1.0),{&M_tilDmsq});

  D.update(U); Dm.update(U); Dtil.update(U);
  std::cout << "# overlap operators set + updated (U=1)." << std::endl;

  // site-area weights (bare local current carries no kappa): A_n = dual_areas[n]
  const int n_sites = static_cast<int>(base.n_sites);
  std::vector<double> w_site(n_sites);
  double sumA=0.0;
  for(int n=0;n<n_sites;n++){ w_site[n]=base.dual_areas[n]; sumA+=w_site[n]; }
  const double Vst = (double)Comp::Nt * sumA;       // spacetime volume for the average

  FermionVector ex, psi, tmp, work, phimm;          // unit source + scratch
  Complex acc_etadag_phi(0,0), acc_etadag_1mD_phi(0,0);
  Complex acc_etadag_phimm(0,0), acc_etadag_1mDdag_phimm(0,0);

  Complex diag00_t0(0,0);    // (D_m^{-1})_xx at (site0,spin0,t=0) -- for the t-invariance check
  auto t0 = std::chrono::steady_clock::now();
  // ONE-TIMESLICE sweep (free t-invariance): x = NS*n + i in [0,Nx) is the t=0 block; n(x)=x/NS.
  // Full-volume trace = Nt * (this one-timeslice trace) -- applied after the loop.
  for(Idx x=0; x<Nx; x++){
    const int n = (int)(x/NS);
    const double A = w_site[n];

    memset(ex.field, 0, Comp::N*CD);
    ex.field[x] = Complex(1.0, 0.0);

    // psi = D_m^{-1} e_x  (op_DmH RHS + op_Dmsq CG), as in the dilute's phi' solve
    op_DmH.from_cpu<N>(tmp.field, ex.field);
    op_Dmsq.solve<N>(psi.field, tmp.field, Comp::TOL_OUTER);
    if(x==0) diag00_t0 = psi.field[0];
    acc_etadag_phi += A * psi.field[x];                       // A_x (D_m^{-1})_xx

    op_oneMinusD.from_cpu<N>(work.field, psi.field);          // work = (1 - D_ov) psi
    acc_etadag_1mD_phi += A * work.field[x];                  // A_x [(1-D_ov)D_m^{-1}]_xx

    if(parity){
      // phimm = tilde D_{m_P}^{-dag} e_x  (op_tilDm RHS + op_tilDmsq CG; normal op -> gives ^{-dag})
      op_tilDm.from_cpu<N>(tmp.field, ex.field);
      op_tilDmsq.solve<N>(phimm.field, tmp.field, Comp::TOL_OUTER);
      acc_etadag_phimm += A * phimm.field[x];                 // A_x (tilde D^{-dag})_xx
      op_DH.from_cpu<N>(work.field, phimm.field);             // work = D_ov^dag phimm
      acc_etadag_1mDdag_phimm += A * (phimm.field[x] - work.field[x]);  // A_x [(1-D_ov^dag) tilde D^{-dag}]_xx
    }

    if((x+1)%16==0 || x+1==Nx){
      const double s = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
      std::cout << "#   col "<<(x+1)<<"/"<<Nx<<" (t=0 block)  ["<<s<<" s]" << std::endl;
    }
  }
  // free-field t-translation invariance: full-volume trace = Nt * (one-timeslice trace)
  acc_etadag_phi          *= (double)Comp::Nt;
  acc_etadag_1mD_phi      *= (double)Comp::Nt;
  acc_etadag_phimm        *= (double)Comp::Nt;
  acc_etadag_1mDdag_phimm *= (double)Comp::Nt;

  // sanity check: the diagonal at t=1 (site0,spin0) must equal t=0 (else t-invariance is broken)
  {
    const Idx xchk = Nx;   // (site0,spin0,t=1)
    memset(ex.field, 0, Comp::N*CD);
    ex.field[xchk] = Complex(1.0, 0.0);
    op_DmH.from_cpu<N>(tmp.field, ex.field);
    op_Dmsq.solve<N>(psi.field, tmp.field, Comp::TOL_OUTER);
    const Complex diag_t1 = psi.field[xchk];
    std::cout << "# t-invariance check: D_m^{-1} diag(site0,spin0)  t0="<<diag00_t0
              << "  t1="<<diag_t1<<"  |diff|="<<std::abs(diag_t1-diag00_t0)
              << "  (should be ~CG tol)" << std::endl;
  }

  // assemble the three bilinears (same per-mass-case logic as jj_corr_dilute)
  const Complex etadag_xi = -acc_etadag_phi;
  Complex xidag_eta, xidag_1mDdag_eta;
  if(parity){
    const Complex inv1pmP = Complex(1.0,0.0)/(Complex(1.0,0.0)+valence_mass);
    xidag_eta        = -inv1pmP*acc_etadag_phimm;
    xidag_1mDdag_eta = -inv1pmP*acc_etadag_1mDdag_phimm;
  } else {
    xidag_eta        =  std::conj(etadag_xi);
    xidag_1mDdag_eta = -std::conj(acc_etadag_1mD_phi);
  }
  const Complex sigma_PS = etadag_xi + xidag_eta;
  const Complex sigma_FS = etadag_xi - xidag_1mDdag_eta;
  std::cout << "# etadag_xi="<<etadag_xi<<"  xidag_eta="<<xidag_eta
            << "  xidag_1mDdag_eta="<<xidag_1mDdag_eta << std::endl;
  std::cout << "# sigma_PS(sum)="<<sigma_PS<<"  -> avg="<<sigma_PS/Vst
            << "  ;  sigma_FS(sum)="<<sigma_FS<<"  -> avg="<<sigma_FS/Vst << std::endl;

  // ---- exact m.i. contact subtraction (condensate_contact_massive_claude.md Sec 10) ----
  // csub_PS: 2 (m0,mF), 1+1/(1+m_P) (mP).   csub_FS: 2 (m0), 2-m_F (mF), 2-(m_P/(1+m_P))^2 (mP).
  // The symmetry-protected residuals should vanish; any leftover that does NOT shrink with npole is a
  // genuine lattice O(m) (term-2) effect, not a finite-Zolotarev sign-function artifact.
  const Complex one(1.0,0.0);
  Complex csub_PS, csub_FS;
  if(std::abs(mass_re)<=1e-15 && std::abs(mass_im)<=1e-15){ csub_PS=2.0; csub_FS=2.0; }          // m0
  else if(parity){ csub_PS = one + one/(one+valence_mass);                                       // mP
                   const Complex r = valence_mass/(one+valence_mass); csub_FS = 2.0 - r*r; }
  else { csub_PS = 2.0; csub_FS = 2.0 - valence_mass; }                                          // mF
  const Complex PSsub = sigma_PS/Vst + csub_PS;
  const Complex FSsub = sigma_FS/Vst + csub_FS;
  std::cout << "# [npole="<<npole<<"] CONTACT-SUBTRACTED:  sigma_PS_sub="<<PSsub
            << "   sigma_FS_sub="<<FSsub
            << "   (csub_PS="<<csub_PS<<", csub_FS="<<csub_FS<<")" << std::endl;

  const std::string esnid = std::string("free")
                          + "_vmRe"+std::to_string(mass_re)+"vmIm"+std::to_string(mass_im);
  const std::string dir_out = "data_"+esnid+"/condensate_deter"+out_tag+"_L"+std::to_string(Comp::N_REFINE)+"/";
  std::filesystem::create_directories(dir_out);
  const std::string h5path = dir_out+"cond.h5";
  const std::string h5tmp  = h5path+".tmp";
  {
    HighFive::File h5(h5tmp, HighFive::File::Truncate);
    h5.createDataSet("N",        std::vector<int>{(int)Comp::N});
    h5.createDataSet("Nt",       std::vector<int>{Comp::Nt});
    h5.createDataSet("n_sites",  std::vector<int>{n_sites});
    h5.createDataSet("N_REFINE", std::vector<int>{Comp::N_REFINE});
    h5.createDataSet("mass_re",  std::vector<double>{mass_re});
    h5.createDataSet("mass_im",  std::vector<double>{mass_im});
    h5.createDataSet("parity",   std::vector<int>{parity?1:0});
    h5.createDataSet("Vst",      std::vector<double>{Vst});     // = Nt * sum_n dual_areas[n]
    write_scalar(h5, "condensate/etadag_xi",        etadag_xi);
    write_scalar(h5, "condensate/xidag_eta",        xidag_eta);
    write_scalar(h5, "condensate/xidag_1mDdag_eta", xidag_1mDdag_eta);
    h5.createDataSet("complete", std::vector<int>{1});
  }
  std::filesystem::rename(h5tmp, h5path);
  std::cout << "# wrote " << h5path << std::endl;

  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].deallocate();
  return 0;
}
