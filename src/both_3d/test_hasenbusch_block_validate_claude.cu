// test_hasenbusch_block_validate_claude.cu
// _claude: BLOCK-vs-SERIAL NSTACK = NF/2 PARITY GATE for the Nf-block Hasenbusch pseudofermion.
//
// Proves that HasenbuschPFBlock<Fermion,Force,NSTACK> (includes/pseudofermion_hasenbusch_block_claude.h)
// reproduces the SERIAL HasenbuschPF<Fermion,Force> (includes/pseudofermion_hasenbusch_claude.h) at
// NSTACK = NF/2 flavors, on ONE thermalized config, driven with the SAME phi. The block path packs the
// Nf/2 = NSTACK stacks into one mrhs block and runs every action / force apply through BlockedMat /
// BlockedForce (block CG + inner Zolotarev multishift, B. Jegerlehner, hep-lat/9612014). For NSTACK>1 the
// columns genuinely INTERACT inside the block solves, so this exercises REAL block packing (e.g. L2 Nf6 ->
// NSTACK=3). The reference is NSTACK independent SERIAL stacks summed: the block action S must match the
// SUM of the per-stack S, and the block MD force -dS_f/dU must match the SUM of the per-stack forces.
// A mismatch localizes to a frame (via the per-frame S / force printout) or to Term A vs Term B (per-frame
// force). NSTACK=1 (NF=2, the default) reduces to the original single-column no-op-block gate.
//
// This driver writes NO gauge configs (parity gate only; it reads one existing config read-only).
//
// Physics: M. Hasenbusch, Phys. Lett. B 519 (2001) 177, hep-lat/0107019 (mass preconditioning),
// generalized to a mass ladder (Hasenbusch-Jansen, hep-lat/0211042). Block path: BlockedMat / BlockedForce
// (B. Jegerlehner, hep-lat/9612014, shifted-system Krylov / block CG).
//
// Usage: ./bin [gsq=8] [Nf=NF] [config_k=0(=last)]
// Build: needs -DLREF= and -DNF= (NSTACK = NF/2; default NF=2 -> NSTACK=1). The runtime Nf argument must
//        equal the compile-time NF. It sums the NSTACK serial stacks to compare against the one block PF.

#include <typeinfo>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <filesystem>
#include <cstdint>
#include <complex>
#include <array>
#include <vector>
#include <map>
#include <memory>
#include <string>
#include <sstream>
#include <cmath>
#include <Eigen/Dense>

using Double = double;
using Idx = std::int32_t;
using Complex = std::complex<double>;
using Face = std::vector<Idx>;
using MS=Eigen::Matrix2cd;
using VD=Eigen::Vector2d;
using VE=Eigen::Vector3d;
using VC=Eigen::VectorXcd;

static constexpr int NS = 2;
static constexpr int DIM = 2;
static constexpr Complex I = Complex(0.0, 1.0);

#define InfoDelta

namespace Comp{
  constexpr bool is_compact=false;
  constexpr int NPARALLEL_DUPDATE=1;
  constexpr int NPARALLEL=NPARALLEL_DUPDATE;
  constexpr int NSTREAMS=NPARALLEL_DUPDATE;
  constexpr int NPARALLEL_GAUGE=16;
  constexpr int NPARALLEL_SORT=16;

#ifndef LREF
#define LREF 1
#endif
  constexpr int N_REFINE=LREF;
  constexpr int NS=2;
  constexpr int Nt=128;

  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE;
  constexpr Idx Nx=NS*N_SITES;
  constexpr Idx N=Nx*Nt;

  const double TOL_INNER=1.0e-9;
  const double TOL_OUTER=1.0e-8;
}

// NSTACK = NF/2 flavors are packed into one block PF; default NF=2 -> NSTACK=1 (original gate).
#ifndef NF
#define NF 2
#endif

const std::string dir = "../../geometry/data/";
#include "../../geometry/geodesic.h"

#include "timer.h"
#include "s2n_simp.h"
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

#include "sparse_matrix_claude.h"
#include "dirac_simp.h"
#include "dirac_ext.h"
#include "sparse_dirac_claude.h"
#include "matpoly_claude.h"
#define GRAD_L4
#include "includes/overlap_wmass_claude.h"
#include "frozen_window_claude.h"
#include "hasenbusch_ladder_claude.h"
#include "pseudofermion_hasenbusch_claude.h"        // SERIAL HasenbuschPF (reference path)
#include "blocked_mat_claude.h"                     // BlockMemPool + BlockedMat (mrhs; matpoly pulled the block kernels)
#include "overlap_force_Nfblock_claude.h"           // BlockedForce (Term A + block Term B); AFTER blocked_mat
#include "pseudofermion_hasenbusch_block_claude.h"  // HasenbuschPFBlock (Nf-BLOCK; the path under test)

// L-infinity (max |component|) force norm over the whole lattice (all links + temporal sites, all t).
// Complements Force::squared_norm() (L2). Templated on the driver types (mirror of the tune driver helper).
template<typename ForceT, typename BaseT>
static double force_inf_norm( const ForceT& f, const BaseT& base ){
  double mx = 0.0;
  for(int t=0; t<Comp::Nt; t++){
    for(Idx il=0; il<base.n_links; il++) mx = std::max( mx, std::abs(f.sp(t,il)) );
    for(Idx ix=0; ix<base.n_sites; ix++) mx = std::max( mx, std::abs(f.tp(t,ix)) );
  }
  return mx;
}

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(6);

  static constexpr int NSTACK = NF/2;   // number of independent Hasenbusch stacks packed into the block PF

  double gsq      = 8.0;
  int    Nf       = NF;
  int    config_k = 0;      // 0 -> last config of the massless ensemble
  if(argc>1) gsq      = atof(argv[1]);
  if(argc>2) Nf       = atoi(argv[2]);
  if(argc>3) config_k = atoi(argv[3]);
  assert(Nf==NF);           // runtime Nf must match compile-time NF -> Nf/2 = NSTACK flavors

  for(int i=0;i<Comp::NSTREAMS;i++) d_MemorySets[i].allocate();

  constexpr Idx N = Comp::N;
  using Base=S2Simp;
  using WilsonDirac=DiracExt<Base, DiracS2Simp>;
  using Gauge=GaugeExt<Base,Comp::Nt,Comp::is_compact>;
  using Force=Gauge;
  using Action=U1WilsonExt<Base>;
  using Fermion=OverlapWMass<WilsonDirac>;
  using Rng=ParallelRngExt<Base, Comp::Nt>;

  Base base(Comp::N_REFINE);
  const double M5 = -1.0, at = 0.2, nu0 = 1.0;
  WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);

  Rng rng(base);
  Gauge U(base);

  // ACTION operator D (n_act poles) + FORCE operator Df (n_force poles), SAME D_W + SAME frozen window;
  // Df on the narrowed window [2*lambda_min, lambda_max] (two-operator split, hmc_hasenbusch_block_claude.cu).
  double lmin_frozen, lmax_frozen;
  frozen_window( Comp::N_REFINE, lmin_frozen, lmax_frozen );

  Fermion D(DW, Complex(0.0,0.0), hasenbusch_naction(Comp::N_REFINE), 0.001);   // action op; window overridden below
  D.set_lambda( lmin_frozen, lmax_frozen );

  const int    nforce    = hasenbusch_nforce( Comp::N_REFINE );
  const double lmin_force = 2.0 * lmin_frozen;
  Fermion Df(DW, Complex(0.0,0.0), nforce, 0.001);                              // force op; narrowed window
  Df.set_lambda( lmin_force, lmax_frozen );

  Action SW( gsq, at, base );

  // ----- massless config dir + a thermalized config (default = last) -----
  const std::string cfgdir = "Nf"+std::to_string(Nf)+"_gsq"+std::to_string(gsq)+"at"+std::to_string(at)
                           + "nu0"+std::to_string(nu0)+"nt"+std::to_string(Comp::Nt)
                           + "L"+std::to_string(Comp::N_REFINE)+"/";
  if(config_k<=0){
    for(int k=1; k<=100000; k++)
      if(std::filesystem::exists(cfgdir+"ckpoint_lat."+std::to_string(k))) config_k = k;
  }
  const std::string str_lat = cfgdir+"ckpoint_lat."+std::to_string(config_k);
  if(!std::filesystem::exists(str_lat)){
    std::cout << "# ERROR: config not found: " << str_lat << std::endl;
    return 1;
  }
  U.read( str_lat );
  D.update( U );
  Df.update( U );

  const std::vector<Complex> masses = hasenbusch_ladder( Comp::N_REFINE );   // frame 0 = 0 (massless)
  const int K = (int)masses.size()-1;

  std::cout << "# ================ Hasenbusch BLOCK-vs-SERIAL PARITY  L=" << Comp::N_REFINE
            << " gsq=" << gsq << " Nf=" << Nf << " (NSTACK=" << NSTACK << ") ================" << std::endl;
  std::cout << "# frozen window: lambda_min=" << lmin_frozen << " lambda_max=" << lmax_frozen
            << " Delta=" << D.Delta() << "  (force n=" << nforce << " lmin_f=" << lmin_force
            << " Delta_f=" << Df.Delta() << ")" << std::endl;
  std::cout << "# ladder {";
  for(size_t i=0;i<masses.size();i++) std::cout << masses[i].real() << (i+1<masses.size()?",":"");
  std::cout << "}  config k=" << config_k << std::endl;

  // ----- build both PFs on the SAME (D, Df, masses) -----
  // SERIAL reference = NSTACK independent stacks (each holds CUDA resources / non-copyable -> unique_ptr).
  std::vector<std::unique_ptr<HasenbuschPF<Fermion,Force>>> spf;
  for(int c=0; c<NSTACK; c++)
    spf.emplace_back( std::make_unique<HasenbuschPF<Fermion,Force>>( D, Df, masses, base ) );
  HasenbuschPFBlock<Fermion,Force,NSTACK> bpf( D, Df, masses, base );  // Nf-block under test

  // (1) SERIAL heatbath: gen EACH stack (independent phi per stack), then re-derive that stack's force eta.
  // gen() already solved the action eta. (2) COPY the SAME per-stack phi into the block PF: column c of the
  // block d_phi[i] (an N*NSTACK block, column-major [c*N + i]) starts at offset c*N, so it takes spf[c]'s phi.
  // Then re-derive the block action eta + force eta on the packed phi.
  for(int c=0; c<NSTACK; c++){
    spf[c]->gen( rng );
    spf[c]->update_eta_force_frames( 0, K );
  }
  for(int i=0; i<=K; i++){
    for(int c=0; c<NSTACK; c++){
      CUDA_CHECK( cudaMemcpy( bpf.d_phi[i] + (size_t)c*N, spf[c]->d_phi[i], N*CD, D2D ) );  // column c -> offset c*N
    }
  }
  bpf.update_eta();
  bpf.update_eta_force_frames( 0, K );

  // ============ (3) ACTION S parity ============
  // Serial reference = SUM of the NSTACK per-stack actions; block packs the same stacks into one S().
  double Ss = 0.0;
  for(int c=0; c<NSTACK; c++) Ss += spf[c]->S();
  const double Sb   = bpf.S();
  const double dS   = std::abs( Ss - Sb );
  const double relS = ( std::abs(Ss) > 0.0 ) ? dS/std::abs(Ss) : dS;
  std::cout << "\n# ---- (3) action S parity ----" << std::endl;
  std::cout << "#   S_serial=" << Ss << "  S_block=" << Sb << "  |diff|=" << dS << "  rel=" << relS << std::endl;
  for(int i=0; i<=K; i++){
    double sfi = 0.0;
    for(int c=0; c<NSTACK; c++) sfi += spf[c]->S_frame(i);
    const double bfi = bpf.S_frame(i);
    const double dfi = std::abs( sfi - bfi );
    const double rfi = ( std::abs(sfi) > 0.0 ) ? dfi/std::abs(sfi) : dfi;
    std::cout << "#     frame " << i << " (m=" << masses[i].real() << "): S_serial=" << sfi
              << " S_block=" << bfi << " |diff|=" << dfi << " rel=" << rfi << std::endl;
  }

  // ============ (4) FORCE parity (whole stack, frames [0,K]) ============
  // Both use the FORCE eta (d_eta_f), already solved above. Serial reference = SUM over the NSTACK stacks;
  // block computes the packed force in one pass. Diff via the Force operators: fd = fS - fB.
  Force fS( base );
  for(int c=0; c<NSTACK; c++){
    Force fc( base );
    spf[c]->get_force_frames( fc, U, 0, K );
    if(c==0) fS  = fc;
    else     fS += fc;
  }
  Force fB( base );
  bpf.get_force_frames( fB, U, 0, K );
  Force fd( base );
  fd  = fS;
  fd += -1.0 * fB;
  const double L2fS  = std::sqrt( fS.squared_norm() );
  const double L2fB  = std::sqrt( fB.squared_norm() );
  const double L2d   = std::sqrt( fd.squared_norm() );
  const double Linfd = force_inf_norm( fd, base );
  const double relF  = ( L2fS > 0.0 ) ? L2d/L2fS : L2d;
  std::cout << "\n# ---- (4) force parity (whole stack) ----" << std::endl;
  std::cout << "#   L2(fS)=" << L2fS << "  L2(fB)=" << L2fB << std::endl;
  std::cout << "#   L2(diff)=" << L2d << "  Linf(diff)=" << Linfd << "  rel=L2(diff)/L2(fS)=" << relF << std::endl;

  // ============ (5) per-frame FORCE parity (localizes a mismatch to a frame / Term A vs Term B) ============
  std::cout << "\n# ---- (5) per-frame force parity ----" << std::endl;
  for(int i=0; i<=K; i++){
    Force fSi( base );
    for(int c=0; c<NSTACK; c++){
      Force fci( base );
      spf[c]->get_force_frames( fci, U, i, i );
      if(c==0) fSi  = fci;
      else     fSi += fci;
    }
    Force fBi( base );
    bpf.get_force_frames( fBi, U, i, i );
    Force fdi( base );
    fdi  = fSi;
    fdi += -1.0 * fBi;
    const double l2s  = std::sqrt( fSi.squared_norm() );
    const double l2d  = std::sqrt( fdi.squared_norm() );
    const double lid  = force_inf_norm( fdi, base );
    const double rel  = ( l2s > 0.0 ) ? l2d/l2s : l2d;
    std::cout << "#     frame " << i << " (m=" << masses[i].real() << "): L2(fS)=" << l2s
              << " L2(diff)=" << l2d << " Linf(diff)=" << lid << " rel=" << rel << std::endl;
  }

  // ============ PASS/FAIL summary (rel-diff < 1e-6 = solver tolerance) ============
  const double tol = 1.0e-6;
  const bool pass_S = relS < tol;
  const bool pass_F = relF < tol;
  std::cout << "\n# ---- summary (tol=" << tol << ") ----" << std::endl;
  std::cout << "#   ACTION : " << (pass_S ? "PASS" : "FAIL") << "  (rel=" << relS << ")" << std::endl;
  std::cout << "#   FORCE  : " << (pass_F ? "PASS" : "FAIL") << "  (rel=" << relF << ")" << std::endl;
  std::cout << "#   OVERALL: " << ((pass_S && pass_F) ? "PASS" : "FAIL") << std::endl;

  std::cout << "\n# ================ block-vs-serial parity done ================" << std::endl;
  return 0;
}
