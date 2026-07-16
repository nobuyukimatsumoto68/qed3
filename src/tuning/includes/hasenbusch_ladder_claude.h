#pragma once

#include <vector>
#include <complex>
#include <cassert>

// hasenbusch_ladder_claude.h (_claude, per NM 2026-07-13; ladders+steps finalized 2026-07-13)
//
// Production Hasenbusch mass ladders + per-stage MD step counts per lattice refinement L, the single
// source of truth for the massless overlap HMC. Ladder entries [1..K] are M_mass COEFFICIENTS (applied via
// OverlapWMass::set_mass_coeff): frame i>0 has m_L = c_i * diag(A_y/Abar). Frame 0 = the PHYSICAL target
// (0 = massless; set via set_mass).
//
//   L=1: 2-stage {0, 0.5}             steps {2,2}     (mult-1)
//   L=2: 2-stage {0, 0.5}             steps {3,3}     (mult-1)
//   L=4: 4-stage {0, 0.16, 0.32, 0.5} steps {2,2,2,4} (heavy frame sub-cycled x2)
// (trajectory length s_tot = tau = 1.2 for L1/L2, 1.0 for L4 = hasenbusch_tau(L). Multi-timescale -> MinimumNorm2ML.)
//
// c_1 near the low-mode tail (nudged UP toward the geo-mean c_1 ~ sqrt(lambda_min lambda_max); below the
// tail is wasteful -- degenerate split + still-singular D_1). L4 4-stage (c = 0.16, 0.32, 0.5).
// c_1: L1 0.5, L2 0.5, L4 0.16. Chosen from the Wilson low-mode distribution + blowup analysis (frozen_window_claude.md).
// Steps: per-stage counts outer(massless)->inner(heavy), nested in MinimumNorm2ML (level-l multiplier =
// steps[l]/steps[l-1]); verified via test_hasenbusch_tune_claude.cu (force norms L2/Linf + per-stage CG + Cost).

// _claude(TUNING sandbox, NM 2026-07-15): STARTING params for the CORRECTED (weaker) gsq range
//   L1 gsq{0.5,1.0,1.5}, L2 gsq{1.0,2.0,3.0}, L4 gsq{2.0,4.0,6.0}, massless. Zolotarev windows +
//   n_action/n_force UNCHANGED (previously frozen). Tune ladder/steps/tau below from dH/acceptance.
//   L1: {0.0, 1.0}      steps {2,2}   tau 1.0  MG 20
//   L2: {0.0, 1.0}      steps {3,3}   tau 1.0  MG 20
//   L4: {0.0, 0.4, 1.0} steps {4,4,4} tau 0.8  MG 100  (SHIP 2026-07-16: 3-frame; the gauge force dominates
//       so a very fine gauge substep MG100 + short tau0.8 give the safety margin, incl. Nf=6. Cheap: gauge
//       substeps have no CG. 5-frame/4-frame ladder experiments were a cost regression, discarded.)
inline std::vector<std::complex<double>> hasenbusch_ladder( const int L ){
  using C = std::complex<double>;
  if( L==1 ) return { C(0.0,0.0), C(1.0,0.0) };
  if( L==2 ) return { C(0.0,0.0), C(1.0,0.0) };
  if( L==4 ) return { C(0.0,0.0), C(0.4,0.0), C(1.0,0.0) };
  assert( false && "hasenbusch_ladder: no ladder defined for this L" );
  return {};
}

// per-stage MD step counts (outer=massless -> inner=heavy), size = ladder size (K+1). Must integer-nest.
inline std::vector<int> hasenbusch_steps( const int L ){
  if( L==1 ) return { 2, 2 };
  if( L==2 ) return { 3, 3 };
  if( L==4 ) return { 4, 4, 4 };
  assert( false && "hasenbusch_steps: no steps defined for this L" );
  return {};
}

// trajectory length tau (= s_tot). L1 -> 1.0, L2 -> 1.0, L4 -> 0.8 (SHIP 2026-07-16: shorter tau for the
// stiff L4, paired with the fine gauge MG100).
inline double hasenbusch_tau( const int L ){
  if( L==1 ) return 1.0;
  if( L==2 ) return 1.0;
  return 0.8;   // L4
}

// gauge substeps MG per innermost(heaviest)-frame MD step. L1/L2 -> 20; L4 -> 100 (SHIP 2026-07-16: the gauge
// force dominates -- a very fine gauge substep controls dH cheaply (no CG per gauge step), safe incl. Nf=6).
inline int hasenbusch_mg( const int L ){
  if( L==1 ) return 20;
  if( L==2 ) return 20;
  return 100;   // L4
}

// _claude (two-operator split-pole HMC, two_operator_force_impl_plan_claude.md): the FORCE overlap
// operator uses a CRUDER Zolotarev sign approx (fewer poles n_f) than the ACTION operator (n=21). The MD
// force is -dS_f/dU built from D_f; heatbath + accept/reject use the accurate n=21 operator. Exact by the
// Metropolis correction (Duane-Kennedy-Pendleton-Roweth 1987). Per-L n_f, tuned by the n_f scan
// (test_hasenbusch_npole_claude.cu, scan {5,7,9,11,21}). Default 11 for all L until the scan fixes them.
inline int hasenbusch_nforce( const int L ){
  if( L==1 ) return 11;
  if( L==2 ) return 11;
  if( L==4 ) return 11;
  assert( false && "hasenbusch_nforce: no n_force defined for this L" );
  return 21;
}

// _claude: ACTION overlap pole count per L (heatbath + accept/reject; the accurate operator). L1/L2 need
// fewer poles (wider window k) than L4 for the same Delta. n=25 (L1/L2) / 31 (L4).
inline int hasenbusch_naction( const int L ){
  if( L==1 ) return 25;
  if( L==2 ) return 25;
  if( L==4 ) return 31;
  assert( false && "hasenbusch_naction: no n_action defined for this L" );
  return 31;
}
