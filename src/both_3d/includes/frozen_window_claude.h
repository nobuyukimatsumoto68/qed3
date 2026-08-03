#pragma once

#include <cassert>

// frozen_window_claude.h (_claude, per NM 2026-07-13)
//
// FROZEN Zolotarev window (lambda_min, lambda_max) per lattice refinement L, the single source of truth
// for the massless overlap runs. Chosen from the Wilson low-mode distribution (eig_wilson_lowmode scan,
// gsq8/12, Nf2/6): lambda_min sits below the observed smallest singular value of D_W, lambda_max above the
// largest, so the sign approximation covers the ensemble. k = lambda_min/lambda_max, n = 21 poles:
//
//   L=1: (0.1,   13)   k = 7.69e-3   Delta = 2.54e-7
//   L=2: (0.06,  8)    k = 7.50e-3   Delta = 2.72e-7
//   L=4: (0.008, 5)    k = 1.60e-3   Delta = 7.08e-6
//
// Delta from zolotarev_delta_claude.cpp (n=21). Use with OverlapWMass::set_lambda(lmin,lmax), which
// rebuilds the Zolotarev fit on k and freezes both edges (no per-config recompute).
//
// NOTE (coverage, measured extrema, worst case Nf2 gsq8): L1 lambda_max=13 < observed 14.6; L2
// lambda_min=0.06 > observed 0.0485; L4 lambda_max=5 < observed 6.5 -- accepted per NM ("looks good").
// Revisit if production configs push outside these edges (watch the "eval below Zolotarev window" warning).

inline void frozen_window( const int L, double& lmin, double& lmax ){
  if( L==1 ){ lmin = 0.1;   lmax = 13.0; }
  else if( L==2 ){ lmin = 0.06;  lmax = 8.0;  }
  else if( L==3 ){ lmin = 0.015; lmax = 8.0;  }  // L3 PROVISIONAL (mirror production/includes; freeze from measured extrema)
  else if( L==4 ){ lmin = 0.008; lmax = 5.0;  }
  else{ assert( false && "frozen_window: no frozen (lambda_min,lambda_max) defined for this L" ); }
}
