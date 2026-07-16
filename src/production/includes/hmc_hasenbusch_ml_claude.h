#pragma once

// hmc_hasenbusch_ml_claude.h (Hasenbusch C6, hasenbusch_massless_impl_plan_claude.md)
//
// MULTI-TIMESCALE (nested Omelyan MN2) integrator, mimicking Grid's
//   Grid/qcd/hmc/integrators/Integrator_algorithm.h  MinimumNorm2::step()
//   (commit 7647576, ~L140). Recursion structure copied verbatim: per level loop over `multiplier`,
//   eps /= 2*multiplier per nesting level, kick pattern lambda,(1-2lambda),lambda with the first/last
//   boundary kicks merged (mm = last_step ? 1 : 2), and the FINEST level (fl = last) does the actual
//   U drift (update_U) while coarser levels recurse. lambda = 0.1931833275037836.
//
// Levels are ordered OUTER (0, coarsest/smallest force/most expensive) -> INNER (fl, gauge, finest).
// For the massless Hasenbusch ladder {0,0.1,0.4} the first try is 3 levels:
//   level 0 = {light-pair: ratio frames 0..K-1}  (few evals),
//   level 1 = {heavy frame K},
//   level 2 = {gauge}  (most evals, does update_U).
// Reduces byte-for-byte to the 2-level MinimumNorm2Hasenbusch when levels = [fermion(0..K), gauge]
// with multipliers {MDsteps-equiv, nsteps_inner} -- the regression gate.

#include <vector>

// ---- polymorphic MD level: computes ITS force group into dS given U ----
template<typename Gauge, typename Force>
struct MDLevel {
  int multiplier;
  unsigned long long cg_iters = 0;   // _claude: CG iterations attributed to THIS level (per-stage cost tuning)
  MDLevel( int m ) : multiplier(m) {}
  virtual ~MDLevel() {}
  virtual void get_force( Force& dS, const Gauge& U ) = 0;   // fermion: update D + eta + force; gauge: gauge force
};

// gauge level (finest): dS = gauge force (no D update needed)
template<typename Action, typename Gauge, typename Force>
struct GaugeLevel : public MDLevel<Gauge,Force> {
  Action* Sg;
  GaugeLevel( Action* Sg_, int m ) : MDLevel<Gauge,Force>(m), Sg(Sg_) {}
  void get_force( Force& dS, const Gauge& U ) override { Sg->get_force( dS, U ); }
};

// fermion level: a frame RANGE [i_lo,i_hi] of every HasenbuschPF stack in pfs. Re-solves only those
// frames' FORCE eta and sums their (Term A + Term B) force. Updates the FORCE operator Df(U) first
// (needed before any fermion force). _claude (two-operator split): the level drives the FORCE op only;
// the ACTION op + its eta are refreshed by HMCHasenbuschML::run just before the accept/reject H().
template<typename Fermion, typename PFptr, typename Gauge, typename Force>
struct FermionGroupLevel : public MDLevel<Gauge,Force> {
  Fermion* fermion_force;
  std::vector<PFptr>& pfs;
  int i_lo, i_hi;
  Force tmp;   // scratch for summing over stacks

  template<typename Base>
  FermionGroupLevel( Fermion* ff, std::vector<PFptr>& p, int lo, int hi, int m, Base& base )
    : MDLevel<Gauge,Force>(m), fermion_force(ff), pfs(p), i_lo(lo), i_hi(hi), tmp(base) {}

  void get_force( Force& dS, const Gauge& U ) override {
    const unsigned long long c0 = get_cg_iters();   // _claude: per-stage CG accounting
    fermion_force->update( U );           // Df(U) current before any fermion force (lambda_max frozen)
    bool first = true;
    for( PFptr pf : pfs ){
      pf->update_eta_force_frames( i_lo, i_hi );
      if(first){ pf->get_force_frames( dS,  U, i_lo, i_hi ); first = false; }
      else     { pf->get_force_frames( tmp, U, i_lo, i_hi ); dS += tmp; }
    }
    this->cg_iters += get_cg_iters() - c0;
  }
};

// ---- MinimumNorm2ML: Grid Integrator_algorithm.h MinimumNorm2::step, recursion mimicked verbatim ----
template<typename Gauge, typename Force>
struct MinimumNorm2ML {
  double trajL;
  int MDsteps;
  std::vector<MDLevel<Gauge,Force>*> levels;   // outer(0) .. inner(fl=gauge)
  const double lambda = 0.1931833275037836;

  MinimumNorm2ML( double trajL_, int MDsteps_, std::vector<MDLevel<Gauge,Force>*> levels_ )
    : trajL(trajL_), MDsteps(MDsteps_), levels(levels_) {}

  void step( Gauge& U, Force& pi, int level, int _first, int _last ) {
    const int fl = (int)levels.size() - 1;
    double eps = trajL / MDsteps * 2.0;
    for(int l=0; l<=level; ++l) eps /= 2.0 * levels[l]->multiplier;

    const int multiplier = levels[level]->multiplier;
    Force dS( U.lattice );
    for(int e=0; e<multiplier; ++e){
      const int first_step = _first && (e==0);
      const int last_step  = _last  && (e==multiplier-1);

      if(first_step){
        levels[level]->get_force( dS, U );
        pi += -lambda*eps * dS;
      }

      if(level==fl) U += 0.5*eps * pi;
      else          step( U, pi, level+1, first_step, 0 );

      levels[level]->get_force( dS, U );
      pi += -(1.0-2.0*lambda)*eps * dS;

      if(level==fl) U += 0.5*eps * pi;
      else          step( U, pi, level+1, 0, last_step );

      const int mm = last_step ? 1 : 2;
      levels[level]->get_force( dS, U );
      pi += -lambda*eps*(double)mm * dS;
    }
  }

  void integrate( Gauge& U, Force& pi ) {
    for(int s=0; s<MDsteps; ++s) step( U, pi, 0, (s==0), (s==MDsteps-1) );
  }
};


// ---- HMC using the ML integrator (levels bake in Sg/fermion/pfs; integrate takes only U,pi) ----
template<typename Rng, typename Action, typename Fermion,
         typename Gauge, typename Force, typename PFptr, typename Integrator>
struct HMCHasenbuschML {
  Rng& rng;
  Action* Sg;
  Fermion* fermion;
  Gauge& U;
  Force& pi;
  std::vector<PFptr>& pfs;
  Integrator* integrator;

  HMCHasenbuschML( Rng& rng_, Action* Sg_, Fermion* fermion_, Gauge& U_, Force& pi_,
                   std::vector<PFptr>& pfs_, Integrator* integrator_ )
    : rng(rng_), Sg(Sg_), fermion(fermion_), U(U_), pi(pi_), pfs(pfs_), integrator(integrator_) {}

  double H() {
    double res = 0.5*pi.squared_norm();
    res += Sg->operator()(U);
    for( PFptr pf : pfs ) res += pf->S();
    return res;
  }

  void run( double& r, double& dH, bool& is_accept, const bool no_reject=false ) {
    pi.gaussian( rng );
    Gauge U0( U );
    for( PFptr pf : pfs ) pf->gen( rng );

    const double h0 = H();
    integrator->integrate( U, pi );
    // _claude (two-operator split): the MD only refreshed the FORCE op/eta. Refresh the ACTION op + its
    // eta at the FINAL U so the accept/reject H() uses the accurate n=21 action (self-contained, no force
    // side effect). h0's action eta came from gen() at U0 (still current there).
    fermion->update( U );
    for( PFptr pf : pfs ) pf->update_eta();
    const double h1 = H();

    dH = h1 - h0;
    r = std::min( 1.0, std::exp(-dH) );
    const double a = rng.uniform();
    if( a < r || no_reject ){ is_accept=true; }
    else { is_accept=false; U = U0; fermion->update( U ); }
  }
};
