#pragma once

// hmc_hasenbusch_claude.h (Hasenbusch C2, hasenbusch_massless_impl_plan_claude.md)
//
// Hasenbusch variants of MinimumNorm2 (integrator.h) and HMC2 (hmc.h). The ONLY change vs the
// originals is that the per-frame fermion force is SELF-CONTAINED: HasenbuschPF::get_force does its
// own set_mass + precalc (Term A + Term B) per frame, so the external
//   fermion->precalc_grad_deviceAsyncLaunch( U, pf->d_eta )
// coupling (which assumes ONE operator + ONE standard force with a scalar pf->d_eta) is DROPPED.
// Everything else -- the 2-level Omelyan MN2 structure (fermion frames on the coarse OUTER timescale,
// gauge on the fine INNER timescale), the kick pattern lambda,(1-2lambda),lambda -- is byte-identical
// to MinimumNorm2 / HMC2. All K+1 Hasenbusch frames ride the single outer fermion kick (Grid-canonical
// 2-level split; the 3-level multi-timescale split is deferred chunk C6).
//
// Originals in integrator.h (MinimumNorm2) and hmc.h (HMC2) are left untouched for A/B and rollback.

#include <random>

// ---- integrator: MinimumNorm2 with a self-contained fermion force (no external precalc_grad) ----
struct MinimumNorm2Hasenbusch {
  const double stot;
  const int nsteps;
  const int nsteps_inner;

  const double tau;
  const double lambda = 0.1931833275037836;

  MinimumNorm2Hasenbusch( const double stot_=1.0, const int nsteps_=10, const int nsteps_inner_=10)
    : stot(stot_)
    , nsteps(nsteps_)
    , nsteps_inner(nsteps_inner_)
    , tau(stot/nsteps)
  {}

  template <typename Action, typename Fermion,
            typename Gauge, typename Force, typename PseudoFermion>
  void integrate( Gauge& U, Force& pi,
                  const Action* Sg, Fermion fermion,
                  std::vector<PseudoFermion>& pfs ) const {
    Force dSg(U.lattice), dSf(U.lattice);

    // 0th outer fermion kick
    for(PseudoFermion pf : pfs){
      pf->update_eta();
      pf->get_force( dSf, U );
#ifdef InfoForce
      dSf.print2log_norm( "# Sf : " );
#endif
      pi += -lambda*tau * dSf;
    }

    const double tau_inner = 0.5*tau/nsteps_inner;
    for(int n=0; n<nsteps; n++) {
      // forward U by 0.5*tau with the gauge force (fine inner timescale)
      Sg->get_force( dSg, U );
#ifdef InfoForce
      dSg.print2log_norm( "# Sg : " );
#endif
      pi += -lambda*tau_inner * dSg;

      for(int n_inner=0; n_inner<nsteps_inner; n_inner++) {
        U += 0.5 * tau_inner * pi;
        Sg->get_force( dSg, U );
#ifdef InfoForce
        dSg.print2log_norm( "# Sg : " );
#endif
        pi += - (1.0 - 2.0*lambda)*tau_inner * dSg;

        U += 0.5 * tau_inner * pi;
        Sg->get_force( dSg, U );
#ifdef InfoForce
        dSg.print2log_norm( "# Sg : " );
#endif
        if(n_inner!=nsteps_inner-1) pi += -2.0*lambda*tau_inner * dSg;
        else pi += -lambda*tau_inner * dSg;
      }
      fermion->update( U );

      // inner (1-2lambda) outer fermion kick
      for(PseudoFermion pf : pfs){
        pf->update_eta();
        pf->get_force( dSf, U );
#ifdef InfoForce
        dSf.print2log_norm( "# Sf : " );
#endif
        pi += -(1.0-2.0*lambda)*tau * dSf;
      }

      // repeat the 0.5*tau gauge step
      Sg->get_force( dSg, U );
#ifdef InfoForce
      dSg.print2log_norm( "# Sg : " );
#endif
      pi += -lambda*tau_inner * dSg;

      for(int n_inner=0; n_inner<nsteps_inner; n_inner++) {
        U += 0.5 * tau_inner * pi;
        Sg->get_force( dSg, U );
#ifdef InfoForce
        dSg.print2log_norm( "# Sg : " );
#endif
        pi += - (1.0 - 2.0*lambda)*tau_inner * dSg;

        U += 0.5 * tau_inner * pi;
        Sg->get_force( dSg, U );
#ifdef InfoForce
        dSg.print2log_norm( "# Sg : " );
#endif
        if(n_inner!=nsteps_inner-1) pi += -2.0*lambda*tau_inner * dSg;
        else pi += -lambda*tau_inner * dSg;
      }
      fermion->update( U );

      // last (2*)lambda outer fermion kick
      for(PseudoFermion pf : pfs){
        pf->update_eta();
        pf->get_force( dSf, U );
#ifdef InfoForce
        dSf.print2log_norm( "# Sf : " );
#endif
        if(n!=nsteps-1) pi += -2.0*lambda*tau * dSf;
        else pi += -lambda*tau * dSf;
      }
    }
  }
};


// ---- HMC: HMC2 with the Hasenbusch heatbath (no external precalc_grad after gen) ----
template <typename Rng, typename Action, typename Fermion,
          typename Gauge, typename Force, typename PseudoFermion,
          typename Integrator>
struct HMCHasenbusch {
  Rng& rng;
  Action* Sg;
  Fermion* fermion;

  Gauge& U;
  Force& pi;
  std::vector<PseudoFermion>& pfs;

  Integrator* integrator;

  static constexpr Idx N = Comp::N;

  HMCHasenbusch()=delete;

  explicit HMCHasenbusch(Rng& rng_,
                         Action* Sg_,
                         Fermion* fermion_,
                         Gauge& U_,
                         Force& pi_,
                         std::vector<PseudoFermion>& pfs_,
                         Integrator* integrator_)
    : rng(rng_)
    , Sg(Sg_)
    , fermion(fermion_)
    , U(U_)
    , pi(pi_)
    , pfs(pfs_)
    , integrator(integrator_)
  {}

  ~HMCHasenbusch(){}

  double H() {
    double res = 0.5*pi.squared_norm();
    res += Sg->operator()(U);
    for(PseudoFermion pf : pfs) res += pf->S();
    return res;
  }

  inline void integrate() {
    integrator->integrate( U, pi, Sg, fermion, pfs );
  }

  void run( double& r,
            double& dH,
            bool& is_accept,
            const bool no_reject = false ) {
    pi.gaussian( rng );

    Gauge U0( U );
    for(PseudoFermion pf : pfs) pf->gen( rng );   // Hasenbusch heatbath (self-contained)

    const double h0 = H();
    integrate();
    const double h1 = H();

    dH = h1-h0;
    r = std::min( 1.0, std::exp(-dH) );
    const double a = rng.uniform();
    if( a < r || no_reject ){
      is_accept=true;
    }
    else {
      is_accept=false;
      U = U0;
      fermion->update( U );
    }
  }
};
