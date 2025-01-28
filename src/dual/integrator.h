#pragma once

// struct IntegratorBase {
//   virtual void operator()( Gauge& U, Force& pi,
//                            const Action& Sg, const Fermion& fermion,
//                            const PseudoFermion& pf ) const = 0;
// };

struct ExplicitLeapfrog {
  const double stot;
  const int nsteps;
  const double tau;

  ExplicitLeapfrog( const double stot_=1.0, const int nsteps_=10)
    : stot(stot_)
    , nsteps(nsteps_)
    , tau(stot/nsteps)
  {}

  template <typename Action, typename Fermion,
            typename Gauge, typename Force, typename PseudoFermion>
  void onestep( Gauge& U, Force& pi,
                Action* Sg, Fermion* fermion,
                PseudoFermion* pf ) const {
    Force dSg(U.lattice), dSf(U.lattice);

    Sg->get_force( dSg, U ); pf->get_force( dSf, U );
#ifdef InfoForce
    dSg.print2log_norm( "# Sg : " );
    dSf.print2log_norm( "# Sf : " );
#endif
    pi += -0.5*tau * ( dSg + dSf );

    U += tau * pi;
    fermion->update( U ); pf->update_eta();
    fermion->precalc_grad_deviceAsyncLaunch( U, pf->d_eta );

    Sg->get_force( dSg, U ); pf->get_force( dSf, U );
#ifdef InfoForce
    dSg.print2log_norm( "# Sg : " );
    dSf.print2log_norm( "# Sf : " );
#endif
    pi += -0.5*tau * ( dSg + dSf );
  }

  template <typename Action, typename Fermion,
            typename Gauge, typename Force, typename PseudoFermion>
  void integrate( Gauge& U, Force& pi,
                  Action* Sg, Fermion* fermion,
                  PseudoFermion* pf ) const {
    for(int n=0; n<nsteps; n++) onestep( U, pi, Sg, fermion, pf );
  }
};




struct ExplicitLeapfrogML {
  const double stot;
  const int nsteps;
  const int nsteps_inner;

  const double tau;
  const double tau_inner;

  ExplicitLeapfrogML( const double stot_=1.0, const int nsteps_=10, const int nsteps_inner_=10)
    : stot(stot_)
    , nsteps(nsteps_)
    , nsteps_inner(nsteps_inner_)
    , tau(stot/nsteps)
    , tau_inner(tau/nsteps_inner)
  {}


  template <typename Action, typename Fermion,
            typename Gauge, typename Force, typename PseudoFermion>
  void integrate( Gauge& U, Force& pi,
                  const Action* Sg, Fermion* fermion,
                  PseudoFermion* pf ) const {
    Force dSg(U.lattice), dSf(U.lattice);

    // 0th
    pf->get_force( dSf, U );
#ifdef InfoForce
    dSf.print2log_norm( "# Sf : " );
#endif
    pi += -0.5*tau * dSf;

    for(int n=0; n<nsteps; n++) {
      { //------------------
        // 0th
        Sg->get_force( dSg, U );
#ifdef InfoForce
        dSg.print2log_norm( "# Sg : " );
#endif
        pi += -0.5*tau_inner * dSg;
        for(int n_inner=0; n_inner<nsteps_inner; n_inner++) {
          U += tau_inner * pi;

          Sg->get_force( dSg, U );
#ifdef InfoForce
          dSg.print2log_norm( "# Sg : " );
#endif
          if(n_inner!=nsteps_inner-1) pi += -1.0*tau_inner * dSg;
          else pi += -0.5*tau_inner * dSg;
        }
      } //------------------
      fermion->update( U ); pf->update_eta();
      fermion->precalc_grad_deviceAsyncLaunch( U, pf->d_eta );

      pf->get_force( dSf, U );
#ifdef InfoForce
      dSf.print2log_norm( "# Sf : " );
#endif
      if(n!=nsteps-1) pi += -1.0*tau * dSf;
      else pi += -0.5*tau * dSf;
    }
  }
};
