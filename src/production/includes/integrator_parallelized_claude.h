#pragma once

// integrator_parallelized_claude.h -- MinimumNorm2Block: a structural copy of MinimumNorm2
// (integrator.h) whose fermion kicks drive the ACTION inversion through a PseudoFermionBlock
// (one block solve_sq -> all NSTACK etas) instead of the serial for(pf : pfs) loop. The FORCE is
// still computed SERIALLY per flavor (fermion->precalc_grad + Force::compute on eta_col(f)) -- the
// blocking of the force is a later step. The gauge sweeps, lambda weights and nsteps_inner cadence
// are byte-identical to MinimumNorm2, so the only numerical change is serial->block on the action solve.
//
// Reference (integrator): 2nd-order minimum-norm (Omelyan/Grid), same as MinimumNorm2.

struct MinimumNorm2Block {
  const double stot;
  const int nsteps;
  const int nsteps_inner;

  const double tau;
  const double lambda = 0.1931833275037836;

  MinimumNorm2Block( const double stot_=1.0, const int nsteps_=10, const int nsteps_inner_=10)
    : stot(stot_)
    , nsteps(nsteps_)
    , nsteps_inner(nsteps_inner_)
    , tau(stot/nsteps)
  {}

  // fermion is the Fermion* (as in HMC2); bpf packs the NSTACK flavors.
  template <typename Action, typename Gauge, typename Force,
            typename Fermion, int NSTACK>
  void integrate( Gauge& U, Force& pi,
                  const Action* Sg, Fermion* fermion,
                  PseudoFermionBlock<Fermion,NSTACK>& bpf ) const {
    Force dSg(U.lattice), dSf(U.lattice);

    // 0th : one block action solve, then serial per-flavor force
    bpf.update_eta();
    for(int f=0; f<NSTACK; f++){
      fermion->precalc_grad_deviceAsyncLaunch( U, bpf.eta_col(f) );
      dSf.compute( U, bpf.eta_col(f), *fermion );
#ifdef InfoForce
      dSf.print2log_norm( "# Sf : " );
#endif
      pi += -lambda*tau * dSf;
    }

    const double tau_inner = 0.5*tau/nsteps_inner;
    for(int n=0; n<nsteps; n++) {
      // forward U by 0.5*tau with gauge force
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
      } //------------------
      fermion->update( U );

      // inner (1-2lambda) fermion kick : block action solve + serial force
      bpf.update_eta();
      for(int f=0; f<NSTACK; f++){
        fermion->precalc_grad_deviceAsyncLaunch( U, bpf.eta_col(f) );
        dSf.compute( U, bpf.eta_col(f), *fermion );
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
      } //------------------
      fermion->update( U );

      // last (2*)lambda fermion kick : block action solve + serial force
      bpf.update_eta();
      for(int f=0; f<NSTACK; f++){
        fermion->precalc_grad_deviceAsyncLaunch( U, bpf.eta_col(f) );
        dSf.compute( U, bpf.eta_col(f), *fermion );
#ifdef InfoForce
        dSf.print2log_norm( "# Sf : " );
#endif
        if(n!=nsteps-1) pi += -2.0*lambda*tau * dSf;
        else pi += -lambda*tau * dSf;
      }
    }
  }

};
