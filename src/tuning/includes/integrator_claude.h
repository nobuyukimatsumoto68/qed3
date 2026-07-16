#pragma once

// Copy of integrator.h with #ifdef PrintH blocks added to
// MinimumNorm2::integrate(). Prints H = S_gauge + S_ferm + S_mom
// at each outer-step boundary (11 points for nsteps=5).
// Format: # PrintH <step> : Sg=... Sf=... Smom=... H=...

// struct IntegratorBase {
//   virtual void operator()( Gauge& U, Force& pi,
//                            const Action& Sg, const Fermion& fermion,
//                            const PseudoFermion& pf ) const = 0;
// };

// struct ExplicitLeapfrog {
//   const double stot;
//   const int nsteps;
//   const double tau;

//   ExplicitLeapfrog( const double stot_=1.0, const int nsteps_=10)
//     : stot(stot_)
//     , nsteps(nsteps_)
//     , tau(stot/nsteps)
//   {}

//   template <typename Action, typename Fermion,
//             typename Gauge, typename Force, typename PseudoFermion>
//   void onestep( Gauge& U, Force& pi,
//                 Action* Sg, Fermion* fermion,
//                 PseudoFermion* pf ) const {
//     Force dSg(U.lattice), dSf(U.lattice);

//     Sg->get_force( dSg, U ); pf->get_force( dSf, U );
// #ifdef InfoForce
//     dSg.print2log_norm( "# Sg : " );
//     dSf.print2log_norm( "# Sf : " );
// #endif
//     pi += -0.5*tau * ( dSg + dSf );

//     U += tau * pi;
//     fermion->update( U ); pf->update_eta();
//     fermion->precalc_grad_deviceAsyncLaunch( U, pf->d_eta );

//     Sg->get_force( dSg, U ); pf->get_force( dSf, U );
// #ifdef InfoForce
//     dSg.print2log_norm( "# Sg : " );
//     dSf.print2log_norm( "# Sf : " );
// #endif
//     pi += -0.5*tau * ( dSg + dSf );
//   }

//   template <typename Action, typename Fermion,
//             typename Gauge, typename Force, typename PseudoFermion>
//   void integrate( Gauge& U, Force& pi,
//                   Action* Sg, Fermion* fermion,
//                   PseudoFermion* pf ) const {
//     for(int n=0; n<nsteps; n++) onestep( U, pi, Sg, fermion, pf );
//   }
// };




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

  template <typename Action, typename Fermion,
            typename Gauge, typename Force, typename PseudoFermion>
  void integrate( Gauge& U, Force& pi,
                  const Action* Sg, Fermion* fermion,
                  std::vector<PseudoFermion*>& pfs ) const {
    Force dSg(U.lattice), dSf(U.lattice);

    // 0th
    for(PseudoFermion* pf : pfs){
      pf->update_eta(); // may be omitted; but for safety
      fermion->precalc_grad_deviceAsyncLaunch( U, pf->d_eta );
      pf->get_force( dSf, U );
#ifdef InfoForce
      dSf.print2log_norm( "# Sf : " );
#endif
      pi += -0.5*tau * dSf;
    }

    for(int n=0; n<nsteps; n++) {
      //------------------
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
      } //------------------
      fermion->update( U );

      for(PseudoFermion* pf : pfs){
        pf->update_eta();
        fermion->precalc_grad_deviceAsyncLaunch( U, pf->d_eta );

        pf->get_force( dSf, U );
#ifdef InfoForce
        dSf.print2log_norm( "# Sf : " );
#endif
        if(n!=nsteps-1) pi += -1.0*tau * dSf;
        else pi += -0.5*tau * dSf;
      }
    }
  }

};




// inf_norm_force: max absolute value over all spatial and temporal components.
template<typename Force>
double inf_norm_force(const Force& f) {
  double val = 0.0;
  for(const auto& row : f.spatial)  for(double v : row) val = std::max(val, std::abs(v));
  for(const auto& row : f.temporal) for(double v : row) val = std::max(val, std::abs(v));
  return val;
}

// https://github.com/paboyle/Grid/blob/4a0aaf07868a39ffafffba7496c75a3c03bc7aab/Grid/qcd/hmc/integrators/Integrator_algorithm.h#L140
struct MinimumNorm2 {
  const double stot;
  int nsteps;
  const int nsteps_inner;

  // tau is recomputed inside integrate() as stot/nsteps so that
  // changing nsteps at runtime takes effect without reconstructing.
  const double lambda = 0.1931833275037836;

  MinimumNorm2( const double stot_=1.0, const int nsteps_=10, const int nsteps_inner_=10)
    : stot(stot_)
    , nsteps(nsteps_)
    , nsteps_inner(nsteps_inner_)
  {}


  template <typename Action, typename Fermion,
            typename Gauge, typename Force, typename PseudoFermion>
  void integrate( Gauge& U, Force& pi,
                  const Action* Sg, Fermion fermion,
                  std::vector<PseudoFermion>& pfs ) const {
    const double tau = stot / nsteps;  // recomputed so runtime nsteps changes take effect
    Force dSg(U.lattice), dSf(U.lattice);

#ifdef PrintH
    int ph_step = -1;  // -1 = true h0 before any kick
    auto print_H = [&](){
      double ph_Sg   = Sg->operator()(U);
      double ph_Sf   = 0.0;
      for(PseudoFermion pf : pfs) ph_Sf += pf->S();
      double ph_Smom = 0.5 * pi.squared_norm();
      std::clog << "# PrintH " << ph_step << " :"
                << " Sg=" << ph_Sg
                << " Sf=" << ph_Sf
                << " Smom=" << ph_Smom
                << " H="   << ph_Sg + ph_Sf + ph_Smom
                << std::endl;
      ph_step++;
    };
#endif

#ifdef PrintEtaPhi
    // print_eta_phi: prints ||eta||_2 / ||phi||_2 for one pseudo-fermion.
    // Large ratio => phi couples to near-zero modes of D^dag D.
    // Format: # EtaPhi <outer> <phase:0=init,1=mid,2=fin> <ipf> : <ratio>
    int ep_outer = -1;
    int ep_phase = 0;
    auto print_eta_phi = [&](auto pf, int ipf){
      CuC np, ne;
      pf->Op_DHD.template dot<Comp::N>(&np, pf->d_phi, pf->d_phi);
      pf->Op_DHD.template dot<Comp::N>(&ne, pf->d_eta, pf->d_eta);
      double ratio = std::sqrt(cuCreal(ne) / cuCreal(np));
      std::clog << "# EtaPhi " << ep_outer << " " << ep_phase << " " << ipf
                << " : " << ratio << std::endl;
    };
#endif

#ifdef PrintHinner
    // print_Hinner: lightweight H monitor inside the inner gauge loop.
    // Skips Sf (no CG); outer=outer step index, half=0/1 for first/second half-sweep.
    auto print_Hinner = [&](int outer, int half, int n_inner){
      double ph_Sg   = Sg->operator()(U);
      double ph_Smom = 0.5 * pi.squared_norm();
      std::clog << "# PrintHinner " << outer << " " << half << " " << n_inner
                << " : Sg=" << ph_Sg
                << " Smom=" << ph_Smom
                << std::endl;
    };
#endif

#ifdef PrintH
    print_H(); // step -1: true h0, before any kick
#endif

    // 0th
#ifdef PrintEtaPhi
    { int _ep_ipf = 0; ep_outer = -1; ep_phase = 0;
#endif
    for(PseudoFermion pf : pfs){
      pf->update_eta(); // may be omitted; but for safety
      fermion->precalc_grad_deviceAsyncLaunch( U, pf->d_eta ); // to omit, Xs in D needs to be flavored
      pf->get_force( dSf, U );
#ifdef InfoForce
      dSf.print2log_norm( "# Sf : " );
      std::clog << "# Sf_inf : " << inf_norm_force(dSf) << std::endl;
#endif
      pi += -lambda*tau * dSf;
#ifdef PrintEtaPhi
      print_eta_phi(pf, _ep_ipf++);
#endif
    }
#ifdef PrintEtaPhi
    }
#endif

#ifdef PrintH
    print_H(); // step 0: after initial lambda kick, before first gauge sweep
#endif


    const double tau_inner = 0.5*tau/nsteps_inner;
    for(int n=0; n<nsteps; n++) {
      // forward U by 0.5*tau with gauge force
      // 0th
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
#ifdef PrintHinner
        print_Hinner(n, 0, n_inner);
#endif
      } //------------------
      fermion->update( U );


      // inner (1-2lambda) pi fermion update
#ifdef PrintEtaPhi
      { int _ep_ipf = 0; ep_outer = n; ep_phase = 1;
#endif
      for(PseudoFermion pf : pfs){
        pf->update_eta();
        fermion->precalc_grad_deviceAsyncLaunch( U, pf->d_eta );

        pf->get_force( dSf, U );
#ifdef InfoForce
        dSf.print2log_norm( "# Sf : " );
        std::clog << "# Sf_inf : " << inf_norm_force(dSf) << std::endl;
#endif
        pi += -(1.0-2.0*lambda)*tau * dSf;
#ifdef PrintEtaPhi
        print_eta_phi(pf, _ep_ipf++);
#endif
      }
#ifdef PrintEtaPhi
      }
#endif

#ifdef PrintH
      print_H(); // step 2n+1: after first gauge sweep + middle fermion kick of outer step n
#endif


      // repeat the 0.5*tau step
      // 0th
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
#ifdef PrintHinner
        print_Hinner(n, 1, n_inner);
#endif
      } //------------------
      fermion->update( U );


      // last (2*)lambda pi fermion update
#ifdef PrintEtaPhi
      { int _ep_ipf = 0; ep_outer = n; ep_phase = 2;
#endif
      for(PseudoFermion pf : pfs){
        pf->update_eta();
        fermion->precalc_grad_deviceAsyncLaunch( U, pf->d_eta );

        pf->get_force( dSf, U );
#ifdef InfoForce
        dSf.print2log_norm( "# Sf : " );
        std::clog << "# Sf_inf : " << inf_norm_force(dSf) << std::endl;
#endif
        if(n!=nsteps-1) pi += -2.0*lambda*tau * dSf;
        else pi += -lambda*tau * dSf;
#ifdef PrintEtaPhi
        print_eta_phi(pf, _ep_ipf++);
#endif
      }
#ifdef PrintEtaPhi
      }
#endif

#ifdef PrintH
      print_H(); // step 2n+2: after second gauge sweep + final fermion kick of outer step n
#endif
    }
  }

};
