#pragma once

#include <random>


/*
  Gauge objects should have:
  Gauge operator+=(const Force& rhs);
*/

/*
  Force objects should have:

  operator+=(const Force& rhs);
  double square() const;
  void rand();
  Force& operator+=(const Force& rhs);
  friend Force operator*(const double a, Force v);
  // friend Force operator*(Force v, const double a);
*/





template <typename Rng, typename Action, typename Fermion,
          typename Gauge, typename Force, typename PseudoFermion>
struct HMC {
  Rng& rng;
  Action& Sg;
  Fermion& fermion; // on device; need .update( Gauge& U )

  Gauge& U;
  Force& pi;
  PseudoFermion& pf;

  const double stot;
  const int nsteps;
  const double tau;

  HMC()=delete;

  explicit HMC(Rng& rng_,
               Action& Sg_,
               Fermion& fermion_,
               Gauge& U_,
               Force& pi_,
               PseudoFermion& pf_,
               const double stot_=1.0, const int nsteps_=10)
    : rng(rng_)
    , Sg(Sg_)
    , fermion(fermion_)
    , U(U_)
    , pi(pi_)
    , pf(pf_)
    , stot(stot_)
    , nsteps(nsteps_)
    , tau(stot/nsteps)
  {}

  double H() {
    double res = 0.0;
    for(const auto elem : pi ) res += elem*elem;
    res *= 0.5;
    res += Sg(U);
    res += pf.S();
    return res;
  }

  void leapfrog_explicit_singlestep() {
    Force dSg(U.lattice), dSf(U.lattice);

    Sg.get_force( dSg, U ); pf.get_force( dSf, U );
#ifdef InfoForce
    dSg.print2log_norm( "# Sg : " );
    dSf.print2log_norm( "# Sf : " );
#endif
    pi += -0.5*tau * ( dSg + dSf );

    U += tau * pi;
    fermion.update( U ); pf.update_eta();

    Sg.get_force( dSg, U ); pf.get_force( dSf, U );
#ifdef InfoForce
    dSg.print2log_norm( "# Sg : " );
    dSf.print2log_norm( "# Sf : " );
#endif
    pi += -0.5*tau * ( dSg + dSf );
  }

  void leapfrog_explicit() {
    for(int n=0; n<nsteps; n++) leapfrog_explicit_singlestep();
  }

  void run( double& r,
            double& dH,
            bool& is_accept,
            const bool no_reject = false ) {
    pi.gaussian( rng );

    Gauge U0( U );
    pf.gen( rng );

    const double h0 = H();
    leapfrog_explicit();
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
    }
  }

};
