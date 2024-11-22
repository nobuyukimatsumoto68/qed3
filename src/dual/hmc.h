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



template <class Force, class Gauge, class Action>
struct HMCPureGauge {
  ParallelRng& rng;
  const Action& S;
  const double stot;
  const int nsteps;
  const double tau;

  HMCPureGauge(ParallelRng& rng, const Action& S_, const double stot_=1.0, const int nsteps_=10)
    : rng(rng)
    , S(S_)
    , stot(stot_)
    , nsteps(nsteps_)
    , tau(stot/nsteps)
  {}

  double H( const Force& pi, const Gauge& W ) const {
    double res = 0.0;
    for(const auto elem : pi ) res += elem*elem;
    res *= 0.5;
    res += S(W);
    return res;
  }

  void leapfrog_explicit_singlestep( Force& pi, Gauge& W ) const {
    pi += -0.5*tau * S.d(W);
    W += tau * pi;
    pi += -0.5*tau * S.d(W);
  }

  void leapfrog_explicit( Force& pi, Gauge& W ) const {
    for(int n=0; n<nsteps; n++) leapfrog_explicit_singlestep(pi,W);
  }

  void run( Gauge& W0,
	    double& r,
	    double& dH,
	    bool& is_accept,
	    const bool no_reject = false ) const {
    Force pi( W0.lattice );
    pi.gaussian( rng );
    Gauge W( W0 );
    const double h0 = H(pi, W);
    leapfrog_explicit( pi, W );
    const double h1 = H(pi, W);

    dH = h1-h0;
    r = std::min( 1.0, std::exp(-dH) );
    const double a = rng.uniform();
    if( a < r || no_reject ){
      W0 = W;
      is_accept=true;
    }
    else is_accept=false;
  }

};





template <class Force, class Gauge, class Action, class Fermion>
struct HMC {
  ParallelRng& rng;
  const Action& S;
  const Fermion& D;
  const double stot;
  const int nsteps;
  const double tau;

  PseudoFermion phi;

  HMC(ParallelRng& rng, const Action& S_, const Fermion& D_,
      const double stot_=1.0, const int nsteps_=10)
    : rng(rng)
    , S(S_)
    , D(D_)
    , stot(stot_)
    , nsteps(nsteps_)
    , tau(stot/nsteps)
    , phi(D)
  {}

  double H( const Force& pi, const Gauge& U ) {
    assert( phi.flag );

    double res = 0.0;
    for(const auto elem : pi ) res += elem*elem;
    res *= 0.5;
    res += S(U);
    res += phi.S();
    return res;
  }

  void leapfrog_explicit_singlestep( Force& pi, Gauge& U ) {
    assert( phi.flag );

    pi += -0.5*tau * ( S.d(U) + phi.dS(U) );

    U += tau * pi;
    phi.calc_eta( U );

    pi += -0.5*tau * ( S.d(U) + phi.dS(U) );
  }

  void leapfrog_explicit( Force& pi, Gauge& U ) {
    for(int n=0; n<nsteps; n++) leapfrog_explicit_singlestep(pi, U );
  }

  void run( Gauge& U0,
	    double& r,
	    double& dH,
	    bool& is_accept,
	    const bool no_reject = false ) {
    Force pi( U0.lattice );
    pi.gaussian( rng );

    Gauge U( U0 );
    phi.gen( U, rng );

    const double h0 = H(pi, U);
    leapfrog_explicit( pi, U );
    const double h1 = H(pi, U);

    dH = h1-h0;
    r = std::min( 1.0, std::exp(-dH) );
    const double a = rng.uniform();
    if( a < r || no_reject ){
      U0 = U;
      is_accept=true;
    }
    else is_accept=false;

    phi.flag=false;
  }

};
