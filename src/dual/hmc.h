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
          typename Gauge, typename Force, typename PseudoFermion,
          typename Integrator>
struct HMC {
  Rng& rng;
  Action* Sg;
  Fermion* fermion; // on device; need .update( Gauge& U )

  Gauge& U;
  Force& pi;
  PseudoFermion* pf;

  Integrator* integrator;

  HMC()=delete;

  explicit HMC(Rng& rng_,
               Action* Sg_,
               Fermion* fermion_,
               Gauge& U_,
               Force& pi_,
               PseudoFermion* pf_,
               Integrator* integrator_)
    : rng(rng_)
    , Sg(Sg_)
    , fermion(fermion_)
    , U(U_)
    , pi(pi_)
    , pf(pf_)
    , integrator(integrator_)
  {}

  ~HMC(){}

  double H() {
    double res = 0.0;
    for(const auto elem : pi ) res += elem*elem;
    res *= 0.5;
    res += Sg->operator()(U);
    res += pf->S();
    return res;
  }

  inline void integrate() {
    integrator->integrate( U, pi, Sg, fermion, pf );
  }

  void run( double& r,
            double& dH,
            bool& is_accept,
            const bool no_reject = false ) {
    pi.gaussian( rng );

    Gauge U0( U );
    pf->gen( rng );

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
      pf->update_eta();
      fermion->precalc_grad_deviceAsyncLaunch( U, pf->d_eta );
    }
  }

};
