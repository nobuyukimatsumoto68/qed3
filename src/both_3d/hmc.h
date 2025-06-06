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

  CuC *d_eta_saved;
  static constexpr Idx N = Comp::N;


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
  {
    CUDA_CHECK(cudaMalloc(&d_eta_saved, N*CD));
  }

  ~HMC(){
    CUDA_CHECK(cudaFree(d_eta_saved));
  }

  double H() {
    double res = pi.squared_norm();
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
    CUDA_CHECK(cudaMemcpy(this->d_eta_saved, pf->d_eta, N*CD, D2D));

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
      // pf->update_eta();
      CUDA_CHECK(cudaMemcpy(pf->d_eta, this->d_eta_saved, N*CD, D2D));
      fermion->precalc_grad_deviceAsyncLaunch( U, pf->d_eta );
    }
  }

};



template <typename Rng, typename Action, typename Fermion,
          typename Gauge, typename Force, typename PseudoFermion,
          typename Integrator>
struct HMC2 {
  Rng& rng;
  Action* Sg;
  Fermion* fermion; // on device; need .update( Gauge& U )

  Gauge& U;
  Force& pi;
  // PseudoFermion* pf;
  std::vector<PseudoFermion*>& pfs;

  Integrator* integrator;

  std::vector<CuC*> d_eta_saveds;
  static constexpr Idx N = Comp::N;


  HMC2()=delete;

  explicit HMC2(Rng& rng_,
                Action* Sg_,
                Fermion* fermion_,
                Gauge& U_,
                Force& pi_,
                std::vector<PseudoFermion*>& pfs_,
                Integrator* integrator_)
    : rng(rng_)
    , Sg(Sg_)
    , fermion(fermion_)
    , U(U_)
    , pi(pi_)
    , pfs(pfs_)
    , integrator(integrator_)
  {
    d_eta_saveds.resize(pfs.size());
    for(int jf=0; jf<pfs.size(); jf++) CUDA_CHECK(cudaMalloc(&d_eta_saveds[jf], N*CD));
  }

  ~HMC2(){
    // CUDA_CHECK(cudaFree(d_eta_saved));
    for(CuC* d_eta_saved : d_eta_saveds ) CUDA_CHECK(cudaFree(d_eta_saved));
  }

  double H() {
    double res = pi.squared_norm();
    res *= 0.5;
    res += Sg->operator()(U);
    for(PseudoFermion* pf : pfs) res += pf->S();
    return res;
  }

  inline void integrate() {
    // integrator->integrate( U, pi, Sg, fermion, pf );
    integrator->integrate( U, pi, Sg, fermion, pfs );
  }

  void run( double& r,
            double& dH,
            bool& is_accept,
            const bool no_reject = false ) {
    pi.gaussian( rng );

    Gauge U0( U );
    for(PseudoFermion* pf : pfs) pf->gen( rng );
    for(int jf=0; jf<pfs.size(); jf++) CUDA_CHECK(cudaMemcpy(this->d_eta_saveds[jf], pfs[jf]->d_eta, N*CD, D2D));

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
      // pf->update_eta();
      for(int jf=0; jf<pfs.size(); jf++) CUDA_CHECK(cudaMemcpy(pfs[jf]->d_eta, this->d_eta_saveds[jf], N*CD, D2D));
      for(PseudoFermion* pf : pfs) fermion->precalc_grad_deviceAsyncLaunch( U, pf->d_eta );
    }
  }

};







template <typename Rng, typename Action,
          typename Gauge, typename Force>
struct HMCPureGauge {
  Rng& rng;
  Action* Sg;

  Gauge& U;
  Force& pi;

  double tmax;
  int nsteps;
  double tau;

  HMCPureGauge()=delete;

  explicit HMCPureGauge(Rng& rng_,
                        Action* Sg_,
                        Gauge& U_,
                        Force& pi_,
                        const double tmax_=1.0,
                        const int nsteps_=10)
    : rng(rng_)
    , Sg(Sg_)
    , U(U_)
    , pi(pi_)
    , tmax(tmax_)
    , nsteps(nsteps_)
    , tau(tmax/nsteps)
  {
  }

  ~HMCPureGauge(){
  }

  double H() {
    double res = 0.5*pi.squared_norm();
    res += Sg->operator()(U);
    return res;
  }

  void onestep( Gauge& U, Force& pi, Action* Sg ) const {
    Force dSg(U.lattice);

    Sg->get_force( dSg, U );
#ifdef InfoForce
    dSg.print2log_norm( "# Sg : " );
#endif
    pi += -0.5*tau * ( dSg );

    U += tau * pi;

    Sg->get_force( dSg, U );
#ifdef InfoForce
    dSg.print2log_norm( "# Sg : " );
#endif
    pi += -0.5*tau * ( dSg );
  }

  void integrate() const {
    for(int n=0; n<nsteps; n++) onestep( U, pi, Sg );
  }

  void run( double& r,
            double& dH,
            bool& is_accept,
            const bool no_reject = false ) {
    pi.gaussian( rng );

    Gauge U0( U );

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
    }
  }

};
