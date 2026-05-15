#pragma once

template <typename Action>
struct Flow {
  Action* Sg;
  // Gauge& U;
  // Force& pi;

  const double tmax;
  const int nsteps;
  const double tau;

  Flow()=delete;
  Flow(Action* Sg_,
       // Gauge& U_,
       // Force& pi_,
       const double tmax_=1.0,
       const int nsteps_=100)
    : Sg(Sg_)
    , tmax(tmax_)
    , nsteps(nsteps_)
    , tau(tmax/nsteps)
  {}

  template<typename Gauge>
  void operator()( Gauge& U ) const {
    Gauge dSg(U.lattice);
    for(int i=0; i<nsteps; i++){
      Sg->get_spatial( dSg, U );
      U += -tau * dSg;
    }
  }


};
