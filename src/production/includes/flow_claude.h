#pragma once

// Wilson gradient-flow smearing. DEFAULT: the flow is driven by the SPATIAL-only action gradient
// (get_spatial) -> it acts within each timeslice (the "on a timeslice" / per-timeslice restriction),
// temporal links are never flowed. With -DFLOW_FULL the flow is driven by the FULL 3D gauge-action
// gradient (get_force = spatial + temporal plaquettes), so temporal links are flowed too and the
// smearing mixes timeslices. Both gradient routines already live in action_ext_claude.h; this header
// only selects which one drives the flow (get_force's spatial block is identical to get_spatial, same
// sign convention, so the swap is consistent).
template <typename Action>
struct Flow {
  Action* Sg;

  const double tmax;
  const int nsteps;
  const double tau;

  Flow()=delete;
  Flow(Action* Sg_,
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
#ifdef FLOW_FULL
      Sg->get_force( dSg, U );    // full 3D Wilson flow (spatial + temporal links)
#else
      Sg->get_spatial( dSg, U );  // spatial-only flow (per-timeslice restriction)
#endif
      U += -tau * dSg;
    }
  }


};
