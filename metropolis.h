#pragma once

template<class Action, class Field>
struct Metropolis {
  const Action& S;
  const double width;

  Metropolis()=delete;
  Metropolis(const Metropolis& other)=delete;

  Metropolis( const Action& S_, const double width_ )
    : S(S_)
    , width(width_)
  {}

  Metropolis& operator=(const Metropolis&)=delete;

  double operator()( Field& U ) const {
    int n_accept = 0;
    int n_tries = 0;

    for (int i=0; i<U.lattice.n_links; i++) {
      Field Up(U);
      Up[i] += U.lattice.rng.RandReal(-1.0, 1.0) * width;

      const double dS = S(Up) - S(U);
      // const double dS = S.local(i, Up) - S.local(i, U);

      const double rate = U.lattice.rng.RandReal();
      if ( rate < std::exp(-dS) ) {
	U[i] = Up[i];
	n_accept++;
      }
      n_tries++;
    }
    return (1.0*n_accept) / (1.0*n_tries);
  }

};
