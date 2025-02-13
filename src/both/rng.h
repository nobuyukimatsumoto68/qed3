#pragma once

#include <random>


struct SingleRng {
  std::mt19937_64 mt;

  std::normal_distribution<double> dist_gaussian;
  double gaussian(){ return dist_gaussian(mt); }

  std::uniform_real_distribution<> dist_01;
  double uniform(){ return dist_01(mt); }

  std::uniform_int_distribution<int> dist_z2;

  SingleRng()
    : dist_z2(1,2)
  {}

  double z2(){ return dist_z2(mt); }

  double RandReal(const double a=0., const double b=1.0) {
    return (b-a) * uniform() + a;
  }

  void SeedRng(const int seed) { mt.seed(seed); }
};


template<typename Lattice>
struct ParallelRng {
  const Lattice& lattice;
  const int seed;
  std::mt19937_64 mt;
  std::vector<std::mt19937_64> mt_link;
  std::vector<std::mt19937_64> mt_site;
  std::uniform_int_distribution<int> dist_z2;

  explicit ParallelRng(const Lattice& lattice_, const int seed=0 )
    : lattice(lattice_)
    , seed(seed)
    , mt_link( lattice.n_links )
    , mt_site( lattice.n_sites )
    , dist_z2(1,2)
  {
    mt.seed(seed);
    for(auto& elem : mt_link) elem.seed( mt() );
    for(auto& elem : mt_site) elem.seed( mt() );
  }

  std::normal_distribution<double> dist_gaussian;
  double gaussian_link( const int il ){ return dist_gaussian(mt_link[il]); }
  double gaussian_site( const int ix ){ return dist_gaussian(mt_site[ix]); }
  double gaussian(){ return dist_gaussian(mt); }

  std::uniform_real_distribution<> dist_01;
  double uniform_link( const int il ){ return dist_01(mt_link[il]); }
  double uniform_site( const int ix ){ return dist_01(mt_site[ix]); }
  double uniform(){ return dist_01(mt); }

  double z2_site( const Idx ix ){ return dist_z2(mt_site[ix]); }

  double RandReal(const double a=0., const double b=1.0) {
    return (b-a) * uniform() + a;
  }

  void reseed(const int seed) {
    mt.seed(seed);
    for(auto& elem : mt_link) elem.seed( mt() );
    for(auto& elem : mt_site) elem.seed( mt() );
  }
};
