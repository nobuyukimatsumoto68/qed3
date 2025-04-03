#pragma once

#include <random>


struct SingleRng {
  std::mt19937_64 mt;

  std::normal_distribution<double> dist_gaussian;
  std::uniform_real_distribution<> dist_01;
  std::uniform_int_distribution<int> dist_z2;

  SingleRng()
    : dist_z2(1,2)
  {}

  double gaussian(){ return dist_gaussian(mt); }
  double uniform(){ return dist_01(mt); }
  double z2(){ return dist_z2(mt); }
  void reseed(const int seed) { mt.seed(seed); }
};


template<typename Lattice>
struct ParallelRng {
  const Lattice& lattice;
  SingleRng master;
  std::vector<SingleRng> links;
  std::vector<SingleRng> sites;

  explicit ParallelRng(const Lattice& lattice_, const int seed=0 )
    : lattice(lattice_)
    , links( lattice.n_links )
    , sites( lattice.n_sites )
  {
    reseed(seed);
  }

  double gaussian_link( const int il ){ return links[il].gaussian(); }
  double gaussian_site( const int ix ){ return sites[ix].gaussian(); }
  double gaussian(){ return master.gaussian(); }

  double uniform_link( const int il ){ return links[il].uniform(); }
  double uniform_site( const int ix ){ return sites[ix].uniform(); }
  double uniform(){ return master.uniform(); }

  double z2_site( const Idx ix ){ return sites[ix].z2(); }

  void reseed(const int seed) {
    master.reseed(seed);
    for(auto& elem : links) elem.reseed( master.mt() );
    for(auto& elem : sites) elem.reseed( master.mt() );
  }
};


template<typename Lattice, int Nt>
struct ParallelRngExt {
  const Lattice& lattice;
  SingleRng master;
  std::vector<std::vector<SingleRng>> links; // [t][il]
  std::vector<std::vector<SingleRng>> sites; // [t][ix]

  explicit ParallelRngExt(const Lattice& lattice_, const int seed=0 )
    : lattice(lattice_)
    , links( Nt, std::vector<SingleRng>(lattice.n_links) )
    , sites( Nt, std::vector<SingleRng>(lattice.n_sites) )
  {
    reseed(seed);
  }

  double gaussian_link( const int s, const int il ){ return links[s][il].gaussian(); }
  double gaussian_site( const int s, const int ix ){ return sites[s][ix].gaussian(); }
  double gaussian(){ return master.gaussian(); }

  double uniform_link( const int s, const int il ){ return links[s][il].uniform(); }
  double uniform_site( const int s, const int ix ){ return sites[s][ix].uniform(); }
  double uniform(){ return master.uniform(); }

  double z2_site( const int s, const Idx ix ){ return sites[s][ix].z2(); }

  void reseed(const int seed) {
    master.reseed(seed);
    for(auto& v : links) for(auto& elem : v) elem.reseed( master.mt() );
    for(auto& v : sites) for(auto& elem : v) elem.reseed( master.mt() );
  }

  void fill_gaussian( std::vector<Complex>& xi ){
    for(int s=0; s<Comp::Nt; s++) {
      for(Idx ix=0; ix<Comp::N_SITES; ix++) {
        for(int a=0; a<Comp::NS; a++) xi[Comp::Nx*s + NS*ix+a] = ( gaussian_site(s,ix)
                                                                   + I*gaussian_site(s,ix) ) / std::sqrt(2.0);
      }
    }
  }
};


