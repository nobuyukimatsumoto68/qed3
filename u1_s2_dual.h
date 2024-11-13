#pragma once

#include <vector>
#include <cmath>


struct U1onS2 {
  using Link = std::array<int,2>; // <int,int>;
  using Face = std::vector<int>;

  Lattice& lattice;
  std::vector<double> field;

  U1onS2()=delete;

  U1onS2(Lattice& lattice_)
    : lattice(lattice_)
    , field(lattice.n_links, 0.0)
  {}

  U1onS2( const U1onS2& other )
    : lattice(other.lattice)
    , field(other.field)
  {}

  U1onS2 & operator=(const U1onS2&) = delete;

  double operator[](const int i) const { return field[i]; }
  double& operator[](const int i) { return field[i]; }

  double operator()(const Link& ell) const { // recommended
    const int il = lattice.map2il.at(ell);
    const int sign = lattice.map2sign.at(ell);
    return sign * field[il] ;
  }

  double plaquette_angle(const int& i_face) const {
    const Face& face = lattice.faces[i_face];
    double sum = 0.0;
    for(int i=0; i<face.size(); i++) {
      const int ix = face[i];
      const int iy = face[(i+1)%face.size()];
      sum += (*this)( Link{ix,iy} );
    }
    return sum;
  }

  double average_plaquette() const {
    double sum = 0.0;
    for(int i=0; i<lattice.n_faces; i++) {
      sum += std::cos(plaquette_angle(i));
    }
    return sum/lattice.n_faces;
  }

  // double average_plaqsq() const {
  //   double sum = 0.0;
  //   for(int i=0; i<lattice.n_faces; i++) {
  //     double val = plaquette_angle(i);
  //     while(val>M_PI) val -= 2.0*M_PI;
  //     while(val<-M_PI) val += 2.0*M_PI;
  //     const double factor = 1.0 / lattice.vps[i];
  //     sum += factor * std::pow( val, 2 );
  //   }
  //   return sum/lattice.n_faces;
  // }

  void dist_plaqsq( double& mean, double& square,
		    const bool is_compact=false ) const {
    mean = 0.0;
    square = 0.0;
    for(int i=0; i<lattice.n_faces; i++) {
      double val = plaquette_angle(i);
      if(is_compact){
	while(val>M_PI) val -= 2.0*M_PI;
	while(val<-M_PI) val += 2.0*M_PI;
      }
      const double factor = 1.0 / lattice.vps[i];

      const double v = factor * std::pow( val, 2 );
      mean += v;
      square += v*v;
    }
    mean /= lattice.n_faces;
    square /= lattice.n_faces;
  }

};




struct U1Wilson {
  const double gR;
  const bool is_compact;

  U1Wilson() = delete;
  U1Wilson(const U1Wilson&) = delete;

  U1Wilson(const double gR_, const bool is_compact_)
    : gR(gR_)
    , is_compact(is_compact_)
  {}

  U1Wilson & operator=(const U1Wilson&) = delete;

  double operator()( const U1onS2& U ) const {
    double res = 0.0;
    for(int i=0; i<U.lattice.n_faces; i++) {
      if(is_compact) res += - 1.0/U.lattice.vps[i] *  std::cos( U.plaquette_angle(i) );
      else res += 0.5/U.lattice.vps[i] * std::pow( U.plaquette_angle(i), 2 );
    }
    res /= gR*gR;
    return res;
  }

  // !! need debug ?
  // double local( const int i_link, const U1onS2& U ) const {
  //   double res = 0.0;
  //   for(int i=0; i<U.lattice.links[i_link].n_faces; i++) {
  //     const int i_face = U.lattice.links[i_link].faces[i];
  //     res += -beta * std::cos(U.plaquette_angle(i_face));
  //   }
  //   return res;
  // }

};


