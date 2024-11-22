#pragma once

#include "s2.h"


struct U1onS2 {
  QfeLatticeS2& lattice;
  const LinkInfo info_ell;
  const FaceInfo info_p;

  using Link = std::array<int,2>; // <int,int>;

  std::vector<double> field;

  U1onS2()=delete;

  U1onS2(QfeLatticeS2& lattice_)
    : lattice(lattice_)
    , info_ell(lattice)
    , info_p(lattice)
    , field(lattice.n_links, 0.0)
  {}

  U1onS2( const U1onS2& other )
    : lattice(other.lattice)
    , info_ell(other.info_ell)
    , info_p(other.info_p)
    , field(other.field)
  {}

  U1onS2 & operator=(const U1onS2&) = delete;

  double operator[](const int i) const { return field[i]; }
  double& operator[](const int i) { return field[i]; }

  double operator()(const Link& ell) const {
    const int il = info_ell.map2il.at(ell);
    const int sign = info_ell.map2sign.at(ell);
    return sign * field[il] ;
  }

  double plaquette_angle(const int& i_face) const {
    const QfeFace& face = lattice.faces[i_face];
    double sum = 0.0;
    for(int i=0; i<face.n_edges; i++) {
      const int ix = face.sites[i];
      const int iy = face.sites[(i+1)%face.n_edges];
      sum += (*this)( Link{ix,iy} );
    }
    return info_p.signs[i_face] * sum;
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
  //     const double factor = 1.0 / info_p.vps[i];
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
      const double factor = 1.0 / info_p.vps[i];

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
      if(is_compact) res += - 1.0/U.info_p.vps[i] *  std::cos( U.plaquette_angle(i) );
      else res += 0.5/U.info_p.vps[i] * std::pow( U.plaquette_angle(i), 2 );
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

