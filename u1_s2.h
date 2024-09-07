#pragma once

#include "s2.h"

struct CompactU1onS2 {
  QfeLatticeS2& lattice;

  std::vector<int> face_signs;
  std::vector<double> field;

  CompactU1onS2()=delete;

  CompactU1onS2(QfeLatticeS2& lattice_, const double beta_)
    : lattice(lattice_)
    , field(lattice.n_links, 0.0)
    , face_signs(lattice.n_faces)
  {
    set_face_signs();
  }

  CompactU1onS2( const CompactU1onS2& other )
    : lattice(other.lattice)
    , field(other.field)
    , face_signs(other.face_signs)
  {}

  CompactU1onS2 & operator=(const CompactU1onS2&) = delete;

  double operator[](const int i) const { return field[i]; }
  double& operator[](const int i) { return field[i]; }

  int face_sign(const int i_face) const {
    const QfeFace& face = lattice.faces[i_face];
    const Vec3 r0 = lattice.r[face.sites[0]];
    const Vec3 r1 = lattice.r[face.sites[1]];
    const Vec3 r2 = lattice.r[face.sites[2]];

    const Vec3 cross = (r1-r0).cross(r2-r0);
    const Vec3 sum = r0+r1+r2;

    const double inner = cross.dot(sum);

    int res = 1;
    if(inner<0) res = -1;
    return res;
  }

  void set_face_signs() { for(int i=0; i<lattice.n_faces; i++) face_signs[i] = face_sign(i); }


  double plaquette_angle(const int& i_face) const {
    const QfeFace& face = lattice.faces[i_face];
    double sum = 0.0;
    // for(const int i_edge : face.edges) sum += field[i_edge];
    for(int i=0; i<face.n_edges; i++) sum += field[face.edges[i]];
    return face_signs[i_face] * sum;
  }


  double average_plaquette() const {
    double sum = 0.0;
    for(int i=0; i<lattice.n_faces; i++) {
      sum += std::cos(plaquette_angle(i));
    }
    return sum/lattice.n_faces;
  }


};


struct U1Wilson {
  double beta;

  U1Wilson() = delete;
  U1Wilson(const U1Wilson&) = delete;

  U1Wilson(const double beta_)
    : beta(beta_)
  {}

  U1Wilson & operator=(const U1Wilson&) = delete;

  double operator()( const CompactU1onS2& U ) const {
    double res = 0.0;
    for(int i=0; i<U.lattice.n_faces; i++) {
      res += -beta * std::cos(U.plaquette_angle(i));
    }
    return res;
  }

  // !! need debug ?
  // double local( const int i_link, const CompactU1onS2& U ) const {
  //   double res = 0.0;
  //   for(int i=0; i<U.lattice.links[i_link].n_faces; i++) {
  //     const int i_face = U.lattice.links[i_link].faces[i];
  //     res += -beta * std::cos(U.plaquette_angle(i_face));
  //   }
  //   return res;
  // }

};


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


