#pragma once


struct U1Wilson {
  using Link = std::array<int,2>; // <int,int>;
  using Face = std::vector<int>;
  using Action=U1Wilson;

  const double beta;

  U1Wilson() = delete;
  U1Wilson(const U1Wilson&) = delete;

  U1Wilson(const double beta_)
    : beta(beta_)
  {}

  Action & operator=(const Action&) = delete;

  template <typename Gauge>
  double operator()( const Gauge& U ) const {
    double res = 0.0;
    for(int i=0; i<U.lattice.n_faces; i++) {
      if constexpr(U.is_compact) res += - 1.0/U.lattice.vps[i] *  std::cos( U.plaquette_angle(i) );
      else res += 0.5/U.lattice.vps[i] * std::pow( U.plaquette_angle(i), 2 );
    }
    res *= beta;
    return res;
  }

  template <typename Force, typename Gauge>
  void get_force( Force& pi, const Gauge& U ) const {
    for(auto& elem : pi) elem=0.0;

    for(int i_face=0; i_face<U.lattice.n_faces; i_face++){
      const Face& face = U.lattice.faces[i_face];
      const double grad = 1.0/U.lattice.vps[i_face] * U.plaquette_angle(face);

      for(int i=0; i<face.size(); i++) {
	const int ix = face[i];
	const int iy = face[(i+1)%face.size()];
	const Link ell{ix, iy};
	pi[ U.lattice.map2il.at(ell) ] += grad * U.lattice.map2sign.at(ell);
      }
    }

    pi *= beta;
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


