#pragma once


struct U1WilsonExt {
  using Link = std::array<int,2>; // <int,int>;
  using Face = std::vector<int>;
  using Action=U1WilsonExt;

  const double beta_t, beta_s;

  U1WilsonExt() = delete;
  U1WilsonExt(const U1WilsonExt&) = delete;

  U1WilsonExt(const double beta_t_,
              const double beta_s_)
    : beta_t(beta_t_)
    , beta_s(beta_s_)
  {}

  Action & operator=(const Action&) = delete;

  template <typename Gauge>
  double operator()( const Gauge& U ) const {
    double res = 0.0;
    // spatial
    for(int s=0; s<U.Nt; s++){
      int i=0;
      for(const Face& face : U.lattice.faces) {
        if constexpr(U.is_compact) res += - 1.0*beta_s*U.lattice.mean_vol/U.lattice.vols[i] *  std::cos( U.plaquette_angle(s, face) );
        else res += 0.5*beta_s*U.lattice.mean_vol/U.lattice.vols[i] * std::pow( U.plaquette_angle(s, face), 2 );
        i++;
      }
    }
    // temporal
    for(int s=0; s<U.Nt; s++){
      for(const Link& link : U.lattice.links) {
        if constexpr(U.is_compact) res += - beta_t*std::cos( U.plaquette_angle(s, link) );
        else res += 0.5*beta_t*std::pow( U.plaquette_angle(s, link), 2 );
      }
    }
    return res;
  }


  template <typename Force, typename Gauge>
  void get_force( Force& pi, const Gauge& U ) const {
    for(Idx i=0; i<pi.spatial.size(); i++) for(Idx j=0; j<pi.spatial[i].size(); j++) pi.spatial[i][j] = 0.0;
    for(Idx i=0; i<pi.temporal.size(); i++) for(Idx j=0; j<pi.temporal[i].size(); j++) pi.temporal[i][j] = 0.0;

    // spatial
    for(int s=0; s<U.Nt; s++){
      for(int i_face=0; i_face<U.lattice.n_faces; i_face++){
        const Face& face = U.lattice.faces[i_face];
        const double grad = beta_s*U.lattice.mean_vol/U.lattice.vols[i_face] * U.plaquette_angle(s, face);

        for(int i=0; i<face.size(); i++) {
          const Idx ix = face[i];
          const Idx iy = face[(i+1)%face.size()];
          const Link ell{ix, iy};
          const Idx il = U.lattice.map2il.at(ell);
          pi.sp( s, il ) += grad * U.lattice.map2sign.at(ell);
        }
      }
    }

    // temporal
    for(int s=0; s<U.Nt; s++){
      for(const Link& link : U.lattice.links) {
        const double grad = beta_t * U.plaquette_angle(s, link);

        pi.sp( s, U.lattice.map2il.at(link) ) += grad;
        pi.tp( s, link[1] ) += grad;
        pi.sp( s+1, U.lattice.map2il.at(link) ) -= grad;
        pi.tp( s, link[0] ) -= grad;
      }
    }
  }

};


