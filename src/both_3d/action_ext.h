#pragma once


struct U1WilsonExt {
  using Link = std::array<int,2>; // <int,int>;
  using Face = std::vector<int>;
  using Action=U1WilsonExt;

  const double beta, at;

  U1WilsonExt() = delete;
  U1WilsonExt(const U1WilsonExt&) = delete;

  U1WilsonExt(const double beta_,
              const double at_)
    : beta(beta_)
    , at(at_)
  {}

  Action & operator=(const Action&) = delete;

  template <typename Gauge>
  double operator()( const Gauge& U ) const {
    const auto& base = U.lattice;

    std::vector<double> tmp(U.Nt, 0.0);
    // spatial
#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL_GAUGE) // collapse(2)
#endif
    for(int s=0; s<U.Nt; s++){
      //for(const Face& face : base.faces) {
      for(Idx i=0; i<base.faces.size(); i++) {
        const Face& face = base.faces[i];
        // Idx i=std::distance( base.faces[0],  );
        // int i=0;
        // const double factor = base.mean_vol/base.vols[i];
        // const double factor = base.mean_vol/base.vols[i];
        // if constexpr(U.is_compact) tmp[s] += - beta_s*factor * ( std::cos( U.plaquette_angle(s, face) ) - 1.0);
        tmp[s] += 0.5*beta*at/base.vols[i]  * std::pow( U.plaquette_angle(s, face), 2 );
        // i++;
      }
    }

    // // temporal
    // if(U.Nt!=1){
#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL_GAUGE)
#endif
    for(int s=0; s<U.Nt ; s++){
      for(const Link& link : base.links) {
        const Idx il = base.map2il.at(link);
        // const double factor = std::pow(base.mean_ell/base.ell[il], 2) * base.link_volume[il]/base.mean_link_volume;
        // if constexpr(U.is_compact) tmp[s] += - beta_t*factor * (std::cos( U.plaquette_angle(s, link) ) - 1.0);
        // else tmp[s] += 0.5*beta_t* factor *std::pow( U.plaquette_angle(s, link), 2 );
        tmp[s] += 0.5*beta/at * base.link_volume[il]/std::pow(base.ell[il],2) *std::pow( U.plaquette_angle(s, link), 2 );
      }
    }
    //}

    double res = 0.0;
    for(int s=0; s<U.Nt; s++){
      res += tmp[s];
    }

    return res;
  }


  template <typename Force, typename Gauge>
  void get_force( Force& pi, const Gauge& U ) const {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL_GAUGE)
#endif
    for(Idx i=0; i<pi.spatial.size(); i++) for(Idx j=0; j<pi.spatial[i].size(); j++) pi.spatial[i][j] = 0.0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL_GAUGE)
#endif
    for(Idx i=0; i<pi.temporal.size(); i++) for(Idx j=0; j<pi.temporal[i].size(); j++) pi.temporal[i][j] = 0.0;

    const auto& base = U.lattice;

    // spatial
#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL_GAUGE) // collapse(2)
#endif
    for(int s=0; s<U.Nt; s++){
      for(int i_face=0; i_face<base.n_faces; i_face++){
        const Face& face = base.faces[i_face];
        double grad;
        // const double factor = base.mean_vol/base.vols[i_face];
        // if constexpr(U.is_compact) grad = beta_s*factor * std::sin( U.plaquette_angle(s, face) );
        // else grad = beta_s*factor * U.plaquette_angle(s, face);
        grad = beta*at/base.vols[i_face] * U.plaquette_angle(s, face);

        for(int i=0; i<face.size(); i++) {
          const Idx ix = face[i];
          const Idx iy = face[(i+1)%face.size()];
          const Link ell{ix, iy};
          const Idx il = base.map2il.at(ell);
          pi.sp( s, il ) += grad * base.map2sign.at(ell);
        }
      }
    }

    // temporal
#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL_GAUGE)
#endif
    for(int s=0; s<U.Nt; s++){
      for(const Link& link : base.links) {
        const Idx il = base.map2il.at(link);
        double grad, grad2;
        // const double factor = std::pow(base.mean_ell/base.ell[il], 2) * base.link_volume[il]/base.mean_link_volume;
        // if constexpr(U.is_compact) {
        //   grad = beta_t *factor* std::sin( U.plaquette_angle(s, link) );
        //   grad2 = beta_t *factor* std::sin( U.plaquette_angle(s-1, link) );
        // }
        // else {
        //   grad = beta_t *factor* U.plaquette_angle(s, link);
        //   grad2 = beta_t *factor* U.plaquette_angle(s-1, link);
        // }
        grad  = beta/at * base.link_volume[il]/std::pow(base.ell[il],2) * U.plaquette_angle(s, link);
        grad2 = beta/at * base.link_volume[il]/std::pow(base.ell[il],2) * U.plaquette_angle(s-1, link);


        pi.sp( s, base.map2il.at(link) ) += grad;
        pi.tp( s, link[1] ) += grad;
        pi.sp( s, base.map2il.at(link) ) -= grad2;
        pi.tp( s, link[0] ) -= grad;
      }
    }
  }

};


