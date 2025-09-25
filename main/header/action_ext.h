#pragma once


template<class Base>
struct U1WilsonExt {
  using Link = std::array<int,2>; // <int,int>;
  using Face = std::vector<int>;
  using Action=U1WilsonExt;

  Base& base;
  const double gsq, at;
  std::vector<double> beta_t;
  std::vector<double> beta_s;

  U1WilsonExt() = delete;
  U1WilsonExt(const U1WilsonExt&) = delete;

  U1WilsonExt(const double gsq_,
              const double at_,
              Base& base_)
    : base(base_)
    , gsq(gsq_)
    , at(at_)
    , beta_t(base.n_links)
    , beta_s(base.n_faces)
  {
    set_beta();
  }

  Action & operator=(const Action&) = delete;

  template <typename Gauge>
  double operator()( const Gauge& U ) const {
    // const auto& base = U.lattice;

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

        // tmp[s] += 0.5*beta*at/base.vols[i]  * std::pow( U.plaquette_angle(s, face), 2 );
        tmp[s] += 0.5*beta_s[i] * std::pow( U.plaquette_angle(s, face), 2 );
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

        // tmp[s] += 0.5*beta/at * base.link_volume[il]/std::pow(base.ell[il],2) *std::pow( U.plaquette_angle(s, link), 2 );
        tmp[s] += 0.5*beta_t[il] * std::pow( U.plaquette_angle(s, link), 2 );
      }
    }
    //}

    double res = 0.0;
    for(int s=0; s<U.Nt; s++){
      res += tmp[s];
    }

    return res;
  }




//   void coo_structure( std::vector<Idx>& is,
//                       std::vector<Idx>& js ) const {
//     const Idx len = 4*lattice.counter_accum.back()*Nt + 8*lattice.n_sites*Nt + 4*lattice.n_sites*Nt;
//     is.resize(len);
//     js.resize(len);

// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL)
// #endif
//     for(int s=0; s<Nt; s++){
//       for(Idx ix=0; ix<lattice.n_sites; ix++){
//         Idx counter = 4*lattice.counter_accum.back()*s + 4*lattice.counter_accum[ix];
//         for(const Idx iy : lattice.nns[ix]){
//           is[counter] = Nx*s+NS*ix; js[counter] = Nx*s+NS*iy; counter++;
//           is[counter] = Nx*s+NS*ix; js[counter] = Nx*s+NS*iy+1; counter++;

//           is[counter] = Nx*s+NS*ix+1; js[counter] = Nx*s+NS*iy; counter++;
//           is[counter] = Nx*s+NS*ix+1; js[counter] = Nx*s+NS*iy+1; counter++;
//         }
//       }
//     }


// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL)
// #endif
//     for(int s=0; s<Nt; s++){
//       for(Idx ix=0; ix<lattice.n_sites; ix++){
//         Idx counter = 4*lattice.counter_accum.back()*Nt + 8*(lattice.n_sites*s + ix);
//         is[counter] = ( Nx*(s+1)+NS*ix )%N; js[counter] = Nx*s+NS*ix; counter++;
//         is[counter] = ( Nx*(s+1)+NS*ix )%N; js[counter] = Nx*s+NS*ix+1; counter++;
//         is[counter] = ( Nx*(s-1)+NS*ix + N )%N; js[counter] = Nx*s+NS*ix; counter++;
//         is[counter] = ( Nx*(s-1)+NS*ix + N )%N; js[counter] = Nx*s+NS*ix+1; counter++;

//         is[counter] = ( Nx*(s+1)+NS*ix+1 )%N; js[counter] = Nx*s+NS*ix; counter++;
//         is[counter] = ( Nx*(s+1)+NS*ix+1 )%N; js[counter] = Nx*s+NS*ix+1; counter++;
//         is[counter] = ( Nx*(s-1)+NS*ix+1 + N )%N; js[counter] = Nx*s+NS*ix; counter++;
//         is[counter] = ( Nx*(s-1)+NS*ix+1 + N )%N; js[counter] = Nx*s+NS*ix+1; counter++;
//       }
//     }

// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL)
// #endif
//     for(int s=0; s<Nt; s++){
//       for(Idx ix=0; ix<lattice.n_sites; ix++){
//         Idx counter = 4*lattice.counter_accum.back()*Nt + 8*lattice.n_sites*Nt + 4*(lattice.n_sites*s + ix);
//         is[counter] = Nx*s+NS*ix; js[counter] = Nx*s+NS*ix; counter++;
//         is[counter] = Nx*s+NS*ix; js[counter] = Nx*s+NS*ix+1; counter++;
//         is[counter] = Nx*s+NS*ix+1; js[counter] = Nx*s+NS*ix; counter++;
//         is[counter] = Nx*s+NS*ix+1; js[counter] = Nx*s+NS*ix+1; counter++;
//       }
//     }
//   }





  template<typename Gauge>
  void coo_format( std::vector<Idx>& is,
                   std::vector<Idx>& js,
                   std::vector<Complex>& vs,
		   const Gauge& U ) const {
    const Idx Nx = Comp::Nx;
    const int Nt = Comp::Nt;

    // const Idx len = 4*lattice.counter_accum.back()*Nt + 8*lattice.n_sites*Nt + 4*lattice.n_sites*Nt;
    // is.resize(len);
    // js.resize(len);
    is.clear();
    js.clear();
    vs.clear();

    // std::cout << "debug." << std::endl;


    // for(Idx i=0; i<v.size(); i++) v[i] = 0.0;

    // spatial
// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL_DUPDATE)  schedule(static)
// #endif
    Idx counter=0;
    for(int s=0; s<Nt; s++){
      for(Idx ib=0; ib<base.faces.size(); ib++) {
        const Face& face = base.faces[ib];

        // tmp[s] += 0.5*beta_s[i] * std::pow( U.plaquette_angle(s, face), 2 );
        for(Idx i=0; i<face.size(); i++) {
          for(Idx j=0; j<face.size(); j++) {
            // std::cout << "debug. pt0" << std::endl;
            const Idx ix = face[i];
            const Idx iy = face[(i+1)%face.size()];
            const Idx jx = face[j];
            const Idx jy = face[(j+1)%face.size()];
            const Link li{ix,iy};
            const Link lj{jx,jy};

            // std::cout << "debug. pt1" << std::endl;

            // (ix,iy), (jx,jy) = 0.5*beta_s[ib] * U.sp( s, BaseLink{ix,iy} ) * U.sp( s, BaseLink{jx,jy} );
            // idx_sp(s, ix,iy), idx_sp(s, jx,jy) = sign(ix,iy) * sign(jx,jy) * 0.5*beta_s[ib];
            is.push_back( U.idx_sp(s, li) );
            // std::cout << "debug. pt2" << std::endl;
            js.push_back( U.idx_sp(s, lj) );
            // std::cout << "debug. pt3" << std::endl;
            vs.push_back( base.map2sign.at(li) * base.map2sign.at(lj) * beta_s[ib] );
            counter++;
            // std::cout << "debug. counter = " << counter << std::endl;
          }}
      } // faces
    } // s


// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL_DUPDATE) schedule(static)
// #endif
    for(int s=0; s<Nt; s++){
      for(const Link& link : base.links) {
        const Idx il = base.map2il.at(link);

        // tmp[s] += 0.5*beta_t[il] * std::pow( U.plaquette_angle(s, link), 2 );
        // sum += U.sp( s, link );
        {
          // sum += 0.5*beta_t[il] * U.sp( s, link ) * U.sp( s, link );
          // idx_sp(s, link),idx_sp(s, link) = 0.5*beta_t[il];
          is.push_back( U.idx_sp(s, il) );
          js.push_back( U.idx_sp(s, il) );
          vs.push_back( beta_t[il] );
          counter++;

          // sum += 0.5*beta_t[il] * U.sp( s, link ) * U.tp( s, link[1] );
          // idx_sp(s, link),idx_tp(s, link[1]) = 0.5*beta_t[il];
          is.push_back( U.idx_sp(s, il) );
          js.push_back( U.idx_tp(s, link[1]) );
          vs.push_back( beta_t[il] );
          counter++;

          // sum -= 0.5*beta_t[il] * U.sp( s, link ) * U.sp( s+1, link );
          // idx_sp(s, link),idx_sp(s+1, link) = -0.5*beta_t[il];
          is.push_back( U.idx_sp(s, il) );
          js.push_back( U.idx_sp(s+1, link) );
          vs.push_back( -beta_t[il] );
          counter++;

          // sum -= 0.5*beta_t[il] * U.sp( s, link ) * U.tp( s, link[0] );
          // idx_sp(s, link),idx_tp(s, link[0]) = -0.5*beta_t[il];
          is.push_back( U.idx_sp(s, il) );
          js.push_back( U.idx_tp(s, link[0]) );
          vs.push_back( -beta_t[il] );
          counter++;
        }

        // sum += U.tp( s, link[1] );
        {
          // sum += 0.5*beta_t[il] * U.tp( s, link[1] ) * U.sp( s, link );
          // idx_tp(s, link[1]),idx_sp(s, link) = 0.5*beta_t[il];
          is.push_back( U.idx_tp(s, link[1]) );
          js.push_back( U.idx_sp(s, il) );
          vs.push_back( beta_t[il] );
          counter++;

          // sum += 0.5*beta_t[il] * U.tp( s, link[1] ) * U.tp( s, link[1] );
          // idx_tp(s, link[1]),idx_sr(s, link[1]) = 0.5*beta_t[il];
          is.push_back( U.idx_tp(s, link[1]) );
          js.push_back( U.idx_tp(s, link[1]) );
          vs.push_back( beta_t[il] );
          counter++;

          // sum -= 0.5*beta_t[il] * U.tp( s, link[1] ) * U.sp( s+1, link );
          // idx_tp(s, link[1]),idx_sp(s+1, link) = -0.5*beta_t[il];
          is.push_back( U.idx_tp(s, link[1]) );
          js.push_back( U.idx_sp(s+1, link) );
          vs.push_back( -beta_t[il] );
          counter++;

          // sum -= 0.5*beta_t[il] * U.tp( s, link[1] ) * U.tp( s, link[0] );
          // idx_tp(s, link[1]),idx_sp(s, link[0]) = -0.5*beta_t[il];
          is.push_back( U.idx_tp(s, link[1]) );
          js.push_back( U.idx_tp(s, link[0]) );
          vs.push_back( -beta_t[il] );
          counter++;
        }

        // sum -= U.sp( s+1, link );
        {
          // sum -= 0.5*beta_t[il] * U.sp( s+1, link ) * U.sp( s, link );
          // idx_sp(s+1, link),idx_sp(s, link) = -0.5*beta_t[il];
          is.push_back( U.idx_sp(s+1, link) );
          js.push_back( U.idx_sp(s, il) );
          vs.push_back( -beta_t[il] );
          counter++;

          // sum -= 0.5*beta_t[il] * U.sp( s+1, link ) * U.tp( s, link[1] );
          // idx_sp(s+1, link),idx_sr(s, link[1]) = -0.5*beta_t[il];
          is.push_back( U.idx_sp(s+1, link) );
          js.push_back( U.idx_tp(s, link[1]) );
          vs.push_back( -beta_t[il] );
          counter++;

          // sum += 0.5*beta_t[il] * U.sp( s+1, link ) * U.sp( s+1, link );
          // idx_sp(s+1, link),idx_sp(s+1, link) = 0.5*beta_t[il];
          is.push_back( U.idx_sp(s+1, link) );
          js.push_back( U.idx_sp(s+1, link) );
          vs.push_back( beta_t[il] );
          counter++;

          // sum += 0.5*beta_t[il] * U.sp( s+1, link ) * U.tp( s, link[0] );
          // idx_sp(s+1, link),idx_sp(s, link[0]) = 0.5*beta_t[il];
          is.push_back( U.idx_sp(s+1, link) );
          js.push_back( U.idx_tp(s, link[0]) );
          vs.push_back( beta_t[il] );
          counter++;
        }

        // sum -= U.tp( s, link[0] );
        {
          // sum -= 0.5*beta_t[il] * U.tp( s, link[0] ) * U.sp( s, link );
          // idx_tp(s, link[0]),idx_sp(s, link) = -0.5*beta_t[il];
          is.push_back( U.idx_tp(s, link[0]) );
          js.push_back( U.idx_sp(s, il) );
          vs.push_back( -beta_t[il] );
          counter++;

          // sum -= 0.5*beta_t[il] * U.tp( s, link[0] ) * U.tp( s, link[1] );
          // idx_tp(s, link[0]),idx_sr(s, link[1]) = -0.5*beta_t[il];
          is.push_back( U.idx_tp(s, link[0]) );
          js.push_back( U.idx_tp(s, link[1]) );
          vs.push_back( -beta_t[il] );
          counter++;

          // sum += 0.5*beta_t[il] * U.tp( s, link[0] ) * U.sp( s+1, link );
          // idx_tp(s, link[0]),idx_sp(s+1, link) = 0.5*beta_t[il];
          is.push_back( U.idx_tp(s, link[0]) );
          js.push_back( U.idx_sp(s+1, link) );
          vs.push_back( beta_t[il] );
          counter++;

          // sum += 0.5*beta_t[il] * U.tp( s, link[0] ) * U.tp( s, link[0] );
          // idx_tp(s, link[0]),idx_sp(s, link[0]) = 0.5*beta_t[il];
          is.push_back( U.idx_tp(s, link[0]) );
          js.push_back( U.idx_tp(s, link[0]) );
          vs.push_back( beta_t[il] );
          counter++;
        }
      }
    }


// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL_DUPDATE) schedule(static)
// #endif
//     for(int s=0; s<Nt; s++){
//       for(Idx ix=0; ix<lattice.n_sites; ix++){
//         double coeff = 0.0;
//         for(const Idx iy : lattice.nns[ix]){
//           const Idx il = lattice.map2il.at(BaseLink{ix,iy});
//           coeff += 0.5 * r*bd.kappa[il];
//         }
//         coeff += r*kappa_t[ix];
//         coeff += M5;

//         Idx counter = 4*lattice.counter_accum.back()*Nt + 8*lattice.n_sites*Nt + 4*(lattice.n_sites*s + ix);
//         const MS tmp2 = coeff * sigma[0];

//         v[counter] = tmp2(0,0); counter++;
//         v[counter] = tmp2(0,1); counter++;

//         v[counter] = tmp2(1,0); counter++;
//         v[counter] = tmp2(1,1); counter++;
//       }
//     }
  }























  template <typename Force, typename Gauge>
  void get_force( Force& pi, const Gauge& U ) const {
    // const auto& base = U.lattice;

#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL_GAUGE)
#endif
    for(Idx i=0; i<pi.spatial.size(); i++) for(Idx j=0; j<pi.spatial[i].size(); j++) pi.spatial[i][j] = 0.0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL_GAUGE)
#endif
    for(Idx i=0; i<pi.temporal.size(); i++) for(Idx j=0; j<pi.temporal[i].size(); j++) pi.temporal[i][j] = 0.0;

    // const auto& base = U.lattice;

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
        grad = beta_s[i_face] * U.plaquette_angle(s, face);

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
        grad  = beta_t[il] * U.plaquette_angle(s, link);
        grad2 = beta_t[il] * U.plaquette_angle(s-1, link);

        pi.sp( s, base.map2il.at(link) ) += grad;
        pi.tp( s, link[1] ) += grad;
        pi.sp( s, base.map2il.at(link) ) -= grad2;
        pi.tp( s, link[0] ) -= grad;
      }
    }
  }

  void set_beta(){
    if(std::abs(gsq)>1.0e-15){
      for(Idx i=0; i<base.faces.size(); i++) {
        beta_s[i] = at/base.vols[i] / gsq;
      }
      for(Idx il=0; il<base.n_links; il++) {
        beta_t[il] = 2.0/at * base.link_volume[il]/std::pow(base.ell[il],2) / gsq;
      }
    }
  }

};


