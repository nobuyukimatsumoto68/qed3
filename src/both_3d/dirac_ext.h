#pragma once

#include "dirac_simp.h"

// template<typename Base>
class DiracExt : private DiracS2Simp{
public:
  using Base=S2Simp;
  // const double M5_t;

  Base& lattice;

  const Idx Nx; // = Comp::N;
  const int Nt; // = Comp::Nt;
  const Idx N;

  std::vector<double> kappa_t;

  DiracExt(Base& lattice_,
           const double m_=0.0,
           const double r_=1.0,
           const double M5_=0.0
           // const double M5_t_=0.0
           )
    : DiracS2Simp(lattice_, m_, r_, M5_)
    , lattice(lattice_)
    , Nx(Comp::Nx)
    , Nt(Comp::Nt)
    , N(Comp::N)
      // , M5_t(M5_t_)
    , kappa_t(lattice.n_sites)
  {
    // set_sigma();
    // set_face_signs();
    // set_ell_and_link_volumes();
    // set_kappa();
    set_kappa_t();
  }


  void coo_structure( std::vector<Idx>& is,
                      std::vector<Idx>& js ) const {
    // const Idx N = Comp::N;
    // const int Nt = Comp::Nt;

    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<lattice.n_sites; ix++){
        // Idx counter = lattice.counter_accum[ix];
        for(const Idx iy : lattice.nns[ix]){
          is.push_back(Nx*s+NS*ix); js.push_back(Nx*s+NS*iy);
          is.push_back(Nx*s+NS*ix); js.push_back(Nx*s+NS*iy+1);

          is.push_back(Nx*s+NS*ix); js.push_back(Nx*s+NS*ix);
          is.push_back(Nx*s+NS*ix); js.push_back(Nx*s+NS*ix+1);

          is.push_back(Nx*s+NS*ix+1); js.push_back(Nx*s+NS*iy);
          is.push_back(Nx*s+NS*ix+1); js.push_back(Nx*s+NS*iy+1);

          is.push_back(Nx*s+NS*ix+1); js.push_back(Nx*s+NS*ix);
          is.push_back(Nx*s+NS*ix+1); js.push_back(Nx*s+NS*ix+1);
        }
      }
    }

    if(Nt==1) return;
    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<lattice.n_sites; ix++){
        is.push_back( ( Nx*(s+1)+NS*ix )%N ); js.push_back( Nx*s+NS*ix );
        is.push_back( ( Nx*(s+1)+NS*ix )%N ); js.push_back( Nx*s+NS*ix+1 );
        is.push_back( ( Nx*(s-1)+NS*ix + N )%N ); js.push_back( Nx*s+NS*ix );
        is.push_back( ( Nx*(s-1)+NS*ix + N )%N ); js.push_back( Nx*s+NS*ix+1 );
        is.push_back( Nx*s+NS*ix ); js.push_back( Nx*s+NS*ix );
        is.push_back( Nx*s+NS*ix ); js.push_back( Nx*s+NS*ix+1 );

        is.push_back( ( Nx*(s+1)+NS*ix+1 )%N ); js.push_back( Nx*s+NS*ix );
        is.push_back( ( Nx*(s+1)+NS*ix+1 )%N ); js.push_back( Nx*s+NS*ix+1 );
        is.push_back( ( Nx*(s-1)+NS*ix+1 + N )%N ); js.push_back( Nx*s+NS*ix );
        is.push_back( ( Nx*(s-1)+NS*ix+1 + N )%N ); js.push_back( Nx*s+NS*ix+1 );
        is.push_back( Nx*s+NS*ix+1 ); js.push_back( Nx*s+NS*ix );
        is.push_back( Nx*s+NS*ix+1 ); js.push_back( Nx*s+NS*ix+1 );
      }
    }
  }



  template<typename Gauge>
  void coo_format( std::vector<Complex>& v,
		   const Gauge& u ) const {
    const Idx Nx = Comp::Nx;
    const int Nt = Comp::Nt;

    // elem.clear();
    // elem.resize( N*Nt );
    for(Idx i=0; i<Nx*Nt; i++) v[i] = 0.0;

// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL) collapse(2)
// #endif
    // Idx counter=0;
    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<lattice.n_sites; ix++){
        Idx counter = lattice.counter_accum.back()*s + lattice.counter_accum[ix];
        for(const Idx iy : lattice.nns[ix]){
          const Idx il = lattice.map2il.at(Link{ix,iy});

          const MS tmp = 0.5 * kappa[il] * ( -r*sigma[0] + gamma(ix, iy) ) * std::exp( I*u.sp(s,Link{ix,iy})) * Omega(ix, iy);
          const MS tmp2 = 0.5 * r*kappa[il] * sigma[0] + M5/lattice.nns[ix].size() * sigma[0];

          // res[NS*ix] += -tmp(0,0)*v[NS*iy] - tmp(0,1)*v[NS*iy+1];
          v[counter] = tmp(0,0); counter++;
          v[counter] = tmp(0,1); counter++;

          // res[NS*ix] += tmp(0,0)*v[NS*ix] + tmp(0,1)*v[NS*ix+1];
          v[counter] = tmp2(0,0); counter++;
          v[counter] = tmp2(0,1); counter++;

          // res[NS*ix+1] += -tmp(1,0)*v[NS*iy] - tmp(1,1)*v[NS*iy+1];
          v[counter] = tmp(1,0); counter++;
          v[counter] = tmp(1,1); counter++;

          // res[NS*ix+1] += tmp(1,0)*v[NS*ix] + tmp(1,1)*v[NS*ix+1];
          v[counter] = tmp2(1,0); counter++;
          v[counter] = tmp2(1,1); counter++;
        }
      }
    }

    if(Nt==1) return;
// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL) collapse(2)
// #endif
    Idx counter = lattice.counter_accum.back()*Nt;
    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<lattice.n_sites; ix++){
        int signP = 1; if(s==Nt-1) signP = -1;
        int signM = 1; if(s==0) signM = -1;

        // Idx counter = 12 * (lattice.n_sites*s + ix);

        // const MS tmpP = 0.5 * signP * ( -sigma[0] + sigma[3] );
        // const MS tmpM = 0.5 * signM * ( -sigma[0] - sigma[3] );
        // const MS tmpD = sigma[0];

        const MS tmpP = 0.5 * signP * kappa_t[ix] * ( -r*sigma[0] + sigma[3] ) * std::exp( I*u.tp(s,ix) );
        const MS tmpM = 0.5 * signM * kappa_t[ix] * ( -r*sigma[0] - sigma[3] ) * std::exp( I*u.tp(s-1,ix) );
        const MS tmpD = r*kappa_t[ix] * sigma[0];

        v[counter] = tmpP(0,0); counter++;
        v[counter] = tmpP(0,1); counter++;
        v[counter] = tmpM(0,0); counter++;
        v[counter] = tmpM(0,1); counter++;
        v[counter] = tmpD(0,0); counter++;
        v[counter] = tmpD(0,1); counter++;

        v[counter] = tmpP(1,0); counter++;
        v[counter] = tmpP(1,1); counter++;
        v[counter] = tmpM(1,0); counter++;
        v[counter] = tmpM(1,1); counter++;
        v[counter] = tmpD(1,0); counter++;
        v[counter] = tmpD(1,1); counter++;
      }
    }
  }

  // void set_kappa() {
  //   for(Idx il=0; il<lattice.n_links; il++) {
  //     kappa[il] = lattice.link_volume[il]/lattice.ell[il]/(lattice.alat/std::sqrt(3.0));
  //   }
  // }

  void set_kappa_t() {
    // const double factor = std::pow(base.mean_ell/base.ell[il], 2) * base.link_volume[il]/base.mean_link_volume;
    for(Idx ix=0; ix<lattice.n_sites; ix++) { // @@@
      kappa_t[ix] = lattice.dual_areas[ix]/lattice.mean_dual_area;
    }
  }



};
