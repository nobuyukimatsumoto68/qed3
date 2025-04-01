#pragma once

#include "dirac_base.h"
// #include "dirac_simp.h"

template<class Base, class BaseDirac>
class DiracExt : public DiracBase {
public:
  Base& lattice;
  BaseDirac bd;

  const Idx Nx; // = Comp::N;
  const int Nt; // = Comp::Nt;
  const Idx N;

  std::vector<double> kappa_t;

  const double m;
  const double r;
  const double M5;

  const double c;

  DiracExt(Base& lattice_,
           // BaseDirac& bd_,
           const double m_=0.0,
           const double r_=1.0,
           const double M5_=0.0,
           const double c_=1.0
           // const double M5_t_=0.0
           )
    : lattice(lattice_)
    , bd(lattice_,m_,r_,M5_)
    , Nx(Comp::Nx)
    , Nt(Comp::Nt)
    , N(Comp::N)
    , m(m_)
    , r(r_)
    , M5(M5_)
    , kappa_t(lattice.n_sites)
    , c(c_)
  {
    set_kappa_t();
  }


  void coo_structure( std::vector<Idx>& is,
                      std::vector<Idx>& js ) const {
    // const Idx N = Comp::N;
    // const int Nt = Comp::Nt;

    const Idx len = 4*lattice.counter_accum.back()*Nt + 8*lattice.n_sites*Nt + 4*lattice.n_sites*Nt;
    // const Idx len = 4*lattice.counter_accum.back();
    is.resize(len);
    js.resize(len);

    // Idx counter=0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL)
#endif
    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<lattice.n_sites; ix++){
        Idx counter = 4*lattice.counter_accum.back()*s + 4*lattice.counter_accum[ix];
        // assert( counter==4*lattice.counter_accum.back()*s + 4*lattice.counter_accum[ix] );
        for(const Idx iy : lattice.nns[ix]){
          is[counter] = Nx*s+NS*ix; js[counter] = Nx*s+NS*iy; counter++;
          is[counter] = Nx*s+NS*ix; js[counter] = Nx*s+NS*iy+1; counter++;

          is[counter] = Nx*s+NS*ix+1; js[counter] = Nx*s+NS*iy; counter++;
          is[counter] = Nx*s+NS*ix+1; js[counter] = Nx*s+NS*iy+1; counter++;
        }
      }
    }

    // return ; // @@@ DEBUG

#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL)
#endif
    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<lattice.n_sites; ix++){
        Idx counter = 4*lattice.counter_accum.back()*Nt + 8*(lattice.n_sites*s + ix);
        // assert( counter==4*lattice.counter_accum.back()*Nt + 8*(lattice.n_sites*s + ix) );
        is[counter] = ( Nx*(s+1)+NS*ix )%N; js[counter] = Nx*s+NS*ix; counter++;
        is[counter] =  ( Nx*(s+1)+NS*ix )%N; js[counter] = Nx*s+NS*ix+1; counter++;
        is[counter] = ( Nx*(s-1)+NS*ix + N )%N; js[counter] = Nx*s+NS*ix; counter++;
        is[counter] =  ( Nx*(s-1)+NS*ix + N )%N; js[counter] = Nx*s+NS*ix+1; counter++;

        is[counter] = ( Nx*(s+1)+NS*ix+1 )%N; js[counter] = Nx*s+NS*ix; counter++;
        is[counter] = ( Nx*(s+1)+NS*ix+1 )%N; js[counter] = Nx*s+NS*ix+1; counter++;
        is[counter] = ( Nx*(s-1)+NS*ix+1 + N )%N; js[counter] = Nx*s+NS*ix; counter++;
        is[counter] = ( Nx*(s-1)+NS*ix+1 + N )%N; js[counter] = Nx*s+NS*ix+1; counter++;
      }
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL)
#endif
    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<lattice.n_sites; ix++){
        Idx counter = 4*lattice.counter_accum.back()*Nt + 8*lattice.n_sites*Nt + 4*(lattice.n_sites*s + ix);
        // assert( counter==4*lattice.counter_accum.back()*Nt + 8*lattice.n_sites*Nt + 4*(lattice.n_sites*s + ix) );
        is[counter] = Nx*s+NS*ix; js[counter] = Nx*s+NS*ix; counter++;
        is[counter] = Nx*s+NS*ix; js[counter] = Nx*s+NS*ix+1; counter++;
        is[counter] = Nx*s+NS*ix+1; js[counter] = Nx*s+NS*ix; counter++;
        is[counter] = Nx*s+NS*ix+1; js[counter] = Nx*s+NS*ix+1; counter++;
      }
    }
  }




  template<typename Gauge>
  void coo_format( std::vector<Complex>& v,
		   const Gauge& u ) const {
    const Idx Nx = Comp::Nx;
    const int Nt = Comp::Nt;

// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL)
// #endif
    for(Idx i=0; i<v.size(); i++) v[i] = 0.0;


// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL)
// #endif
    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<lattice.n_sites; ix++){
        Idx counter = 4*lattice.counter_accum.back()*s + 4*lattice.counter_accum[ix];
        for(const Idx iy : lattice.nns[ix]){
          const Idx il = lattice.map2il.at(Link{ix,iy});

          const MS tmp = 0.5 * bd.kappa[il] * ( -r * sigma[0] + bd.gamma(ix, iy) ) * std::exp( I*u.sp(s,Link{ix,iy})) * bd.Omega(ix, iy);

          // res[NS*ix] += -tmp(0,0)*v[NS*iy] - tmp(0,1)*v[NS*iy+1];
          v[counter] = tmp(0,0); counter++;
          v[counter] = tmp(0,1); counter++;

          // res[NS*ix+1] += -tmp(1,0)*v[NS*iy] - tmp(1,1)*v[NS*iy+1];
          v[counter] = tmp(1,0); counter++;
          v[counter] = tmp(1,1); counter++;
        }
      }
    }


// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL)
// #endif
    for(int s=0; s<Nt; s++){
      int signP = 1;
      int signM = 1;
      if(s==Nt-1) signP = -1;
      if(s==0) signM = -1;

      for(Idx ix=0; ix<lattice.n_sites; ix++){
        Idx counter = 4*lattice.counter_accum.back()*Nt + 8*(lattice.n_sites*s + ix);

        const MS tmpP = 0.5 * signP * kappa_t[ix] * ( -r*sigma[0] + sigma[3] ) * std::exp( I*u.tp(s,ix) );
        const MS tmpM = 0.5 * signM * kappa_t[ix] * ( -r*sigma[0] - sigma[3] ) * std::exp( I*u.tp(s-1,ix) );

        v[counter] = tmpP(0,0); counter++;
        v[counter] = tmpP(0,1); counter++;
        v[counter] = tmpM(0,0); counter++;
        v[counter] = tmpM(0,1); counter++;

        v[counter] = tmpP(1,0); counter++;
        v[counter] = tmpP(1,1); counter++;
        v[counter] = tmpM(1,0); counter++;
        v[counter] = tmpM(1,1); counter++;
      }
    }


// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL)
// #endif
    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<lattice.n_sites; ix++){
        double coeff = 0.0;
        for(const Idx iy : lattice.nns[ix]){
          const Idx il = lattice.map2il.at(Link{ix,iy});
          coeff += 0.5 * r*bd.kappa[il];
        }
        coeff += r*kappa_t[ix];
        coeff += M5;

        Idx counter = 4*lattice.counter_accum.back()*Nt + 8*lattice.n_sites*Nt + 4*(lattice.n_sites*s + ix);
        const MS tmp2 = coeff * sigma[0];

        v[counter] = tmp2(0,0); counter++;
        v[counter] = tmp2(0,1); counter++;

        v[counter] = tmp2(1,0); counter++;
        v[counter] = tmp2(1,1); counter++;
      }
    }
  }


  void volume_matrix( std::vector<COOEntry>& elem, const double pow ) const {
    elem.clear();

    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<lattice.n_sites; ix++){
        const Idx i = Nx*s+NS*ix;
        const double v = lattice.dual_areas[ix]/lattice.mean_dual_area;

        elem.push_back( COOEntry( std::pow(v, pow), i, i ) );
        elem.push_back( COOEntry( std::pow(v, pow), i+1, i+1 ) );
      }}
  }





  // template<typename Gauge>
  // void d_coo_format( std::vector<COOEntry>& elem,
  //       	     const Gauge& u,
  //       	     const Link& el ) const {
  //   const Idx ix0 = el[0];
  //   const Idx iy0 = el[1];

  //   elem.clear();
  //   {
  //     // pos
  //     const Idx ix = ix0;
  //     for(int jj=0; jj<lattice.nns[ix].size(); jj++){
  //       const Idx iy = lattice.nns[ix][jj];
  //       if(iy!=iy0) continue;
  //       const Idx il = lattice.map2il.at(Link{ix,iy});
  //       const MS tmp = 0.5 * kappa[il] * ( -r *sigma[0] + gamma(ix, iy) ) * I*std::exp( I* u(Link{ix,iy})) * Omega(ix, iy);

  //       // res[NS*ix] += -tmp(0,0)*v[NS*iy] - tmp(0,1)*v[NS*iy+1];
  //       elem.push_back(COOEntry(tmp(0,0),NS*ix,NS*iy));
  //       elem.push_back(COOEntry(tmp(0,1),NS*ix,NS*iy+1));

  //       // res[NS*ix+1] += -tmp(1,0)*v[NS*iy] - tmp(1,1)*v[NS*iy+1];
  //       elem.push_back(COOEntry(tmp(1,0),NS*ix+1,NS*iy));
  //       elem.push_back(COOEntry(tmp(1,1),NS*ix+1,NS*iy+1));
  //     }
  //   }

  //   {
  //     // neg
  //     const Idx iy = iy0;
  //     for(int jj=0; jj<lattice.nns[iy].size(); jj++){
  //       const Idx ix = lattice.nns[iy0][jj];
  //       if(ix!=ix0) continue;
  //       const Idx il = lattice.map2il.at(Link{ix,iy});
  //       const MS tmp = -0.5 * kappa[il] * ( -r *sigma[0] + gamma(iy, ix) ) * I*std::exp( I* u(Link{iy,ix})) * Omega(iy, ix);

  //       // res[NS*iy] += -tmp(0,0)*v[NS*ix] - tmp(0,1)*v[NS*ix+1];
  //       elem.push_back(COOEntry(tmp(0,0),NS*iy,NS*ix));
  //       elem.push_back(COOEntry(tmp(0,1),NS*iy,NS*ix+1));

  //       // res[NS*iy+1] += -tmp(1,0)*v[NS*ix] - tmp(1,1)*v[NS*ix+1];
  //       elem.push_back(COOEntry(tmp(1,0),NS*iy+1,NS*ix));
  //       elem.push_back(COOEntry(tmp(1,1),NS*iy+1,NS*ix+1));
  //     }
  //   }
  // }




  void set_kappa_t() {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL)
#endif
    for(Idx ix=0; ix<lattice.n_sites; ix++) {
      if(Nt!=1) kappa_t[ix] = c * lattice.dual_areas[ix]/lattice.mean_dual_area;
      else kappa_t[ix] = 0.0;
    }
  }



};
