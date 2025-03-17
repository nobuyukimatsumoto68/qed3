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

    Idx counter=0;
// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL)
// #endif
    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<lattice.n_sites; ix++){
        // Idx counter = 4*lattice.counter_accum.back()*s + 4*lattice.counter_accum[ix];
        for(const Idx iy : lattice.nns[ix]){
          is[counter] = Nx*s+NS*ix; js[counter] = Nx*s+NS*iy; counter++;
          is[counter] = Nx*s+NS*ix; js[counter] = Nx*s+NS*iy+1; counter++;
          // is.push_back(Nx*s+NS*ix); js.push_back(Nx*s+NS*iy);
          // is.push_back(Nx*s+NS*ix); js.push_back(Nx*s+NS*iy+1);

          // is.push_back(Nx*s+NS*ix); js.push_back(Nx*s+NS*ix);
          // is.push_back(Nx*s+NS*ix); js.push_back(Nx*s+NS*ix+1);

          is[counter] = Nx*s+NS*ix+1; js[counter] = Nx*s+NS*iy; counter++;
          is[counter] = Nx*s+NS*ix+1; js[counter] = Nx*s+NS*iy+1; counter++;
          // is.push_back(Nx*s+NS*ix+1); js.push_back(Nx*s+NS*iy);
          // is.push_back(Nx*s+NS*ix+1); js.push_back(Nx*s+NS*iy+1);

          // is.push_back(Nx*s+NS*ix+1); js.push_back(Nx*s+NS*ix);
          // is.push_back(Nx*s+NS*ix+1); js.push_back(Nx*s+NS*ix+1);
        }
      }
    }

    // return ; // @@@ DEBUG

// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL)
// #endif
    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<lattice.n_sites; ix++){
        // Idx counter = 4*lattice.counter_accum.back()*Nt + 8*(lattice.n_sites*s + ix);
        is[counter] = ( Nx*(s+1)+NS*ix )%N; js[counter] = Nx*s+NS*ix; counter++;
        is[counter] =  ( Nx*(s+1)+NS*ix )%N; js[counter] = Nx*s+NS*ix+1; counter++;
        is[counter] = ( Nx*(s-1)+NS*ix + N )%N; js[counter] = Nx*s+NS*ix; counter++;
        is[counter] =  ( Nx*(s-1)+NS*ix + N )%N; js[counter] = Nx*s+NS*ix+1; counter++;
        // is.push_back( Nx*s+NS*ix ); js.push_back( Nx*s+NS*ix );
        // is.push_back( Nx*s+NS*ix ); js.push_back( Nx*s+NS*ix+1 );

        is[counter] = ( Nx*(s+1)+NS*ix+1 )%N; js[counter] = Nx*s+NS*ix; counter++;
        is[counter] = ( Nx*(s+1)+NS*ix+1 )%N; js[counter] = Nx*s+NS*ix+1; counter++;
        is[counter] = ( Nx*(s-1)+NS*ix+1 + N )%N; js[counter] = Nx*s+NS*ix; counter++;
        is[counter] = ( Nx*(s-1)+NS*ix+1 + N )%N; js[counter] = Nx*s+NS*ix+1; counter++;
        // is.push_back( Nx*s+NS*ix+1 ); js.push_back( Nx*s+NS*ix );
        // is.push_back( Nx*s+NS*ix+1 ); js.push_back( Nx*s+NS*ix+1 );
      }
    }

// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL)
// #endif
    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<lattice.n_sites; ix++){
        // Idx counter = 4*lattice.counter_accum.back()*Nt + 8*lattice.n_sites*Nt + 4*(lattice.n_sites*s + ix);
        is[counter] = Nx*s+NS*ix; js[counter] = Nx*s+NS*ix; counter++;
        is[counter] = Nx*s+NS*ix; js[counter] = Nx*s+NS*ix+1; counter++;
        is[counter] = Nx*s+NS*ix+1; js[counter] = Nx*s+NS*ix; counter++;
        is[counter] = Nx*s+NS*ix+1; js[counter] = Nx*s+NS*ix+1; counter++;
      }
    }
  }




  // void coo_structure( std::vector<Idx>& is,
  //                     std::vector<Idx>& js ) const {
  //   // const Idx N = Comp::N;
  //   // const int Nt = Comp::Nt;

  //   for(int s=0; s<Nt; s++){
  //     for(Idx ix=0; ix<lattice.n_sites; ix++){
  //       // Idx counter = lattice.counter_accum[ix];
  //       for(const Idx iy : lattice.nns[ix]){
  //         is.push_back(Nx*s+NS*ix); js.push_back(Nx*s+NS*iy);
  //         is.push_back(Nx*s+NS*ix); js.push_back(Nx*s+NS*iy+1);

  //         is.push_back(Nx*s+NS*ix); js.push_back(Nx*s+NS*ix);
  //         is.push_back(Nx*s+NS*ix); js.push_back(Nx*s+NS*ix+1);

  //         is.push_back(Nx*s+NS*ix+1); js.push_back(Nx*s+NS*iy);
  //         is.push_back(Nx*s+NS*ix+1); js.push_back(Nx*s+NS*iy+1);

  //         is.push_back(Nx*s+NS*ix+1); js.push_back(Nx*s+NS*ix);
  //         is.push_back(Nx*s+NS*ix+1); js.push_back(Nx*s+NS*ix+1);
  //       }
  //     }
  //   }

  //   if(Nt==1) return;
  //   for(int s=0; s<Nt; s++){
  //     for(Idx ix=0; ix<lattice.n_sites; ix++){
  //       is.push_back( ( Nx*(s+1)+NS*ix )%N ); js.push_back( Nx*s+NS*ix );
  //       is.push_back( ( Nx*(s+1)+NS*ix )%N ); js.push_back( Nx*s+NS*ix+1 );
  //       is.push_back( ( Nx*(s-1)+NS*ix + N )%N ); js.push_back( Nx*s+NS*ix );
  //       is.push_back( ( Nx*(s-1)+NS*ix + N )%N ); js.push_back( Nx*s+NS*ix+1 );
  //       is.push_back( Nx*s+NS*ix ); js.push_back( Nx*s+NS*ix );
  //       is.push_back( Nx*s+NS*ix ); js.push_back( Nx*s+NS*ix+1 );

  //       is.push_back( ( Nx*(s+1)+NS*ix+1 )%N ); js.push_back( Nx*s+NS*ix );
  //       is.push_back( ( Nx*(s+1)+NS*ix+1 )%N ); js.push_back( Nx*s+NS*ix+1 );
  //       is.push_back( ( Nx*(s-1)+NS*ix+1 + N )%N ); js.push_back( Nx*s+NS*ix );
  //       is.push_back( ( Nx*(s-1)+NS*ix+1 + N )%N ); js.push_back( Nx*s+NS*ix+1 );
  //       is.push_back( Nx*s+NS*ix+1 ); js.push_back( Nx*s+NS*ix );
  //       is.push_back( Nx*s+NS*ix+1 ); js.push_back( Nx*s+NS*ix+1 );
  //     }
  //   }
  // }






//   template<typename Gauge>
//   void coo_format( std::vector<Complex>& v,
// 		   const Gauge& u ) const {
//     const Idx Nx = Comp::Nx;
//     const int Nt = Comp::Nt;

//     for(Idx i=0; i<Nx*Nt; i++) v[i] = 0.0;

//     // Idx counter=0;
// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL) collapse(2)
// #endif
//     for(int s=0; s<Nt; s++){
//       for(Idx ix=0; ix<lattice.n_sites; ix++){
//         Idx counter = lattice.counter_accum.back()*s + lattice.counter_accum[ix];
//         for(const Idx iy : lattice.nns[ix]){
//           // std::cout << "debug. iy = " << iy << std::endl;
//           const Idx il = lattice.map2il.at(Link{ix,iy});

//           const MS tmp = 0.5 * bd.kappa[il] * ( -r * sigma[0] + bd.gamma(ix, iy) ) * std::exp( I*u.sp(s,Link{ix,iy})) * bd.Omega(ix, iy);
//           // std::cout << "debug. tmp = " << tmp << std::endl;
//           const MS tmp2 = 0.5 * r*bd.kappa[il] * sigma[0] + M5/lattice.nns[ix].size() * sigma[0];

//           // res[NS*ix] += -tmp(0,0)*v[NS*iy] - tmp(0,1)*v[NS*iy+1];
//           v[counter] = tmp(0,0); counter++;
//           v[counter] = tmp(0,1); counter++;

//           // res[NS*ix] += tmp(0,0)*v[NS*ix] + tmp(0,1)*v[NS*ix+1];
//           v[counter] = tmp2(0,0); counter++;
//           v[counter] = tmp2(0,1); counter++;

//           // res[NS*ix+1] += -tmp(1,0)*v[NS*iy] - tmp(1,1)*v[NS*iy+1];
//           v[counter] = tmp(1,0); counter++;
//           v[counter] = tmp(1,1); counter++;

//           // res[NS*ix+1] += tmp(1,0)*v[NS*ix] + tmp(1,1)*v[NS*ix+1];
//           v[counter] = tmp2(1,0); counter++;
//           v[counter] = tmp2(1,1); counter++;
//         }
//       }
//     }

//     if(Nt==1) return;
//     // Idx counter = lattice.counter_accum.back()*Nt;

// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL) collapse(2)
// #endif
//     for(int s=0; s<Nt; s++){
//       for(Idx ix=0; ix<lattice.n_sites; ix++){
//         int signP = 1; if(s==Nt-1) signP = -1;
//         int signM = 1; if(s==0) signM = -1;

//         Idx counter = lattice.counter_accum.back()*Nt + 12 * (lattice.n_sites*s + ix);

//         // const MS tmpP = 0.5 * signP * ( -sigma[0] + sigma[3] );
//         // const MS tmpM = 0.5 * signM * ( -sigma[0] - sigma[3] );
//         // const MS tmpD = sigma[0];

//         const MS tmpP = 0.5 * signP * kappa_t[ix] * ( -r*sigma[0] + sigma[3] ) * std::exp( I*u.tp(s,ix) );
//         const MS tmpM = 0.5 * signM * kappa_t[ix] * ( -r*sigma[0] - sigma[3] ) * std::exp( I*u.tp(s-1,ix) );
//         const MS tmpD = r*kappa_t[ix] * sigma[0];

//         v[counter] = tmpP(0,0); counter++;
//         v[counter] = tmpP(0,1); counter++;
//         v[counter] = tmpM(0,0); counter++;
//         v[counter] = tmpM(0,1); counter++;
//         v[counter] = tmpD(0,0); counter++;
//         v[counter] = tmpD(0,1); counter++;

//         v[counter] = tmpP(1,0); counter++;
//         v[counter] = tmpP(1,1); counter++;
//         v[counter] = tmpM(1,0); counter++;
//         v[counter] = tmpM(1,1); counter++;
//         v[counter] = tmpD(1,0); counter++;
//         v[counter] = tmpD(1,1); counter++;
//       }
//     }
//   }




  template<typename Gauge>
  void coo_format( std::vector<Complex>& v,
		   const Gauge& u ) const {
    const Idx Nx = Comp::Nx;
    const int Nt = Comp::Nt;

// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL)
// #endif
    // for(Idx i=0; i<Nx*Nt; i++) v[i] = 0.0;
    for(Idx i=0; i<v.size(); i++) v[i] = 0.0;

    // Idx counter=0;
// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL)
// #endif
    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<lattice.n_sites; ix++){
        Idx counter = 4*lattice.counter_accum.back()*s + 4*lattice.counter_accum[ix];
        // assert( counter == 4*lattice.counter_accum.back()*s + 4*lattice.counter_accum[ix]);
        for(const Idx iy : lattice.nns[ix]){
          // std::cout << "debug. iy = " << iy << std::endl;
          const Idx il = lattice.map2il.at(Link{ix,iy});

          const MS tmp = 0.5 * bd.kappa[il] * ( -r * sigma[0] + bd.gamma(ix, iy) ) * std::exp( I*u.sp(s,Link{ix,iy})) * bd.Omega(ix, iy);
          // std::cout << "debug. tmp = " << tmp << std::endl;
          // const MS tmp2 = 0.5 * r*bd.kappa[il] * sigma[0] + M5/lattice.nns[ix].size() * sigma[0];

          // res[NS*ix] += -tmp(0,0)*v[NS*iy] - tmp(0,1)*v[NS*iy+1];
          v[counter] = tmp(0,0); counter++;
          v[counter] = tmp(0,1); counter++;

          // // res[NS*ix] += tmp(0,0)*v[NS*ix] + tmp(0,1)*v[NS*ix+1];
          // v[counter] = tmp2(0,0); counter++;
          // v[counter] = tmp2(0,1); counter++;

          // res[NS*ix+1] += -tmp(1,0)*v[NS*iy] - tmp(1,1)*v[NS*iy+1];
          v[counter] = tmp(1,0); counter++;
          v[counter] = tmp(1,1); counter++;

          // // res[NS*ix+1] += tmp(1,0)*v[NS*ix] + tmp(1,1)*v[NS*ix+1];
          // v[counter] = tmp2(1,0); counter++;
          // v[counter] = tmp2(1,1); counter++;
        }
      }
    }

    // return ; // @@@ DEBUG


// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL)
// #endif
    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<lattice.n_sites; ix++){
        int signP = 1; if(s==Nt-1) signP = -1;
        int signM = 1; if(s==0) signM = -1;

        Idx counter = 4*lattice.counter_accum.back()*Nt + 8*(lattice.n_sites*s + ix);
        // assert( counter == 4*lattice.counter_accum.back()*Nt + 8*(lattice.n_sites*s + ix) );

        // const MS tmpP = 0.5 * signP * ( -sigma[0] + sigma[3] );
        // const MS tmpM = 0.5 * signM * ( -sigma[0] - sigma[3] );
        // const MS tmpD = sigma[0];

        const MS tmpP = 0.5 * signP * kappa_t[ix] * ( -r*sigma[0] + sigma[3] ) * std::exp( I*u.tp(s,ix) );
        const MS tmpM = 0.5 * signM * kappa_t[ix] * ( -r*sigma[0] - sigma[3] ) * std::exp( I*u.tp(s-1,ix) );
        // const MS tmpD = r*kappa_t[ix] * sigma[0];

        v[counter] = tmpP(0,0); counter++;
        v[counter] = tmpP(0,1); counter++;
        v[counter] = tmpM(0,0); counter++;
        v[counter] = tmpM(0,1); counter++;
        // v[counter] = tmpD(0,0); counter++;
        // v[counter] = tmpD(0,1); counter++;

        v[counter] = tmpP(1,0); counter++;
        v[counter] = tmpP(1,1); counter++;
        v[counter] = tmpM(1,0); counter++;
        v[counter] = tmpM(1,1); counter++;
        // v[counter] = tmpD(1,0); counter++;
        // v[counter] = tmpD(1,1); counter++;
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
        // assert( counter == 4*lattice.counter_accum.back()*Nt + 8*lattice.n_sites*Nt + 4*(lattice.n_sites*s + ix) );
        // const Idx counter = 4*lattice.counter_accum.back()*Nt + 4*s;
        const MS tmp2 = coeff * sigma[0];

        v[counter] = tmp2(0,0); counter++;
        v[counter] = tmp2(0,1); counter++;

        v[counter] = tmp2(1,0); counter++;
        v[counter] = tmp2(1,1); counter++;
      }
    }

  }

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
