#pragma once

#include <cmath>


template<typename F, typename Grad, typename Lattice>
struct PseudoFermion {
  // using Complex = std::complex<double>;
  using T = CuC;

  using Link = std::array<Idx,2>; // <int,int>;
  using Face = std::vector<Idx>;
  // using Fermion = Overlap;

  static constexpr Complex I = Complex(0.0, 1.0);
  const int NS=2;

  MatPoly& Op_DHD;
  F& f_DH;
  Grad& f_mgrad_DHD;

  CuC *d_phi, *d_eta;
  static constexpr Idx N = Comp::N;

  Lattice& lattice;

  PseudoFermion()=delete;

  explicit PseudoFermion( MatPoly& Op_DHD_,
                          F& f_DH_,
                          Grad& f_mgrad_DHD_,
                          Lattice& lattice_)
    : Op_DHD(Op_DHD_)
    , f_DH(f_DH_)
    , f_mgrad_DHD(f_mgrad_DHD_)
    , lattice(lattice_)
  {
    CUDA_CHECK(cudaMalloc(&d_phi, N*CD));
    CUDA_CHECK(cudaMalloc(&d_eta, N*CD));
  }

  ~PseudoFermion(){
    CUDA_CHECK(cudaFree(d_phi));
    CUDA_CHECK(cudaFree(d_eta));
  }

  template<class Rng>
  void gen( Rng& rng ) {
    std::vector<Complex> xi(N, 0.0);

    for(Idx ix=0; ix<Comp::N_SITES; ix++) {
      for(int a=0; a<Comp::NS; a++) xi[NS*ix+a] = ( rng.gaussian_site(ix)
                                                    + I*rng.gaussian_site(ix) ) / std::sqrt(2.0);
    }

    CuC *d_xi;
    CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
    CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));
    f_DH( d_phi, d_xi );
    CUDA_CHECK(cudaFree(d_xi));

    update_eta();
  }

  inline void update_eta() {
    Op_DHD.solve<N>( d_eta, d_phi, Comp::TOL_OUTER );
  } // outer CG

  double S() const {
    CuC tmp;
    Op_DHD.dot<N>(&tmp, d_phi, d_eta);
    return real(tmp);
  }


  template<class Gauge>
  inline double get_force( const Gauge& U, const Link& ell ) const {
    return f_mgrad_DHD( ell, U, d_eta );
  }


  template<typename Gauge>
  void get_force( Gauge& pi, const Gauge& u ) const {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL2)
#endif
    for(int ell=0; ell<lattice.n_links; ell++) pi[ell] = get_force( u, lattice.links[ell] );
  }


};


