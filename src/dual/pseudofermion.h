#pragma once

// #include <vector>
#include <cmath>


// template<class Gauge, class Force, class Rng, typename F, typename Grad>
template<typename F, typename Grad>
struct PseudoFermion {
  using Complex = std::complex<double>;
  using T = CuC;

  using Link = std::array<int,2>; // <int,int>;
  using Face = std::vector<int>;
  // using Gauge = U1onS2;
  // using Force = U1onS2;
  // using Rng = ParallelRng;
  using Fermion = Overlap;

  static constexpr Complex I = Complex(0.0, 1.0);
  const int NS=2;

  // const Dirac1fonS2& D;
  // const CGCUDA cg;
  // Fermion* fermion;
  MatPoly& Op_DHD;
  F& f_DH;
  Grad& f_mgrad_DHD;

  CuC *d_phi, *d_eta;
  static constexpr Idx N = Comp::N;

  PseudoFermion()=delete;

  explicit PseudoFermion( // Fermion* fermion_,
                         MatPoly& Op_DHD_,
                         F& f_DH_,
                         Grad& f_mgrad_DHD_)
    : // fermion(fermion_)
    // ,
    Op_DHD(Op_DHD_)
    , f_DH(f_DH_)
    , f_mgrad_DHD(f_mgrad_DHD_)
  {
    CUDA_CHECK(cudaMalloc(&d_phi, N*CD));
    CUDA_CHECK(cudaMalloc(&d_eta, N*CD));
  }

  ~PseudoFermion(){
    CUDA_CHECK(cudaFree(d_phi));
    CUDA_CHECK(cudaFree(d_eta));
  }

  // Complex operator[](const int i) const { return phi[i]; }
  // Complex& operator[](const int i) { return phi[i]; }

  // Complex operator()(const int ix, const int a) const { return phi[NS*ix+a]; }
  // Complex& operator()(const int ix, const int a) { return phi[NS*ix+a]; }


  // void multD( std::vector<Complex>& Dxi, const std::vector<Complex>& xi, const Gauge& U ) const {
  //   assert( Dxi.size()==D.lattice.n_sites*NS );
  //   assert( xi.size()==D.lattice.n_sites*NS );

  //   Complex D_coo[cg.sparse.len], D_csr[cg.sparse.len];
  //   D.coo_format(D_coo, U);
  //   cg.sparse.coo2csr( D_csr, D_coo );
  //   cg.sparse.mult<Complex>( Dxi.data(), xi.data(), D_csr );
  // }

  // void multDH( std::vector<Complex>& DHxi, const std::vector<Complex>& xi, const Gauge& U ) const {
  //   assert( DHxi.size()==D.lattice.n_sites*NS );
  //   assert( xi.size()==D.lattice.n_sites*NS );

  //   Complex D_coo[cg.sparse.len], D_csrH[cg.sparse.len];
  //   D.coo_format(D_coo, U);
  //   cg.sparse.coo2csrH( D_csrH, D_coo );
  //   cg.sparse.multT<Complex>( DHxi.data(), xi.data(), D_csrH );    
  // }

  // void multDHD( std::vector<Complex>& DHDxi, const std::vector<Complex>& xi, const Gauge& U ) const {
  //   assert( DHDxi.size()==D.lattice.n_sites*NS );
  //   assert( xi.size()==D.lattice.n_sites*NS );

  //   Complex D_coo[cg.sparse.len], D_csr[cg.sparse.len], D_csrH[cg.sparse.len];
  //   D.coo_format(D_coo, U);
  //   cg.sparse.coo2csr_csrH( D_csr, D_csrH, D_coo );

  //   std::vector<Complex> tmp( D.lattice.n_sites*NS );
  //   cg.sparse.mult<Complex>( tmp.data(), xi.data(), D_csr );
  //   cg.sparse.multT<Complex>( DHDxi.data(), tmp.data(), D_csrH );    
  // }

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

  // auto begin(){ return phi.begin(); }
  // auto end(){ return phi.end(); }
  // auto begin() const { return phi.begin(); }
  // auto end() const { return phi.end(); }

  // Complex dot( const std::vector<Complex>& eta1, const std::vector<Complex>& xi) const {
  //   assert( eta1.size()==xi.size() );
  //   Complex res = 0.0;
  //   for(int i=0; i<eta1.size(); i++) res += std::conj(eta1[i]) * xi[i];
  //   return res;
  // }

  // Complex dot( const std::vector<Complex>& eta1 ) const {
  //   return dot(this->phi, eta1);
  // }

  double S() const {
    CuC tmp;
    Op_DHD.dot<N>(&tmp, d_phi, d_eta);
    return real(tmp);
  }


  template<class Gauge>
  inline double get_force( const Gauge& U, const Link& ell ) const {
    // const int N = D.lattice.n_sites*NS;
    // std::vector<Complex> dD;
    // std::vector<int> is;
    // std::vector<int> js;
    // D.d_coo_format( dD, is, js, U, ell );
    // std::vector<Complex> dD_eta(N);
    // cg.sparse.multcoo( dD_eta, eta, dD, is, js );
    // std::vector<Complex> DH_dD_eta(N);
    // multDH( DH_dD_eta, dD_eta, U );
    return f_mgrad_DHD( ell, U, d_eta );
  }


  template<typename Gauge, typename Force>
  void get_force( Force& pi, const Gauge& U ) const {
    // Force pi( U.lattice ); // 0 initialized
// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL)
// #endif
    for(int ell=0; ell<U.lattice.n_links; ell++) pi[ell] = get_force( U, U.lattice.links[ell] );
    // return pi;
  }


};


