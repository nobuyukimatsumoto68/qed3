#pragma once

#include <vector>
#include <cmath>

#include "cg_cuda.h"

struct PseudoFermion {
  using Complex = std::complex<double>;
  using Link = std::array<int,2>; // <int,int>;
  using Face = std::vector<int>;
  using Gauge = U1onS2;
  using Force=U1onS2;
  using Rng = ParallelRng;

  static constexpr Complex I = Complex(0.0, 1.0);
  const int NS=2;

  const Dirac1fonS2& D;
  const CGCUDA cg;

  std::vector<Complex> phi;
  std::vector<Complex> eta;

  PseudoFermion()=delete;

  PseudoFermion(const Dirac1fonS2& D_)
    : D(D_)
    , cg(D)
    , phi(D.lattice.n_sites*NS, 0.0)
    , eta(D.lattice.n_sites*NS, 0.0)
  {}


  Complex operator[](const int i) const { return phi[i]; }
  Complex& operator[](const int i) { return phi[i]; }

  Complex operator()(const int ix, const int a) const { return phi[NS*ix+a]; }
  Complex& operator()(const int ix, const int a) { return phi[NS*ix+a]; }


  void multD( std::vector<Complex>& Dxi, const std::vector<Complex>& xi, const Gauge& U ) const {
    assert( Dxi.size()==D.lattice.n_sites*NS );
    assert( xi.size()==D.lattice.n_sites*NS );

    Complex D_coo[cg.sparse.len], D_csr[cg.sparse.len];
    D.coo_format(D_coo, U);
    cg.sparse.coo2csr( D_csr, D_coo );
    cg.sparse.mult<Complex>( Dxi.data(), xi.data(), D_csr );
  }


  void multDH( std::vector<Complex>& DHxi, const std::vector<Complex>& xi, const Gauge& U ) const {
    assert( DHxi.size()==D.lattice.n_sites*NS );
    assert( xi.size()==D.lattice.n_sites*NS );

    Complex D_coo[cg.sparse.len], D_csrH[cg.sparse.len];
    D.coo_format(D_coo, U);
    cg.sparse.coo2csrH( D_csrH, D_coo );
    cg.sparse.multT<Complex>( DHxi.data(), xi.data(), D_csrH );    
  }


  void multDHD( std::vector<Complex>& DHDxi, const std::vector<Complex>& xi, const Gauge& U ) const {
    assert( DHDxi.size()==D.lattice.n_sites*NS );
    assert( xi.size()==D.lattice.n_sites*NS );

    Complex D_coo[cg.sparse.len], D_csr[cg.sparse.len], D_csrH[cg.sparse.len];
    D.coo_format(D_coo, U);
    cg.sparse.coo2csr_csrH( D_csr, D_csrH, D_coo );

    std::vector<Complex> tmp( D.lattice.n_sites*NS );
    cg.sparse.mult<Complex>( tmp.data(), xi.data(), D_csr );
    cg.sparse.multT<Complex>( DHDxi.data(), tmp.data(), D_csrH );    
  }


  void gen( const Gauge& U, Rng& rng ) {
    const int N = D.lattice.n_sites*NS;
    std::vector<Complex> xi(N, 0.0);

    for(int ix=0; ix<D.lattice.n_sites; ix++) for(int a=0; a<NS; a++) xi[NS*ix+a] = ( rng.gaussian_site(ix) + I*rng.gaussian_site(ix) ) / std::sqrt(2.0);

    multDH( phi, xi, U );

    update_eta(U);
  }


  void update_eta( const Gauge& U ) { cg( eta.data(), phi.data(), U ); }


  auto begin(){ return phi.begin(); }
  auto end(){ return phi.end(); }
  auto begin() const { return phi.begin(); }
  auto end() const { return phi.end(); }


  Complex dot( const std::vector<Complex>& eta1, const std::vector<Complex>& xi) const {
    assert( eta1.size()==xi.size() );
    Complex res = 0.0;
    for(int i=0; i<eta1.size(); i++) res += std::conj(eta1[i]) * xi[i];
    return res;
  }

  Complex dot( const std::vector<Complex>& eta1 ) const {
    return dot(this->phi, eta1);
  }

  double S() const { return dot( eta ).real(); }


  double get_force( const Gauge& U, const Link& ell ) const {
    const int N = D.lattice.n_sites*NS;

    std::vector<Complex> dD;
    std::vector<int> is;
    std::vector<int> js;
    D.d_coo_format( dD, is, js, U, ell );

    std::vector<Complex> dD_eta(N);
    cg.sparse.multcoo( dD_eta, eta, dD, is, js );

    std::vector<Complex> DH_dD_eta(N);
    multDH( DH_dD_eta, dD_eta, U );

    return -2.0 * dot( eta, DH_dD_eta ).real();
  }


  Force dS( const Gauge& U ) const {
    Force pi( U.lattice ); // 0 initialized
#ifdef _OPENMP
#pragma omp parallel for num_threads(nparallel)
#endif
    for(int ell=0; ell<U.lattice.n_links; ell++) pi[ell] = get_force( U, U.lattice.links[ell] );
    return pi;
  }


};


