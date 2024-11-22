#pragma once

#include <vector>
#include <cmath>

#include "cg_cuda.h"


struct PseudoFermion {
  using Link = std::array<int,2>; // <int,int>;
  using Face = std::vector<int>;
  using Gauge = U1onS2;
  using Rng = ParallelRng;

  using Complex = std::complex<double>;
  static constexpr Complex I = Complex(0.0, 1.0);
  const int NS=2;

  // using MS=Eigen::Matrix2cd;
  // using Spinor=Eigen::Vector2cd;
  // using VE=Eigen::Vector3d;
  // using VC=Eigen::VectorXcd;

  // const Lattice& lattice;
  const Dirac1fonS2& D;
  const CGCUDA cg;

  std::vector<Complex> phi;

  PseudoFermion()=delete;

  PseudoFermion(const Dirac1fonS2& D_)
    : D(D_)
    , cg(D)
    , phi(D.lattice.n_sites*NS, 0.0)
  {}

  PseudoFermion(const Dirac1fonS2& D_, const Gauge& U, Rng& rng)
    : D(D_)
    , cg(D)
    , phi(D.lattice.n_sites*NS, 0.0)
  {
    gen( U, rng );
  }

  // PseudoFermion( const PseudoFermion& other )
  //   : D(other.D)
  //   , cg(other.cg)
  //   , phi(other.phi)
  // {}

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

    // Complex D_coo[cg.sparse.len], D_csrH[cg.sparse.len];
    // D.coo_format(D_coo, U);
    // cg.sparse.coo2csrH( D_csrH, D_coo );
    // cg.sparse.multT<Complex>( phi.data(), xi.data(), D_csrH );
    multDH( phi, xi, U );
  }

  std::vector<Complex> get_eta( const Gauge& U ) const {
    const int N = D.lattice.n_sites*NS;
    std::vector<Complex> eta(N, 0.0);
    cg( eta.data(), phi.data(), U );
    return eta;
  };

  using Force=U1onS2;


  auto begin(){ return phi.begin(); }
  auto end(){ return phi.end(); }
  auto begin() const { return phi.begin(); }
  auto end() const { return phi.end(); }


  Complex dot( const std::vector<Complex> eta, const std::vector<Complex> phi) const {
    assert( eta.size()==phi.size() );
    Complex res = 0.0;
    for(int i=0; i<eta.size(); i++) res += std::conj(eta[i]) * phi[i];
    return res;
  }
  

  double get_force( const Gauge& U, const Link& ell ) const {
    const int N = D.lattice.n_sites*NS;
    std::vector<Complex> eta = get_eta(U);

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


  Force get_force( const Gauge& U ) const {
    Force pi( U.lattice ); // 0 initialized
    for(int ell=0; ell<U.lattice.n_links; ell++) pi[ell] = get_force( U, U.lattice.links[ell] );
    return pi;
  }


  // PseudoFermion & operator=(const PseudoFermion&) = delete;
  // PseudoFermion& operator=(const PseudoFermion& other){
  //   if (this == &other) return *this;
  //   assert(&D==&other.D);
  //   eta = other.eta;
  //   phi = other.phi;
  //   return *this;
  // }




  // double operator()(const Link& ell) const { // recommended
  //   const int il = lattice.map2il.at(ell);
  //   const int sign = lattice.map2sign.at(ell);
  //   return sign * field[il];
  // }

  // PseudoFermion& operator+=(const PseudoFermion& rhs){
  //   for(int i=0; i<field.size(); i++) field[i] += rhs.field[i];
  //   return *this;
  // }

  // PseudoFermion& operator*=(const double rhs){
  //   for(int i=0; i<field.size(); i++) field[i] *= rhs;
  //   return *this;
  // }

  // PseudoFermion& operator/=(const double rhs){
  //   for(int i=0; i<field.size(); i++) field[i] /= rhs;
  //   return *this;
  // }

  // friend PseudoFermion operator*(PseudoFermion v, const double a) {
  //   for(int i=0; i<v.field.size(); i++) v.field[i] *= a;
  //   return v;
  // }

  // friend PseudoFermion operator*(const double a, PseudoFermion v) {
  //   for(int i=0; i<v.field.size(); i++) v.field[i] *= a;
  //   return v;
  // }

};


