#pragma once

double Ylm_real( const int ell, const int em, const double theta, const double phi ){
  const double Pell = std::assoc_legendre( ell, std::abs(em), std::cos(theta) );

  double trig;
  if( em>0 ) trig = std::cos(std::abs(em)*phi);
  else if( em<0 ) trig = std::sin(std::abs(em)*phi);
  else trig = 1.0/std::sqrt(2.0);

  double factor = (2.0*ell+1)/(4.0*M_PI);
  factor *= std::tgamma( ell-std::abs(em)+1 );
  factor /= std::tgamma( ell+std::abs(em)+1 );
  factor = std::pow(-1.0,std::abs(em)) * std::sqrt(2.0) * std::sqrt(factor);

  return factor * Pell * trig;
}

double Ylm_real( const int ell, const int em, const VE r ){
  assert( std::abs( r.norm()-1.0 )<1.0e-14 );
  const double theta = std::acos(r[2]);
  const double phi = std::atan2(r[1], r[0]);
  return Ylm_real( ell, em, theta, phi );
}


Complex Ylm( const int ell, const int em, const double theta, const double phi ){
  const double Pell = std::assoc_legendre( ell, std::abs(em), std::cos(theta) );

  Complex trig = std::exp( I*(em*phi) );

  double factor = (2.0*ell+1)/(4.0*M_PI);
  factor *= std::tgamma( ell-std::abs(em)+1 );
  factor /= std::tgamma( ell+std::abs(em)+1 );
  factor = std::pow(-1.0,std::abs(em)) * std::sqrt(factor);

  return factor * Pell * trig;
}

Complex Ylm( const int ell, const int em, const VE r ){
  assert( std::abs( r.norm()-1.0 )<1.0e-14 );
  const double theta = std::acos(r[2]);
  const double phi = std::atan2(r[1], r[0]);
  return Ylm( ell, em, theta, phi );
}



// template<typename Lattice>
struct FermionVector {
  // const Lattice& lattice;
  // std::vector<Complex> field;
  Complex* field;
  // MatPoly& Op;

  // static constexpr Complex I = Complex(0.0, 1.0);
  // static constexpr int NS = Comp::NS;

  int Nt;

  // Rng& rng;

  // cublasHandle_t handle;
  // cudaStream_t stream;
  // const bool is_external;

  // Fermionvector() = delete;

  FermionVector() // Rng& rng_)
    : Nt(Comp::Nt)
      // , Op(Op_)
      // , rng(rng_)
      // , field(Comp::N, 0.0)
  {
    CUDA_CHECK( cudaMallocHost( &field, Comp::N*CD ) );
    memset(field, 0, Comp::N*CD);
  }
  ~FermionVector(){
    CUDA_CHECK(cudaFreeHost(field));
  }


  // explicit FermionVector(const Lattice& lattice_,
  //                        const int Nt_,
  //                        Rng& rng_)
  //   : lattice(lattice_)
  //   , Nt(Nt_)
  //     // , Op(Op_)
  //   , rng(rng_)
  //   , field(Comp::Nx*Nt, 0.0)
  // {}

  // GaugeForce& operator=(const GaugeForce& other){
  //   if (this == &other) return *this;

  //   assert(&lattice==&other.lattice);
  //   field = other.field;
  //   return *this;
  // }

  // auto begin(){ return field.begin(); }
  // auto end(){ return field.end(); }
  // auto begin() const { return field.begin(); }
  // auto end() const { return field.end(); }

  // GaugeForce & operator=(const GaugeForce&) = delete;

  // Complex operator()(const Idx ix, const int i) const { return field[NS*ix+i]; }
  // Complex& operator()(const Idx ix, const int i) { return field[NS*ix+i]; }

  Complex operator()(const int s, const Idx ix, const int i) const { return field[Comp::Nx*s+NS*ix+i]; }
  Complex& operator()(const int s, const Idx ix, const int i) { return field[Comp::Nx*s+NS*ix+i]; }

  // void set_pt_source(const Idx ix, const int i) {
  //   // for(auto& elem : field) elem = 0.0;
  //   memset(field, 0, Comp::N*CD);
  //   // field(ix, i) = rng.z2_site( ix ) + I*rng.z2_site( ix );
  //   // field(ix, i) /= std::sqrt(2.0);
  //   (*this)(0, ix, i) = 1.0;
  // }

  void set_pt_source(const int s, const Idx ix, const int i) {
    // for(auto& elem : field) elem = 0.0;
    memset(field, 0, Comp::N*CD);
    // field(ix, i) = rng.z2_site( ix ) + I*rng.z2_site( ix );
    // field(ix, i) /= std::sqrt(2.0);
    (*this)(s, ix, i) = 1.0;
  }

  template <typename Rng>
  void set_random_gauge(Rng& rng, const double width=1.0) {
    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<rng.lattice.n_sites; ix++){
        (*this)(s,ix,0) = width*rng.gaussian_site(s,ix);
        (*this)(s,ix,1) = (*this)(s,ix,0);
      }
    }
  }

  template <typename Rng>
  void fill_z2_source(Rng& rng) {
    memset(field, 0, Comp::N*CD);
    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<rng.lattice.n_sites; ix++){
        for(int i=0; i<2; i++){
          (*this)(s, ix, i) = rng.CZ2_site( s, ix );
        }}}
  }

  template <typename Rng>
  void fill_z2_wall_source(Rng& rng, const int s) {
    memset(field, 0, Comp::N*CD);
    for(Idx ix=0; ix<rng.lattice.n_sites; ix++){
      for(int i=0; i<2; i++){
        (*this)(s, ix, i) = rng.CZ2_site( s, ix );
      }}
  }

  void gauge_trsf(const FermionVector& gauge, const double sign=1.0) {
    for(Idx i=0; i<Comp::N; i++) field[i] *= std::exp( sign*I*gauge.field[i] );
  }

  VC slice( const int s ) const {
    VC res = VC::Zero( Comp::Nx );
    for(Idx ix=0; ix<Comp::N_SITES; ix++){
      for(int i=0; i<2; i++){
        res[NS*ix+i] = (*this)(s, ix, i);
      }}
    return res;
  }


  void mult_sigma3() {
    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<Comp::N_SITES; ix++){
        (*this)(s, ix, 1) *= -1.0;
      }}
  }

  template<typename Lattice>
  void mult_Y(const int ell, const int em, const Lattice& lattice ) {
    for(Idx ix=0; ix<Comp::N_SITES; ix++){
      const VE r = lattice.sites[ix];
      const double theta = std::acos(r[2]);
      const double phi = std::atan2(r[1], r[0]);
      const double ylm = Ylm_real(ell, em, theta, phi);
      for(int s=0; s<Nt; s++){
        (*this)(s, ix, 0) *= ylm;
        (*this)(s, ix, 1) *= ylm;
      }
    }
  }

  // void set_random() {
  //   for(Idx ix=0; ix<lattice.n_sites; ix++){
  //     for(int i=0; i<NS; i++){
  //       (*this)(ix, i) = rng.z2_site( ix ) + I*rng.z2_site( ix );
  //       (*this)(ix, i) /= std::sqrt(2.0*lattice.n_sites);
  //     }}
  // }


};















template<class Gauge>
struct GaugeVector {
  using BaseLink = std::array<int,2>; // <int,int>;

  const int Nt;
  const Idx Ng;
  std::vector<Complex> field;
  const Gauge U;

  GaugeVector(const Gauge& U_) // Rng& rng_)
    : Nt(Comp::Nt)
    , Ng( Comp::Nt*(Comp::N_LINKS + Comp::N_SITES) )
    , field(Ng)
    , U(U_)
  {
    set_zero();
  }

  auto begin(){ return field.begin(); }
  auto end(){ return field.end(); }
  auto begin() const { return field.begin(); }
  auto end() const { return field.end(); }

  // Complex operator()(const Idx ix, const int i) const { return field[NS*ix+i]; }
  // Complex& operator()(const Idx ix, const int i) { return field[NS*ix+i]; }

  // Complex operator()(const int s, const Idx ix, const int i) const { return field[Comp::Nx*s+NS*ix+i]; }
  // Complex& operator()(const int s, const Idx ix, const int i) { return field[Comp::Nx*s+NS*ix+i]; }

  Complex sp(const int s, const Idx& il) const { return field[ U.idx_sp(s,il) ]; }
  Complex& sp(const int s, const Idx& il) { return field[ U.idx_sp(s,il) ]; }
  Complex sp(const int s, const BaseLink& ell) const { return field[ U.idx_sp(s,ell) ]; }
  Complex& sp(const int s, const BaseLink& ell) { return field[ U.idx_sp(s,ell) ]; }

  Complex tp(const int s, const Idx& ix) const { return field[ U.idx_tp(s,ix) ]; }
  Complex& tp(const int s, const Idx& ix) { return field[ U.idx_tp(s,ix) ]; }

  void set_zero(){
    for(auto& elem : field) elem = 0.0;
  }

  Complex dot( const GaugeVector& other ) const {
    Complex res = 0.0;
    for(Idx i=0; i<Ng; i++) res += std::conj(this->field[i]) * other.field[i];
    return res;
  }




  // void set_pt_source(const Idx ix, const int i) {
  //   set_zero();
  //   (*this)(ix, i) = 1.0;
  // }
  // void set_pt_source(const int s, const Idx ix, const int i) {
  //   set_zero();
  //   (*this)(s, ix, i) = 1.0;
  // }



};
