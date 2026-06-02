#pragma once

double Ylm_real( const int ell, const int em, const double theta, const double phi ){
  const double Pell = std::assoc_legendre( ell, std::abs(em), std::cos(theta) ); // no Condon-Shortley phase

  double trig;
  if( em>0 ) trig = std::cos(std::abs(em)*phi);
  else if( em<0 ) trig = std::sin(std::abs(em)*phi);
  else trig = 1.0/std::sqrt(2.0);

  double factor = (2.0*ell+1)/(4.0*M_PI);
  factor *= std::tgamma( ell-std::abs(em)+1 );
  factor /= std::tgamma( ell+std::abs(em)+1 );
  // factor = std::pow(-1.0,std::abs(em)) * std::sqrt(2.0) * std::sqrt(factor);
  factor = std::sqrt(2.0) * std::sqrt(factor);

  return factor * Pell * trig;
}

double Ylm_real( const int ell, const int em, const VE r ){
  assert( std::abs( r.norm()-1.0 )<1.0e-14 );
  const double theta = std::acos(r[2]);
  const double phi = std::atan2(r[1], r[0]);
  return Ylm_real( ell, em, theta, phi );
}


// Complex Ylm( const int ell, const int em, const double theta, const double phi ){
//   const double Pell = std::assoc_legendre( ell, std::abs(em), std::cos(theta) );

//   Complex trig = std::exp( I*(em*phi) );

//   double factor = (2.0*ell+1)/(4.0*M_PI);
//   factor *= std::tgamma( ell-std::abs(em)+1 );
//   factor /= std::tgamma( ell+std::abs(em)+1 );
//   factor = std::pow(-1.0,std::abs(em)) * std::sqrt(factor);

//   return factor * Pell * trig;
// }

// Complex Ylm( const int ell, const int em, const VE r ){
//   assert( std::abs( r.norm()-1.0 )<1.0e-14 );
//   const double theta = std::acos(r[2]);
//   const double phi = std::atan2(r[1], r[0]);
//   return Ylm( ell, em, theta, phi );
// }



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

  FermionVector& operator=(const FermionVector& other) {
    if (this == &other) return *this;
    memcpy( this->data(), other.data(), size()*CD );
    return *this;
  }

  Idx size() const { return Comp::N; }
  Complex* data() { return field; }
  Complex* data() const { return field; }


  Complex operator()(const int s, const Idx ix, const int i) const { return field[Comp::Nx*s+NS*ix+i]; }
  Complex& operator()(const int s, const Idx ix, const int i) { return field[Comp::Nx*s+NS*ix+i]; }

  // }

  void set_pt_source(const int s, const Idx ix, const int i) {
    memset(field, 0, Comp::N*CD);
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


  void fill_one() {
    memset(field, 0, Comp::N*CD);
    for(Idx i=0; i<Comp::N; i++) field[i] = 1;
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
  void fill_z2_pt_source(Rng& rng, const int s, const Idx ix) {
    memset(field, 0, Comp::N*CD);
    for(int i=0; i<2; i++){
      (*this)(s, ix, i) = rng.CZ2_site( s, ix );
    }
  }


  template <typename Rng>
  void fill_z2_wall_source(Rng& rng, const int s) {
    memset(field, 0, Comp::N*CD);
    for(Idx ix=0; ix<rng.lattice.n_sites; ix++){
      for(int i=0; i<2; i++){
        (*this)(s, ix, i) = rng.CZ2_site( s, ix );
      }}
  }

  template <typename Rng>
  void time_spin_dilution(Rng& rng, const int t_s, const int t_block, const int spin) {
    memset(field, 0, Comp::N*CD);
    const int interval = Comp::Nt / t_block;
    for(int t=t_s; t<Comp::Nt; t+=interval){
      for(Idx ix=0; ix<Comp::N_SITES; ix++){
        (*this)(t, ix, spin) = rng.CZ2_site(t, ix);
      }
    }
  }

  void accumulate_loop(std::vector<Complex>& summed_trace_per_timeslice,
                       const FermionVector& phi,
                       const int t_s, const int t_block, const int spin) const {
    const double c = (spin == 0) ? +1.0 : -1.0;
    const int interval = Comp::Nt / t_block;
    for(int t=t_s; t<Comp::Nt; t+=interval){
      Complex loop_t(0.0, 0.0);
      for(Idx ix=0; ix<Comp::N_SITES; ix++){
        loop_t += std::conj((*this)(t, ix, spin)) * c * phi(t, ix, spin);
      }
      summed_trace_per_timeslice[t] += loop_t / static_cast<double>(Comp::N_SITES);
    }
  }

  void accumulate_loop_gamma(std::vector<Complex>& L,
                             const FermionVector& Gamma_phi,
                             const int t_s, const int t_block, const int spin) const {
    const int interval = Comp::Nt / t_block;
    for(int t=t_s; t<Comp::Nt; t+=interval){
      Complex loop_t(0.0, 0.0);
      for(Idx ix=0; ix<Comp::N_SITES; ix++){
        loop_t += std::conj((*this)(t, ix, spin)) * Gamma_phi(t, ix, spin);
      }
      L[t] += loop_t / static_cast<double>(Comp::N_SITES);
    }
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

  void mult_sigma1() {
    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<Comp::N_SITES; ix++){
        auto tmp0 = (*this)(s, ix, 0);
        auto tmp1 = (*this)(s, ix, 1);
        (*this)(s, ix, 0) = tmp1;
        (*this)(s, ix, 1) = tmp0;
        // std::swap( (*this)(s, ix, 0), (*this)(s, ix, 1) );
      }}
  }

  void mult_sigma2() {
    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<Comp::N_SITES; ix++){
        // std::swap( (*this)(s, ix, 0), (*this)(s, ix, 1) );
        // (*this)(s, ix, 0) *= -I;
        // (*this)(s, ix, 1) *= I;
        auto tmp0 = (*this)(s, ix, 0);
        auto tmp1 = (*this)(s, ix, 1);
        (*this)(s, ix, 0) = -I*tmp1;
        (*this)(s, ix, 1) = I*tmp0;
      }}
  }

  void mult_sigma3() {
    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<Comp::N_SITES; ix++){
        (*this)(s, ix, 1) *= -1.0;
      }}
  }

  void mult_sigma( const int a) {
    if(a==1) mult_sigma1();
    else if(a==2) mult_sigma2();
    else if(a==3) mult_sigma3();
    else if(a==0) ;
    else assert(false);
  }

  template<typename Lattice>
  void mult_Ylm_real(const int ell, const int em, const Lattice& lattice ) {
    for(Idx ix=0; ix<Comp::N_SITES; ix++){
      const VE r = lattice.sites[ix];
      const double area = lattice.dual_areas[ix];
      // const double theta = std::acos(r[2]);
      // const double phi = std::atan2(r[1], r[0]);
      const double ylm = Ylm_real(ell, em, r);
      for(int s=0; s<Nt; s++){
        (*this)(s, ix, 0) *= area * ylm;
        (*this)(s, ix, 1) *= area * ylm;
      }
    }
  }

  template<typename Lattice>
  void mult_Ylm_real_nomeasure(const int ell, const int em, const Lattice& lattice ) {
    for(Idx ix=0; ix<Comp::N_SITES; ix++){
      const VE r = lattice.sites[ix];
      // const double area = lattice.dual_areas[ix];
      // const double theta = std::acos(r[2]);
      // const double phi = std::atan2(r[1], r[0]);
      const double ylm = Ylm_real(ell, em, r);
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






struct FermionMatrix {
  std::vector<FermionVector> eta;

  FermionMatrix() // Rng& rng_)
    : eta(Comp::NS)
  {}

  // copy constructor
  FermionMatrix( const FermionMatrix& other )
    : eta(Comp::NS)
  {
    for(int spin=0; spin<2; spin++){
      memcpy( eta[spin].data(), other.eta[spin].data(), eta[spin].size()*CD );
    }
  }


  ~FermionMatrix(){}

  FermionVector operator[](const int j) const { return eta[j]; }
  FermionVector& operator[](const int j) { return eta[j]; }

  inline Complex operator()(const int s, const Idx ix,
                            const int i, const int j
                            ) const { return eta[j](s, ix, i); }
  inline Complex& operator()(const int s, const Idx ix,
                             const int i, const int j) { return eta[j](s, ix, i); }

  void mult_sigma( const int a) {
    for(int spin=0; spin<2; spin++) eta[spin].mult_sigma(a);
  }

  template<typename Lattice>
  void mult_Ylm(const int ell, const int em, const Lattice& lattice ) {
    for(int spin=0; spin<2; spin++) eta[spin].mult_Ylm_real(ell, em, lattice);
  }

  MS get_spinmatrix( const int s, const Idx ix ) const {
    MS res = MS::Zero();
    for(int i=0; i<Comp::NS; i++){
      for(int j=0; j<Comp::NS; j++){
        res(i, j) = (*this)(s, ix, i, j);
      }}
    return res;
  }


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
