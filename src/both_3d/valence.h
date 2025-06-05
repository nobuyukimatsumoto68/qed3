#pragma once

// template<typename Lattice>
struct FermionVector {
  // const Lattice& lattice;
  std::vector<Complex> field;
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
    , field(Comp::N, 0.0)
  {}

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

  auto begin(){ return field.begin(); }
  auto end(){ return field.end(); }
  auto begin() const { return field.begin(); }
  auto end() const { return field.end(); }

  // GaugeForce & operator=(const GaugeForce&) = delete;

  Complex operator()(const Idx ix, const int i) const { return field[NS*ix+i]; }
  Complex& operator()(const Idx ix, const int i) { return field[NS*ix+i]; }

  Complex operator()(const int s, const Idx ix, const int i) const { return field[Comp::Nx*s+NS*ix+i]; }
  Complex& operator()(const int s, const Idx ix, const int i) { return field[Comp::Nx*s+NS*ix+i]; }

  void set_pt_source(const Idx ix, const int i) {
    for(auto& elem : field) elem = 0.0;
    // field(ix, i) = rng.z2_site( ix ) + I*rng.z2_site( ix );
    // field(ix, i) /= std::sqrt(2.0);
    (*this)(ix, i) = 1.0;
  }

  void set_pt_source(const int s, const Idx ix, const int i) {
    for(auto& elem : field) elem = 0.0;
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

  void gauge_trsf(const FermionVector& gauge, const double sign=1.0) {
    for(Idx i=0; i<field.size(); i++) field[i] *= std::exp( sign*I*gauge.field[i] );
  }


  // void set_random() {
  //   for(Idx ix=0; ix<lattice.n_sites; ix++){
  //     for(int i=0; i<NS; i++){
  //       (*this)(ix, i) = rng.z2_site( ix ) + I*rng.z2_site( ix );
  //       (*this)(ix, i) /= std::sqrt(2.0*lattice.n_sites);
  //     }}
  // }


};
