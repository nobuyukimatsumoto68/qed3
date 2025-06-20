#pragma once

template<typename Lattice, typename Rng>
struct FermionVector {
  const Lattice& lattice;
  std::vector<Complex> field;
  // MatPoly& Op;

  static constexpr Complex I = Complex(0.0, 1.0);
  static constexpr int NS = Comp::NS;

  Rng& rng;

  // cublasHandle_t handle;
  // cudaStream_t stream;
  // const bool is_external;

  // Fermionvector() = delete;

  explicit FermionVector(const Lattice& lattice_,
                         // MatPoly& Op_,
                         Rng& rng_)
    : lattice(lattice_)
    // , Op(Op_)
    , rng(rng_)
    , field(lattice.n_sites*NS, 0.0)
  {}

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

  void set_pt_source(const Idx ix, const int i) {
    for(auto& elem : field) elem = 0.0;
    // field(ix, i) = rng.z2_site( ix ) + I*rng.z2_site( ix );
    // field(ix, i) /= std::sqrt(2.0);
    (*this)(ix, i) = 1.0;
  }

  void set_random() {
    for(Idx ix=0; ix<lattice.n_sites; ix++){
      for(int i=0; i<NS; i++){
        (*this)(ix, i) = rng.z2_site( ix ) + I*rng.z2_site( ix );
        (*this)(ix, i) /= std::sqrt(2.0*lattice.n_sites);
      }}
  }


};
