#pragma once

template<typename Lattice>
struct GaugeForce {
  using Link = std::array<int,2>; // <int,int>;
  using Face = std::vector<int>;

  const Lattice& lattice;
  std::vector<double> field;

  GaugeForce() = delete;
  explicit GaugeForce(const Lattice& lattice_)
    : lattice(lattice_)
    , field(lattice.n_links, 0.0)
  {}

  // GaugeForce( const U1onS2& other )
  //   : lattice(other.lattice)
  //   , field(other.field)
  // {}

  GaugeForce& operator=(const GaugeForce& other){
    if (this == &other) return *this;

    assert(&lattice==&other.lattice);
    field = other.field;
    return *this;
  }

  auto begin(){ return field.begin(); }
  auto end(){ return field.end(); }
  auto begin() const { return field.begin(); }
  auto end() const { return field.end(); }

  // GaugeForce & operator=(const GaugeForce&) = delete;

  double operator[](const int i) const { return field[i]; }
  double& operator[](const int i) { return field[i]; }

  double operator()(const Link& ell) const { // recommended
    const int il = lattice.map2il.at(ell);
    const int sign = lattice.map2sign.at(ell);
    return sign * field[il] ;
  }

  template <typename Rng>
  void gaussian(Rng& rng, const double width=1.0) {
    for(int il=0; il<lattice.n_links; il++) field[il] = width*rng.gaussian_link(il);
  }


};
