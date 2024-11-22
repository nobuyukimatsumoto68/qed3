#pragma once

struct GaugeForce {
  using Link = std::array<int,2>; // <int,int>;
  using Face = std::vector<int>;

  Lattice& lattice;
  std::vector<double> field;

  GaugeForce() = delete;
  GaugeForce(const Lattice& lattice_)
    : lattice(lattice_)
    , field(lattice.n_links, 0.0)
  {}

  GaugeForce( const U1onS2& other )
    : lattice(other.lattice)
    , field(other.field)
  {}

  GaugeForce & operator=(const GaugeForce&) = delete;

  double operator[](const int i) const { return field[i]; }
  double& operator[](const int i) { return field[i]; }

  double operator()(const Link& ell) const { // recommended
    const int il = lattice.map2il.at(ell);
    const int sign = lattice.map2sign.at(ell);
    return sign * field[il] ;
  }

};
