#pragma once

#include <vector>
#include <cmath>


template<typename Lattice>
struct U1onS2 {
  using Link = std::array<int,2>; // <int,int>;
  using Face = std::vector<int>;
  // using Lattice = S2Trivalent;

  Lattice& lattice;
  std::vector<double> field;

  U1onS2()=delete;

  explicit U1onS2(Lattice& lattice_)
    : lattice(lattice_)
    , field(lattice.n_links, 0.0)
  {}

  U1onS2( const U1onS2& other )
    : lattice(other.lattice)
    , field(other.field)
  {}

  U1onS2& operator=(const U1onS2& other){
    if (this == &other) return *this;

    assert(&lattice==&other.lattice);
    field = other.field;
    return *this;
  }

  double operator[](const int i) const { return field[i]; }
  double& operator[](const int i) { return field[i]; }

  auto begin(){ return field.begin(); }
  auto end(){ return field.end(); }
  auto begin() const { return field.begin(); }
  auto end() const { return field.end(); }

  double operator()(const Link& ell) const { // recommended
    const int il = lattice.map2il.at(ell);
    const int sign = lattice.map2sign.at(ell);
    return sign * field[il];
  }

  template <typename Force>
  U1onS2& operator+=(const Force& rhs){
    for(int i=0; i<field.size(); i++) field[i] += rhs.field[i];
    return *this;
  }
  template <typename Force>
  friend Force operator+(U1onS2 v, const Force& w) { v += w; return v; }

  U1onS2& operator*=(const double rhs){
    for(int i=0; i<field.size(); i++) field[i] *= rhs;
    return *this;
  }

  U1onS2& operator/=(const double rhs){
    for(int i=0; i<field.size(); i++) field[i] /= rhs;
    return *this;
  }

  friend U1onS2 operator*(U1onS2 v, const double a) {
    for(int i=0; i<v.field.size(); i++) v.field[i] *= a;
    return v;
  }

  friend U1onS2 operator*(const double a, U1onS2 v) {
    for(int i=0; i<v.field.size(); i++) v.field[i] *= a;
    return v;
  }

  double norm() const {
    double sum = 0.0;
    for(const auto elem : field) sum += elem*elem;
    return std::sqrt( sum );
  }

  void print2log_norm(const std::string comment=" : ") const {
    std::clog << comment << norm() << std::endl;
  }

  double plaquette_angle(const Face& face) const {
    double sum = 0.0;
    for(int i=0; i<face.size(); i++) {
      const int ix = face[i];
      const int iy = face[(i+1)%face.size()];
      sum += (*this)( Link{ix,iy} );
    }
    return sum;
  }

  double plaquette_angle(const int i_face) const {
    return plaquette_angle(lattice.faces[i_face]);
  }

  // double average_plaquette() const {
  //   double sum = 0.0;
  //   for(int i=0; i<lattice.n_faces; i++) {
  //     sum += std::cos(plaquette_angle(i));
  //   }
  //   return sum/lattice.n_faces;
  // }

  // void dist_plaqsq( double& mean, double& square) const {
  //   mean = 0.0;
  //   square = 0.0;
  //   for(int i=0; i<lattice.n_faces; i++) {
  //     double val = plaquette_angle(i);
  //     if constexpr(is_compact){
  //       while(val>M_PI) val -= 2.0*M_PI;
  //       while(val<-M_PI) val += 2.0*M_PI;
  //     }
  //     const double factor = 1.0 / lattice.vps[i];

  //     const double v = factor * std::pow( val, 2 );
  //     mean += v;
  //     square += v*v;
  //   }
  //   mean /= lattice.n_faces;
  //   square /= lattice.n_faces;
  // }

  template <typename Rng>
  void gaussian(Rng& rng, const double width=1.0) {
    for(int il=0; il<lattice.n_links; il++) field[il] = width*rng.gaussian_link(il);
  }

};




  // double average_plaqsq() const {
  //   double sum = 0.0;
  //   for(int i=0; i<lattice.n_faces; i++) {
  //     double val = plaquette_angle(i);
  //     while(val>M_PI) val -= 2.0*M_PI;
  //     while(val<-M_PI) val += 2.0*M_PI;
  //     const double factor = 1.0 / lattice.vps[i];
  //     sum += factor * std::pow( val, 2 );
  //   }
  //   return sum/lattice.n_faces;
  // }
