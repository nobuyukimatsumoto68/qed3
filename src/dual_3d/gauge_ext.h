#pragma once

template<typename BaseLattice, bool is_compact=false>
struct GaugeExt {
  using Link = std::array<int,2>; // <int,int>;
  using Face = std::vector<int>;
  using Lattice = S2Trivalent;

  Lattice& lattice;

  std::vector<std::vector<double>> spatial;
  std::vector<std::vector<double>> temporal;

  using Gauge = GaugeExt;

  Gauge()=delete;

  Gauge(Lattice& lattice_,
        const int Nt_)
    : lattice(lattice_)
    , spatial(lattice.Nt,
              std::vector<double>(lattice.base.n_links, 0.0))
    , temporal(lattice.Nt,
               std::vector<double>(lattice.base.n_sites, 0.0))
  {}

  Gauge( const Gauge& other )
    : lattice(other.lattice)
    , spatial(other.spatial)
    , temporal(other.temporal)
  {}


  Gauge& operator=(const Gauge& other){
    if (this == &other) return *this;

    assert(&lattice==&other.lattice);
    spatial = other.spatial;
    temporal = other.temporal;
    return *this;
  }


  double operator()(const int s, const Link& ell) const { // recommended
    const int il = lattice.base.map2il.at(ell);
    const int sign = lattice.base.map2sign.at(ell);
    return sign * spatial[s][il];
  }


  inline double operator()(const int s, const Idx& ix) const { // recommended
    return temporal[s][ix];
  }


  template <typename Force>
  Gauge& operator+=(const Force& rhs){
    for(Idx i=0; i<spatial.size(); i++) spatial[i] += rhs.spatial[i];
    for(Idx i=0; i<temporal.size(); i++) temporal[i] += rhs.temporal[i];
    return *this;
  }
  template <typename Force>
  friend Force operator+(Gauge v, const Force& w) { v += w; return v; }

  Gauge& operator*=(const double rhs){
    for(Idx i=0; i<spatial.size(); i++) spatial[i] *= rhs;
    for(Idx i=0; i<temporal.size(); i++) temporal[i] *= rhs;
    return *this;
  }

  Gauge& operator/=(const double rhs){
    for(Idx i=0; i<spatial.size(); i++) spatial[i] /= rhs;
    for(Idx i=0; i<temporal.size(); i++) temporal[i] /= rhs;
    return *this;
  }

  friend Gauge operator*(Gauge v, const double a) {
    for(Idx i=0; i<spatial.size(); i++) v.spatial[i] *= rhs;
    for(Idx i=0; i<temporal.size(); i++) v.temporal[i] *= rhs;
    return v;
  }

  friend Gauge operator*(const double a, Gauge v) {
    for(Idx i=0; i<spatial.size(); i++) v.spatial[i] *= rhs;
    for(Idx i=0; i<temporal.size(); i++) v.temporal[i] *= rhs;
    return v;
  }


  double plaquette_angle(const int s, const Face& face) const {
    double sum = 0.0;
    for(Idx i=0; i<face.size(); i++) {
      const Idx ix = face[i];
      const Idx iy = face[(i+1)%face.size()];
      sum += (*this)( s, Link{ix,iy} );
    }
    return sum;
  }

  double plaquette_angle(const int s, const Link& link) const {
    double sum = 0.0;

    sum += (*this)( s, link );
    sum += (*this)( s, link[1] );
    sum -= (*this)( s+1, link );
    sum -= (*this)( s, link[0] );

    return sum;
  }





};
