#pragma once

template<typename Lattice, int Nt_, bool is_compact_=false>
struct GaugeExt {
  using BaseLink = std::array<int,2>; // <int,int>;
  using BaseFace = std::vector<int>;

  static constexpr bool is_compact = is_compact_;
  static constexpr int Nt = Nt_;

  Lattice& lattice;

  std::vector<std::vector<double>> spatial; // [t][il]
  std::vector<std::vector<double>> temporal; // [t][il]

  using Gauge = GaugeExt;

  GaugeExt()=delete;

  GaugeExt(Lattice& lattice_)
    : lattice(lattice_)
    , spatial(Nt,
              std::vector<double>(lattice.n_links, 0.0))
    , temporal(Nt,
               std::vector<double>(lattice.n_sites, 0.0))
  {}

  GaugeExt( const GaugeExt& other )
    : lattice(other.lattice)
    , spatial(other.spatial)
    , temporal(other.temporal)
  {}


  GaugeExt& operator=(const GaugeExt& other){
    if (this == &other) return *this;

    assert(&lattice==&other.lattice);
    spatial = other.spatial;
    temporal = other.temporal;
    return *this;
  }


  double sp(const int s, const BaseLink& ell) const { // recommended
    const int il = lattice.map2il.at(ell);
    const int sign = lattice.map2sign.at(ell);
    return sign * spatial[(s+Nt)%Nt][il];
  }

  inline double sp(const int s, const Idx& il) const { // recommended
    return spatial[(s+Nt)%Nt][il];
  }


  inline double tp(const int s, const Idx& ix) const { // recommended
    return temporal[(s+Nt)%Nt][ix];
  }

  double& sp(const int s, const BaseLink& ell) { // recommended
    const int il = lattice.map2il.at(ell);
    const int sign = lattice.map2sign.at(ell);
    assert(sign==1);
    return spatial[(s+Nt)%Nt][il];
  }

  inline double& sp(const int s, const Idx& il) { // recommended
    return spatial[(s+Nt)%Nt][il];
  }

  inline double& tp(const int s, const Idx& ix) { // recommended
    return temporal[(s+Nt)%Nt][ix];
  }


  template <typename Force>
  GaugeExt& operator+=(const Force& rhs){
    for(Idx i=0; i<spatial.size(); i++) for(Idx j=0; j<spatial[i].size(); j++) spatial[i][j] += rhs.spatial[i][j];
    for(Idx i=0; i<temporal.size(); i++) for(Idx j=0; j<temporal[i].size(); j++) temporal[i][j] += rhs.temporal[i][j];
    return *this;
  }

  template <typename Force>
  friend Force operator+(GaugeExt v, const Force& w) { v += w; return v; }

  GaugeExt& operator*=(const double rhs){
    for(Idx i=0; i<spatial.size(); i++) for(Idx j=0; j<spatial[i].size(); j++) spatial[i][j] *= rhs;
    for(Idx i=0; i<temporal.size(); i++) for(Idx j=0; j<temporal[i].size(); j++) temporal[i][j] *= rhs;
    return *this;
  }

  GaugeExt& operator/=(const double rhs){
    for(Idx i=0; i<spatial.size(); i++) for(Idx j=0; j<spatial[i].size(); j++) spatial[i][j] /= rhs;
    for(Idx i=0; i<temporal.size(); i++) for(Idx j=0; j<temporal[i].size(); j++) temporal[i][j] /= rhs;
    return *this;
  }

  friend GaugeExt operator*(GaugeExt v, const double a) {
    v *= a;
    return v;
  }

  friend GaugeExt operator*(const double a, GaugeExt v) {
    v *= a;
    return v;
  }


  double squared_norm() const {
    double sum = 0.0;
    for(Idx i=0; i<spatial.size(); i++) for(Idx j=0; j<spatial[i].size(); j++) sum += spatial[i][j]*spatial[i][j];
    for(Idx i=0; i<temporal.size(); i++) for(Idx j=0; j<temporal[i].size(); j++) sum += temporal[i][j]*temporal[i][j];
    return sum;
  }

  inline double norm() const {
    return std::sqrt( squared_norm() );
  }

  void print2log_norm(const std::string comment=" : ") const {
    std::clog << comment << norm() << std::endl;
  }

  double plaquette_angle(const int s, const BaseFace& face) const {
    double sum = 0.0;
    for(Idx i=0; i<face.size(); i++) {
      const Idx ix = face[i];
      const Idx iy = face[(i+1)%face.size()];
      sum += sp( s, BaseLink{ix,iy} );
    }
    return sum;
  }

  double plaquette_angle(const int s, const BaseLink& link) const {
    double sum = 0.0;

    sum += sp( s, link );
    sum += tp( s, link[1] );
    sum -= sp( s+1, link );
    sum -= tp( s, link[0] );

    return sum;
  }

  template <typename Rng>
  void gaussian(Rng& rng, const double width=1.0) {
    for(Idx i=0; i<spatial.size(); i++) for(Idx j=0; j<spatial[i].size(); j++) spatial[i][j] = width*rng.gaussian_link(i,j);
    for(Idx i=0; i<temporal.size(); i++) for(Idx j=0; j<temporal[i].size(); j++) temporal[i][j] = width*rng.gaussian_site(i,j);
  }

  // double Mod(double a, double b=2.0*M_PI){
  //   a += 100*b;
  //   // int p = int(std::floor(a / b));
  //   // double r = a - p*b;
  //   return fmod(a, b);
  // }


  // void project() {
  //   for(Idx i=0; i<spatial.size(); i++) for(Idx j=0; j<spatial[i].size(); j++) spatial[i][j] = Mod( spatial[i][j], 2.0*M_PI );
  //   for(Idx i=0; i<temporal.size(); i++) for(Idx j=0; j<temporal[i].size(); j++) temporal[i][j] = Mod( temporal[i][j], 2.0*M_PI );
  // }


  template <typename Rng>
  void random_gauge_trsf(Rng& rng, const double width=1.0) {

    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<lattice.n_sites; ix++){
        const double phi = width*rng.gaussian_site(s,ix);
        for(const Idx iy : lattice.nns[ix]){
          const Idx il = lattice.map2il.at(BaseLink{ix,iy});
          const int sign = lattice.map2sign.at(BaseLink{ix,iy});
          sp(s,il) += sign*phi;
        }
        tp(s,ix) += phi;
        tp(s-1,ix) -= phi;
      }
    }

    // project();
  }


  void gauge_trsf(const FermionVector& gauge) {

    for(int s=0; s<Nt; s++){
      for(Idx ix=0; ix<lattice.n_sites; ix++){
        const double phi = gauge(s,ix,0).real();
        for(const Idx iy : lattice.nns[ix]){
          const Idx il = lattice.map2il.at(BaseLink{ix,iy});
          const int sign = lattice.map2sign.at(BaseLink{ix,iy});
          sp(s,il) += sign*phi;
        }
        tp(s,ix) += phi;
        tp(s-1,ix) -= phi;
      }
    }

    // project();
  }


  template<typename Gauge, typename ComplexType, typename Fermion>
  void compute( const Gauge& u, const ComplexType* d_eta, const Fermion& D ){

// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL_GRAD) schedule(static)
// #endif
    for(int s=0; s<Nt; s++)
      for(Idx ell=0; ell<lattice.n_links; ell++)
        sp(s,ell) = D.grad_deviceAsyncLaunch( std::pair<int, BaseLink>(s,lattice.links[ell]),
                                              u, d_eta );

// #ifdef _OPENMP
// #pragma omp parallel for num_threads(Comp::NPARALLEL_GRAD) schedule(static)
// #endif
    for(int s=0; s<Nt; s++)
      for(Idx ix=0; ix<lattice.n_sites; ix++)
        tp(s,ix) = D.grad_deviceAsyncLaunch( std::pair<int, Idx>(s,ix),
                                             u, d_eta );
  }


  void ckpoint( const std::string& str ) const {
    std::ofstream of( str, std::ios::out | std::ios::binary | std::ios::trunc);
    if(!of) assert(false);
    double tmp = 0.0;

    for(Idx i=0; i<spatial.size(); i++){
      for(Idx j=0; j<spatial[i].size(); j++){
        // std::cout << "debug. i,j = " << i << j << std::endl;
        tmp = spatial[i][j];
        of.write( (char*) &tmp, sizeof(double) );
      }
    }
    for(Idx i=0; i<temporal.size(); i++){
      for(Idx j=0; j<temporal[i].size(); j++){
        // std::cout << "debug. i,j = " << i << j << std::endl;
        tmp = temporal[i][j];
        of.write( (char*) &tmp, sizeof(double) );
      }
    }

    of.close();
  }

  void read( const std::string& str ) {
    std::ifstream ifs( str, std::ios::in | std::ios::binary );
    if(!ifs) assert(false);
    double tmp;

    // for(Idx i=0; i<Lx*Ly; ++i){
    //   ifs.read((char*) &tmp, sizeof(double) );
    //   (*this)[i] = tmp;
    // }
    for(Idx i=0; i<spatial.size(); i++){
      for(Idx j=0; j<spatial[i].size(); j++){
        // std::cout << "debug. i,j = " << i << j << std::endl;
        ifs.read( (char*) &tmp, sizeof(double) );
        // std::cout << tmp << std::endl;
        spatial[i][j] = tmp;
      }
    }
    for(Idx i=0; i<temporal.size(); i++){
      for(Idx j=0; j<temporal[i].size(); j++){
        ifs.read( (char*) &tmp, sizeof(double) );
        // std::cout << tmp << std::endl;
        temporal[i][j] = tmp;
      }
    }

    ifs.close();
  }

};
