#pragma once

#include <vector>
#include <cmath>

#include "cg_cuda.h"


struct PseudoFermion {
  using Link = std::array<int,2>; // <int,int>;
  using Face = std::vector<int>;

  using Complex = std::complex<double>;
  static constexpr Complex I = Complex(0.0, 1.0);

  using MS=Eigen::Matrix2cd;
  using VD=Eigen::Vector2d;
  using VE=Eigen::Vector3d;
  using VC=Eigen::VectorXcd;

  // const Lattice& lattice;
  const Dirac1fonS2& D;
  const CGCUDA cg;

  std::vector<VD> eta;
  std::vector<VD> phi;

  PseudoFermion()=delete;

  PseudoFermion(const Dirac1fonS2& D_)
    : D(D_)
    , cg(D)
    , eta(D.lattice.n_sites, VD::Zero())
    , phi(D.lattice.n_sites, VD::Zero())
  {}

  PseudoFermion( const PseudoFermion& other )
    : D(other.D)
    , cg(cg)
    , eta(other.eta)
    , phi(other.phi)
  {}

  // PseudoFermion & operator=(const PseudoFermion&) = delete;
  // PseudoFermion& operator=(const PseudoFermion& other){
  //   if (this == &other) return *this;
  //   assert(&D==&other.D);
  //   eta = other.eta;
  //   phi = other.phi;
  //   return *this;
  // }

  // double operator[](const int i) const { return field[i]; }
  // double& operator[](const int i) { return field[i]; }

  // auto begin(){ return field.begin(); }
  // auto end(){ return field.end(); }
  // auto begin() const { return field.begin(); }
  // auto end() const { return field.end(); }

  








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


