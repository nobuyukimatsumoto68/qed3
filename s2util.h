#pragma once

#include "s2.h"



struct LinkInfo {
  QfeLatticeS2& lattice;

  using Link = std::array<int,2>; // <int,int>;
  const std::map<const Link, const int> map2il;
  const std::map<const Link, const int> map2sign;

  LinkInfo()=delete;

  LinkInfo(QfeLatticeS2& lattice_)
    : lattice(lattice_)
    , map2il(get_map2il())
    , map2sign(get_map2sign())
  {
  }

  LinkInfo(const LinkInfo& other)
    : lattice(other.lattice)
    , map2il(other.map2il)
    , map2sign(other.map2sign)
  {}

  LinkInfo & operator=(const LinkInfo&) = delete;

  std::map<const Link, const int> get_map2il() const {
    std::map<const Link, const int> map2il_;

    for(int il=0; il<lattice.n_links; il++) {
      const auto link = lattice.links[il];
      const int i = link.sites[0];
      const int j = link.sites[1];

      map2il_.insert( { Link{i,j}, il } );
      map2il_.insert( { Link{j,i}, il } );
    }

    return map2il_;
  }

  std::map<const Link, const int> get_map2sign() const {
    std::map<const Link, const int> map2sign_;

    for(int il=0; il<lattice.n_links; il++) {
      const auto link = lattice.links[il];
      const int i = link.sites[0];
      const int j = link.sites[1];

      map2sign_.insert( { Link{i,j}, +1 } );
      map2sign_.insert( { Link{j,i}, -1 } );
    }

    return map2sign_;
  }

};





struct FaceInfo {
  QfeLatticeS2& lattice;

  std::vector<int> signs; // index: ia (Evan's label for faces)
  std::vector<double> vps; // evan's face label

  int Np;

  FaceInfo()=delete;

  FaceInfo(QfeLatticeS2& lattice_)
    : lattice(lattice_)
  {
    set_signs();
    set_vps();
    assert( vps.size()==signs.size() );
    Np = vps.size();
  }

  FaceInfo(const FaceInfo& other)
    : lattice(other.lattice)
    , signs(other.signs)
    , vps(other.vps)
    , Np(other.Np)
  {
    assert( vps.size()==signs.size() );
  }

  FaceInfo & operator=(const FaceInfo&) = delete;

  int sign(const int i_face) const {
    const QfeFace& face = lattice.faces[i_face];
    const Vec3 r0 = lattice.r[face.sites[0]];
    const Vec3 r1 = lattice.r[face.sites[1]];
    const Vec3 r2 = lattice.r[face.sites[2]];

    const Vec3 cross = (r1-r0).cross(r2-r0);
    const Vec3 sum = r0+r1+r2;

    const double inner = cross.dot(sum);

    int res = 1;
    if(inner<0) res = -1;
    return res;
  }

  void set_signs() { for(int i=0; i<lattice.n_faces; i++) signs.push_back( sign(i) ); }

  double vp(const int i_face) const {
    const QfeFace& face = lattice.faces[i_face];
    const Vec3 r0 = lattice.r[face.sites[0]];
    const Vec3 r1 = lattice.r[face.sites[1]];
    const Vec3 r2 = lattice.r[face.sites[2]];

    // const Vec3 p = circumcenter(r0, r1, r2);
    // assert( std::abs( (p-r0).norm() - (p-r1).norm() )<1.0e-14 );
    // assert( std::abs( (p-r0).norm() - (p-r2).norm() )<1.0e-14 );
    // const double a = (r0-p).norm();
    // const double b = (r1-p).norm();
    // const double c = (r0-r1).norm();
    // const double s = 0.5*(a+b+c);
    // const double tmp = s * (s-a) * (s-b) * (s-c);
    // return std::sqrt( tmp );

    const double a = std::acos(r0.dot(r1));
    const double b = std::acos(r1.dot(r2));
    const double c = std::acos(r2.dot(r0));
    const double s = 0.5*(a+b+c);

    double tantan = std::tan(0.5*s);
    tantan *= std::tan(0.5*(s-a));
    tantan *= std::tan(0.5*(s-b));
    tantan *= std::tan(0.5*(s-c));
    return 4.0 * std::atan( std::sqrt(tantan) );
  }

  void set_vps() {
    for(int i=0; i<lattice.n_faces; i++) vps.push_back( vp(i) );

    double sum = 0.0;
    for(auto elem : vps) sum += elem;
    assert( std::abs(sum-4.0*M_PI)<1.0e-10 );
  }

  Vec3 circumcenter(const Vec3& r0, const Vec3& r1, const Vec3& r2) const {
    const Vec3 r10 = r1 - r0;
    const Vec3 r20 = r2 - r0;

    const Vec3 tmp1 = r10.squaredNorm() * r20 - r20.squaredNorm() * r10;
    const Vec3 cross = r10.cross(r20);
    const Vec3 numer = tmp1.cross(cross);
    const double denom = 2.0*cross.squaredNorm();

    return numer/denom + r0;
  }



};
