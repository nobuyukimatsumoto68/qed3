#pragma once

#include <Eigen/Dense>


// template<class Lattice>
class DiracBase {
public:
  // using MS=Eigen::Matrix2cd;
  // using VD=Eigen::Vector2d;
  // using VE=Eigen::Vector3d;
  // using VC=Eigen::VectorXcd;
  // using Link = std::array<Idx,2>; // <int,int>;
  // Lattice& lattice;

  std::array<MS, 4> sigma;
  // std::vector<double> kappa; // evan's link label

  static constexpr int NS = 2;
  static constexpr int DIM = 2;
  static constexpr Complex I = Complex(0.0, 1.0);

  // const double m;
  // const double r;
  // const double M5;

  DiracBase()
  {
    set_sigma();
    // set_kappa();
  }

  void set_sigma(){
    assert(NS==2);
    sigma[0] << 1,0,0,1;
    sigma[1] << 0,1,1,0;
    sigma[2] << 0,-I,I,0;
    sigma[3] << 1,0,0,-1;
  }
  // virtual void set_kappa() = 0;

};
