#pragma once

#include <Eigen/Dense>


// template<class Lattice>
class DiracBase {
public:
  std::array<MS, 4> sigma;

  static constexpr int NS = 2;
  static constexpr int DIM = 2;
  static constexpr Complex I = Complex(0.0, 1.0);

  DiracBase()
  {
    set_sigma();
  }

  void set_sigma(){
    assert(NS==2);
    sigma[0] << 1,0,0,1;
    sigma[1] << 0,1,1,0;
    sigma[2] << 0,-I,I,0;
    sigma[3] << 1,0,0,-1;
  }

};
