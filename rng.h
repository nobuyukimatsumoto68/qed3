#pragma once

#include <random>

struct Rng {
  std::mt19937_64 mt;

  std::normal_distribution<double> dist_gaussian;
  double gaussian(){ return dist_gaussian(mt); }

  std::uniform_real_distribution<> dist_01;
  double uniform(){ return dist_01(mt); }

  double RandReal(const double a=0., const double b=1.0) {
    return (b-a) * uniform() + a;
  }

  void SeedRng(const int seed) { mt.seed(seed); }
};
