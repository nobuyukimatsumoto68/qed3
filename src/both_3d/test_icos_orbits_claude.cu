// test_icos_orbits_claude.cu
// Build the native Ih orbit table on the qed3 s2n_simp lattice and verify:
//   (i) the group closes to |Ih|=120,
//   (ii) every element permutes the lattice sites (max_match_dist ~ 0),
//   (iii) identity maps each site to itself,
//   (iv) orbit-size histogram is sensible.
// Host-only (g++ + Eigen). Run from src/both_3d/.  Usage: ./a.out [L=1]

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <cstdint>
#include <Eigen/Dense>

using Idx = std::int32_t;
using VE = Eigen::Vector3d;

namespace Comp { constexpr int NPARALLEL = 1; }

const std::string dir = "../../geometry/data/";

#include "s2n_simp.h"
#include "icos_orbits_claude.h"

int main(int argc, char* argv[]){
  const int L = (argc>1) ? std::atoi(argv[1]) : 1;

  S2Simp base(L);
  IcosOrbits<S2Simp> orb(base);

  std::cout << "===== icos orbit test L=" << L << " =====" << std::endl;
  std::cout << "n_sites=" << base.n_sites
            << "  |group|=" << orb.group.size()
            << "  max_match_dist=" << orb.max_match_dist << std::endl;

  bool id_ok=true;
  for(Idx i=0; i<base.n_sites; i++) if(orb.orbits[i][0]!=i) id_ok=false;
  std::cout << "identity maps sites to themselves: " << (id_ok?"YES":"NO") << std::endl;

  std::map<Idx,int> hist;
  for(Idx i=0; i<base.n_sites; i++) hist[ orb.orbit_size(i) ]++;
  std::cout << "orbit-size histogram (distinct images : #sites):" << std::endl;
  for(const auto& kv : hist) std::cout << "  " << kv.first << " : " << kv.second << std::endl;

  return 0;
}
