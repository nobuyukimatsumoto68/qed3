// three_tri_check_claude.cu
// Validate the three-triangle shape geometry: instance count, orbit count/sizes at L=1,2,4,
// and that each instance's two shared links cancel (boundary loop = 6 distinct outer links).
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <set>
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
double Ylm_real(const int ell, const int em, const VE r){ return 0.0; } // stub (op() unused here)
#include "icos_orbits_claude.h"
#include "wilson_shapes_claude.h"

int main(){
  for(const int L : {1,2,4}){
    S2Simp base(L);
    IcosOrbits<S2Simp> orb(base);
    WilsonShapes<S2Simp> shp(base, orb);
    auto cand = shp.three_triangles();
    auto orbs = shp.orbits_from(cand);
    // check link cancellation on the first instance: 3 faces, 9 directed edges, 2 shared pairs cancel
    const auto& in0 = cand[0];
    std::map<std::array<Idx,2>,int> edge;
    for(const Idx f : in0.faces){
      const auto& fc = base.faces[f];
      for(int e=0;e<3;e++){
        Idx a=fc[e], b=fc[(e+1)%3];
        std::array<Idx,2> k{std::min(a,b),std::max(a,b)};
        edge[k]++;
      }
    }
    int shared=0, outer=0;
    for(auto& kv: edge){ if(kv.second==2) shared++; else outer++; }
    std::cout << "L=" << L << "  n_faces=" << base.n_faces
              << "  three_tri instances=" << cand.size()
              << "  orbits=" << orbs.size() << "  [";
    for(size_t o=0;o<orbs.size();o++) std::cout << orbs[o].size() << (o+1<orbs.size()?",":"");
    std::cout << "]  first-inst shared-links=" << shared << " outer-links=" << outer << std::endl;
  }
  return 0;
}
