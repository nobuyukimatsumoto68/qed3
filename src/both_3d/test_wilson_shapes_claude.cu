// test_wilson_shapes_claude.cu
// Build the triangle (face) and rectangle (link) shape orbits on the qed3 lattice and report
// orbit counts / sizes, plus the face_signs distribution (relevant to the rectangle 2-face sum).
// Host-only (g++ + Eigen). Run from src/both_3d/.  Usage: ./a.out [L=1]

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
double Ylm_real(const int ell, const int em, const VE r);  // decl only (op methods unused in this test)
#include "icos_orbits_claude.h"
#include "wilson_shapes_claude.h"

int main(int argc, char* argv[]){
  const int L = (argc>1) ? std::atoi(argv[1]) : 1;

  S2Simp base(L);
  IcosOrbits<S2Simp> orb(base);
  WilsonShapes<S2Simp> shp(base, orb);

  std::cout << "===== wilson shapes L=" << L << " =====" << std::endl;
  std::cout << "n_faces=" << base.n_faces << "  n_links=" << base.n_links << std::endl;

  auto report = [&](const char* name, const std::vector<std::vector<WilsonShapes<S2Simp>::Instance>>& orbs){
    std::cout << name << " orbits: " << orbs.size() << "  [";
    for(size_t o=0; o<orbs.size(); o++) std::cout << orbs[o].size() << (o+1<orbs.size()?",":"");
    std::cout << "]" << std::endl;
  };
  report("triangle       ", shp.orbits_from( shp.triangles() ));
  report("rectangle      ", shp.orbits_from( shp.rectangles() ));
  report("twisted-rect   ", shp.orbits_from( shp.twisted_rectangles() ));
  report("figure-8       ", shp.orbits_from( shp.figure8s() ));
  report("twisted-fig-8  ", shp.orbits_from( shp.twisted_figure8s() ));

  int np=0, nm=0;
  for(const int s : base.face_signs){ if(s>0) np++; else nm++; }
  std::cout << "face_signs: +1 x " << np << " , -1 x " << nm << std::endl;

  return 0;
}
