// Standalone geometry tool: print the orbit-index -> shape-type map for the 5-shape basis at a given
// refinement L. Reuses the CUDA-free lattice (S2Simp) + IcosOrbits + WilsonShapes headers; the shape
// ENUMERATION needs only the geometry (no gauge). Build: g++ -O2 -std=c++17 -Iincludes -I<eigen>.
//   ./shape_orbit_map_claude.o <L>
#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <string>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <cassert>
#include <Eigen/Dense>

using Idx = std::int32_t;
using VE = Eigen::Vector3d;
using Link = std::array<Idx,2>;
using Face = std::vector<Idx>;

const std::string dir = "../../geometry/data/";   // s2n_simp.h reads geometry from here

#include "s2n_simp.h"
#include "icos_orbits_claude.h"
#include "wilson_shapes_claude.h"

int main(int argc, char** argv)
{
  const int n_refine = (argc>1) ? std::atoi(argv[1]) : 2;
  using Base = S2Simp;
  Base base(n_refine);
  IcosOrbits<Base> orb(base);
  WilsonShapes<Base> shp(base, orb);
  using Inst = typename WilsonShapes<Base>::Instance;

  std::vector<std::vector<Inst>> all[7] = {
    shp.orbits_from( shp.triangles() ),
    shp.orbits_from( shp.rectangles() ),
    shp.orbits_from( shp.figure8s() ),
    shp.orbits_from( shp.three_triangles() ),
    shp.orbits_from( shp.four_triangles() ),
    shp.orbits_from( shp.trios() ),
    shp.orbits_from( shp.site_contours() ),
  };
  const char* SHNAME[7] = {"triangle","rect","fig8","three-tri","star","trio","five-six"};
  int off = 0;
  std::cout << "L=" << n_refine << "  orbit-index -> shape-type map:\n";
  for(int is=0; is<7; is++){
    const int n = (int)all[is].size();
    std::cout << "  " << SHNAME[is] << ": orbits [" << off << ".." << off+n-1 << "]  (" << n << " orbits)\n";
    off += n;
  }
  std::cout << "  total n_orbits = " << off << "\n";
  return 0;
}
