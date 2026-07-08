// face_sign_check_claude.cu
// Decisive check: are the stored face orientations mixed (so face_signs matters) or all +1
// (so the shapes triangle == the density plaquette operator)? Counts +1 vs -1 at L=1,2,4.
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <map>
#include <set>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <cstdint>
#include <string>
#include <Eigen/Dense>

using Idx = std::int32_t;
using VE = Eigen::Vector3d;

namespace Comp { constexpr int NPARALLEL = 1; }

const std::string dir = "../../geometry/data/";

#include "s2n_simp.h"

int main(int argc, char* argv[]){
  for(const int L : {1,2,4}){
    S2Simp base(L);
    int nplus = 0;
    int nminus = 0;
    for(const int s : base.face_signs){
      if(s>0) nplus++;
      else nminus++;
    }
    std::cout << "L=" << L << "  n_faces=" << base.n_faces
              << "  face_signs: +1=" << nplus << "  -1=" << nminus << std::endl;
  }
  return 0;
}
