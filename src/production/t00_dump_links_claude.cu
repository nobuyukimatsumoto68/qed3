// t00_dump_links_claude.cu
// Dump the PRIMAL spatial link table of the S2Simp lattice that the distillation V / DiracS2Simp are built on,
// so the offline python T_00 vertex can index a gauge config's spatial-link phase sp[t,il] with the EXACT
// runtime ordering (the dats/*links*.dat are DUAL links; the primal il ordering is runtime-built).
//
// Config-INDEPENDENT (pure geometry): for each DIRECTED spatial link (ix,iy) it writes
//     ix  iy  il  sign  kappa
// where il = base.map2il[{ix,iy}], sign = base.map2sign[{ix,iy}] (canonical-orientation sign), kappa=D2.kappa[il].
// The offline vertex is then  W(t)_{ix,iy} = 0.5 kappa[il] (e^a sigma_a)(ix,iy) Omega(ix,iy) exp(i*sign*sp[t,il])
// (r=0 naive e.sigma energy density); gamma/Omega come from omega_n<L>.dat / alpha_n<L>.dat (same files
// SpinStructureSimp reads, so they match D_ov exactly).  u = sign*sp matches gauge_ext.h U.sp(t,ell).
//
// Build/run: bash tmp_claude.sh  (reads log t00_dump_links_claude.log, output primal_links_n<L>_claude.dat).

#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <Eigen/Dense>

using Idx = std::int32_t;
using Double = double;
using Complex = std::complex<double>;
using Link = std::array<Idx, 2>;
using MS = Eigen::Matrix2cd;    // 2x2 complex (spinor)
using VE = Eigen::Vector3d;     // 3D Euclidean vector (sites)
using VC = Eigen::VectorXcd;    // complex vectors

const std::string dir = "../../geometry/data/";

#ifndef N_REFINE_CLI
#define N_REFINE_CLI 1
#endif
namespace Comp {
  constexpr int N_REFINE = N_REFINE_CLI;
  constexpr int NS = 2;
}

#include "s2n_simp.h"

int main() {
  S2Simp base(Comp::N_REFINE);
  // kappa[il] = 2 link_volume[il] / ell[il] / mean_ell   (matches DiracS2Simp::set_kappa, dirac_simp.h:360)

  const std::string out_name = "primal_links_n" + std::to_string(Comp::N_REFINE) + "_claude.dat";
  std::ofstream out(out_name);
  out << "# ix iy il sign kappa   (directed spatial links; u_xy(t)=sign*theta_sp[t,il], U=exp(i u))\n";
  out << "# n_sites=" << base.n_sites << " n_links=" << base.n_links << "\n";
  out << std::setprecision(17);
  Idx ndir = 0;
  for (Idx ix = 0; ix < base.n_sites; ix++) {
    for (Idx iy : base.nns[ix]) {
      const Idx il = base.map2il.at(Link{ix, iy});
      const int sgn = base.map2sign.at(Link{ix, iy});
      const double kap = 2.0 * base.link_volume[il] / base.ell[il] / base.mean_ell;
      out << ix << " " << iy << " " << il << " " << sgn << " " << kap << "\n";
      ndir++;
    }
  }
  out.close();
  // per-site dual areas + mean_ell (kappa_t[ix] = dual_areas[ix]/mean_ell/at, dirac_ext.h:452) for the exact free D_W
  const std::string out2_name = "primal_sites_n" + std::to_string(Comp::N_REFINE) + "_claude.dat";
  std::ofstream out2(out2_name);
  out2 << "# ix dual_area   (mean_ell on the next header line)\n";
  out2 << "# mean_ell=" << std::setprecision(17) << base.mean_ell << "\n";
  for (Idx ix = 0; ix < base.n_sites; ix++) {
    out2 << ix << " " << std::setprecision(17) << base.dual_areas[ix] << "\n";
  }
  out2.close();
  std::cout << "# wrote "
 << out_name << "  (n_sites=" << base.n_sites
            << " n_links=" << base.n_links << " directed=" << ndir << ")\n";
  return 0;
}
