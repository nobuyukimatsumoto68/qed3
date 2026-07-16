#pragma once

// icos_orbits_claude.h
//
// Native full-icosahedral-group (Ih, |Ih|=120) orbit table on the qed3 s2n_simp lattice.
// Reuses only the MATH from s2ising/groups.h: build the C2 (edge) and C3 (face) generators
// from the base icosahedron, close the group under matrix multiplication (NO multiplication-
// table file needed), and apply the 120 rotation matrices to the lattice's own `base.sites`
// (unit VE) -> orbits[ix][ig] = the site index that ix maps to under group element ig, all in
// the qed3 labelling. No s2ising link / no cross-basis index map required.
//
// A Wilson loop (ordered site sequence) maps under g by  site -> orbits[site][ig]; its icos
// orbit is the dedup of the 120 images. Ref: s2ising/groups.h (Rotation, Orbits).
//
// Assumes VE = Eigen::Vector3d and Idx are already defined (as in s2n_simp.h). Include AFTER
// s2n_simp.h.

#include <vector>
#include <array>
#include <cmath>
#include <cassert>
#include <Eigen/Dense>

template<class Base>
struct IcosOrbits {
  static constexpr int NIh = 120;

  const Base& base;
  std::vector<Eigen::Matrix3d> group;      // 120 rotation matrices (SO(3) + inversion)
  std::vector<std::vector<Idx>> orbits;    // orbits[ix][ig] = site index of group[ig]*sites[ix]
  double max_match_dist = 0.0;             // largest |g x - nearest site| (permutation-check residual)

  IcosOrbits(const Base& base_)
    : base(base_)
  {
    build_group();
    build_orbits();
  }

  // Rodrigues rotation about a (not-necessarily-unit) axis by `angle`.
  static Eigen::Matrix3d rot_axis_angle(const VE& axis, const double angle){
    const VE a = axis.normalized();
    Eigen::Matrix3d K;
    K <<   0.0,   -a.z(),  a.y(),
         a.z(),     0.0,  -a.x(),
        -a.y(),   a.x(),    0.0;
    return Eigen::Matrix3d::Identity() + std::sin(angle)*K + (1.0-std::cos(angle))*(K*K);
  }

  void build_group(){
    // base icosahedron vertices = the 12 sites with 5 neighbours (all 12 sites at L=1).
    std::vector<Idx> bv;
    for(Idx i=0; i<base.n_sites; i++){
      if((int)base.nns[i].size()==5) bv.push_back(i);
    }
    assert(bv.size()==12 && "expected 12 five-fold base icosahedron vertices");

    // C2 axis: midpoint direction of the nearest base-vertex pair (an icosahedron edge).
    double best = -2.0;
    Idx ei=bv[0], ej=bv[1];
    for(size_t a=0; a<bv.size(); a++){
      for(size_t b=a+1; b<bv.size(); b++){
        const double d = base.sites[bv[a]].dot(base.sites[bv[b]]);
        if(d>best){ best=d; ei=bv[a]; ej=bv[b]; }
      }
    }
    const Eigen::Matrix3d ms = rot_axis_angle( base.sites[ei]+base.sites[ej], M_PI );

    // C3 axis: face-centre direction (edge ei-ej plus the base vertex nearest to both).
    double best2 = -4.0;
    Idx ek=bv[0];
    for(const Idx k : bv){
      if(k==ei || k==ej) continue;
      const double d = base.sites[k].dot(base.sites[ei]) + base.sites[k].dot(base.sites[ej]);
      if(d>best2){ best2=d; ek=k; }
    }
    const Eigen::Matrix3d mt = rot_axis_angle( base.sites[ei]+base.sites[ej]+base.sites[ek], 2.0*M_PI/3.0 );

    const Eigen::Matrix3d minv = -Eigen::Matrix3d::Identity();

    // closure: BFS from identity under the generators {ms, mt, -1}
    const std::array<Eigen::Matrix3d,3> gens = { ms, mt, minv };
    group.clear();
    group.push_back( Eigen::Matrix3d::Identity() );
    size_t head=0;
    while(head < group.size()){
      const Eigen::Matrix3d g = group[head++];
      for(const Eigen::Matrix3d& s : gens){
        const Eigen::Matrix3d gs = s*g;
        bool found=false;
        for(const Eigen::Matrix3d& h : group){
          if((h-gs).norm()<1.0e-9){ found=true; break; }
        }
        if(!found){
          group.push_back(gs);
          assert((int)group.size()<=NIh && "group closure exceeded 120 -- bad generators");
        }
      }
    }
    assert((int)group.size()==NIh && "group closure did not reach 120");
  }

  void build_orbits(){
    orbits.assign(base.n_sites, std::vector<Idx>(NIh, -1));
    max_match_dist = 0.0;
    for(Idx ix=0; ix<base.n_sites; ix++){
      const VE x = base.sites[ix];
      for(int ig=0; ig<NIh; ig++){
        const VE gx = group[ig]*x;
        Idx best=-1;
        double bestd=1.0e18;
        for(Idx iy=0; iy<base.n_sites; iy++){
          const double d = (gx-base.sites[iy]).squaredNorm();
          if(d<bestd){ bestd=d; best=iy; }
        }
        if(std::sqrt(bestd) > max_match_dist) max_match_dist = std::sqrt(bestd);
        orbits[ix][ig] = best;
      }
    }
    // permutation property: every g*site must coincide with a lattice site
    assert(max_match_dist < 1.0e-8 && "Ih element does not permute the sites (orientation mismatch?)");
  }

  const std::vector<Idx>& operator[](const Idx ix) const { return orbits[ix]; }

  // distinct-image count of a site's orbit (= |Ih| / |stabilizer|)
  Idx orbit_size(const Idx ix) const {
    std::vector<Idx> im = orbits[ix];
    std::sort(im.begin(), im.end());
    im.erase(std::unique(im.begin(), im.end()), im.end());
    return (Idx)im.size();
  }
};
