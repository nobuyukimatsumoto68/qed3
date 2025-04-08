#pragma once

#include <array>
#include <vector>
#include <Eigen/Dense>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>


struct S2Trivalent {
  using Link = std::array<Idx,2>; // <int,int>;
  using Face = std::vector<Idx>;

  std::map<const Link, const Idx> map2il;
  std::map<const Link, const int> map2sign;

  Idx n_sites;
  Idx n_links; // same
  Idx n_faces;

  std::vector<Idx> counter_accum;
  std::vector<double> ell;
  double mean_ell;

  std::vector<VE> sites;
  std::vector<std::vector<Idx> > nns;
  std::vector<Link> links;
  std::vector<Face> faces;
  std::vector<int> face_signs;

  std::vector<double> vols; // hex vol
  double mean_vol;

  std::vector<double> link_volume; // triangular vol
  double mean_link_volume;

  std::vector<double> dual_areas;
  double mean_dual_area;

  std::vector<VE> dual_sites;
  std::vector<Link> dual_links; // dual to _link


  S2Trivalent(const int n_refine)
  {
    const S2Simp simp(n_refine);

    sites = simp.dual_sites;
    n_sites = sites.size();
    nns = simp.dual_nns;

    {
      Idx counter=0;
      for(Idx ix=0; ix<n_sites; ix++){
        counter_accum.push_back(counter);
        for(Idx iy : nns[ix]) counter++;
      }
      counter_accum.push_back(counter);
    }


    links = simp.dual_links;
    n_links = links.size();

    faces = simp.dual_faces;
    n_faces = faces.size();

    vols = simp.dual_areas;
    mean_vol = simp.mean_dual_area;

    dual_sites = simp.sites;
    dual_links = simp.links;

    set_ell_link_volume();

    dual_areas = simp.vols;
    mean_dual_area = simp.mean_vol;

    face_signs = simp.dual_face_signs;

    {
      for(Idx il=0; il<n_links; il++) {
	const Link link = links[il];
	const Idx i = link[0];
	const Idx j = link[1];

	map2il.insert( { Link{i,j}, il } );
	map2il.insert( { Link{j,i}, il } );
      }

      for(Idx il=0; il<n_links; il++) {
	const auto link = links[il];
	const Idx i = link[0];
	const Idx j = link[1];

	map2sign.insert( { Link{i,j}, +1 } );
	map2sign.insert( { Link{j,i}, -1 } );
      }
    }

  }


  void set_ell_link_volume(){
    ell.resize( n_links );
    link_volume.resize( n_links );
    for(Idx il=0; il<n_links; il++) {
      const Idx ix = links[il][0];
      const Idx iy = links[il][1];

      const Idx iA = dual_links[il][0];
      const Idx iB = dual_links[il][1];

      double ellA=0.0, ellB=0.0;
      double areaA=0.0, areaB=0.0;

      const VE x = sites[ix];
      const VE y = sites[iy];
      {
        const VE p = dual_sites[iA];

        double a_ = std::acos( x.dot(p) /(x.norm()* p.norm()) );
	double b_ = std::acos( y.dot(p) /(y.norm()* p.norm()) );
	double c_ = std::acos( x.dot(y)/(x.norm()*y.norm()) ); // ell

	double s_ = 0.5*(a_+b_+c_);
	double tmp = std::tan(0.5*s_) * std::tan(0.5*(s_-a_)) * std::tan(0.5*(s_-b_)) * std::tan(0.5*(s_-c_));
	double area_ = 4.0*std::atan( std::sqrt( tmp ) );

	ellA = c_;
	areaA = area_;
      }
      {
        const VE p = dual_sites[iB];

        double a_ = std::acos( x.dot(p) /(x.norm()* p.norm()) );
	double b_ = std::acos( y.dot(p) /(y.norm()* p.norm()) );
	double c_ = std::acos( x.dot(y)/(x.norm()*y.norm()) ); // ell

	double s_ = 0.5*(a_+b_+c_);
	double tmp = std::tan(0.5*s_) * std::tan(0.5*(s_-a_)) * std::tan(0.5*(s_-b_)) * std::tan(0.5*(s_-c_));
	double area_ = 4.0*std::atan( std::sqrt( tmp ) );

	ellB = c_;
	areaB = area_;
      }

      assert( std::abs(ellA-ellB)<1.0e-14 );
      ell[il] = ellA;
      link_volume[il] = (areaA + areaB);
    }

    mean_link_volume = 0.0;
    for(const double elem : link_volume) mean_link_volume+=elem;
    assert( std::abs(mean_link_volume-4.0*M_PI)<1.0e-12 );
    mean_link_volume /= link_volume.size();

    mean_ell = 0.0;
    for(const double elem : ell) mean_ell+=elem;
    mean_ell /= ell.size();
  }

  inline int nn(const Idx ix) const {
    return 3;
  }


};
