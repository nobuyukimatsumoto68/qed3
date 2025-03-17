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


// std::string dir = "/mnt/hdd_barracuda/qed3/dats/";
// std::string dir = "/mnt/hdd_barracuda/qed3/dats_save/";


struct S2Trivalent {
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

  std::vector<double> vols; // hex vol
  double mean_vol;

  std::vector<double> link_volume; // triangular vol
  double mean_link_volume;

  std::vector<double> dual_areas;
  double mean_dual_area;

  std::vector<VE> simp_sites;
  std::vector<Link> simp_links; // dual to _link

  // const int n_refine;


  S2Trivalent(const int n_refine)
  {
    const S2Simp simp(n_refine);

    sites = simp.dual_sites;
    n_sites = sites.size();
    // std::cout << "n_sites = " << n_sites;

    nns = simp.dual_nns;

    {
      Idx counter=0;
      for(Idx ix=0; ix<n_sites; ix++){
        // Idx counter_tmp=0;
        counter_accum.push_back(counter);
        for(int jj=0; jj<nns[ix].size(); jj++) counter++;
        // counter += counter_tmp*8;
      }
      counter_accum.push_back(counter);
    }


    links = simp.dual_links;
    n_links = links.size();

    faces = simp.dual_faces;
    n_faces = faces.size();

    vols = simp.dual_areas;
    mean_vol = simp.mean_dual_area;

    simp_sites = simp.sites;
    simp_links = simp.links;

    set_ell_link_volume();

    dual_areas = simp.vols;
    mean_dual_area = simp.mean_vol;


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

      const Idx iA = simp_links[il][0];
      const Idx iB = simp_links[il][1];

      double ellA=0.0, ellB=0.0;
      double areaA=0.0, areaB=0.0;

      const VE x = sites[ix];
      const VE y = sites[iy];
      {
        const VE p = simp_sites[iA];

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
        const VE p = simp_sites[iB];

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
