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

    {
      std::cout << "# reading nns" << std::endl;
      std::ifstream file(dir+"nns_dual_n"+std::to_string(n_refine)+"_singlepatch.dat");

      std::string str;
      while (std::getline(file, str)){
        std::istringstream iss(str);
        Idx v1, v2, v3;
        iss >> v1;
        iss >> v2;
        iss >> v3;
        nns.push_back( std::vector<Idx>{v1,v2,v3} );
      }
    }

    {
      Idx counter=0;
      for(Idx ix=0; ix<n_sites; ix++){
        // Idx counter_tmp=0;
        counter_accum.push_back(counter);
        for(int jj=0; jj<nns[ix].size(); jj++) counter+=8;
        // counter += counter_tmp*8;
      }
      counter_accum.push_back(counter);
    }


    links = simp.dual_links;
    n_links = links.size();

    {
      std::cout << "# reading faces" << std::endl;
      std::ifstream file(dir+"face_dual_n"+std::to_string(n_refine)+".dat");
      assert(file.is_open());
      std::string str;
      while (std::getline(file, str)){
        std::istringstream iss(str);
        Idx v;
        std::vector<Idx> face;
        while( iss >> v ) face.push_back( v );
        faces.push_back( face );
      }
      n_faces = faces.size();
    }

    vols = simp.dual_areas;
    mean_vol = simp.mean_dual_area;

    simp_sites = simp.sites;
    simp_links = simp.links;

    set_ell_link_volume();

    dual_areas = simp.vols;
    mean_dual_area = simp.mean_vol;

    // sites = simp.dual_sites;
    // {
    //   std::cout << "# reading simplicial points" << std::endl;
    //   std::ifstream file(dir+"pts_n"+std::to_string(n_refine)+"_singlepatch.dat");

    //   std::string str;
    //   while (std::getline(file, str)){
    //     std::istringstream iss(str);
    //     double v1, v2, v3;
    //     iss >> v1;
    //     iss >> v2;
    //     iss >> v3;
    //     simp_sites.push_back( VE(v1, v2, v3) );
    //   }
    // }
    // {
    //   std::cout << "# reading dual points" << std::endl;
    //   std::ifstream file(dir+"pts_dual_n"+std::to_string(n_refine)+"_singlepatch.dat");

    //   std::string str;
    //   while (std::getline(file, str)){
    //     std::istringstream iss(str);
    //     double v1, v2, v3;
    //     iss >> v1;
    //     iss >> v2;
    //     iss >> v3;
    //     sites.push_back( VE(v1, v2, v3) );
    //   }
    // }
    // {
    //   std::cout << "# reading links" << std::endl;
    //   std::ifstream file(dir+"dual_links_n"+std::to_string(n_refine)+"_singlepatch.dat");

    //   std::string str;
    //   while (std::getline(file, str)){
    //     std::istringstream iss(str);
    //     Idx v1, v2;
    //     iss >> v1;
    //     iss >> v2;
    //     links.push_back( Link{v1,v2} );
    //   }
    // }

    // n_sites = sites.size();
    // n_links = links.size();

    // {
    //   std::cout << "# reading vols" << std::endl;
    //   std::ifstream file(dir+"dualtriangleareas_n"+std::to_string(n_refine)+"_singlepatch.dat");
    //   assert(file.is_open());
    //   std::string str;
    //   while (std::getline(file, str)){
    //     std::istringstream iss(str);
    //     double v1;
    //     iss >> v1;
    //     dual_areas.push_back(v1);
    //   }

    //   assert( dual_areas.size()==n_sites );
    //   mean_dual_area = 0.0;
    //   for(const double elem : dual_area) mean_dual_area+=elem;
    //   assert( std::abs(mean_dual_area-4.0*M_PI)<1.0e-12 );
    //   mean_dual_area /= dual_area.size();
    // }
    // {
    //   Idx counter=0;
    //   mean_vol=0.;
    //   for(double elem : vols){
    //     mean_vol += elem;
    //     counter++;
    //   }
    //   mean_vol /= counter;
    // }

    // // alat = std::sqrt( 8.0*M_PI/std::sqrt(3.0)/n_sites );
    // // alat = std::sqrt( mean_vol*4.0/std::sqrt(3.0) );


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

    // {
    //   for(Idx i=0; i<n_faces; i++) vps.push_back(vp(i));

    //   double sum = 0.0;
    //   for(auto elem : vps) sum += elem;
    //   assert( std::abs(sum-4.0*M_PI)<1.0e-10 );
    // }

    // set_simp_links();
  }



  // void set_simp_links() {
  //   const double threshold=0.42188 * 2.0 / n_refine;

  //   for(const Link& link : links){ // trivalent
  //     const VE x1 = sites[link[0]];
  //     const VE x2 = sites[link[1]];

  //     // Idx ip1, ip2;
  //     std::vector<Idx> tmp;
  //     for(Idx ip=0; ip<simp_sites.size(); ip++){
  //       const VE x0 = simp_sites[ip];
  //       const double d01 = (x0-x1).norm();
  //       const double d02 = (x0-x2).norm();
  //       if(d01<threshold && d02<threshold) tmp.push_back(ip);
  //     }
  //     assert( tmp.size()==2 );

  //     const Idx min = std::min(tmp[0], tmp[1]);
  //     const Idx max = std::max(tmp[0], tmp[1]);
  //     simp_links.push_back( Link{min,max} );
  //   }
  // }


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

        // double a_ = (x-p).norm();
	// double b_ = (y-p).norm();
	// double c_ = (x-y).norm(); // ell

	// double s_ = 0.5*(a_+b_+c_);
	// double tmp = s_ * (s_-a_) * (s_-b_) * (s_-c_);
	// double area_ = std::sqrt( tmp );

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

        // double a_ = (x-p).norm();
	// double b_ = (y-p).norm();
	// double c_ = (x-y).norm(); // ell

	// double s_ = 0.5*(a_+b_+c_);
	// double tmp = s_ * (s_-a_) * (s_-b_) * (s_-c_);
	// double area_ = std::sqrt( tmp );

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


  // double vp(const Idx i_face) const {
  //   double res = 0.0;

  //   const auto face = faces[i_face];
  //   const Idx n = face.size();

  //   const VE r0 = sites[face[0]];
  //   for(Idx j=1; j<n-1; j++){
  //     const Idx k=j+1;

  //     const VE r1 = sites[face[j]];
  //     const VE r2 = sites[face[k]];

  //     const double a = std::acos(r0.dot(r1));
  //     const double b = std::acos(r1.dot(r2));
  //     const double c = std::acos(r2.dot(r0));
  //     const double s = 0.5*(a+b+c);

  //     double tantan = std::tan(0.5*s);
  //     tantan *= std::tan(0.5*(s-a));
  //     tantan *= std::tan(0.5*(s-b));
  //     tantan *= std::tan(0.5*(s-c));

  //     const double tmp = 4.0 * std::atan( std::sqrt(tantan) );
  //     res += tmp;
  //   }

  //   return res;
  // }

  inline int nn(const Idx ix) const {
    return 3;
  }


};
