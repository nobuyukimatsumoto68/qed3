#pragma once

#include "s2.h"

struct S2Simp {
  using Link = std::array<Idx,2>; // <Idx,Idx>;
  using Face = std::vector<Idx>;
  using MS=Eigen::Matrix2cd;
  using VD=Eigen::Vector2d;
  using VE=Eigen::Vector3d;
  using VC=Eigen::VectorXcd;

  std::map<const Link, const Idx> map2il;
  std::map<const Link, const int> map2sign;

  Idx n_sites;
  Idx n_links;
  Idx n_faces;

  std::vector<Idx> counter_accum;
  std::vector<double> ell; // evan's link label
  double mean_ell;
  double alat; // from volume

  std::vector<VE> sites;
  std::vector<std::vector<Idx>> nns;
  std::vector<Link> links;
  std::vector<Face> faces;

  std::vector<double> vols; // triangular vol
  double mean_vol;

  // std::vector<double> ell; // evan's link label
  std::vector<double> ellstarA; // evan's link label
  std::vector<double> ellstarB; // evan's link label
  std::vector<double> link_volume; // evan's link label
  double mean_link_volume;


  std::vector<VE> dual_sites;
  std::vector<Link> dual_links;

  S2Simp(const int n_refine)
  {
    {
      std::cout << "# reading simplicial points" << std::endl;
      std::ifstream file(dir+"pts_n"+std::to_string(n_refine)+"_singlepatch.dat");

      std::string str;
      while (std::getline(file, str)){
	std::istringstream iss(str);
	double v1, v2, v3;
	iss >> v1;
	iss >> v2;
	iss >> v3;
	sites.push_back( VE(v1, v2, v3) );
      }
    }
    n_sites = sites.size();

    {
      std::cout << "# reading dual points" << std::endl;
      std::ifstream file(dir+"pts_dual_n"+std::to_string(n_refine)+"_singlepatch.dat");

      std::string str;
      while (std::getline(file, str)){
	std::istringstream iss(str);
	double v1, v2, v3;
	iss >> v1;
	iss >> v2;
	iss >> v3;
	dual_sites.push_back( VE(v1, v2, v3) );
      }
    }


    {
      std::cout << "# reading links" << std::endl;
      std::ifstream file(dir+"links_n"+std::to_string(n_refine)+"_singlepatch.dat");

      std::string str;
      while (std::getline(file, str)){
	std::istringstream iss(str);
	Idx v1, v2;
	iss >> v1;
	iss >> v2;
	links.push_back( Link{v1,v2} );
      }
    }
    n_links = links.size();

    {
      std::cout << "# reading dual links" << std::endl;
      std::ifstream file(dir+"dual_links_n"+std::to_string(n_refine)+"_singlepatch.dat");

      std::string str;
      while (std::getline(file, str)){
	std::istringstream iss(str);
	Idx v1, v2;
	iss >> v1;
	iss >> v2;
	dual_links.push_back( Link{v1,v2} );
      }
    }

    {
      std::ifstream file(dir+"face_n"+std::to_string(n_refine)+"_singlepatch.dat");

      std::string str;
      // std::cout << "debug file = " << std::endl;
      while (std::getline(file, str)){
        std::istringstream iss(str);
        // std::cout << "debug = " << std::endl;
        Idx v;
        Face face;
        while( iss >> v ) {
          // std::cout << "face = " << v << std::endl;
          face.push_back( v );
        }
        faces.push_back( face );
      }
      n_faces = faces.size();
    }

    {
      std::cout << "# reading nns" << std::endl;
      std::ifstream file(dir+"nns_n"+std::to_string(n_refine)+"_singlepatch.dat");

      std::string str;
      while (std::getline(file, str)){
	std::istringstream iss(str);
	Idx v;
	std::vector<Idx> nn;
	while( iss >> v ) nn.push_back( v );
	nns.push_back( nn );
      }
    }

    // {
    //   std::cout << "# reading vols" << std::endl;
    //   std::ifstream file(dir+"dualtriangleareas_n"+std::to_string(n_refine)+"_singlepatch.dat");
    //   assert(file.is_open());
    //   std::string str;
    //   while (std::getline(file, str)){
    //     std::istringstream iss(str);
    //     double v1;
    //     iss >> v1;
    //     vol.push_back(v1);
    //   }

    //   assert( vol.size()==n_sites );
    // }


    {
      Idx counter=0;
      for(Idx ix=0; ix<n_sites; ix++){
        // Idx counter_tmp=0;
        counter_accum.push_back(counter);
        for(int jj=0; jj<nns[ix].size(); jj++) counter+=8;
        // counter += counter_tmp*8;
      }
    }

    {
      // for(Idx ix=0; ix<n_sites; ix++){
      //   for(int jj=0; jj<nns[ix].size(); jj++){
      //     const Idx iy =nns[ix][jj];
      //     if(iy>ix) continue;

      //     const int il = sites[ix].links[jj];
      //     links[il] = Link{ix,iy}; // ix<iy
      //   }
      // }

      for(Idx il=0; il<n_links; il++) {
	const Link link = links[il];
	const Idx i = link[0];
	const Idx j = link[1];

	map2il.insert( { Link{i,j}, il } );
	map2il.insert( { Link{j,i}, il } );
      }

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


    {
      mean_vol = 0.0;
      Idx counter=0;
      std::cout << "# reading vols" << std::endl;
      std::ifstream file(dir+"dualtriangleareas_n"+std::to_string(n_refine)+"_singlepatch.dat");
      assert(file.is_open());
      std::string str;
      while (std::getline(file, str)){
	std::istringstream iss(str);
	double v1;
	iss >> v1;
        vols.push_back(v1);
        mean_vol += v1;
	counter++;
      }
      assert( vols.size()==n_faces );
      assert( std::abs(mean_vol-4.0*M_PI)<1.0e-12 );

      mean_vol /= counter;
      alat = std::sqrt( mean_vol*4.0/std::sqrt(3.0) );
    }
    // alat = std::sqrt( 8.0*M_PI/std::sqrt(3.0)/n_sites );
    // alat = std::sqrt( 8.0*M_PI/std::sqrt(3.0)/this->n_faces );

    set_ell_ellstar_linkvols();
    // set_ell();
  }

  // inline int nn(const Idx ix) const {
  //   return sites[ix].nn;
  // }



  // void set_ell() {
  //   ell.resize( n_links );
  //   mean_ell = 0.0;
  //   for(Idx il=0; il<n_links; il++) {
  //     const Link link = links[il];
  //     const Idx ix = links[il][0];
  //     const Idx iy = links[il][1];

  //     const VE x = sites[ix];
  //     const VE y = sites[iy];
  //     ell[il] = std::acos( x.dot(y)/(x.norm()*y.norm()) ); // ell
  //     mean_ell += ell[il];
  //   }
  //   mean_ell /= n_links;
  // }

  void set_ell_ellstar_linkvols(){
    ell.resize( n_links );
    ellstarA.resize( n_links );
    ellstarB.resize( n_links );
    link_volume.resize( n_links );

    for(Idx il=0; il<n_links; il++) {
      const Link link = links[il];

      const VE x = sites[ links[il][0] ];
      const VE y = sites[ links[il][1] ];

      const Idx iA = *std::min_element( dual_links[il].begin(), dual_links[il].end() );
      const Idx iB = *std::max_element( dual_links[il].begin(), dual_links[il].end() );

      double ellA=0.0, ellB=0.0;
      double areaA=0.0, areaB=0.0;
      {
        const VE p = dual_sites[iA];

        double a_ = std::acos( x.dot(p) /(x.norm()* p.norm()) );
        double b_ = std::acos( y.dot(p) /(y.norm()* p.norm()) );
        double c_ = std::acos( x.dot(y)/(x.norm()*y.norm()) ); // ell

        double s_ = 0.5*(a_+b_+c_);
        double tmp = std::tan(0.5*s_) * std::tan(0.5*(s_-a_)) * std::tan(0.5*(s_-b_)) * std::tan(0.5*(s_-c_));
        double area_ = 4.0*std::atan( std::sqrt( tmp ) );
        // double area_ = a_+b_+c_ - M_PI;

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
        // double area_ = a_+b_+c_ - M_PI;

        ellB = c_;
        areaB = area_;
      }

      assert( std::abs(ellA-ellB)<1.0e-14 );
      ell[il] = ellA;
      ellstarA[il] = areaA/ell[il];
      ellstarB[il] = areaB/ell[il];
      link_volume[il] = areaA + areaB;
    }

    mean_link_volume = 0.0;
    for(const double elem : link_volume) mean_link_volume+=elem;
    // std::cout << "debug. sum = " << sum << std::endl;
    assert( std::abs(mean_link_volume-4.0*M_PI)<1.0e-12 );
    mean_link_volume /= link_volume.size();

    mean_ell = 0.0;
    for(const double elem : ell) mean_ell+=elem;
    mean_ell /= ell.size();
  }








};
