#pragma once

// #include "s2.h"

struct S2Simp {
  using Link = std::array<Idx,2>; // <int,int>;
  using Face = std::vector<Idx>;

  std::map<const Link, const Idx> map2il;
  std::map<const Link, const int> map2sign;

  Idx n_sites;
  Idx n_links;
  Idx n_faces;

  std::vector<Idx> counter_accum;
  std::vector<double> ell; // evan's link label
  double mean_ell;

  std::vector<VE> sites;
  std::vector<std::vector<Idx>> nns;
  std::vector<Link> links;
  std::vector<Face> faces;
  std::vector<int> face_signs;

  std::vector<double> vols; // triangular vol
  double mean_vol;

  std::vector<double> link_volume; // evan's link label
  double mean_link_volume;

  std::vector<VE> dual_sites;
  std::vector<std::vector<Idx>> dual_nns;
  std::vector<Link> dual_links;
  std::vector<Face> dual_faces;
  std::vector<int> dual_face_signs;

  std::vector<double> dual_areas;
  double mean_dual_area;

  double alat;

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
      if(!file) assert(false);

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
      std::cout << "# reading nns" << std::endl;
      std::ifstream file(dir+"nns_dual_n"+std::to_string(n_refine)+"_singlepatch.dat");
      if(!file) assert(false);

      std::string str;
      while (std::getline(file, str)){
        std::istringstream iss(str);
        Idx v1, v2, v3;
        iss >> v1;
        iss >> v2;
        iss >> v3;
        dual_nns.push_back( std::vector<Idx>{v1,v2,v3} );
      }
    }


    {
      std::cout << "# reading links" << std::endl;
      std::ifstream file(dir+"links_n"+std::to_string(n_refine)+"_singlepatch.dat");
      if(!file) assert(false);

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
    // std::cout << "debug. n_links = " << n_links << std::endl;

    {
      std::cout << "# reading dual links" << std::endl;
      std::ifstream file(dir+"dual_links_n"+std::to_string(n_refine)+"_singlepatch.dat");
      if(!file) assert(false);

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
      if(!file) assert(false);

      std::string str;
      while (std::getline(file, str)){
        std::istringstream iss(str);
        Idx v;
        Face face;
        while( iss >> v ) {
          face.push_back( v );
        }
        faces.push_back( face );
      }
      n_faces = faces.size();
    }

    {
      std::cout << "# reading nns" << std::endl;
      std::ifstream file(dir+"nns_n"+std::to_string(n_refine)+"_singlepatch.dat");
      if(!file) assert(false);

      std::string str;
      while (std::getline(file, str)){
	std::istringstream iss(str);
	Idx v;
	std::vector<Idx> nn;
	while( iss >> v ) nn.push_back( v );
	nns.push_back( nn );
      }
    }

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
        dual_faces.push_back( face );
      }
    }

    // std::cout << "# debug1" << std::endl;
    {
      Idx counter=0;
      std::cout << "# n_sites = " << n_sites << std::endl;
      for(Idx ix=0; ix<n_sites; ix++){
        std::cout << "# ix = " << ix << std::endl;
        counter_accum.push_back(counter);
        std::cout << "# ix = " << ix << std::endl;
        std::cout << "# nns[ix] = " << nns[ix][0] << std::endl;
        for(Idx iy : nns[ix]) counter++;
      }
      counter_accum.push_back(counter);
    }
    // std::cout << "# debug2" << std::endl;
    {
      for(Idx il=0; il<n_links; il++) {
	const Link link = links[il];
	const Idx i = link[0];
	const Idx j = link[1];

	map2il.insert( { Link{i,j}, il } );
	map2il.insert( { Link{j,i}, il } );
      }
      // std::cout << "# debug3" << std::endl;
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

    // std::cout << "# debug4" << std::endl;
    {
      mean_vol = 0.0;
      Idx counter=0;
      // std::cout << "# reading vols" << std::endl;
      // std::ifstream file(dir+"dualtriangleareas_n"+std::to_string(n_refine)+"_singlepatch.dat");
      // assert(file.is_open());
      // std::string str;
      // while (std::getline(file, str)){
      // std::istringstream iss(str);
      // double v1;
      // iss >> v1;
      for(const Face& face : faces ){
        Idx i0=face[0];
        Idx i1=face[1];
        Idx i2=face[2];

        const VE x0 = sites[ i0 ];
        const VE x1 = sites[ i1 ];
        const VE x2 = sites[ i2 ];

        double a_ = std::acos( x0.dot(x1) );
        double b_ = std::acos( x1.dot(x2) );
        double c_ = std::acos( x2.dot(x0) );

        double s_ = 0.5*(a_+b_+c_);
        double tmp = std::tan(0.5*s_) * std::tan(0.5*(s_-a_)) * std::tan(0.5*(s_-b_)) * std::tan(0.5*(s_-c_));
        double area = 4.0*std::atan( std::sqrt( tmp ) );

        vols.push_back( area );
        mean_vol += area;
        counter++;
      }
      assert( vols.size()==n_faces );
      assert( std::abs(mean_vol-4.0*M_PI)<1.0e-10 );

      mean_vol /= counter;
      alat = std::sqrt( mean_vol*4.0/std::sqrt(3.0) );
    }

    // std::cout << "# debug5" << std::endl;

    set_ell_ellstar_linkvols();
    // std::cout << "# debug6" << std::endl;
    set_dual_areas();
    // std::cout << "# debug7" << std::endl;
    set_facesigns();
  }

  void set_facesigns(){
    {
      face_signs.clear();
      for(const Face& face : faces){
        const VE x0 = sites[face[0]];
        const VE x1 = sites[face[1]];
        const VE x2 = sites[face[2]];

        const double dot = ((x2-x1).cross(x0-x1)).dot(x1);

        const int sign = dot>0 ? 1 : -1;
        face_signs.push_back(sign);
      }
    }
    {
      dual_face_signs.clear();
      for(const Face& face : dual_faces){
        const VE x0 = sites[face[0]];
        const VE x1 = sites[face[1]];
        const VE x2 = sites[face[2]];

        const double dot = ((x2-x1).cross(x0-x1)).dot(x1);

        const int sign = dot>0 ? 1 : -1;
        dual_face_signs.push_back(sign);
      }
    }
  }


  void set_ell_ellstar_linkvols(){
    ell.resize( n_links );
    link_volume.resize( n_links );

    for(Idx il=0; il<n_links; il++) {
      // std::cout << "il = " << il << std::endl;
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
      link_volume[il] = areaA + areaB;
      // std::cout << "areas = " << areaA << " " << areaB << std::endl;
    }

    mean_link_volume = 0.0;
    for(const double elem : link_volume) {
      mean_link_volume+=elem;
    }
    // std::cout << "debug. mean_link_volume = " << 
    // mean_link_volume << std::endl;
    assert( std::abs(mean_link_volume-4.0*M_PI)<1.0e-10 );
    mean_link_volume /= link_volume.size();

    mean_ell = 0.0;
    for(const double elem : ell) mean_ell+=elem;
    mean_ell /= ell.size();
  }

  void set_dual_areas(){
    dual_areas.clear();
    dual_areas.resize(n_sites);

#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL)
#endif
    for(Idx ip=0; ip<n_sites; ip++){
      double area=0.0;

      const VE p = sites[ip];
      const Face face = dual_faces[ip];

      for(int i=0; i<face.size(); i++){
        const Idx ix = face[i];
        const Idx iy = face[(i+1)%face.size()];

        const VE x = dual_sites[ ix ];
        const VE y = dual_sites[ iy ];

        double a_ = std::acos( x.dot(p) / (x.norm()* p.norm()) );
        double b_ = std::acos( y.dot(p) / (y.norm()* p.norm()) );
        double c_ = std::acos( x.dot(y) / (x.norm()*y.norm()) ); // ell

        double s_ = 0.5*(a_+b_+c_);
        double tmp = std::tan(0.5*s_) * std::tan(0.5*(s_-a_)) * std::tan(0.5*(s_-b_)) * std::tan(0.5*(s_-c_));
        area += 4.0*std::atan( std::sqrt( tmp ) );
      }

      dual_areas[ip] = area;
    } // end for ix

    mean_dual_area = 0.0;
    for(const double elem : dual_areas) mean_dual_area+=elem;
    assert( std::abs(mean_dual_area-4.0*M_PI)<1.0e-10 );
    mean_dual_area /= dual_areas.size();
  }






};

    // for(Idx il=0; il<n_links; il++) {
    //   const Link link = links[il];

    //   const VE x = sites[ links[il][0] ];
    //   const VE y = sites[ links[il][1] ];

    //   const Idx iA = *std::min_element( dual_links[il].begin(), dual_links[il].end() );
    //   const Idx iB = *std::max_element( dual_links[il].begin(), dual_links[il].end() );

    //   double ellA=0.0, ellB=0.0;
    //   double areaA=0.0, areaB=0.0;
    //   {
    //     const VE p = dual_sites[iA];

    //     double a_ = std::acos( x.dot(p) /(x.norm()* p.norm()) );
    //     double b_ = std::acos( y.dot(p) /(y.norm()* p.norm()) );
    //     double c_ = std::acos( x.dot(y)/(x.norm()*y.norm()) ); // ell

    //     double s_ = 0.5*(a_+b_+c_);
    //     double tmp = std::tan(0.5*s_) * std::tan(0.5*(s_-a_)) * std::tan(0.5*(s_-b_)) * std::tan(0.5*(s_-c_));
    //     double area_ = 4.0*std::atan( std::sqrt( tmp ) );
    //     // double area_ = a_+b_+c_ - M_PI;

    //     ellA = c_;
    //     areaA = area_;
    //   }
    //   {
    //     const VE p = dual_sites[iB];

    //     double a_ = std::acos( x.dot(p) /(x.norm()* p.norm()) );
    //     double b_ = std::acos( y.dot(p) /(y.norm()* p.norm()) );
    //     double c_ = std::acos( x.dot(y)/(x.norm()*y.norm()) ); // ell

    //     double s_ = 0.5*(a_+b_+c_);
    //     double tmp = std::tan(0.5*s_) * std::tan(0.5*(s_-a_)) * std::tan(0.5*(s_-b_)) * std::tan(0.5*(s_-c_));
    //     double area_ = 4.0*std::atan( std::sqrt( tmp ) );
    //     // double area_ = a_+b_+c_ - M_PI;

    //     ellB = c_;
    //     areaB = area_;
    //   }
