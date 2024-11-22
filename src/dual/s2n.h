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



struct Lattice {
  using Link = std::array<int,2>; // <int,int>;
  using Face = std::vector<int>;
  using MS=Eigen::Matrix2cd;
  using VD=Eigen::Vector2d;
  using VE=Eigen::Vector3d;
  using VC=Eigen::VectorXcd;

  const int n_refine;

  std::vector<VE> simp_sites;
  std::vector<VE> sites;
  std::vector<std::vector<int> > nns;
  std::vector<Link> links;

  std::vector<std::array<double, 3>> u; // site, M=A,B,C
  std::vector<std::array<VD, 3>> v; // site, M=A,B,C
  std::vector<double> vol; // site

  std::vector<Face> faces;
  std::vector<double> vps; // face

  std::map<const Link, const int> map2il;
  std::map<const Link, const int> map2sign;

  double mean_vol;

  int n_sites;
  int n_links;
  int n_faces;

  Lattice(const int n_refine)
    : n_refine(n_refine)
  {
    {
      std::cout << "# reading simplicial points" << std::endl;
      std::ifstream file("../../dats/pts_n"+std::to_string(n_refine)+"_singlepatch.dat");

      std::string str;
      while (std::getline(file, str)){
	std::istringstream iss(str);
	double v1, v2, v3;
	iss >> v1;
	iss >> v2;
	iss >> v3;
	simp_sites.push_back( VE(v1, v2, v3) );
      }
    }
    {
      std::cout << "# reading dual points" << std::endl;
      std::ifstream file("../../dats/pts_dual_n"+std::to_string(n_refine)+"_singlepatch.dat");

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
    {
      std::cout << "# reading nns" << std::endl;
      std::ifstream file("../../dats/nns_dual_n"+std::to_string(n_refine)+"_singlepatch.dat");

      std::string str;
      while (std::getline(file, str)){
	std::istringstream iss(str);
	int v1, v2, v3;
	iss >> v1;
	iss >> v2;
	iss >> v3;
	nns.push_back( std::vector<int>{v1,v2,v3} );
      }
    }
    {
      std::cout << "# reading links" << std::endl;
      std::ifstream file("../../dats/dual_links_n"+std::to_string(n_refine)+"_singlepatch.dat");

      std::string str;
      while (std::getline(file, str)){
	std::istringstream iss(str);
	int v1, v2;
	iss >> v1;
	iss >> v2;
	links.push_back( Link{v1,v2} );
      }
    }

    n_sites = sites.size();
    n_links = links.size();

    {
      std::cout << "# reading vs" << std::endl;
      std::ifstream file("../../dats/vs_n"+std::to_string(n_refine)+"_singlepatch.dat");

      std::string str;
      while (std::getline(file, str)){
	std::istringstream iss(str);
	double v1, v2, v3, v4, v5, v6;
	iss >> v1;
	iss >> v2;
	const VD vv1(v1,v2);
	iss >> v3;
	iss >> v4;
	const VD vv2(v3,v4);
	iss >> v5;
	iss >> v6;
	const VD vv3(v5,v6);
	v.push_back( std::array<VD,3>{vv1,vv2,vv3} );
      }
      assert( v.size()==n_sites );
    }

    {
      std::cout << "# reading us" << std::endl;
      std::ifstream file("../../dats/us_n"+std::to_string(n_refine)+"_singlepatch.dat");

      std::string str;
      while (std::getline(file, str)){
	std::istringstream iss(str);
	double v1, v2, v3;
	iss >> v1;
	iss >> v2;
	iss >> v3;
	u.push_back( std::array<double,3>{v1,v2,v3} );
      }

      assert( u.size()==n_sites );
    }

    {
      std::cout << "# reading vols" << std::endl;
      std::ifstream file("../../dats/dualtriangleareas_n"+std::to_string(n_refine)+"_singlepatch.dat");
      assert(file.is_open());
      std::string str;
      while (std::getline(file, str)){
	std::istringstream iss(str);
	double v1;
	iss >> v1;
	vol.push_back(v1);
      }

      assert( vol.size()==n_sites );
    }
    {
      int counter=0;
      for(double elem : vol){
	mean_vol += elem;
	counter++;
      }
      mean_vol /= counter;
    }

    {
      std::cout << "# reading faces" << std::endl;
      std::ifstream file("../../dats/face_dual_n"+std::to_string(n_refine)+".dat");
      assert(file.is_open());
      std::string str;
      while (std::getline(file, str)){
	std::istringstream iss(str);
	int v;
	std::vector<int> face;
	while( iss >> v ) face.push_back( v );
	faces.push_back( face );
      }
      n_faces = faces.size();
    }

    {
      for(int il=0; il<n_links; il++) {
	const Link link = links[il];
	const int i = link[0];
	const int j = link[1];

	map2il.insert( { Link{i,j}, il } );
	map2il.insert( { Link{j,i}, il } );
      }

      for(int il=0; il<n_links; il++) {
	const auto link = links[il];
	const int i = link[0];
	const int j = link[1];

	map2sign.insert( { Link{i,j}, +1 } );
	map2sign.insert( { Link{j,i}, -1 } );
      }
    }

    {
      for(int i=0; i<n_faces; i++) vps.push_back(vp(i));

      double sum = 0.0;
      for(auto elem : vps) sum += elem;
      assert( std::abs(sum-4.0*M_PI)<1.0e-10 );
    }  
  }


  double vp(const int i_face) const {
    double res = 0.0;

    const auto face = faces[i_face];
    const int n = face.size();

    const VE r0 = sites[face[0]];
    for(int j=1; j<n-1; j++){
      const int k=j+1;

      const VE r1 = sites[face[j]];
      const VE r2 = sites[face[k]];

      const double a = std::acos(r0.dot(r1));
      const double b = std::acos(r1.dot(r2));
      const double c = std::acos(r2.dot(r0));
      const double s = 0.5*(a+b+c);

      double tantan = std::tan(0.5*s);
      tantan *= std::tan(0.5*(s-a));
      tantan *= std::tan(0.5*(s-b));
      tantan *= std::tan(0.5*(s-c));
      
      const double tmp = 4.0 * std::atan( std::sqrt(tantan) );
      res += tmp;
    }

    return res;
  }

};
