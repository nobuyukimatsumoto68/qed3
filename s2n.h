#pragma once

#include <cmath>
#include <complex>
#include <vector>
#include <array>
#include <map>

#include <string>
#include <sstream>

#include <Eigen/Dense>

using Complex = std::complex<double>;
const Complex I = Complex(0.0, 1.0);

using Link = std::array<int,2>; // <int,int>;

constexpr int NS = 2;
using MS=Eigen::Matrix2cd;

constexpr int DIM = 2;
using VD=Eigen::Vector2d;

constexpr int EDIM = 3;
using VE=Eigen::Vector3d;


struct Lattice {
  const int n_refine;

  std::vector<VE> simp_sites;
  std::vector<VE> sites;
  std::vector<std::vector<int> > nns;
  std::vector<Link> links;

  // std::vector<std::array<VD, 3>> exM; // site, M=A,B,C
  std::vector<std::array<double, 3>> u; // site, M=A,B,C
  std::vector<std::array<VD, 3>> v; // site, M=A,B,C
  std::vector<double> vol; // site

  int n_sites;
  int n_links;

  Lattice(const int n_refine)
    : n_refine(n_refine)
  {
    {
      std::cout << "reading simplicial points" << std::endl;
      std::ifstream file("pts_n"+std::to_string(n_refine)+".dat");

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
      std::cout << "reading dual points" << std::endl;
      std::ifstream file("pts_dual_n"+std::to_string(n_refine)+".dat");

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
      std::cout << "reading nns" << std::endl;
      std::ifstream file("nns_dual_n"+std::to_string(n_refine)+".dat");

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
      std::cout << "reading links" << std::endl;
      std::ifstream file("dual_links_n"+std::to_string(n_refine)+".dat");

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
      std::cout << "reading vs" << std::endl;
      std::ifstream file("vs_n"+std::to_string(n_refine)+".dat");

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
      // std::cout << "set. size=" << exM.size() << std::endl;
      assert( v.size()==n_sites );
      // std::cout << "assert done. " << std::endl;
    }

    {
      std::cout << "reading us" << std::endl;
      std::ifstream file("us_n"+std::to_string(n_refine)+".dat");

      std::string str;
      while (std::getline(file, str)){
	std::istringstream iss(str);
	double v1, v2, v3;
	iss >> v1;
	iss >> v2;
	iss >> v3;
	// const VD e3(v5,v6);
	u.push_back( std::array<double,3>{v1,v2,v3} );
      }

      // std::cout << "set. size=" << exM.size() << std::endl;
      assert( u.size()==n_sites );
      // std::cout << "assert done. " << std::endl;
    }

    {
      std::cout << "reading vols" << std::endl;
      // std::ifstream file("dualtriangleareas_n"+std::to_string(n_refine)+".dat");
      std::ifstream file("dualtriangleareas_n"+std::to_string(n_refine)+".dat");
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

    // {
    //   std::cout << "reading es" << std::endl;
    //   std::ifstream file("estars_n"+std::to_string(n_refine)+".dat");

    //   std::string str;
    //   std::string file_contents;
    //   while (std::getline(file, str)){
    // 	std::istringstream iss(str);
    // 	double v1, v2, v3, v4, v5, v6;
    // 	iss >> v1;
    // 	iss >> v2;
    // 	const VD e1(v1,v2);
    // 	iss >> v3;
    // 	iss >> v4;
    // 	const VD e2(v3,v4);
    // 	iss >> v5;
    // 	iss >> v6;
    // 	const VD e3(v5,v6);
    // 	exM.push_back( std::array<VD,3>{e1,e2,e3} );
    // 	std::cout << v6 << std::endl;
    //   }

    //   std::cout << "set. size=" << exM.size() << std::endl;
    //   assert( exM.size()==n_sites );
    //   std::cout << "assert done. " << std::endl;
    // }



  }




};
