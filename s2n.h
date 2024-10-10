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

  int n_sites;
  int n_links;

  Lattice(const int n_refine)
    : n_refine(n_refine)
  {
    {
      std::cout << "reading simplicial points" << std::endl;
      std::ifstream file("pts_n"+std::to_string(n_refine)+".dat");

      std::string str;
      std::string file_contents;
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
      std::string file_contents;
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
      std::string file_contents;
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
      std::string file_contents;
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
  }




};
