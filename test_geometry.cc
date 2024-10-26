#include <iostream>

#include "s2.h"

int main(){

  const int q=5; // icosahedron
  const int n_refine=2;

  QfeLatticeS2 lattice(q, n_refine);

  // std::cout << "checking typedefs:" << std::endl;
  // std::cout << "Vec3 = " << typeid(Vec3).name() << std::endl;
  // std::cout << std::endl << std::endl;

  // // ###############

  // std::cout << "checking member variables of QfeLatticeS2:" << std::endl;
  // std::cout << std::endl << std::endl;

  // std::cout << "int q = " << lattice.q << std::endl;
  // std::cout << std::endl << std::endl;

  {
    int counter=0;
    std::cout << "std::vector<Vec3> r = " << std::endl;
    std::cout << "(vertex coordinates)" << std::endl;
    for(auto elem : lattice.r) {
      std::cout << counter << " : ";
      std::cout << elem.transpose() << std::endl;
      counter++;
    }
    std::cout << std::endl
	      << "(counter = " << counter << " elements)" << std::endl
	      << std::endl;
  }

  // {
  //   int counter=0;
  //   std::cout << "std::vector<int> antipode = " << std::endl;
  //   std::cout << "(antipode of each site (0 by default))" << std::endl;
  //   for(auto elem : lattice.antipode) {
  //     std::cout << counter << " : ";
  //     std::cout << elem << std::endl;
  //     counter++;
  //   }
  //   std::cout << std::endl
  // 	      << "(counter = " << counter << " elements)" << std::endl
  // 	      << std::endl;
  // }

  // {
  //   int counter=0;
  //   std::cout << "std::vector<int> site_orbit = " << std::endl;
  //   std::cout << "(orbit id for each site)" << std::endl;
  //   for(auto elem : lattice.site_orbit) {
  //     std::cout << counter << " : ";
  //     std::cout << elem << std::endl;
  //     counter++;
  //   }
  //   std::cout << std::endl
  // 	      << "(counter = " << counter << " elements)" << std::endl
  // 	      << std::endl;
  // }

  // {
  //   int counter=0;
  //   std::cout << "std::vector<int> face_orbit = " << std::endl;
  //   std::cout << "(orbit id for each face)" << std::endl;
  //   for(auto elem : lattice.face_orbit) {
  //     std::cout << counter << " : ";
  //     std::cout << elem << std::endl;
  //     counter++;
  //   }
  //   std::cout << std::endl
  // 	      << "(counter = " << counter << " elements)" << std::endl
  // 	      << std::endl;
  // }

  // {
  //   int counter=0;
  //   std::cout << "std::vector<Vec3> orbit_xi = " << std::endl;
  //   std::cout << "(barycentric coordinates for each orbit)" << std::endl;
  //   for(auto elem : lattice.orbit_xi) {
  //     std::cout << counter << " : ";
  //     std::cout << elem.transpose() << std::endl;
  //     counter++;
  //   }
  //   std::cout << std::endl
  // 	      << "(counter = " << counter << " elements)" << std::endl
  // 	      << std::endl;
  // }

  // std::cout << "Vec3 first_face_r[4] = " << std::endl;
  // std::cout << "(coordinates of first face vertices)" << std::endl;
  // for(int i=0; i<4; i++) std::cout << i << " : "
  // 				   << lattice.first_face_r[i].transpose() << std::endl;
  // std::cout << std::endl << std::endl;

  // {
  //   int counter=0;
  //   std::cout << "std::vector<GrpElemO3> G = " << std::endl;
  //   std::cout << "(symmetry group elements)" << std::endl;
  //   for(auto elem : lattice.G) {
  //     std::cout << counter << " : ";
  //     std::cout << elem.q << std::endl;
  //     counter++;
  //   }
  //   std::cout << std::endl
  // 	      << "(counter = " << counter << " elements)" << std::endl
  // 	      << std::endl;
  // }
      
  // {
  //   int counter=0;
  //   std::cout << "std::vector<int> site_g = " << std::endl;
  //   std::cout << "(group element for each site)" << std::endl;
  //   for(auto elem : lattice.site_g) {
  //     std::cout << counter << " : ";
  //     std::cout << elem << std::endl;
  //     counter++;
  //   }
  //   std::cout << std::endl
  // 	      << "(counter = " << counter << " elements)" << std::endl
  // 	      << std::endl;
  // }

  // std::cout << "std::vector<std::vector<Complex>> ylm = " << std::endl;
  // std::cout << "(spherical harmonics)" << std::endl;
  // {
  //   int counter1=0;
  //   for(auto vector : lattice.ylm) {
  //     std::cout << counter1 << " : ";
  //     int counter2=0;
  //     for(auto elem : vector) {
  // 	std::cout << elem << " ";
  // 	counter2++;
  //     }
  //     std::cout << "(counter2 = " << counter2 << " elements)" << std::endl;

  //     counter1++;
  //     std::cout << std::endl;
  //   }
  //   std::cout << "(counter1 = " << counter1 << " elements)" << std::endl
  // 	      << std::endl;
  // }


  // std::cout << "inherited from QfeLattice:" << std::endl;
  // std::cout << std::endl << std::endl;

  // std::cout << "int n_sites = " << lattice.n_sites << std::endl;
  // std::cout << std::endl << std::endl;

  // std::cout << "int n_links = " << lattice.n_links << std::endl;
  // std::cout << std::endl << std::endl;

  // std::cout << "int n_faces = " << lattice.n_faces << std::endl;
  // std::cout << std::endl << std::endl;

  // std::cout << "int n_cells = " << lattice.n_cells << std::endl;
  // std::cout << std::endl << std::endl;

  // std::cout << "double vol = " << lattice.vol << std::endl;
  // std::cout << std::endl << std::endl;


  // std::cout << "std::vector<QfeSite> sites = " << std::endl;
  {
    int counter=0;
    for(auto elem : lattice.sites) {
      // std::cout << counter << " : "
      // 		<< elem.nn << std::endl;
      // std::cout << "wt = " << elem.wt << std::endl;
      // std::cout << "nn = " << elem.nn << std::endl;

      // std::cout << "links = ";
      // for(int i=0; i<elem.nn; i++) std::cout << elem.links[i] << " ";
      // std::cout << std::endl;

      // std::cout << "neighbors = ";
      for(int i=0; i<elem.nn; i++) std::cout << elem.neighbors[i] << " ";
      std::cout << std::endl;

      // std::cout << "id = " << elem.id << std::endl;

      counter++;
      // std::cout << std::endl;
    }
    // std::cout << "(counter = " << counter << " elements)" << std::endl
    // 	      << std::endl;
  }

  // std::cout << "std::vector<QfeLink> links = " << std::endl;
  // {
  //   int counter=0;
  //   for(auto elem : lattice.links) {
  //     std::cout << counter << " : " << std::endl;
  //     std::cout << "wt = " << elem.wt << std::endl;

  //     std::cout << "sites = ";
  //     for(int i=0; i<2; i++) std::cout << elem.sites[i] << " ";
  //     std::cout << std::endl;

  //     std::cout << "n_faces = " << elem.n_faces << std::endl;

  //     std::cout << "faces = ";
  //     for(int i=0; i<elem.n_faces; i++) std::cout << elem.faces[i] << " ";
  //     std::cout << std::endl;

  //     counter++;
  //     std::cout << std::endl;
  //   }
  //   std::cout << "(counter = " << counter << " elements)" << std::endl
  // 	      << std::endl;
  // }

  // std::cout << "std::vector<QfeFace> faces = " << std::endl;
  // {
  //   int counter=0;
  //   for(auto elem : lattice.faces) {
  //     std::cout << counter << " : " << std::endl;
  //     std::cout << "wt = " << elem.wt << std::endl;
  //     std::cout << "n_edges = " << elem.n_edges << std::endl;
  //     std::cout << "n_cells = " << elem.n_cells << std::endl;

  //     std::cout << "cells = ";
  //     for(int i=0; i<elem.n_cells; i++) std::cout << elem.cells[i] << " ";
  //     std::cout << std::endl;

  //     std::cout << "edges = ";
  //     for(int i=0; i<elem.n_edges; i++) std::cout << elem.edges[i] << " ";
  //     std::cout << std::endl;

  //     std::cout << "sites = ";
  //     for(int i=0; i<elem.n_edges; i++) std::cout << elem.sites[i] << " ";
  //     std::cout << std::endl;

  //     counter++;
  //     std::cout << std::endl;
  //   }
  //   std::cout << "(counter = " << counter << " elements)" << std::endl
  // 	      << std::endl;
  // }

  // std::cout << "std::vector<QfeCell> cells = " << std::endl;
  // {
  //   int counter=0;
  //   for(auto elem : lattice.cells) {
  //     std::cout << counter << " : " << std::endl;
  //     std::cout << "n_faces = " << elem.n_faces << std::endl;

  //     std::cout << "faces = ";
  //     for(int i=0; i<elem.n_faces; i++) std::cout << elem.faces[i] << " ";
  //     std::cout << std::endl;

  //     std::cout << "sites = ";
  //     for(int i=0; i<elem.n_faces; i++) std::cout << elem.sites[i] << " ";
  //     std::cout << std::endl;

  //     counter++;
  //     std::cout << std::endl;
  //   }
  //   std::cout << "(counter = " << counter << " elements)" << std::endl
  // 	      << std::endl;
  // }


  // std::vector<QfeFace> faces;
  // std::vector<QfeCell> cells;

  // // symmetrically distinct sites
  // std::vector<int> distinct_n_sites;  // number of sites for each distinct id
  // std::vector<int> distinct_first;    // representative site for distinct group
  // int n_distinct;

  // QfeRng rng;

  return 0;
}
