#include <iostream>
#include <cstdlib>
#include <ctime>
// #include <boost/math/special_functions/bessel.hpp>

// #include "statistics.h"
#include "dirac_s2.h"


#include <algorithm>
#include <iterator>


int main(int argc, char* argv[]){

  // geometry
  const int q=5; // icosahedron
  const int n_refine=1; // no refinement

  // ----------------------------------

  Eigen::IOFormat fmt(15, 0,
		      "\t",
		      "\n",
		      "",
		      "",
		      "",
		      "");

  // ----------------------------------

  QfeLatticeS2 lattice(q, n_refine);
  Dirac1fonS2 D(lattice);

  {
    const int iN = 0;
    const int iS = 6;
    assert( std::abs(lattice.r[iN][2]-1.0)<1.0e-15 );
    assert( std::abs(lattice.r[iS][2]+1.0)<1.0e-15 );

    double alpha0_N = 0.0;
    double alpha0_S = 3.0*M_PI/2.0 - 0.5 * 0.628319;

    { // N
      D.alpha[iN].resize(lattice.sites[iN].nn);
      for(int jj=0; jj<lattice.sites[iN].nn; jj++){
	const int il = D.link_oriented[iN][jj];
	D.omega[il] = 0.0;

	D.alpha[iN][jj] = double_mod( alpha0_N + 2.0*M_PI/5.0 * jj );
      }
    } // end N
    { // S
      D.alpha[iS].resize(lattice.sites[iS].nn);
      for(int jj=0; jj<lattice.sites[iS].nn; jj++){
	const int il = D.link_oriented[iS][jj];
	D.omega[il] = 0.0;

	D.alpha[iS][jj] = double_mod( alpha0_S + 2.0*M_PI/5.0 * jj );
      }
    } // end S
    { // nn of N
      for(int jj=0; jj<lattice.sites[iN].nn; jj++){
	const int il = D.link_oriented[iN][jj];
	const int ix = D.nn_oriented[iN][jj];
	D.alpha[ix].resize(lattice.sites[ix].nn);

	int kk0=0; // find kk0 that corresponds to il
	for(; kk0<lattice.sites[ix].nn; kk0++) if(D.link_oriented[ix][kk0]==il) break;
	assert(kk0!=lattice.sites[ix].nn);

	D.alpha[ix][kk0] = D.alpha[iN][jj] - D.omega[il];

	for(int dk=1; dk<lattice.sites[ix].nn; dk++){
	  const int kk = (kk0+dk)%lattice.sites[ix].nn;
	  D.alpha[ix][kk] = double_mod( D.alpha[ix][kk0] + 2.0*M_PI/5.0 * jj );
	}
      }
    } // end nn of N
    { // nn of S
      for(int jj=0; jj<lattice.sites[iS].nn; jj++){
	const int il = D.link_oriented[iS][jj];
	const int ix = D.nn_oriented[iS][jj];
	D.alpha[ix].resize(lattice.sites[ix].nn);

	int kk0=0; // find kk0 that corresponds to il
	for(; kk0<lattice.sites[ix].nn; kk0++) if(D.link_oriented[ix][kk0]==il) break;
	assert(kk0!=lattice.sites[ix].nn);
	D.alpha[ix][kk0] = D.alpha[iS][jj] - D.omega[il];

	for(int dk=1; dk<lattice.sites[ix].nn; dk++){
	  const int kk = (kk0+dk)%lattice.sites[ix].nn;
	  D.alpha[ix][kk] = double_mod( D.alpha[ix][kk0] + 2.0*M_PI/5.0 * jj );
	}
      }
    } // end nn of S
    { // set omega
      for(int il=0; il<lattice.n_links; il++){
	const QfeLink link = lattice.links[il];
	const int ix = std::min(link.sites[0], link.sites[1]);
	const int iy = std::max(link.sites[0], link.sites[1]);

	const auto itx_ell = std::find(D.link_oriented[ix].begin(),
				       D.link_oriented[ix].end(),
				       il);
	const auto ity_ell = std::find(D.link_oriented[iy].begin(),
				       D.link_oriented[iy].end(),
				       il);
	const int ixl = std::distance(D.link_oriented[ix].begin(), itx_ell);
	const int iyl = std::distance(D.link_oriented[iy].begin(), ity_ell);

	D.omega[il] = double_mod( D.alpha[ix][ixl] - D.alpha[iy][iyl] - M_PI );
      }
    }
    
  }




  // {
  //   for(auto elem : D.omega) {
  //     std::cout << std::setw(25) << elem << " ";
  //   }
  //   std::cout << std::endl;
  //   // std::cout << projection(vec).format(fmt).transpose() << std::endl;
  // }


  {
    for(int ia=0; ia<lattice.n_faces; ia++){
      double delta_alpha_sum = 0.0;
      double omega_sum = 0.0;
      for(int i=0; i<3; i++){
	int il = lattice.faces[ia].edges[i];
	omega_sum += D.omega[il];

	// ----------------------------

	const QfeLink link = lattice.links[il];
	const int ix = std::min(link.sites[0], link.sites[1]);
	const int iy = std::max(link.sites[0], link.sites[1]);

	const auto itx_ell = std::find(D.link_oriented[ix].begin(),
				       D.link_oriented[ix].end(),
				       il);
	const auto ity_ell = std::find(D.link_oriented[iy].begin(),
				       D.link_oriented[iy].end(),
				       il);
	const int ixl = std::distance(D.link_oriented[ix].begin(), itx_ell);
	const int iyl = std::distance(D.link_oriented[iy].begin(), ity_ell);

	double alpha_x = D.alpha[ix][ixl];
	double alpha_y = D.alpha[iy][iyl] + M_PI;
	double diff = alpha_x - (alpha_y+D.omega[il]);
	std::cout << "diff = " << double_mod(diff) << std::endl;

	// ----------------------------

      }
      std::cout << "ia = " << ia
		<< ", sum (directed) = " << double_mod( D.face_signs[ia]*omega_sum ) << std::endl;
      
      // ------------------------

    }
  }



  // std::cout << "# x y z" << std::endl;
  // std::cout << std::scientific << std::setprecision(15);

  // {
  //   auto vec = lattice.r[ix];
  //   for(auto elem : vec) {
  //     std::cout << std::setw(25) << elem << " ";
  //   }
  //   std::cout << std::endl;
  //   std::cout << projection(vec).format(fmt).transpose() << std::endl;
  // }
  // {
  //   auto vec = lattice.r[iy];
  //   for(auto elem : vec) {
  //     std::cout << std::setw(25) << elem << " ";
  //   }
  //   std::cout << std::endl;
  //   std::cout << projection(vec).format(fmt).transpose() << std::endl;
  // }


  return 0;
}
