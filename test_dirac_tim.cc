#include <iostream>
#include <fstream>
#include <cstdlib>

#include "dirac_s2.h"
#include "dirac_s2_tim.h"

int main(int argc, char* argv[]){

  // geometry
  const int q=5; // icosahedron
  const int n_refine=2;

  QfeLatticeS2 lattice(q, n_refine);
  Dirac1fonS2 D0(lattice, n_refine);

  Dirac1fonS2Tim D;

  // ----------------------------------

  auto mat = D.matrix_form_alternative(D0);
  // std::cout << mat << std::endl;

  auto ev = mat.eigenvalues();
  std::cout << "# spec = " << std::endl;
  for(int i=0; i<ev.size(); i++){
    std::cout << ev[i].real() << " " << ev[i].imag() << std::endl;
  }

  // ----------------------------------

  // {
  //   std::vector<int> Evan2Tim;
  //   for(int ix=0; ix<lattice.n_sites; ix++){
  //     for(int iy=0; iy<lattice.n_sites; iy++){
  // 	if(D.Tim2Evan[iy] == ix) {
  // 	  Evan2Tim.push_back(iy);
  // 	  break;
  // 	}
  //     }
  //   }

  //   std::cout << "Tim2Evan" << std::endl;
  //   for(auto elem : D.Tim2Evan) std::cout << elem << " ";
  //   std::cout << std::endl;
  //   std::cout << "Evan2Tim" << std::endl;
  //   for(auto elem : Evan2Tim) std::cout << elem << " ";
  //   std::cout << std::endl;

  //   {
  //     using Link = std::array<int,2>; // <int,int>;
  //     double TOL=1.0e-6;
  //     {
  // 	for(int ix=0; ix<lattice.n_sites; ix++){
  // 	  for(int jj=0; jj<lattice.sites[ix].nn; jj++){
  // 	    const int iy = lattice.sites[ix].neighbors[jj];
	  
  // 	    const double alpha1 = D0.alpha.at(Link{ix,iy});
  // 	    double alpha2 = D0.alpha.at(Link{iy,ix});
  // 	    double omega12 = D0.omega.at(Link{ix,iy});

  // 	    double diff = (alpha2 + M_PI + omega12) - alpha1;
  // 	    assert( std::abs(Mod(diff))<TOL );
  // 	  }}
  //     }
  //     {
  // 	for(int ix=0; ix<lattice.n_sites; ix++){
  // 	  for(int jj=0; jj<lattice.sites[ix].nn; jj++){
  // 	    const int iy = lattice.sites[ix].neighbors[jj];

  // 	    const int ixTim = Evan2Tim[ix];
  // 	    const int iyTim = Evan2Tim[iy];
  // 	    const int jjTim = D.find_jj( ixTim, iyTim );
  // 	    const int kkTim = D.find_jj( iyTim, ixTim );

  // 	    const double alpha1 = D.alpha[ixTim][jjTim];
  // 	    const double alpha2 = D.alpha[iyTim][kkTim];
  // 	    double omega12 = D.omega[ixTim][jjTim];
  // 	    double omega21 = D.omega[iyTim][kkTim];

  // 	    double diff1 = (alpha2 + M_PI + omega12) - alpha1;
  // 	    double diff2 = omega12 + omega21;
  // 	    assert( std::abs(Mod(diff1))<TOL );
  // 	    assert( std::abs(Mod(diff2))<TOL );
  // 	  }}
  //     }
  //     {
  // 	for(int ia=0; ia<lattice.n_faces; ia++){
  // 	  double sum = 0.0;
  // 	  double sum2 = 0.0;

  // 	  for(int i=0; i<3; i++){
  // 	    int ix = lattice.faces[ia].sites[i];
  // 	    int iy = lattice.faces[ia].sites[(i+1)%3];

  // 	    const int ixTim = Evan2Tim[ix];
  // 	    const int iyTim = Evan2Tim[iy];
  // 	    const int jjTim = D.find_jj( ixTim, iyTim );
  // 	    const int kkTim = D.find_jj( iyTim, ixTim );

  // 	    const double alpha1 = D.alpha[ixTim][jjTim];
  // 	    const double alpha2 = D.alpha[iyTim][kkTim];
  // 	    double omega12 = D.omega[ixTim][jjTim];
  // 	    double omega21 = D.omega[iyTim][kkTim];

  // 	    sum -= omega12;
  // 	    sum += alpha1;
  // 	    sum -= alpha2 + M_PI;

  // 	    sum2 -= omega12;
  // 	  }

  // 	  while(sum2>M_PI) sum2 -= 2.0*M_PI;
  // 	  while(sum2<-M_PI) sum2 += 2.0*M_PI;
  // 	  std::cout << "sum2(1) : " << D0.face_signs[ia] * sum2 << std::endl;
  // 	  assert( std::abs(Mod(-std::abs(Mod(sum)))) < TOL );
  // 	}
  // 	{
  // 	  for(int ia=0; ia<lattice.n_faces; ia++){
  // 	    double sum = 0.0;
  // 	    double sum2 = 0.0;

  // 	    for(int i=0; i<3; i++){
  // 	      int ix = lattice.faces[ia].sites[i];
  // 	      int iy = lattice.faces[ia].sites[(i+1)%3];
  // 	      sum -= D0.omega.at(Link{ix,iy});
  // 	      sum += D0.alpha.at(Link{ix,iy});
  // 	      sum -= D0.alpha.at(Link{iy,ix}) + M_PI;

  // 	      sum2 -= D0.omega.at(Link{ix,iy});
  // 	    }

  // 	    while(sum2>M_PI) sum2 -= 2.0*M_PI;
  // 	    while(sum2<-M_PI) sum2 += 2.0*M_PI;
  // 	    std::cout << "sum2(2) : " << D0.face_signs[ia] * sum2 << std::endl;
  // 	    assert( std::abs(Mod(-std::abs(Mod(sum)))) < TOL );
  // 	  }
  // 	}
  //     }
  //     {
  // 	for(int ix=0; ix<lattice.n_sites; ix++){
  // 	  const int ixTim = Evan2Tim[ix];

  // 	  const int iy0 = lattice.sites[ix].neighbors[0];
  // 	  const double alpha0 = D0.alpha.at(Link{ix,iy0});

  // 	  const int iy0Tim = Evan2Tim[iy0];
  // 	  const int jj0Tim = D.find_jj( ixTim, iy0Tim );
  // 	  const double alpha0Tim = D.alpha[ixTim][jj0Tim];

  // 	  for(int jj=1; jj<lattice.sites[ix].nn; jj++){
  // 	    const int iy1 = lattice.sites[ix].neighbors[jj];
  // 	    const double alpha1 = D0.alpha.at(Link{ix,iy1});

  // 	    const int iy1Tim = Evan2Tim[iy1];
  // 	    const int jj1Tim = D.find_jj( ixTim, iy1Tim );
  // 	    const double alpha1Tim = D.alpha[ixTim][jj1Tim];

  // 	    std::cout << alpha1 - alpha0 << std::endl;
  // 	    std::cout << alpha1Tim - alpha0Tim << std::endl;
  // 	    double diff = alpha1 - alpha0 - (alpha1Tim - alpha0Tim);
  // 	    double diff2 = std::abs(Mod(-std::abs(Mod(diff))));
  // 	    std::cout << diff2 << std::endl;
  // 	    assert( diff2<1.0e-5 );
  // 	    std::cout << std::endl;
  // 	  }
  // 	}
  //     }

  //   }

  // }

  // ----------------------------------

  return 0;
}



  // for(int ix=0; ix<lattice.n_sites; ix++){
  //   for(int jj=0; jj<lattice.sites[ix].nn; jj++){
  //     const int iy = lattice.sites[ix].neighbors[jj];
  //     auto mat1 = ( D.sigma[0] - D.gamma(ix, iy) ) * D.Omega(ix, iy);
  //     auto mat2 = D.Omega(ix, iy) * ( D.sigma[0] - D.gamma(iy, ix, M_PI) );
  //     std::cout << mat1-mat2 << std::endl;
  //   }}

  // ----------------------------------
