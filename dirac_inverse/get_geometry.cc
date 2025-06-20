#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include <algorithm>

#include <stdfloat>
// #include "geodesic2.h"
using Double = std::float64_t;
# include "geodesic.h"

#include "s2.h"

using namespace Geodesic;

using Idx = std::int32_t;
// using Complex = std::complex<double>;

using Link = std::array<Idx,2>; // <Idx,Idx>;
using Face = std::vector<Idx>;

// using MS=Eigen::Matrix2cd;
// using VD=Eigen::Vector2d;
// using VE=Eigen::Vector3d;
// using VC=Eigen::VectorXcd;
// using MS=Eigen::Matrix<Complex, 2, 2>;
// using VD=Eigen::Matrix<Double, 2, 1>;
// using VE=Eigen::Matrix<Double, 3, 1>;
using VD=V2;
using VE=V3;
// using VC=Eigen::Matrix<Complex, Eigen::Dynamic, 1>;


std::string dir = "/mnt/hdd_barracuda/qed3/dats/";

int main(int argc, char* argv[]){

  const int q=5; // icosahedron
  int n_refine=1;
  if(argc==2) n_refine = atoi(argv[1]);

  QfeLatticeS2 lattice(q, n_refine);

  // --------------------------
  // simplicial sites
  std::vector<VE> simp_sites;
  {
    for(auto& vec : lattice.r) {
      // VE site(std::vector<Double>({vec[0], vec[1], vec[2]}));
      VE site({vec[0], vec[1], vec[2]});
      simp_sites.push_back( site );
    }
  }

  // --------------------------
  // simp face & dual sites
  std::vector<Face> simp_faces;
  {
    int counter=0;
    for(auto& elem : lattice.faces) {
      Face face;
      for(int i=0; i<3; i++) face.push_back(elem.sites[i]);
      simp_faces.push_back( face );
      // for(int i=0; i<3; i++) std::cout << face[i] << " ";
      // std::cout << std::endl;
    }
  }
  // std::cout << simp_faces.size() << std::endl;
  // assert( simp_faces.size()==simp_sites.size() );

  std::vector<VE> dual_sites;
  std::vector<std::vector<Idx>> dual_faces_unsorted(simp_sites.size());
  {
    // for( auto& elem : simp_faces ){
    Idx ix=0; // dual site label=simp face label

    // std::cout << "debug. check1" << std::endl;
    for( const Face& face : simp_faces ){
      // const QfeFace& face = lattice.faces[iA];
      Vec3 r0, r1, r2; // r0,1: link
      // if(face.dual_sites[0]==link.dual_sites[0] && face.dual_sites[1]==link.dual_sites[1]){
      r0 = simp_sites[face[0]];
      r1 = simp_sites[face[1]];
      r2 = simp_sites[face[2]];

      dual_faces_unsorted[face[0]].push_back(ix);
      dual_faces_unsorted[face[1]].push_back(ix);
      dual_faces_unsorted[face[2]].push_back(ix);

      //
      const Vec3 p = circumcenter(r0, r1, r2).transpose();
      assert( std::abs( (p-r0).norm() - (p-r1).norm() )<1.0e-14 );
      assert( std::abs( (p-r0).norm() - (p-r2).norm() )<1.0e-14 );
      Vec3 tmp = p/p.norm();
      assert( std::abs( tmp.norm()-1.0 ) < 1.0e-14 );

      // Eigen::Matrix<Idx,1,3> face = Eigen::Map<Eigen::Matrix<Idx,1,3> >(elem.data());
      // VE tmp = VE::Zero();
      // for(int i=0; i<3; i++) tmp += simp_sites[face[i]];
      // tmp /= tmp.norm();
      dual_sites.push_back(tmp);
      ix++;
      // std::cout << "debug. check2" << std::endl;
      // std::cout << tmp.transpose() << std::endl;
    }
  }
  const int n_dual_sites = dual_sites.size();
  std::cout << "n_dual_sites = " << n_dual_sites << std::endl;

  // --------------------------
  // dual links
  std::vector<Link> dual_links;
  const double threshold = 0.42188 * 2.0 / n_refine;
  {
    for(Idx i=0; i<n_dual_sites; i++){
      Pt x( dual_sites[i] );
      for(Idx j=i+1; j<n_dual_sites; j++){
        Pt y( dual_sites[j] );
        const double ell = geodesicLength( x, y );
        // std::cout << "ell = " << ell << std::endl;
        if(ell<threshold) {
          // std::cout << i << " " << j << " " << ell << std::endl;
          Idx min, max;
          if(i<j) { min=i; max=j; }
          else { min=j; max=i; }
          dual_links.push_back( Link{min,max} );
        }
      }
    }
  }
  const int n_dual_links = dual_links.size();
  std::cout << "n_dual_links = " << n_dual_links << std::endl;

  // --------------------------
  // nearest neighbor
  std::vector<std::vector<Idx>> dual_nns;
  {
    for(Idx i=0; i<n_dual_sites; i++){
      Pt x( dual_sites[i] );
      std::vector<Idx> nn;
      for(Idx j=0; j<n_dual_sites; j++){
        if(i==j) continue;
        Pt y( dual_sites[j] );
        const double ell = geodesicLength( x, y );
        if(ell<threshold) nn.push_back(j);
      }
      assert( nn.size()==3 ); // this is dual nn
      dual_nns.push_back( nn );
    }
  }

  // --------------------------
  // dual faces
  std::vector<std::vector<Idx>> dual_faces;
  {
    // std::vector<std::vector<Idx>> dual_faces_unsorted;
    // for(const VE& cen : simp_dual_sites){
    //   std::vector<Idx> dual_face_unsorted;
    //   for(Idx ix=0; ix<n_dual_sites; ix++){
    //     const VE x = dual_sites[ix];
    //     if(geodesicLength(Pt(cen), Pt(x))<threshold) dual_face_unsorted.push_back(ix);
    //   }
    //   dual_faces_unsorted.push_back( dual_face_unsorted );
    // }

    const int nf = dual_faces_unsorted.size();

    std::vector<std::vector<Link>> facelinks;
    for(const auto& dual_face_unsorted : dual_faces_unsorted){

      std::vector<Link> facelink;
      // std::cout << "debug. dual_face_unsorted.size() = " << dual_face_unsorted.size() << std::endl;
      for(Idx i=0; i<dual_face_unsorted.size(); i++){ const Idx ix = dual_face_unsorted[i];
        for(Idx j=i+1; j<dual_face_unsorted.size(); j++){ const Idx jx = dual_face_unsorted[j];
          Link link({ix,jx});
          if(std::find(dual_links.begin(), dual_links.end(), link) != dual_links.end()) {
            facelink.push_back(link);
            continue;
          }
          link = Link{jx,ix};
          if(std::find(dual_links.begin(), dual_links.end(), link) != dual_links.end()) {
            facelink.push_back(link);
            continue;
          }
          // assert(false);
        }
      }
      facelinks.push_back(facelink);
      // std::cout << "debug. facelink.size() = " << facelink.size() << std::endl;
    }

    std::vector<std::vector<Link>> list_facelinkordered;
    std::vector<std::vector<int>> list_facelinksign;
    std::vector<std::vector<Link>> list_facelinkorderedsigned;

    int counter=0;
    // facelink
    for(auto& facelink : facelinks){
      std::vector<Link> facelinkordered;
      std::vector<int> facelinksign;
      std::vector<Link> facelinkorderedsigned;

      Idx il=0;
      // std::cout << "facelink.size() = "<< facelink.size() << std::endl;
      while(true){
        // std::cout << "il = " << il << std::endl;
        Link fl = facelink[il];
        facelinkordered.push_back(fl);
        facelink.erase(facelink.begin()+il);

        VE x0 = simp_sites[counter];
        VE x1 = dual_sites[fl[0]];
        VE x2 = dual_sites[fl[1]];

        int sign = 1;
        VE x10 = x1-x0;
        VE x21 = x2-x1;
        VE tmp = x10.cross(x21);
        Double tmp2 = tmp.dot(x0);
        if(tmp2 < 0.0) sign = -1;
        facelinksign.push_back( sign );

        Link sgnd({fl[0], fl[1]});
        if(sign<0) sgnd = Link{fl[1], fl[0]};
        facelinkorderedsigned.push_back( sgnd );

        if(facelink.size()==0) break;
        const Idx next = sgnd[1];
        for(Idx jl=0; jl<facelink.size(); jl++){
          Link tmp = facelink[jl];
          if(std::find(tmp.begin(), tmp.end(), next) != tmp.end()) il = jl;
        }
      }
      counter++;

      list_facelinkordered.push_back(facelinkordered);
      list_facelinksign.push_back(facelinksign);
      list_facelinkorderedsigned.push_back(facelinkorderedsigned);
    }

    for(const auto& dual_links : list_facelinkorderedsigned){
      std::vector<Idx> face;
      for(Idx il=0; il<dual_links.size(); il++) {
        // std::cout << dual_links[il][0] << " ";
        face.push_back( dual_links[il][0] );
      }
      // std::cout << std::endl;
      dual_faces.push_back(face);
    }
  }


  {
    std::ofstream ofs(dir+"pts_n"+std::to_string(n_refine)+"_singlepatch.dat");
    ofs << std::scientific << std::setprecision(25);
    for(const auto& vec : simp_sites) {
      for(const auto& elem : vec) {
        ofs << std::setw(50) << elem << " ";
      }
      ofs << std::endl;
    }
  }
  {
    std::ofstream ofs(dir+"pts_dual_n"+std::to_string(n_refine)+"_singlepatch.dat");
    ofs << std::scientific << std::setprecision(25);
    for(const auto& vec : dual_sites) {
      for(const auto& elem : vec) {
        ofs << std::setw(50) << elem << " ";
      }
      ofs << std::endl;
    }
  }
  {
    std::ofstream ofs(dir+"nns_dual_n"+std::to_string(n_refine)+"_singlepatch.dat");
    ofs << std::scientific << std::setprecision(25);
    for(const auto& vec : dual_nns) {
      for(const auto& elem : vec) {
        ofs << std::setw(50) << elem << " ";
      }
      ofs << std::endl;
    }
  }
  {
    std::ofstream ofs(dir+"dual_links_n"+std::to_string(n_refine)+"_singlepatch.dat");
    ofs << std::scientific << std::setprecision(25);
    for(const auto& vec : dual_links) {
      for(const auto& elem : vec) {
        ofs << std::setw(50) << elem << " ";
      }
      ofs << std::endl;
    }
  }

  // vs, us, dualtriangleareas
  {
    std::ofstream ofs(dir+"face_n"+std::to_string(n_refine)+"_singlepatch.dat");
    ofs << std::scientific << std::setprecision(25);
    for(const auto& vec : simp_faces) {
      for(const auto& elem : vec) {
        ofs << std::setw(50) << elem << " ";
      }
      ofs << std::endl;
    }
  }
  {
    std::ofstream ofs(dir+"face_dual_n"+std::to_string(n_refine)+".dat");
    ofs << std::scientific << std::setprecision(25);
    for(const auto& vec : dual_faces) {
      for(const auto& elem : vec) {
        ofs << std::setw(50) << elem << " ";
      }
      ofs << std::endl;
    }
  }
  // {
  //   std::ofstream ofs(dir+"facesign_dual_n"+std::to_string(n_refine)+".dat");
  //   ofs << std::scientific << std::setprecision(25);
  //   for(const auto& vec : faces) {
  //     for(const auto& elem : vec) {
  //       ofs << std::setw(50) << elem << " ";
  //     }
  //     ofs << std::endl;
  //   }
  // }



  return 0;
}


  // {
  //   int counter=0;
  //   std::ofstream ofs("nearest_neighbor_table.dat");
  //   ofs << "# ix iy il" << std::endl;
  //   ofs << std::scientific << std::setprecision(25);
  //   for(int ix=0; ix<lattice.n_dual_sites; ix++) {
  //     const auto x = lattice.dual_sites[ix];
  //     for(int iw=0; iw<x.nn; iw++){
  // 	ofs << std::setw(10) << ix << " ";
  // 	ofs << std::setw(10) << x.neighbors[iw] << " ";
  // 	ofs << std::setw(10) << x.dual_links[iw] << " ";
  // 	ofs << std::endl;
  //     }
  //   }
  // }
