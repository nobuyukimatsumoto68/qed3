#include <iostream>
#include <iomanip>
#include <fstream>

#include <algorithm>

#include "geodesic.h"
#include "s2.h"

using Idx = std::int32_t;
// using Complex = std::complex<double>;

using Link = std::array<Idx,2>; // <Idx,Idx>;
using Face = std::vector<Idx>;

// using MS=Eigen::Matrix2cd;
// using VD=Eigen::Vector2d;
// using VE=Eigen::Vector3d;
// using VC=Eigen::VectorXcd;
// using MS=Eigen::Matrix<Complex, 2, 2>;
using VD=Eigen::Matrix<Double, 2, 1>;
using VE=Eigen::Matrix<Double, 3, 1>;
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
      VE site;
      site << vec[0], vec[1], vec[2];
      simp_sites.push_back( site );
    }
  }

  // --------------------------
  // dual sites
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
  std::cout << simp_faces.size() << std::endl;
  // assert( simp_faces.size()==simp_sites.size() );

  std::vector<VE> sites;
  {
    for( auto& elem : simp_faces ){
      Eigen::Matrix<Idx,1,3> face = Eigen::Map<Eigen::Matrix<Idx,1,3> >(elem.data());
      VE tmp = VE::Zero();
      for(int i=0; i<3; i++) tmp += simp_sites[face[i]];
      tmp /= tmp.norm();
      sites.push_back(tmp);
      // std::cout << tmp.transpose() << std::endl;
    }
  }
  const int n_sites = sites.size();

  // --------------------------
  // dual links
  std::vector<Link> links;
  const double threshold = 0.42188 * 2.0 / n_refine;
  {
    for(Idx i=0; i<n_sites; i++){
      Pt x( sites[i] );
      for(Idx j=i+1; j<n_sites; j++){
        Pt y( sites[j] );
        const double ell = geodesicLength( x, y );
        if(ell<threshold) {
          // std::cout << i << " " << j << " " << ell << std::endl;
          Idx min, max;
          if(i<j) { min=i; max=j; }
          else { min=j; max=i; }
          links.push_back( Link{min,max} );
        }
      }
    }
  }
  const int n_links = links.size();

  // --------------------------
  // nearest neighbor
  std::vector<std::vector<Idx>> nns;
  {
    for(Idx i=0; i<n_sites; i++){
      Pt x( sites[i] );
      std::vector<Idx> nn;
      for(Idx j=0; j<n_sites; j++){
        if(i==j) continue;
        Pt y( sites[j] );
        const double ell = geodesicLength( x, y );
        if(ell<threshold) nn.push_back(j);
      }
      nns.push_back( nn );
    }
  }

  // --------------------------
  // dual faces
  std::vector<std::vector<Idx>> faces;
  {
    std::vector<std::vector<Idx>> dual_faces_unsorted;
    for(const VE& cen : simp_sites){
      std::vector<Idx> dual_face_unsorted;
      for(Idx ix=0; ix<n_sites; ix++){
        const VE x = sites[ix];
        if(geodesicLength(Pt(cen), Pt(x))<threshold) dual_face_unsorted.push_back(ix);
      }
      dual_faces_unsorted.push_back( dual_face_unsorted );
    }

    const int nf = dual_faces_unsorted.size();

    std::vector<std::vector<Link>> facelinks;
    for(const auto& dual_face_unsorted : dual_faces_unsorted){

      std::vector<Link> facelink;
      for(int i=0; i<dual_face_unsorted.size(); i++){ const Idx ix = dual_face_unsorted[i];
        for(int j=i+1; j<dual_face_unsorted.size(); j++){ const Idx jx = dual_face_unsorted[j];
          Link link({ix,jx});
          if(std::find(links.begin(), links.end(), link) != links.end()) facelink.push_back(link);
        }
      }
      facelinks.push_back(facelink);
      // std::cout << facelink.size() << std::endl;
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

      int il=0;
      while(true){
        Link fl = facelink[il];
        facelinkordered.push_back(fl);
        facelink.erase(facelink.begin()+il);

        VE x0 = simp_sites[counter];
        VE x1 = sites[fl[0]];
        VE x2 = sites[fl[1]];

        int sign = 1;
        if((x1-x0).cross(x2-x1).dot(x0) < 0) sign = -1;
        facelinksign.push_back( sign );

        Link sgnd({fl[0], fl[1]});
        if(sign<0) sgnd = Link{fl[1], fl[0]};
        facelinkorderedsigned.push_back( sgnd );

        if(facelink.size()==0) break;
        const Idx next = sgnd[1];
        for(int jl=0; jl<facelink.size(); jl++){
          Link tmp = facelink[jl];
          if(std::find(tmp.begin(), tmp.end(), next) != tmp.end()) il = jl;
        }
      }
      counter++;

      list_facelinkordered.push_back(facelinkordered);
      list_facelinksign.push_back(facelinksign);
      list_facelinkorderedsigned.push_back(facelinkorderedsigned);
    }

    for(const auto& links : list_facelinkorderedsigned){
      std::vector<Idx> face;
      for(int il=0; il<links.size(); il++) {
        std::cout << links[il][0] << " ";
        face.push_back( links[il][0] );
      }
      std::cout << std::endl;
      faces.push_back(face);
    }
  }


  {
    std::ofstream ofs(dir+"pts_n"+std::to_string(n_refine)+"_singlepatch.dat");
    ofs << std::scientific << std::setprecision(15);
    for(const auto& vec : simp_sites) {
      for(const auto& elem : vec) {
        ofs << std::setw(25) << elem << " ";
      }
      ofs << std::endl;
    }
  }
  {
    std::ofstream ofs(dir+"pts_dual_n"+std::to_string(n_refine)+"_singlepatch.dat");
    ofs << std::scientific << std::setprecision(15);
    for(const auto& vec : sites) {
      for(const auto& elem : vec) {
        ofs << std::setw(25) << elem << " ";
      }
      ofs << std::endl;
    }
  }
  {
    std::ofstream ofs(dir+"nns_dual_n"+std::to_string(n_refine)+"_singlepatch.dat");
    ofs << std::scientific << std::setprecision(15);
    for(const auto& vec : nns) {
      for(const auto& elem : vec) {
        ofs << std::setw(25) << elem << " ";
      }
      ofs << std::endl;
    }
  }
  {
    std::ofstream ofs(dir+"dual_links_n"+std::to_string(n_refine)+"_singlepatch.dat");
    ofs << std::scientific << std::setprecision(15);
    for(const auto& vec : links) {
      for(const auto& elem : vec) {
        ofs << std::setw(25) << elem << " ";
      }
      ofs << std::endl;
    }
  }

  // vs, us, dualtriangleareas
  {
    std::ofstream ofs(dir+"face_n"+std::to_string(n_refine)+"_singlepatch.dat");
    ofs << std::scientific << std::setprecision(15);
    for(const auto& vec : simp_faces) {
      for(const auto& elem : vec) {
        ofs << std::setw(25) << elem << " ";
      }
      ofs << std::endl;
    }
  }
  {
    std::ofstream ofs(dir+"face_dual_n"+std::to_string(n_refine)+"_singlepatch.dat");
    ofs << std::scientific << std::setprecision(15);
    for(const auto& vec : faces) {
      for(const auto& elem : vec) {
        ofs << std::setw(25) << elem << " ";
      }
      ofs << std::endl;
    }
  }
  // {
  //   std::ofstream ofs(dir+"facesign_dual_n"+std::to_string(n_refine)+".dat");
  //   ofs << std::scientific << std::setprecision(15);
  //   for(const auto& vec : faces) {
  //     for(const auto& elem : vec) {
  //       ofs << std::setw(25) << elem << " ";
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
  //   ofs << std::scientific << std::setprecision(15);
  //   for(int ix=0; ix<lattice.n_sites; ix++) {
  //     const auto x = lattice.sites[ix];
  //     for(int iw=0; iw<x.nn; iw++){
  // 	ofs << std::setw(10) << ix << " ";
  // 	ofs << std::setw(10) << x.neighbors[iw] << " ";
  // 	ofs << std::setw(10) << x.links[iw] << " ";
  // 	ofs << std::endl;
  //     }
  //   }
  // }
