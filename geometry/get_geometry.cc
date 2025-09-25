#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include <algorithm>
#include <filesystem>

#include <stdfloat>
using Double = std::float64_t;
# include "geodesic.h"

#include "s2.h"

using namespace Geodesic;

using Idx = std::int32_t;

using Link = std::array<Idx,2>; // <Idx,Idx>;
using Face = std::vector<Idx>;

using VD=V2;
using VE=V3;


std::string dir = "./data/";

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
    }
  }

  std::vector<VE> dual_sites;
  std::vector<std::vector<Idx>> dual_faces_unsorted(simp_sites.size());
  {
    Idx ix=0; // dual site label=simp face label

    for( const Face& face : simp_faces ){
      Vec3 r0, r1, r2; // r0,1: link
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

      dual_sites.push_back(tmp);
      ix++;
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
        if(ell<threshold) {
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
  // simp links
  std::vector<Link> simp_links;
  {
    const double threshold=0.42188 * 2.0 / n_refine;

    for(const Link& link : dual_links){ // trivalent
      const VE x1 = dual_sites[link[0]];
      const VE x2 = dual_sites[link[1]];

      std::vector<Idx> tmp;
      for(Idx ip=0; ip<simp_sites.size(); ip++){
        const VE x0 = simp_sites[ip];
        const double d01 = (x0-x1).norm();
        const double d02 = (x0-x2).norm();
        if(d01<threshold && d02<threshold) tmp.push_back(ip);
      }
      assert( tmp.size()==2 );

      const Idx min = std::min(tmp[0], tmp[1]);
      const Idx max = std::max(tmp[0], tmp[1]);
      simp_links.push_back( Link{min,max} );
    }
  }
  assert( simp_links.size()==dual_links.size() );

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
      dual_nns.push_back( nn );
    }
  }

  // --------------------------
  // SIMP nearest neighbor
  std::vector<std::vector<Idx>> simp_nns;
  {
    for(Idx i=0; i<lattice.n_sites; i++){
      std::vector<Idx> tmp;
      for(int jj=0; jj<lattice.sites[i].nn; jj++) tmp.push_back( lattice.sites[i].neighbors[jj] );
      simp_nns.push_back(tmp);
    }
  }


  // --------------------------
  // dual faces
  std::vector<std::vector<Idx>> dual_faces;
  {
    const int nf = dual_faces_unsorted.size();

    std::vector<std::vector<Link>> facelinks;
    for(const auto& dual_face_unsorted : dual_faces_unsorted){

      std::vector<Link> facelink;
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
        }
      }
      facelinks.push_back(facelink);
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
      while(true){
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
        face.push_back( dual_links[il][0] );
      }
      dual_faces.push_back(face);
    }
  }

  std::filesystem::create_directory(dir);

  {
    std::ofstream ofs(dir+"pts_n"+std::to_string(n_refine)+".dat");
    ofs << std::scientific << std::setprecision(25);
    for(const auto& vec : simp_sites) {
      for(const auto& elem : vec) {
        ofs << std::setw(50) << elem << " ";
      }
      ofs << std::endl;
    }
  }
  {
    std::ofstream ofs(dir+"pts_dual_n"+std::to_string(n_refine)+".dat");
    ofs << std::scientific << std::setprecision(25);
    for(const auto& vec : dual_sites) {
      for(const auto& elem : vec) {
        ofs << std::setw(50) << elem << " ";
      }
      ofs << std::endl;
    }
  }
  {
    std::ofstream ofs(dir+"nns_n"+std::to_string(n_refine)+".dat");
    ofs << std::scientific << std::setprecision(25);
    for(const auto& vec : simp_nns) {
      for(const auto& elem : vec) {
        ofs << std::setw(50) << elem << " ";
      }
      ofs << std::endl;
    }
  }
  {
    std::ofstream ofs(dir+"nns_dual_n"+std::to_string(n_refine)+".dat");
    ofs << std::scientific << std::setprecision(25);
    for(const auto& vec : dual_nns) {
      for(const auto& elem : vec) {
        ofs << std::setw(50) << elem << " ";
      }
      ofs << std::endl;
    }
  }
  {
    std::ofstream ofs(dir+"dual_links_n"+std::to_string(n_refine)+".dat");
    ofs << std::scientific << std::setprecision(25);
    for(const auto& vec : dual_links) {
      for(const auto& elem : vec) {
        ofs << std::setw(50) << elem << " ";
      }
      ofs << std::endl;
    }
  }
  {
    std::ofstream ofs(dir+"links_n"+std::to_string(n_refine)+".dat");
    ofs << std::scientific << std::setprecision(25);
    for(const auto& vec : simp_links) {
      for(const auto& elem : vec) {
        ofs << std::setw(50) << elem << " ";
      }
      ofs << std::endl;
    }
  }

  // vs, us, dualtriangleareas
  {
    std::ofstream ofs(dir+"face_n"+std::to_string(n_refine)+".dat");
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



  return 0;
}


