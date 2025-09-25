#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <map>
#include <limits>
#include <stdfloat>

#include <Eigen/Dense>

using Double = std::float64_t;
#include "geodesic.h"
#include "integral.h"


using namespace Geodesic;

using Idx = std::int32_t;

using Link = std::array<Idx,2>; // <Idx,Idx>;
using Face = std::vector<Idx>;
using VD=Eigen::Matrix<Double, 2, 1>;
using VE=Eigen::Matrix<Double, 3, 1>;



Double TOL=1.0e-14;

constexpr Double TOLLOOSE=1.0e-6;

const int limit = 4000; // 1000;
Double epsabs = 1.0e-15; // 0.;
Double epsrel = 1.0e-13; // TOLLOOSE;
int key = 5;


std::string dir = "./data/";


int main(int argc, char* argv[]){
  int n_refine=1;
  if(argc==2) n_refine = atoi(argv[1]);

  std::vector<VE> simp_sites;
  {
    std::ifstream file(dir+"pts_n"+std::to_string(n_refine)+".dat");

    std::string str;
    // skip the header line
    // https://stackoverflow.com/questions/16068997/how-to-skip-the-header-line-in-a-text-file-and-read-back-the-rest-of-the-data-to
    // file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    while (std::getline(file, str)){
      std::istringstream iss(str);
      Double v1, v2, v3;
      iss >> v1;
      iss >> v2;
      iss >> v3;
      simp_sites.push_back( VE(v1, v2, v3) );
    }
  }

  std::vector<VE> sites;
  {
    std::ifstream file(dir+"pts_dual_n"+std::to_string(n_refine)+".dat");

    std::string str;
    // skip the header line
    // https://stackoverflow.com/questions/16068997/how-to-skip-the-header-line-in-a-text-file-and-read-back-the-rest-of-the-data-to
    // file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    while (std::getline(file, str)){
      std::istringstream iss(str);
      Double v1, v2, v3;
      iss >> v1;
      iss >> v2;
      iss >> v3;
      sites.push_back( VE(v1, v2, v3) );
    }
  }

  std::vector<Link> links;
  {
    std::ifstream file(dir+"dual_links_n"+std::to_string(n_refine)+".dat");
    std::string str;
    while (std::getline(file, str)){
      std::istringstream iss(str);
      Idx v1, v2;
      iss >> v1;
      iss >> v2;
      links.push_back( Link({v1, v2}) );
    }
  }

  std::vector<std::vector<Idx>> nns;
  {
    std::ifstream file(dir+"nns_dual_n"+std::to_string(n_refine)+".dat");
    std::string str;
    while (std::getline(file, str)){
      std::istringstream iss(str);
      Idx v1, v2, v3;
      iss >> v1;
      iss >> v2;
      iss >> v3;
      nns.push_back( std::vector<Idx>{v1,v2,v3} );
    }
  }

  std::vector<Face> faces;
  {
    std::ifstream file(dir+"face_dual_n"+std::to_string(n_refine)+".dat");
    assert(file.is_open());
    std::string str;
    while (std::getline(file, str)){
      std::istringstream iss(str);
      Idx v;
      std::vector<Idx> face;
      while( iss >> v ) face.push_back( v );
      faces.push_back( face );
    }
  }


  std::vector<Double> omegas;
  std::vector<Sol> sols;
  for(const Link& link : links){

    Pt x1( sites[link[0]] );
    Pt x2( sites[link[1]] );

    Sol sol = SolveGeodesics( x1, x2 );
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (limit);
    double result, error;

    if(!sol.is_split){
      F f1 = [&](Double s){ return sol.Dphi(s) * std::cos( sol.theta(s) ); };
      gsl_function F;
      F.function = &unwrap;
      F.params = &f1;
      gsl_integration_qag(&F, 0., sol.ell, epsabs, epsrel, limit,
                          key, w, &result, &error);
    }
    else{
      double result1, result2;
      double error1, error2;
      {
        F f1 = [&](Double s){ return sol.Dphi(s, 0) * std::cos( sol.theta(s, 0) ); };
        gsl_function F;
        F.function = &unwrap;
        F.params = &f1;
        gsl_integration_qag(&F, 0.0, sol.sE, epsabs, epsrel, limit,
                            key, w, &result1, &error1);
      }
      {
        F f1 = [&](Double s){ return sol.Dphi(s, 1) * std::cos( sol.theta(s, 1) ); };
        gsl_function F;
        F.function = &unwrap;
        F.params = &f1;
        gsl_integration_qag(&F, sol.sE, sol.ell, epsabs, epsrel, limit,
                            key, w, &result2, &error2);
      }
      result = result1+result2;
      error=error1+error2;
    }
    omegas.push_back( result );
    sols.push_back( sol );
  }

  {
    const Double phi0 = _M_PI/48.0;
    for(const auto& site : sites) assert( std::abs(Mod(Pt(site).xi[1]) - phi0)>TOLLOOSE );
    std::vector<Idx> section;
    for(Idx il=0; il<links.size(); il++){
      Link link = links[il];
      const Double phi1 = Pt(sites[link[0]]).xi[1];
      const Double phi2 = Pt(sites[link[1]]).xi[1];
      Double dphi1 = phi1-phi0;
      Double dphi2 = phi2-phi0;
      Double dphi = phi1-phi2;

      if(dphi1 >= _M_PI) dphi1 -= 2.0*_M_PI;
      else if(dphi1 < -_M_PI) dphi1 += 2.0*_M_PI;
      if(dphi2 >= _M_PI) dphi2 -= 2.0*_M_PI;
      else if(dphi2 < -_M_PI) dphi2 += 2.0*_M_PI;
      if(dphi >= _M_PI) dphi -= 2.0*_M_PI;
      else if(dphi < -_M_PI) dphi += 2.0*_M_PI;

      if( std::abs(dphi1) + std::abs(dphi2) - std::abs(dphi) < TOLLOOSE ){
        section.push_back( il );
      }
    }
    for(Idx il : section) omegas[il] -= 2.0*_M_PI;
  }

  std::map<const Link, Double> omega;
  for(Idx il=0; il<links.size(); il++){
    Link link = links[il];
    omega.insert( { link, omegas[il] } );
    omega.insert( { Link{link[1], link[0]}, -omegas[il] } );
  }


  std::vector<Double> alpha0s, alpha1s;
  for( auto& sol : sols){
    Double alpha0, alpha1;
    if(!sol.is_split){
      const Double theta0 = sol.theta(0);
      const Double dtheta_ds0 = sol.Dtheta(0);
      const Double dphi_ds0 = sol.Dphi(0);
      assert( std::abs(dtheta_ds0*dtheta_ds0 + std::sin(theta0)*std::sin(theta0)*dphi_ds0*dphi_ds0-1.0)<TOLLOOSE );

      VD e0( dtheta_ds0, std::sin(theta0)*dphi_ds0 );
      alpha0 = _acos( e0[0] );
      if(e0[1]<0) alpha0 *= -1.0;
      assert( std::abs(std::cos(alpha0)-e0[0]) + std::abs(std::sin(alpha0)-e0[1]) < TOLLOOSE );

      const Double theta1 = sol.theta(sol.ell);
      const Double dtheta_ds1 = sol.Dtheta(sol.ell);
      const Double dphi_ds1 = sol.Dphi(sol.ell);
      assert( std::abs(dtheta_ds1*dtheta_ds1 + std::sin(theta1)*std::sin(theta1)*dphi_ds1*dphi_ds1-1.0)<TOLLOOSE );

      VD e1( dtheta_ds1, std::sin(theta1)*dphi_ds1 );
      alpha1 = _acos( e1[0] );
      if(e1[1]<0) alpha1 *= -1.0;
      assert( std::abs(std::cos(alpha1)-e1[0]) + std::abs(std::sin(alpha1)-e1[1]) < TOLLOOSE );
    }
    else{
      const Double theta0 = sol.theta(0, 0);
      const Double dtheta_ds0 = sol.Dtheta(0, 0);
      const Double dphi_ds0 = sol.Dphi(0, 0);
      assert( std::abs(dtheta_ds0*dtheta_ds0 + std::sin(theta0)*std::sin(theta0)*dphi_ds0*dphi_ds0-1.0)<TOLLOOSE );

      VD e0( dtheta_ds0, std::sin(theta0)*dphi_ds0 );
      alpha0 = _acos( e0[0] );
      if(e0[1]<0) alpha0 *= -1.0;
      assert( std::abs(std::cos(alpha0)-e0[0]) + std::abs(std::sin(alpha0)-e0[1]) < TOLLOOSE );

      const Double theta1 = sol.theta(sol.ell, 1);
      const Double dtheta_ds1 = sol.Dtheta(sol.ell, 1);
      const Double dphi_ds1 = sol.Dphi(sol.ell, 1);
      assert( std::abs(dtheta_ds1*dtheta_ds1 + std::sin(theta1)*std::sin(theta1)*dphi_ds1*dphi_ds1-1.0)<TOLLOOSE );

      VD e1( dtheta_ds1, std::sin(theta1)*dphi_ds1 );
      alpha1 = _acos( e1[0] );
      if(e1[1]<0) alpha1 *= -1.0;
      assert( std::abs(std::cos(alpha1)-e1[0]) + std::abs(std::sin(alpha1)-e1[1]) < TOLLOOSE );
    }
    alpha0s.push_back( alpha0 );
    alpha1s.push_back( Mod(alpha1 + _M_PI) );
  }

  std::map<const Link, Double> alpha;
  for(Idx il=0; il<links.size(); il++){
    Link link = links[il];
    alpha.insert( { link, alpha0s[il] } );
    alpha.insert( { Link{link[1],link[0]}, alpha1s[il] } );
  }

  for(Idx ix=0; ix<sites.size(); ix++){
    const V3 e0x = Pt( sites[ix] ).e0();
    const V3 e1x = Pt( sites[ix] ).e1();
    for(Idx iy : nns[ix]){
      const V3 diff = sites[iy]-sites[ix];
      const Double c0 = e0x.dot(diff);
      const Double c1 = e1x.dot(diff);
      Double al = std::atan(c1/c0);
      const int br = decide_branch( al-alpha.at(Link({ix, iy})) );
      al -= _M_PI*br;
      alpha[Link({ix, iy})] = al;
    }
  }


  {
    std::clog << "# checking spin structure" << std::endl;
    for(Idx ix=0; ix<sites.size(); ix++){
      for(Idx iy : nns[ix]){
        const Double alpha1 = alpha.at(Link{ix,iy});
        Double alpha2 = alpha.at(Link{iy,ix});
        Double omega12 = omega.at(Link{ix,iy});

        Double diff = (alpha2 + _M_PI + omega12) - alpha1;
        assert( isModdable(diff, TOL) );

        Double om = alpha1 - (alpha2 + _M_PI);
        const int br = decide_branch( om-omega12 );
        om -= _M_PI*br;
        omega[Link({ix, iy})] = om;
      }}
  }

  {
    std::clog << "# checking deficits" << std::endl;
    int counter=0;
    for(auto& face : faces){
      int sign = 1;
      {
        VE x0 = simp_sites[counter];
        VE x1 = sites[face[0]];
        VE x2 = sites[face[1]];
        if((x1-x0).cross(x2-x1).dot(x0) < 0) sign = -1;
      }

      Double sum = 0.0;
      for(int i=0; i<face.size(); i++){
        const Idx ix = face[i];
        const Idx j = (i+1)%face.size();
        const Idx jx = face[j];
        sum += omega.at( Link{ix,jx} );
      }
      sum *= sign;
      Double mod = Mod(sum, 4.0*_M_PI);
      if(mod>2.0*_M_PI) mod -= 4.0*_M_PI;
      assert( -1.5 * 4.0*_M_PI/face.size() < mod && mod < 0.0 );
      counter++;
    }
  }



  {
    std::ofstream ofs(dir+"omega_dual_n"+std::to_string(n_refine)+".dat");
    ofs << std::scientific << std::setprecision(15);
    for(Idx il=0; il<links.size(); il++){
      Link link = links[il];
      ofs << std::setw(25) << link[0] << " ";
      ofs << std::setw(25) << link[1] << " ";
      ofs << std::setw(25) << omegas[il] << std::endl;
    }
  }

  {
    std::ofstream ofs(dir+"alpha_dual_n"+std::to_string(n_refine)+".dat");
    ofs << std::scientific << std::setprecision(15);
    for(Idx il=0; il<links.size(); il++){
      Link link = links[il];
      ofs << std::setw(25) << link[0] << " ";
      ofs << std::setw(25) << link[1] << " ";
      ofs << std::setw(25) << alpha0s[il] << std::endl;

      ofs << std::setw(25) << link[1] << " ";
      ofs << std::setw(25) << link[0] << " ";
      ofs << std::setw(25) << alpha1s[il] << std::endl;
    }
  }

  // dualtriangleareas
  std::vector<Face> simp_faces;
  {
    std::ifstream file(dir+"face_n"+std::to_string(n_refine)+".dat");

    std::string str;
    while (std::getline(file, str)){
      std::istringstream iss(str);
      Idx v;
      Face face;
      while( iss >> v ) {
        face.push_back( v );
      }
      simp_faces.push_back( face );
    }
  }
  std::vector<Double> triangleareas;
  for(const Face& face : simp_faces){
    Double area = sphericalarea( Pt(simp_sites[face[0]]),
                                 Pt(simp_sites[face[1]]),
                                 Pt(simp_sites[face[2]]));
    triangleareas.push_back( area );
  }
  assert( simp_faces.size()==triangleareas.size() );
  std::cout << sites.size() << std::endl;
  std::cout << simp_faces.size() << std::endl;
  std::cout << triangleareas.size() << std::endl;
  assert( sites.size()==triangleareas.size() );

  {
    std::ofstream ofs(dir+"dualtriangleareas_n"+std::to_string(n_refine)+".dat");
    ofs << std::scientific << std::setprecision(15);
    for(Double elem : triangleareas) ofs << std::setw(25) << elem << std::endl;
  }



  return 0;
}

