#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <map>
#include <limits>

#include <Eigen/Dense>

// #include "geodesic.h"
// #include "integral.h"

using Idx = std::int32_t;
using Complex = std::complex<double>;
const Complex I = Complex(0.0, 1.0);

#include "geodesic.h"
#include "integral.h"


using Link = std::array<Idx,2>; // <Idx,Idx>;
using Face = std::vector<Idx>;
using MS=Eigen::Matrix2cd;
using VD=Eigen::Vector2d;
using VE=Eigen::Vector3d;
using VC=Eigen::VectorXcd;


const int n_refine=12;

int main(int argc, char* argv[]){

  std::vector<VE> simp_sites;
  {
    std::ifstream file("dats/pts_n"+std::to_string(n_refine)+"_singlepatch.dat");

    std::string str;
    // skip the header line
    // https://stackoverflow.com/questions/16068997/how-to-skip-the-header-line-in-a-text-file-and-read-back-the-rest-of-the-data-to
    // file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    while (std::getline(file, str)){
      std::istringstream iss(str);
      double v1, v2, v3;
      iss >> v1;
      iss >> v2;
      iss >> v3;
      simp_sites.push_back( VE(v1, v2, v3) );
    }
  }

  std::vector<VE> sites;
  {
    std::ifstream file("dats/pts_dual_n"+std::to_string(n_refine)+"_singlepatch.dat");

    std::string str;
    // skip the header line
    // https://stackoverflow.com/questions/16068997/how-to-skip-the-header-line-in-a-text-file-and-read-back-the-rest-of-the-data-to
    // file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    while (std::getline(file, str)){
      std::istringstream iss(str);
      double v1, v2, v3;
      iss >> v1;
      iss >> v2;
      iss >> v3;
      sites.push_back( VE(v1, v2, v3) );
    }
  }

  std::vector<Link> links;
  {
    std::ifstream file("dats/dual_links_n"+std::to_string(n_refine)+"_singlepatch.dat");

    std::string str;
    // skip the header line
    // https://stackoverflow.com/questions/16068997/how-to-skip-the-header-line-in-a-text-file-and-read-back-the-rest-of-the-data-to
    // file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
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
    std::ifstream file("dats/nns_dual_n"+std::to_string(n_refine)+"_singlepatch.dat");
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
    std::ifstream file("dats/face_dual_n"+std::to_string(n_refine)+".dat");
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


  std::vector<double> omegas;
  std::vector<Sol> sols;
  // int counter=0;
  for(const Link& link : links){
    // {
    // Link link = links[74];
    // std::cout << "counter = " << counter << std::endl;
    // counter++;
    // Link link = links[2];

    Pt x1( sites[link[0]] );
    Pt x2( sites[link[1]] );

    // std::cout << "x1.x = " << x1.x.transpose() << std::endl;
    // std::cout << "x1.xi = " << x1.xi.transpose() << std::endl;
    // std::cout << "x2.x = " << x2.x.transpose() << std::endl;
    // std::cout << "x2.xi = " << x2.xi.transpose() << std::endl;
    // {
    //   I2 sign1, sign2;
    //   getSign( sign1, sign2, x1, x2 );
    //   std::cout << "sign1 = " << sign1.transpose() << std::endl;
    //   std::cout << "sign2 = " << sign2.transpose() << std::endl;
    // }

    Sol sol = SolveGeodesics( x1, x2 );
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    double result, error;
    if(!sol.is_split){
      F f1 = [&](double s){ return sol.Dphi(s) * std::cos( sol.theta(s) ); };
      gsl_function F;
      F.function = &unwrap;
      F.params = &f1;
      gsl_integration_qag(&F, 0., sol.ell, 0, 1e-7, 1000,
                          3, w, &result, &error);
    }
    else{
      double result1, result2;
      double error1, error2;
      {
        F f1 = [&](double s){ return sol.Dphi(s, 0) * std::cos( sol.theta(s, 0) ); };
        gsl_function F;
        F.function = &unwrap;
        F.params = &f1;
        gsl_integration_qag(&F, 0.0, sol.sE, 0, 1e-7, 1000,
                            3, w, &result1, &error1);
      }
      {
        F f1 = [&](double s){ return sol.Dphi(s, 1) * std::cos( sol.theta(s, 1) ); };
        gsl_function F;
        F.function = &unwrap;
        F.params = &f1;
        gsl_integration_qag(&F, sol.sE, sol.ell, 0, 1e-7, 1000,
                            3, w, &result2, &error2);
      }
      result = result1+result2;
      error=error1+error2;
    }
    omegas.push_back( result );
    // printf ("result          = % .18f\n", result);
    // printf ("estimated error = % .18f\n", error);
    sols.push_back( sol );
  }

  {
    const double phi0 = M_PI/48.0;
    for(const auto site : sites) assert( std::abs(Mod(Pt(site).xi[1]) - phi0)>TOLLOOSE );
    std::vector<Idx> section;
    for(Idx il=0; il<links.size(); il++){
      Link link = links[il];
      const double phi1 = Pt(sites[link[0]]).xi[1];
      const double phi2 = Pt(sites[link[1]]).xi[1];
      //if( std::abs(Mod2(phi1-phi0)) + std::abs(Mod2(phi0-phi2)) <= std::abs(Mod2(phi1-phi2)) ){
      // if( std::abs(std::acos(std::cos(phi1-phi0))) + std::abs(std::acos(std::cos(phi0-phi2))) <= std::abs(std::acos(std::cos(phi1-phi2))) ){
      double dphi1 = phi1-phi0;
      double dphi2 = phi2-phi0;
      double dphi = phi1-phi2;

      if(dphi1 >= M_PI) dphi1 -= 2.0*M_PI;
      else if(dphi1 < -M_PI) dphi1 += 2.0*M_PI;
      if(dphi2 >= M_PI) dphi2 -= 2.0*M_PI;
      else if(dphi2 < -M_PI) dphi2 += 2.0*M_PI;
      if(dphi >= M_PI) dphi -= 2.0*M_PI;
      else if(dphi < -M_PI) dphi += 2.0*M_PI;

      // std::cout << phi1 << " " << phi2 << " " << phi0 << std::endl;
      // std::cout << std::abs(dphi1) << " " << std::abs(dphi2) << " " << std::abs(dphi) << std::endl;
      // if( phi2 < phi0 && phi0 < phi2 || phi2 < phi0 && phi0 < phi1 ){
      if( std::abs(dphi1) + std::abs(dphi2) - std::abs(dphi) < TOLLOOSE ){
        section.push_back( il );
        // std::cout << "yes" << std::endl;
      }
    }
    for(Idx il : section) omegas[il] -= 2.0*M_PI;
    // for(auto elem : section) std::cout << elem << " ";
    // std::cout << std::endl;
    // std::cout << section.size() << std::endl;
    // std::cout << links.size() << std::endl;
  }

  std::map<const Link, const double> omega;
  for(Idx il=0; il<links.size(); il++){
    Link link = links[il];
    omega.insert( { link, omegas[il] } );
    omega.insert( { Link{link[1], link[0]}, -omegas[il] } );
  }


  std::vector<double> alpha0s, alpha1s;
  for( auto& sol : sols){
    double alpha0, alpha1;
    if(!sol.is_split){
      const double theta0 = sol.theta(0);
      const double dtheta_ds0 = sol.Dtheta(0);
      const double dphi_ds0 = sol.Dphi(0);
      assert( std::abs(dtheta_ds0*dtheta_ds0 + std::sin(theta0)*std::sin(theta0)*dphi_ds0*dphi_ds0-1.0)<TOLLOOSE );

      VD e0( dtheta_ds0, std::sin(theta0)*dphi_ds0 );
      alpha0 = std::acos( e0[0] );
      if(e0[1]<0) alpha0 *= -1.0;
      assert( std::abs(std::cos(alpha0)-e0[0]) + std::abs(std::sin(alpha0)-e0[1]) < TOLLOOSE );

      const double theta1 = sol.theta(sol.ell);
      const double dtheta_ds1 = sol.Dtheta(sol.ell);
      const double dphi_ds1 = sol.Dphi(sol.ell);
      assert( std::abs(dtheta_ds1*dtheta_ds1 + std::sin(theta1)*std::sin(theta1)*dphi_ds1*dphi_ds1-1.0)<TOLLOOSE );

      VD e1( dtheta_ds1, std::sin(theta1)*dphi_ds1 );
      alpha1 = std::acos( e1[0] );
      if(e1[1]<0) alpha1 *= -1.0;
      assert( std::abs(std::cos(alpha1)-e1[0]) + std::abs(std::sin(alpha1)-e1[1]) < TOLLOOSE );
    }
    else{
      const double theta0 = sol.theta(0, 0);
      const double dtheta_ds0 = sol.Dtheta(0, 0);
      const double dphi_ds0 = sol.Dphi(0, 0);
      assert( std::abs(dtheta_ds0*dtheta_ds0 + std::sin(theta0)*std::sin(theta0)*dphi_ds0*dphi_ds0-1.0)<TOLLOOSE );

      VD e0( dtheta_ds0, std::sin(theta0)*dphi_ds0 );
      alpha0 = std::acos( e0[0] );
      if(e0[1]<0) alpha0 *= -1.0;
      assert( std::abs(std::cos(alpha0)-e0[0]) + std::abs(std::sin(alpha0)-e0[1]) < TOLLOOSE );

      const double theta1 = sol.theta(sol.ell, 1);
      const double dtheta_ds1 = sol.Dtheta(sol.ell, 1);
      const double dphi_ds1 = sol.Dphi(sol.ell, 1);
      assert( std::abs(dtheta_ds1*dtheta_ds1 + std::sin(theta1)*std::sin(theta1)*dphi_ds1*dphi_ds1-1.0)<TOLLOOSE );

      VD e1( dtheta_ds1, std::sin(theta1)*dphi_ds1 );
      alpha1 = std::acos( e1[0] );
      if(e1[1]<0) alpha1 *= -1.0;
      assert( std::abs(std::cos(alpha1)-e1[0]) + std::abs(std::sin(alpha1)-e1[1]) < TOLLOOSE );
    }
    alpha0s.push_back( alpha0 );
    alpha1s.push_back( alpha1 );
  }

  std::map<const Link, const double> alpha;
  for(Idx il=0; il<links.size(); il++){
    Link link = links[il];
    // std::cout << link[0] << " " << link[1] << std::endl;
    alpha.insert( { link, alpha0s[il] } );
    alpha.insert( { Link{link[1],link[0]}, alpha1s[il] } );
  }

  double TOL=1.0e-6;
  {
    std::cout << "# checking spin structure" << std::endl;
    for(Idx ix=0; ix<sites.size(); ix++){
      for(Idx iy : nns[ix]){
        // std::cout << ix << " " << iy << std::endl;
        const double alpha1 = alpha.at(Link{ix,iy});
        double alpha2 = alpha.at(Link{iy,ix});
        double omega12 = omega.at(Link{ix,iy});

        // double diff = (alpha2 + M_PI + omega12) - alpha1;
        // assert( std::abs(Mod(diff))<TOL );
      }}
  }
  {
    std::cout << "# checking deficits" << std::endl;
    int counter=0;
    for(auto face : faces){
      int sign = 1;
      {
        VE x0 = simp_sites[counter];
        VE x1 = sites[face[0]];
        VE x2 = sites[face[1]];
        if((x1-x0).cross(x2-x1).dot(x0) < 0) sign = -1;
      }
      // std::cout << "sign = " << sign << std::endl;

      double sum = 0.0;
      for(int i=0; i<face.size(); i++){
        // std::cout << "face.size() = " << face.size() << std::endl;
        // std::cout << "i = " << i << std::endl;
        const Idx ix = face[i];
        // std::cout << "ix = " << ix << std::endl;
        const Idx j = (i+1)%face.size();
        // std::cout << "j = " << j << std::endl;
        const Idx jx = face[j];
        // std::cout << "jx = " << jx << std::endl;
        sum += omega.at( Link{ix,jx} );
      }
      sum *= sign;
      // std::cout << "sum = " << sum << std::endl;
      double mod = Mod(sum, 4.0*M_PI);
      // std::cout << "sum (mod4pi) = " << mod << std::endl;
      if(mod>2.0*M_PI) mod -= 4.0*M_PI;
      std::cout << "sum (mod4pi, repr) = " << mod << std::endl;
      counter++;
    }
  }



  return 0;
}

