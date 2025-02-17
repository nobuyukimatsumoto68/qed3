#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <map>
#include <limits>
#include <stdfloat>

#include <Eigen/Dense>

// #include "geodesic.h"
// #include "integral.h"




using Double = std::float64_t;
// using Double = std::float128_t;
#include "geodesic.h"
#include "integral.h"


// #include "dirac_s2.h"
// #include "s2.h"


using namespace Geodesic;

using Idx = std::int32_t;
// using Complex = std::complex<double>;
// const Complex I = Complex(0.0, 1.0);



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



Double TOL=1.0e-14;

constexpr Double TOLLOOSE=1.0e-4;

const int limit = 4000; // 1000;
Double epsabs = 1.0e-13; // 0.;
Double epsrel = 1.0e-13; // TOLLOOSE;
int key = 5;


std::string dir = "/mnt/hdd_barracuda/qed3/dats/";


// we only need omega, alpha

int main(int argc, char* argv[]){
  int n_refine=1;
  if(argc==2) n_refine = atoi(argv[1]);

  //  QfeLatticeS2 lattice(5, n_refine);


  // Double tht = _M_PI/n_refine;
  Double tht = 1.0/n_refine;
  // Double tht2 = -_M_PI/11.0;
  Eigen::Matrix<Double, 3, 3> rot;
  rot << std::cos(tht), 0.0, std::sin(tht),
    0.0, 1.0, 0.0,
    -std::sin(tht), 0.0, std::cos(tht);
  // Eigen::Matrix<Double, 3, 3> rot2;
  // rot2 <<
  //   std::cos(tht2), std::sin(tht2), 0.0,
  //   -std::sin(tht2), std::cos(tht2), 0.0,
  //   0.0, 0.0, 1.0;

  std::vector<VE> sites;
  {
    std::ifstream file(dir+"pts_n"+std::to_string(n_refine)+"_singlepatch.dat");

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
    std::ifstream file(dir+"links_n"+std::to_string(n_refine)+"_singlepatch.dat");
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
    std::ifstream file(dir+"nns_n"+std::to_string(n_refine)+"_singlepatch.dat");
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
    std::ifstream file(dir+"face_n"+std::to_string(n_refine)+".dat");
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





  // std::vector<VE> sites;
  for(Idx ix=0; ix<sites.size(); ix++) {
    VE tmp = rot*sites[ix];
    sites[ix] = tmp;
    // // VE x = lattice.r[ix];
    // sites.push_back(x);
  }




  std::vector<Double> omegas;
  std::vector<Sol> sols;
  int counter=0;
  // std::vector<Link> links;
  for(const Link& link : links){
    // {
    // Link link = links[74];
    // std::cout << "counter = " << counter << std::endl;
    counter++;
    // Link link = links[2];
    Idx ix1 = link[0];
    Idx ix2 = link[1];

    Pt x1( sites[ix1] );
    Pt x2( sites[ix2] );

    std::cout << "counter = " << counter << std::endl;
    std::cout << "x1.x = " << x1.x << std::endl;
    // std::cout << "x1.xi = " << x1.xi.transpose() << std::endl;
    std::cout << "x2.x = " << x2.x << std::endl;
    // std::cout << "x2.xi = " << x2.xi.transpose() << std::endl;
    // {
    //   I2 sign1, sign2;
    //   getSign( sign1, sign2, x1, x2 );
    //   std::cout << "sign1 = " << sign1.transpose() << std::endl;
    //   std::cout << "sign2 = " << sign2.transpose() << std::endl;
    // }

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
    // printf ("result          = % .18f\n", result);
    // printf ("estimated error = % .18f\n", error);
    sols.push_back( sol );
  }

  {
    // const Double phi0 = 0.3;
    const Double phi0 = _M_PI/48.0;
    // const Double phi0 = 0.0;
    // for(Idx ix=0; ix<lattice.n_sites; ix++) assert( std::abs(Mod(Pt(sites[ix]).xi[1]) - phi0)>TOLLOOSE );
    std::vector<Idx> section;
    for(Idx il=0; il<links.size(); il++){
      const Link& link = links[il];
      Idx ix1 = link[0];
      Idx ix2 = link[1];

      const Double phi1 = Pt(sites[ix1]).xi[1];
      const Double phi2 = Pt(sites[ix2]).xi[1];
      Double dphi1 = phi1-phi0;
      Double dphi2 = phi2-phi0;
      Double dphi = phi1-phi2;

      if(dphi1 >= _M_PI) dphi1 -= 2.0*_M_PI;
      else if(dphi1 < -_M_PI) dphi1 += 2.0*_M_PI;
      if(dphi2 >= _M_PI) dphi2 -= 2.0*_M_PI;
      else if(dphi2 < -_M_PI) dphi2 += 2.0*_M_PI;
      if(dphi >= _M_PI) dphi -= 2.0*_M_PI;
      else if(dphi < -_M_PI) dphi += 2.0*_M_PI;

      // std::cout << phi1 << " " << phi2 << " " << phi0 << std::endl;
      // std::cout << std::abs(dphi1) << " " << std::abs(dphi2) << " " << std::abs(dphi) << std::endl;
      //if( phi2 < phi0 && phi0 < phi2 || phi2 < phi0 && phi0 < phi1 ){
      if( std::abs(dphi1) + std::abs(dphi2) - std::abs(dphi) < TOLLOOSE ){
        section.push_back( il );
        // std::cout << "yes" << std::endl;
      }
    }
    for(Idx il : section) omegas[il] -= 2.0*_M_PI;
    // for(auto elem : section) std::cout << elem << " ";
    // std::cout << std::endl;
    // std::cout << section.size() << std::endl;
    // std::cout << links.size() << std::endl;
  }

  std::map<const Link, Double> omega;
  for(Idx il=0; il<links.size(); il++){
    // Link link = links[il];
    const Link& link = links[il];
    Idx ix1 = link[0];
    Idx ix2 = link[1];

    omega.insert( { Link{ix1, ix2},  omegas[il] } );
    omega.insert( { Link{ix2, ix1}, -omegas[il] } );
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
      Double theta0, dtheta_ds0, dphi_ds0;
      std::cout << "sol.sE = " << sol.sE << std::endl;
      if( !isModdable(sol.sE) ){
        theta0 = sol.theta(0, 0);
        dtheta_ds0 = sol.Dtheta(0, 0);
        dphi_ds0 = sol.Dphi(0, 0);
        std::cout << "debug. " <<
          dtheta_ds0 << " " << theta0 << " " << dphi_ds0 << std::endl <<
          std::abs(dtheta_ds0*dtheta_ds0 + std::sin(theta0)*std::sin(theta0)*dphi_ds0*dphi_ds0-1.0)
                  << std::endl;
        assert( std::abs(dtheta_ds0*dtheta_ds0 + std::sin(theta0)*std::sin(theta0)*dphi_ds0*dphi_ds0-1.0)<TOLLOOSE );
      }
      else{
        theta0 = sol.theta(sol.sE, 1);
        dtheta_ds0 = sol.Dtheta(sol.sE, 1);
        dphi_ds0 = sol.Dphi(sol.sE, 1);
        assert( std::abs(dtheta_ds0*dtheta_ds0 + std::sin(theta0)*std::sin(theta0)*dphi_ds0*dphi_ds0-1.0)<TOLLOOSE );
      }

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
      if(std::abs(e1[1]) > TOLLOOSE && e1[1]<0) alpha1 *= -1.0;
      std::cout << "debug. " << std::abs(std::cos(alpha1)-e1[0]) + std::abs(std::sin(alpha1)-e1[1])
                << std::endl;
      assert( std::abs(std::cos(alpha1)-e1[0]) + std::abs(std::sin(alpha1)-e1[1]) < TOLLOOSE );
    }
    alpha0s.push_back( alpha0 );
    alpha1s.push_back( Mod(alpha1 + _M_PI) );
  }

  std::map<const Link, Double> alpha;
  for(Idx il=0; il<links.size(); il++){
    Link link = links[il];
    // const Link link = lattice.links[il];
    Idx ix1 = link[0];
    Idx ix2 = link[1];

    // omega.insert( { Link{ix1, ix2},  omegas[il] } );
    // omega.insert( { Link{ix1, ix2}, -omegas[il] } );
    alpha.insert( { Link{ix1, ix2}, alpha0s[il] } );
    alpha.insert( { Link{ix2, ix1}, alpha1s[il] } );
  }

  // for(Idx ix=0; ix<sites.size(); ix++){
  // for(Idx ix=0; ix<lattice.n_sites; ix++){
  //   V3 e0x = Pt( sites[ix] ).e0();
  //   V3 e1x = Pt( sites[ix] ).e1();

  //   const auto x = lattice.sites[ix];
  //   for(int iw=0; iw<x.nn; iw++){
  //     Idx iy = x.neighbors[iw];
  //     // for(Idx iy : nns[ix]){
  //     const V3 diff = sites[iy] - sites[ix];
  //     const Double c0 = e0x.dot(diff);
  //     const Double c1 = e1x.dot(diff);
  //     Double al = std::atan(c1/c0);

  //     if(isModdable(Pt(sites[ix]).xi[0])){
  //       al = Pt( sites[iy] ).xi[1];
  //     }
  //     else{
  //       std::cout << "c0 = " << c0 << std::endl;
  //       std::cout << "c1 = " << c1 << std::endl;
  //       std::cout << "alpha = " << al << std::endl;
  //       std::cout << "check = " << alpha.at(Link({ix, iy})) << std::endl;
  //       std::cout << "diff = " << al-alpha.at(Link({ix, iy})) << std::endl;
  //       std::cout << "diff = " << isModdable(al-alpha.at(Link({ix, iy})), _M_PI) << std::endl;
  //       const int br = decide_branch( al-alpha.at(Link({ix, iy})) );
  //       // // std::cout << "branch = " << br << std::endl;
  //       al -= _M_PI*br;
  //       // std::cout << "alpha = " << al << std::endl;
  //       // std::cout << "diff = " << al-alpha.at(Link({ix, iy})) << std::endl;
  //       alpha[Link({ix, iy})] = al;
  //       // std::cout << "diff = " << al-alpha.at(Link({ix, iy})) << std::endl;
  //       // std::cout << "diff = " << isModdable(al-alpha.at(Link({ix, iy}))) << std::endl;
  //     }
  //   }
  // }


  {
    std::clog << "# checking spin structure" << std::endl;
    for(Idx ix=0; ix<sites.size(); ix++){
      // const VE x = sites[ix];
      for(int iw=0; iw<nns[ix].size(); iw++){
        Idx iy = nns[ix][iw];
        // for(Idx iy : nns[ix]){
        // std::cout << ix << " " << iy << std::endl;
        const Double alpha1 = alpha.at(Link{ix,iy});
        Double alpha2 = alpha.at(Link{iy,ix});
        Double omega12 = omega.at(Link{ix,iy});

        Double diff = (alpha2 + _M_PI + omega12) - alpha1;
        // std::cout << std::abs(Mod(diff)) << std::endl;
        // std::cout << isModdable(diff, TOL) << std::endl;
        assert( isModdable(diff, TOL) );

        Double om = alpha1 - (alpha2 + _M_PI);
        const int br = decide_branch( om-omega12 );
        om -= _M_PI*br;
        omega[Link({ix, iy})] = om;
      }}
  }

  // // make omegas exact
  // // {
  // //   for(Idx il=0; il<links.size(); il++){
  // //     Link link = links[il];
  // //     Idx ix = link[0];
  // //     Idx iy = link[1];

  // //     const Double alpha1 = alpha.at(Link{ix,iy});
  // //     Double alpha2 = alpha.at(Link{iy,ix});
  // //     // Double omega12 = omega.at(Link{ix,iy});
  // //     // Double diff = (alpha2 + _M_PI + omega12) - alpha1;
  // //     Double omega12 = alpha1 - alpha2 - _M_PI;
  // //     // std::cout << std::abs(Mod(diff)) << std::endl;
  // //     // std::cout << isModdable(diff, TOL) << std::endl;
  // //     // assert( isModdable(diff, TOL) );
  // //     // }
  // //     omegas[il] = omega12;
  // //   }
  // // }

  // // omega.clear();
  // // for(Idx il=0; il<links.size(); il++){
  // //   Link link = links[il];
  // //   omega.insert( { link, omegas[il] } );
  // //   omega.insert( { Link{link[1], link[0]}, -omegas[il] } );
  // // }

  // // {
  // //   std::clog << "# checking spin structure (2)" << std::endl;
  // //   for(Idx ix=0; ix<sites.size(); ix++){
  // //     for(Idx iy : nns[ix]){
  // //       // std::cout << ix << " " << iy << std::endl;
  // //       const Double alpha1 = alpha.at(Link{ix,iy});
  // //       Double alpha2 = alpha.at(Link{iy,ix});
  // //       Double omega12 = omega.at(Link{ix,iy});

  // //       Double diff = (alpha2 + _M_PI + omega12) - alpha1;
  // //       // std::cout << std::abs(Mod(diff)) << std::endl;
  // //       // std::cout << isModdable(diff, TOL) << std::endl;
  // //       assert( isModdable(diff, TOL) );
  // //     }}
  // // }

  {
    std::clog << "# checking deficits @@@ NEED DEBUG" << std::endl;
    int counter=0;
    for(auto& face : faces){
      // for(int ia=0; ia<lattice.n_faces; ia++){
      int sign = 1;
      {
        VE x0 = sites[ face[0] ];
        VE x1 = sites[ face[1] ];
        VE x2 = sites[ face[2] ];
        // if((x1-x0).cross(x2-x1).dot(x0) < 0) sign = -1;
        VE sum = x0+x1+x2;
        if((x1-x0).cross(x2-x0).dot(sum) < 0) sign = -1;
      }
      std::cout << "sign = " << sign << std::endl;

      Double sum = 0.0;
      //for(int i=0; i<lattice.n_faces; i++){
      for(int i=0; i<3; i++){
        // std::cout << "face.size() = " << face.size() << std::endl;
        // std::cout << "i = " << i << std::endl;
        // const Idx ix = face[i];
        // // std::cout << "ix = " << ix << std::endl;
        // const Idx j = (i+1)%face.size();
        // // std::cout << "j = " << j << std::endl;
        // const Idx jx = face[j];
        Idx ix = face[i];
        Idx jx = face[(i+1)%3];

        // std::cout << "jx = " << jx << std::endl;
        sum += omega.at( Link{ix,jx} );
      }
      sum *= sign;
      std::cout << "sum = " << sum << std::endl;
      Double mod = Mod(sum, 4.0*_M_PI);
      std::cout << "sum (mod4pi) = " << mod << std::endl;
      if(mod>2.0*_M_PI) mod -= 4.0*_M_PI;
      std::clog << "# sum (mod4pi, repr) = " << mod << std::endl;
      assert( (-1.5 * 4.0*_M_PI/faces.size() < mod && mod < 0.0) );
      // if(! (-1.5 * 4.0*_M_PI/lattice.n_faces < mod && mod < 0.0) ){
      //   std::cout << "@@@@@@@@@@ ia = " << ia << std::endl;
      // };
      counter++;
    }
  }


  {
    std::ofstream ofs(dir+"omega_n"+std::to_string(n_refine)+"_singlepatch.dat");
    ofs << std::scientific << std::setprecision(15);
    for(Idx il=0; il<links.size(); il++){
      Link link = links[il];
      // const auto& link = lattice.links[il];
      Idx ix1 = link[0];
      Idx ix2 = link[1];

      ofs << std::setw(25) << ix1 << " ";
      ofs << std::setw(25) << ix2 << " ";
      ofs << std::setw(25) << omegas[il] << std::endl;
    }
  }

  {
    std::ofstream ofs(dir+"alpha_n"+std::to_string(n_refine)+"_singlepatch.dat");
    ofs << std::scientific << std::setprecision(15);
    for(Idx il=0; il<links.size(); il++){
      Link link = links[il];
      // const auto& link = lattice.links[il];
      Idx ix1 = link[0];
      Idx ix2 = link[1];

      ofs << std::setw(25) << ix1 << " ";
      ofs << std::setw(25) << ix2 << " ";
      ofs << std::setw(25) << alpha0s[il] << std::endl;

      ofs << std::setw(25) << ix2 << " ";
      ofs << std::setw(25) << ix1 << " ";
      ofs << std::setw(25) << alpha1s[il] << std::endl;
    }
  }

  // // dualtriangleareas
  // std::vector<Face> simp_faces;
  // {
  //   std::ifstream file(dir+"face_n"+std::to_string(n_refine)+"_singlepatch.dat");

  //   std::string str;
  //   // std::cout << "debug file = " << std::endl;
  //   while (std::getline(file, str)){
  //     std::istringstream iss(str);
  //     // std::cout << "debug = " << std::endl;
  //     Idx v;
  //     Face face;
  //     while( iss >> v ) {
  //       // std::cout << "face = " << v << std::endl;
  //       face.push_back( v );
  //     }
  //     simp_faces.push_back( face );
  //   }
  // }
  // std::vector<Double> triangleareas;
  // for(const Face& face : simp_faces){
  //   Double area = sphericalarea( Pt(simp_sites[face[0]]),
  //                                Pt(simp_sites[face[1]]),
  //                                Pt(simp_sites[face[2]]));
  //   triangleareas.push_back( area );
  // }
  // assert( simp_faces.size()==triangleareas.size() );
  // std::cout << sites.size() << std::endl;
  // std::cout << simp_faces.size() << std::endl;
  // std::cout << triangleareas.size() << std::endl;
  // assert( sites.size()==triangleareas.size() );

  // {
  //   std::ofstream ofs(dir+"dualtriangleareas_n"+std::to_string(n_refine)+"_singlepatch.dat");
  //   ofs << std::scientific << std::setprecision(15);
  //   for(Double elem : triangleareas) ofs << std::setw(25) << elem << std::endl;
  // }

  // Double vbar = 0.0;
  // for(const Double& vx : triangleareas) vbar += vx;
  // vbar /= triangleareas.size();
  // const Double alat = std::sqrt( 4.0/std::sqrt(3.0) * vbar );

  // // u, v
  // {
  //   std::ofstream ofs(dir+"vs_n"+std::to_string(n_refine)+"_singlepatch.dat");
  //   ofs << std::scientific << std::setprecision(15);

  //   for(Idx ix=0; ix<nns.size(); ix++){
  //     Double len1 = geodesicLength( Pt(sites[ix]), Pt(sites[nns[ix][0]]) );
  //     Double len2 = geodesicLength( Pt(sites[ix]), Pt(sites[nns[ix][1]]) );
  //     Double len3 = geodesicLength( Pt(sites[ix]), Pt(sites[nns[ix][2]]) );

  //     Double al1 = alpha.at( Link{ix, nns[ix][0]} );
  //     Double al2 = alpha.at( Link{ix, nns[ix][1]} );
  //     Double al3 = alpha.at( Link{ix, nns[ix][2]} );

  //     VD e1, e2, e3;
  //     e1 << std::cos(al1), std::sin(al1);
  //     e2 << std::cos(al2), std::sin(al2);
  //     e3 << std::cos(al3), std::sin(al3);

  //     const VD ell1 = len1 * e1;
  //     const VD ell2 = len2 * e2;
  //     const VD ell3 = len3 * e3;

  //     const VD d13 = ell1 - ell3;
  //     const VD d23 = ell2 - ell3;

  //     Eigen::Matrix4d mat;
  //     mat <<
  //       d13[0], d23[0], 0., 0.,
  //       d13[1], d23[1], 0., 0.,
  //       0., 0., d13[0], d23[0],
  //       0., 0., d13[1], d23[1];

  //     Eigen::Vector4d b;
  //     b << 1.0, 0.0, 0.0, 1.0;
  //     b *= alat/std::sqrt(3.0);

  //     Eigen::Vector4d vv = mat.inverse() * b;

  //     VD v1, v2;
  //     v1 << vv[0], vv[2];
  //     v2 << vv[1], vv[3];
  //     VD v3 = -v1-v2;

  //     ofs << std::setw(25) << v1[0] << " ";
  //     ofs << std::setw(25) << v1[1] << " ";
  //     ofs << std::setw(25) << v2[0] << " ";
  //     ofs << std::setw(25) << v2[1] << " ";
  //     ofs << std::setw(25) << v3[0] << " ";
  //     ofs << std::setw(25) << v3[1] << std::endl;
  //   }
  // }



  // {
  //   std::ofstream ofs(dir+"us_n"+std::to_string(n_refine)+"_singlepatch.dat");
  //   ofs << std::scientific << std::setprecision(15);

  //   for(Idx ix=0; ix<nns.size(); ix++){
  //     Double len1 = geodesicLength( Pt(sites[ix]), Pt(sites[nns[ix][0]]) );
  //     Double len2 = geodesicLength( Pt(sites[ix]), Pt(sites[nns[ix][1]]) );
  //     Double len3 = geodesicLength( Pt(sites[ix]), Pt(sites[nns[ix][2]]) );

  //     Double al1 = alpha.at( Link{ix, nns[ix][0]} );
  //     Double al2 = alpha.at( Link{ix, nns[ix][1]} );
  //     Double al3 = alpha.at( Link{ix, nns[ix][2]} );

  //     VD e1, e2, e3;
  //     e1 << std::cos(al1), std::sin(al1);
  //     e2 << std::cos(al2), std::sin(al2);
  //     e3 << std::cos(al3), std::sin(al3);

  //     const VD ell1 = len1 * e1;
  //     const VD ell2 = len2 * e2;
  //     const VD ell3 = len3 * e3;

  //     Eigen::Matrix2d mat;
  //     mat <<
  //       ell2[0], ell3[0],
  //       ell2[1], ell3[1];

  //     Eigen::Vector2d b;
  //     b << -ell1[0], -ell1[1];

  //     Eigen::Vector2d vv = mat.inverse() * b;

  //     VE u;
  //     u << 1.0, vv[0], vv[1];
  //     u *= 3.0/u.mean();

  //     ofs << std::setw(25) << u[0] << " ";
  //     ofs << std::setw(25) << u[1] << " ";
  //     ofs << std::setw(25) << u[2] << std::endl;
  //   }
  // }




  return 0;
}

