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

#include "s2.h"


using namespace Geodesic;

using Idx = std::int32_t;

using Link = std::array<Idx,2>; // <Idx,Idx>;
using Face = std::vector<Idx>;
using VD=Eigen::Matrix<Double, 2, 1>;
using VE=Eigen::Matrix<Double, 3, 1>;



Double TOL=1.0e-14;

constexpr Double TOLLOOSE=1.0e-6;

const int limit = 4000; // 1000;
Double epsabs = 1.0e-13; // 0.;
Double epsrel = 1.0e-13; // TOLLOOSE;
int key = 5;

std::string dir = "./data/";

int main(int argc, char* argv[]){
  int n_refine=1;
  if(argc==2) n_refine = atoi(argv[1]);

  QfeLatticeS2 lattice(5, n_refine);


  Double tht = 1.0/n_refine;
  Eigen::Matrix<Double, 3, 3> rot;
  rot << std::cos(tht), 0.0, std::sin(tht),
    0.0, 1.0, 0.0,
    -std::sin(tht), 0.0, std::cos(tht);

  std::vector<VE> sites;
  for(Idx ix=0; ix<lattice.n_sites; ix++) {
    VE x = rot*lattice.r[ix];
    sites.push_back(x);
  }




  std::vector<Double> omegas;
  std::vector<Sol> sols;
  int counter=0;
  for(const auto& link : lattice.links){
    counter++;
    Idx ix1 = link.sites[0];
    Idx ix2 = link.sites[1];

    Pt x1( sites[ix1] );
    Pt x2( sites[ix2] );

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
    for(Idx ix=0; ix<lattice.n_sites; ix++) assert( std::abs(Mod(Pt(sites[ix]).xi[1]) - phi0)>TOLLOOSE );
    std::vector<Idx> section;
    for(Idx il=0; il<lattice.n_links; il++){
      const auto& link = lattice.links[il];
      Idx ix1 = link.sites[0];
      Idx ix2 = link.sites[1];

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

      if( std::abs(dphi1) + std::abs(dphi2) - std::abs(dphi) < TOLLOOSE ){
        section.push_back( il );
      }
    }
    for(Idx il : section) omegas[il] -= 2.0*_M_PI;
  }

  std::map<const Link, Double> omega;
  for(Idx il=0; il<lattice.n_links; il++){
    const auto& link = lattice.links[il];
    Idx ix1 = link.sites[0];
    Idx ix2 = link.sites[1];

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
      if( !isModdable(sol.sE) ){
        theta0 = sol.theta(0, 0);
        dtheta_ds0 = sol.Dtheta(0, 0);
        dphi_ds0 = sol.Dphi(0, 0);
        // std::cout << "debug. " <<
        //   dtheta_ds0 << " " << theta0 << " " << dphi_ds0 << std::endl <<
        //   std::abs(dtheta_ds0*dtheta_ds0 + std::sin(theta0)*std::sin(theta0)*dphi_ds0*dphi_ds0-1.0)
        //           << std::endl;
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
      // std::cout << "debug. " << std::abs(std::cos(alpha1)-e1[0]) + std::abs(std::sin(alpha1)-e1[1])
      //           << std::endl;
      assert( std::abs(std::cos(alpha1)-e1[0]) + std::abs(std::sin(alpha1)-e1[1]) < TOLLOOSE );
    }
    alpha0s.push_back( alpha0 );
    alpha1s.push_back( Mod(alpha1 + _M_PI) );
  }

  std::map<const Link, Double> alpha;
  for(Idx il=0; il<lattice.n_links; il++){
    const auto& link = lattice.links[il];
    Idx ix1 = link.sites[0];
    Idx ix2 = link.sites[1];

    alpha.insert( { Link{ix1, ix2}, alpha0s[il] } );
    alpha.insert( { Link{ix2, ix1}, alpha1s[il] } );
  }

  {
    std::clog << "# checking spin structure" << std::endl;
    for(Idx ix=0; ix<sites.size(); ix++){
      const auto x = lattice.sites[ix];
      for(int iw=0; iw<x.nn; iw++){
        Idx iy = x.neighbors[iw];
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
    for(int ia=0; ia<lattice.n_faces; ia++){
      int sign = 1;
      {
        VE x0 = sites[ lattice.faces[ia].sites[0] ];
        VE x1 = sites[ lattice.faces[ia].sites[1] ];
        VE x2 = sites[ lattice.faces[ia].sites[2] ];
        VE sum = x0+x1+x2;
        if((x1-x0).cross(x2-x0).dot(sum) < 0) sign = -1;
      }

      Double sum = 0.0;
      for(int i=0; i<3; i++){
        Idx ix = lattice.faces[ia].sites[i];
        Idx jx = lattice.faces[ia].sites[(i+1)%3];

        sum += omega.at( Link{ix,jx} );
      }
      sum *= sign;
      // std::cout << "sum = " << sum << std::endl;
      Double mod = Mod(sum, 4.0*_M_PI);
      // std::cout << "sum (mod4pi) = " << mod << std::endl;
      if(mod>2.0*_M_PI) mod -= 4.0*_M_PI;
      // std::clog << "# sum (mod4pi, repr) = " << mod << std::endl;
      assert( (-1.5 * 4.0*_M_PI/lattice.n_faces < mod && mod < 0.0) );
      counter++;
    }
  }



  {
    std::ofstream ofs(dir+"omega_n"+std::to_string(n_refine)+".dat");
    ofs << std::scientific << std::setprecision(15);
    for(Idx il=0; il<lattice.n_links; il++){
      const auto& link = lattice.links[il];
      Idx ix1 = link.sites[0];
      Idx ix2 = link.sites[1];

      ofs << std::setw(25) << ix1 << " ";
      ofs << std::setw(25) << ix2 << " ";
      ofs << std::setw(25) << omegas[il] << std::endl;
    }
  }

  {
    std::ofstream ofs(dir+"alpha_n"+std::to_string(n_refine)+".dat");
    ofs << std::scientific << std::setprecision(15);
    for(Idx il=0; il<lattice.n_links; il++){
      const auto& link = lattice.links[il];
      Idx ix1 = link.sites[0];
      Idx ix2 = link.sites[1];

      ofs << std::setw(25) << ix1 << " ";
      ofs << std::setw(25) << ix2 << " ";
      ofs << std::setw(25) << alpha0s[il] << std::endl;

      ofs << std::setw(25) << ix2 << " ";
      ofs << std::setw(25) << ix1 << " ";
      ofs << std::setw(25) << alpha1s[il] << std::endl;
    }
  }

  return 0;
}

