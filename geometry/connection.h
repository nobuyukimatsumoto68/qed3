#pragma once

#include <gsl/gsl_integration.h>

using SpinMatrix = Eigen::Matrix2cd;


using Func=std::function<double(const double)>;
double unwrap(double x, void *p) {
  auto fp = static_cast<Func*>(p);
  return (*fp)(x);
}

constexpr Complex I = Complex(0.0, 1.0);

struct DualSpinStructure : public Sphere {
  const RefinedIcosahedronDual& dual;

  std::array<SpinMatrix, 4> sigma;
  // std::vector<double> kappa; // evan's link label

  void set_sigma(){
    sigma[0] << 1,0,0,1;
    sigma[1] << 0,1,1,0;
    sigma[2] << 0,-I,I,0;
    sigma[3] << 1,0,0,-1;
  }

  std::vector<double> alphas;
  std::vector<double> omegas;

  const double TOL;
  const double phi0;

  DualSpinStructure(const RefinedIcosahedronDual& dual_,
                    const double TOL_=4.0e-12, // okay up to L~160
                    const double phi0_=M_PI/48.0
                    )
    : dual( dual_ )
    , TOL(TOL_)
    , phi0(phi0_)
  {
    for(const V3& r : dual.vertices){
      const double phi_ = phi(r);
      assert( std::norm(phi0-phi_)>1.0e-14 );
    }

    FillAlphas();
    FillOmegas();
    check_omega();

    set_sigma();
  }

  void check_omega() const {
    for(Idx i=0; i<dual.simplicial.NVertices(); i++){
      const DualFace face = dual.faces[i];
      double omega_sum = 0.0;
      for(auto elem : face) {
        omega_sum += omegas[ dual.directedlinkidx( elem.first, elem.second ) ];
      }
      while(omega_sum>2.0*M_PI) omega_sum -= 4.0*M_PI;
      while(omega_sum<=-2.0*M_PI) omega_sum += 4.0*M_PI;

      assert( std::abs(omega_sum+dual.vols[i])<TOL );
    }
  }

  // DirectedLink idx2DirectedLinkReverse( const Idx il ) const {
  //   const int df = il%3;
  //   const Idx iff = il/3;
  //   const FaceCoords f1 = dual.idx2FaceCoords(iff);

  //   FaceCoords f2;
  //   dual.shift(f2, f1, df);

  //   return std::make_pair(dual.idx(f2), df);
  // }

  void FillAlphas(){
    alphas.clear();
    alphas.resize( dual.NDirectedLinks() );

    for(Idx if1=0; if1<dual.NVertices(); if1++){
      const FaceCoords f1 = dual.idx2FaceCoords(if1); //{ n[0],n[1],n[2],type };
      const V3 rf1 = dual.vertices[if1];

      for(int df=0; df<3; df++){
        FaceCoords f2;
        dual.shift( f2, f1, df);
        const Idx if2 = dual.idx( f2 );
        const V3 rf2 = dual.vertices[if2];

        alphas[dual.directedlinkidx(if1, df)] = alpha( rf1, rf2, 0.0 );
      }
    }
  }



  const int limit = 10000; // 1000;
  double epsabs = 1.0e-15; // 0.;
  double epsrel = 1.0e-13; // TOLLOOSE;
  int key = 6;



  void FillOmegas(){
    omegas.clear();
    omegas.resize( dual.NDirectedLinks() );

    for(Idx if1=0; if1<dual.NVertices(); if1++){
      const FaceCoords f1 = dual.idx2FaceCoords(if1); //{ n[0],n[1],n[2],type };
      const V3 rf1 = dual.vertices[if1];

      for(int df=0; df<3; df++){
        FaceCoords f2;
        dual.shift( f2, f1, df);
        const Idx if2 = dual.idx( f2 );
        const V3 rf2 = dual.vertices[if2];

        double omega12;
        {
          gsl_integration_workspace * w = gsl_integration_workspace_alloc (limit);
          double result, error;

          Func f1 = [&](const double stilde){
            const V3 p_ = p(rf1, rf2, stilde);
            const double dphi_dstilde_ = dphi_dstilde(rf1, rf2, stilde);
            const double theta_ = theta( p_ );

            return -dphi_dstilde_ * std::cos( theta_ );
          };
          gsl_function F;
          F.function = &unwrap;
          F.params = &f1;
          gsl_integration_qag(&F, 0., 1., epsabs, epsrel, limit,
                              key, w, &result, &error);

          omega12 = result;
          gsl_integration_workspace_free (w);
        }

        {
          const double alpha_1 = alpha( rf1,rf2,0.0);
          const double alpha_2 = alpha( rf1,rf2,1.0);

          double diff = alpha_2 - omega12 - alpha_1;
          // std::cout << "debug. diff = " << diff << std::endl;
          if(diff>M_PI) diff -= 2.0*M_PI;
          else if(diff<=-M_PI) diff += 2.0*M_PI;
          if(std::abs(diff)>TOL){
            std::cout << "debug. diff = " << std::abs(diff) << std::endl;
            assert( std::abs(diff)<TOL );
          }
        }

        const double phi1 = phi(rf1);
        const double phi2 = phi(rf2);

        double dphi1 = phi1-phi0;
        double dphi2 = phi2-phi0;
        double dphi = phi1-phi2;

        if(dphi1 >= M_PI) dphi1 -= 2.0*M_PI;
        else if(dphi1 < -M_PI) dphi1 += 2.0*M_PI;
        if(dphi2 >= M_PI) dphi2 -= 2.0*M_PI;
        else if(dphi2 < -M_PI) dphi2 += 2.0*M_PI;
        if(dphi >= M_PI) dphi -= 2.0*M_PI;
        else if(dphi < -M_PI) dphi += 2.0*M_PI;
        if( std::abs(dphi1) + std::abs(dphi2) - std::abs(dphi) < 1.0e-14 ) omega12 -= 2.0*M_PI;

        omegas[dual.directedlinkidx(if1, df)] = omega12;
      }
    }
  }

  SpinMatrix Omega( const FaceCoords f, const int df ) const {
    const double omega = omegas[dual.directedlinkidx(dual.idx(f), df)];
    return std::cos(0.5*omega) * sigma[0] + I*std::sin(0.5*omega)*sigma[3];
  }

  SpinMatrix P( const FaceCoords f, const int df ) const {
    const double alpha = alphas[dual.directedlinkidx(dual.idx(f), df)];
    const SpinMatrix tmp = 0.5*(sigma[0] - std::cos(alpha)*sigma[1] - std::sin(alpha)*sigma[2]);
    return tmp * Omega(f,df);
  }


};
