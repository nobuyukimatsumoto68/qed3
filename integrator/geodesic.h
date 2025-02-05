#include <iostream>
#include <cmath>

#include <Eigen/Dense>


// constexpr double TOLLOOSE=1.0e-10;
constexpr double TOLLOOSE=1.0e-6;
constexpr double EPSNUMDER=1.0e-6;


double Mod(double a, double b=2.0*M_PI){
  int p = int(std::floor(a / b));
  // if(a<0) p += 1;
  double r = a - p*b;
  return r;
}

double Mod2(double a){
  double tmp = Mod(a);
  if(tmp>M_PI) tmp -= 2.0*M_PI;
  return tmp;
}

template<typename V>
V Mod2(const V& a){
  double b=2.0*M_PI;

  V res = a;
  for(auto& elem : res){
    // int p = int(std::floor((elem+2.0*b) / b));
    // int p = int(std::floor((elem+2.0*b) / b));
    // // if(a<0) p += 1;
    // double r = elem - (p+2)*b;
    // // if(r>M_PI) r -= 2.0*M_PI;
    elem = Mod(elem);
  }
  return res;
}


int sgn(const double a){
  int res = 1;
  if(a<0) res *= -1;
  return res;
}


bool isModdable( const double v, const double q=2.0*M_PI, const double TOL=TOLLOOSE ){
  const double tmp1 = Mod(std::abs(v),q);
  const double tmp2 = Mod(-std::abs(v),q);
  return tmp1<TOL || tmp2<TOL;
}


using I2=Eigen::Vector2i;

constexpr int DIM = 2;
using V2=Eigen::Vector2d;

constexpr int EDIM = 3;
using V3=Eigen::Vector3d;

V3 embedding3D( const V2& xi ){
  return V3(std::sin(xi[0])*std::cos(xi[1]),
	    std::sin(xi[0])*std::sin(xi[1]),
	    std::cos(xi[0]));
}

V2 projectionS2( const V3& x ){
  const double r = x.norm();
  const double theta = std::acos(x(2)/r);
  // std::cout << "debug. ismoddable = " << isModdable(theta, M_PI) << std::endl;
  double phi = 0.0;
  if(!isModdable(theta, M_PI)){
    // std::cout << "debug. arg = " << x(0)/(r*std::sin(theta)) << std::endl;
    double arg = x(0)/(r*std::sin(theta));
    if(arg<=-1.0) phi=M_PI;
    else if(arg>=1.0) phi=0.0;
    else {
      phi = std::acos( x(0)/(r*std::sin(theta)) );
      if(x(1)<0.0) phi *= -1.0;
    }
  }
  return V2(theta, phi);
}


struct Pt{
  const V3 x;
  const V2 xi;
  const bool is_singular;

  Pt() = delete;

  Pt( const V3& x )
    : x(x)
    , xi(projectionS2(x))
    , is_singular( isModdable(xi(0), M_PI) )
  {}

  Pt( const V2& xi )
    : x(embedding3D(xi))
    , xi(xi)
    , is_singular( isModdable(xi(0), M_PI) )
  {}

  Pt(const double x1,
     const double x2,
     const double x3)
    : x(V3(x1,x2,x3))
    , xi(projectionS2(x))
    , is_singular( isModdable(xi(0), M_PI) )
  {}

  Pt( const double xi1,
      const double xi2 )
    : x(embedding3D(V2(xi1,xi2)))
    , xi(V2(xi1,xi2))
    , is_singular( isModdable(xi(0), M_PI) )
  {}


  double theta() const { return xi(0); }
  double phi() const {
    assert( !is_singular );
    return xi(1);
  }

  V3 e0() const {
    return V3(std::cos(xi[0])*std::cos(xi[1]),
              std::cos(xi[0])*std::sin(xi[1]),
              -std::sin(xi[0]));
  }

  V3 e1() const {
    return V3(-std::sin(xi[0])*std::sin(xi[1]),
              std::sin(xi[0])*std::cos(xi[1]),
              0.0);
  }

  double operator[](const int mu) const { return x(mu); }
};


double geodesicLength(const Pt& x1, const Pt& x2){
  const double inner = x1.x.dot( x2.x );
  return std::acos( inner );
}


using F=std::function<double(const double)>;
using VF=std::array<F, 2*DIM>;

struct Sol{
  const bool is_split;
  const double ell;
  const double sE;
  const Pt p1;
  const Pt p2;

  std::vector<VF> solutions;

  Sol() = delete;

  Sol(const VF sol,
      const double ell,
      const Pt p1,
      const Pt p2)
    : is_split(false)
    , ell(ell)
    , sE(ell)
    , p1(p1)
    , p2(p2)
  {
    solutions.push_back(sol);
  }

  Sol(const VF sol1,
      const VF sol2,
      const double ell,
      const double sE,
      const Pt p1,
      const Pt p2)
    : is_split(true)
    , ell(ell)
    , sE(sE)
    , p1(p1)
    , p2(p2)
  {
    solutions.push_back(sol1);
    solutions.push_back(sol2);
  }

  double theta(const double s, const int level=0){
    assert(level<solutions.size());
    return solutions[level][0](s);
  }
  double phi(const double s, const int level=0){
    assert(level<solutions.size());
    return solutions[level][1](s);
  }
  double Dtheta(const double s, const int level=0){
    assert(level<solutions.size());
    return solutions[level][2](s);
  }
  double Dphi(const double s, const int level=0){
    assert(level<solutions.size());
    return solutions[level][3](s);
  }
};



void getSign(I2& sign1, I2& sign2,
             const Pt& x1, const Pt& x2, const double eps=EPSNUMDER){
  // const Pt& x1, const Pt& x2, const double eps=10.0*TOLLOOSE){
  const V3 p = x1.x;
  const V3 q = x2.x;

  V2 deriv1, deriv2;
  // const V3 p1 = embedding3D(projectionS2(p + TOLLOOSE*(q-p)));
  // const V3 p2 = embedding3D(projectionS2(p + (TOLLOOSE+eps)*(q-p)));
  const V3 p2 = embedding3D(projectionS2(p + eps*(q-p)));
  // std::cout << "debug. p1.xi = " << Pt(p1).xi.transpose() << std::endl;
  // std::cout << "debug. p2.xi = " << Pt(p2).xi.transpose() << std::endl;
  // const V3 diff1 = p2-p1;
  const V3 diff1 = p2-p;
  deriv1 << diff1.dot( x1.e0() ), diff1.dot( x1.e1() );
  // deriv1 << diff1.dot( Pt(p1).e0() ), diff1.dot( Pt(p1).e1() );
  // std::cout << deriv1.transpose() << std::endl;

  // const V3 q2 = embedding3D(projectionS2(q - eps*(q-p)));
  // const V3 q1 = embedding3D(projectionS2(q - TOLLOOSE*(q-p)));
  // const V3 q2 = embedding3D(projectionS2(q - (TOLLOOSE+eps)*(q-p)));
  // std::cout << "debug." << projectionS2( q ) << std::endl;
  // std::cout << "debug." << projectionS2(q - eps*(q-p)) << std::endl;
  const V3 q2 = embedding3D(projectionS2(q - eps*(q-p)));
  // std::cout << "debug. q1.xi = " << Pt(q1).xi.transpose() << std::endl;
  // std::cout << "debug. q2.xi = " << Pt(q2).xi.transpose() << std::endl;

  // const V3 diff2 = -q2+q;
  // const V3 diff2 = -q2+q1;
  // const V3 diff2 = q1-q2;
  const V3 diff2 = q-q2;
  deriv2 << diff2.dot( x2.e0() ), diff2.dot( x2.e1() );
  // deriv2 << diff2.dot( Pt(q2).e0() ), diff2.dot( Pt(q2).e1() );

  // int sgn2 = 0;
  // const double dphi = x2.xi[1] - x1.xi[1];
  // if( TOLLOOSE<dphi && dphi<M_PI ) sgn2 = 1;
  // else if( -M_PI<dphi && dphi<-TOLLOOSE ) sgn2 = -1;
  // else if( M_PI<dphi ) sgn2 = -1;
  // else if( dphi<-M_PI ) sgn2 = 1;
  // else if(std::abs(dphi)<M_PI) sgn2 = 0;
  // else assert(false);
  // deriv1[1] = sgn2;
  // deriv2[1] = sgn2;

  // std::cout << "x2.e0 = " << x2.e0().transpose() << std::endl;
  // std::cout << "x2.e1 = " << x2.e1().transpose() << std::endl;
  // std::cout << deriv2.transpose() << std::endl;

  // const V2 deriv1F = ( projectionS2(p + eps*(q-p)) - projectionS2(p) )/eps;
  // const V2 deriv1B = -( projectionS2(p - eps*(q-p)) - projectionS2(p) )/eps;
  // V2 deriv1;
  // if(deriv1F.norm() < deriv1B.norm()) deriv1 = deriv1F;
  // else deriv1 = deriv1B;

  // const V2 deriv2F = ( projectionS2(q + eps*(q-p)) - projectionS2(q) )/eps;
  // const V2 deriv2B = -( projectionS2(q - eps*(q-p)) - projectionS2(q) )/eps;
  // V2 deriv2;
  // if(deriv2F.norm() < deriv2B.norm()) deriv2 = deriv2F;
  // else deriv2 = deriv2B;

  sign1 = deriv1.array().sign().matrix().cast<int>();
  sign2 = deriv2.array().sign().matrix().cast<int>();
}


double getCosPhi0(const Pt& x1, const Pt& x2, const int sign1, const int sign2){
  const double t12 = std::tan( x1.theta() )/std::tan( x2.theta() );
  const double numer = sign2 * ( t12*std::cos( x1.phi() ) - sign1*std::cos( x2.phi() ) );
  const double denom = std::sqrt( t12*t12 + 1.0 - 2.0*sign1*t12*std::cos( x1.phi() - x2.phi() ) );
  return numer/denom;
}


double getPhi0(const Pt& x1, const Pt& x2, const int sign1, const int sign2){
  return std::acos( getCosPhi0(x1, x2, sign1, sign2) );
}

void getPhi0WithRatioCheck(double& phi0, bool& is_sol,
			   const Pt& x1, const Pt& x2,
			   const int sign1, const int sign2,
			   const double TOL=TOLLOOSE){
  phi0 = getPhi0(x1, x2, sign1, sign2);
  const double ratio = sign1*std::tan(x1.theta())*std::sin(x1.phi()-phi0) / ( std::tan(x2.theta())*std::sin(x2.phi()-phi0) );
  is_sol = std::abs(std::abs(ratio)-1.0) < TOL;
}


Sol SolveGeodesicsConstPhi( const Pt& x1, const Pt& x2 ){
  assert( isModdable( x1.phi()-x2.phi() ) );

  const double ell = std::abs( x2.theta()-x1.theta() );

  F theta = [=](const double s){ return x1.theta() + (x2.theta()-x1.theta()) * s/ell; };
  F phi = [=](const double s){ return x1.phi(); };
  F Dtheta = [=](const double s){ return (x2.theta()-x1.theta()) * 1.0/ell; };
  F Dphi = [=](const double s){ return 0.0; };

  return Sol( VF{theta, phi, Dtheta, Dphi}, ell, x1, x2 );
}


Sol SolveGeodesicsDeltaPhiEqPi( const Pt& x1, const Pt& x2 ){
  assert( (!isModdable( x1.phi()-x2.phi(), 2.0*M_PI )) && isModdable( x1.phi()-x2.phi(), M_PI ) );

  const double ell = geodesicLength( x1, x2 );

  const double sN = geodesicLength( x1, Pt(0., 0., 1.) );
  const double sS = geodesicLength( x2, Pt(0., 0., -1.) );

  double sE, thetaE;
  if(sN<sS){
    sE = sN;
    thetaE = 0.0;
  }
  else{
    sE = sS;
    thetaE = M_PI;
  }

  F theta1 = [=](const double s){ return x1.theta() + (thetaE-x1.theta()) * s/sE; };
  F theta2 = [=](const double s){ return thetaE + (x2.theta()-thetaE) * (s-sE)/(ell-sE); };
  F phi1 = [=](const double s){ return x1.phi(); };
  F phi2 = [=](const double s){ return x2.phi(); };

  F Dtheta1 = [=](const double s){ return (thetaE-x1.theta()) * 1.0/sE; };
  F Dtheta2 = [=](const double s){ return (x2.theta()-thetaE) * 1.0/(ell-sE); };
  F Dphi1 = [=](const double s){ return 0.0; };
  F Dphi2 = [=](const double s){ return 0.0; };

  return Sol( VF{theta1, phi1, Dtheta1, Dphi1},
	      VF{theta2, phi2, Dtheta2, Dphi2},
	      ell, sE, x1, x2 );
}


Sol SolveGeodesicsEndPtPiHalf( const Pt& x1, const Pt& x2 ){
  assert( isModdable( x1.theta()-0.5*M_PI ) || isModdable( x2.theta()-0.5*M_PI ) );

  const double ell = geodesicLength( x1, x2 );

  double phi0, phip, thetap, s0;
  if( isModdable( x1.theta()-0.5*M_PI ) ){
    phi0 = x1.phi();
    phip = x2.phi();
    thetap = x2.theta();
    s0 = 0.0;
  }
  else{
    phi0 = x2.phi();
    phip = x1.phi();
    thetap = x1.theta();
    s0 = ell;
  }

  const int sign_phi = sgn( Mod( x2.phi()-x1.phi()+M_PI ) - M_PI );
  const int sign_theta = sgn( x2.theta()-x1.theta() );

  const double tmp = std::tan(thetap)*std::tan(thetap) * std::sin(phip-phi0)*std::sin(phip-phi0);
  const double absk = std::sqrt( 1.0/( 1.0 + 1.0/tmp ) );
  
  F theta = [=](const double s){ return std::acos( -sign_theta * std::sqrt(1.0-absk*absk) * std::sin(s-s0) ); };
  F phi = [=](const double s){ return phi0 + std::atan( sign_phi * absk * std::tan(s-s0) ); };

  F Dtheta = [=](const double s){
    const double r = sign_theta * std::sqrt(1.0-absk*absk);
    return r * std::cos(s-s0) / std::sqrt( 1.0 - r*r*std::sin(s-s0)*std::sin(s-s0) );
  };
  F Dphi = [=](const double s){
    const double r = sign_phi * absk;
    return r / (std::cos(s-s0)*std::cos(s-s0)) / ( 1.0 + r*r * std::tan(s-s0)*std::tan(s-s0) );
  };

  return Sol( VF{theta, phi, Dtheta, Dphi}, ell, x1, x2 );
}


Sol SolveGeodesicsMonotonic( const Pt& x1, const Pt& x2 ){
  I2 sign1, sign2;
  getSign( sign1, sign2, x1, x2 );
  assert( sign1(0) == sign2(0) );
  assert( sign1(1) == sign2(1) );

  const int sign_theta = sign1(0);
  const int sign_phi = sign1(1);

  const double ell = geodesicLength( x1, x2 );

  std::vector<double> phi0s;
  for(const int sign1 : {-1,1}){
    for(const int sign2 : {-1,1}){
      double phi0;
      bool is_sol;
      getPhi0WithRatioCheck( phi0, is_sol, x1, x2, sign1, sign2 );
      if(is_sol) phi0s.push_back(phi0);
    }}

  F theta, phi, Dtheta, Dphi;
  bool is_ok = false;
  for(const double phi0 : phi0s){
    const double tmp = std::tan(x1.theta())*std::tan(x1.theta()) * std::sin(x1.phi()-phi0)*std::sin(x1.phi()-phi0);
    const double absk = std::sqrt( 1.0/( 1.0 + 1.0/tmp ) );

    for(int br=-4; br<=4; br++){
      const double s0 = -std::atan( sign_phi * std::tan(x1.phi()-phi0)/absk ) + br*M_PI;

      theta = [=](const double s){ return std::acos( -sign_theta * std::sqrt(1.0-absk*absk) * std::sin(s-s0) ); };
      Dtheta = [=](const double s){
	const double r = sign_theta * std::sqrt(1.0-absk*absk);
	return r * std::cos(s-s0) / std::sqrt( 1.0 - r*r*std::sin(s-s0)*std::sin(s-s0) );
      };
      F phi_tmp = [=](const double s){ return phi0 + std::atan( sign_phi * absk * std::tan(s-s0) ); };

      const double diff1 = phi_tmp(0.0) - x1.phi();
      const double diff2 = theta(ell) - x2.theta();
      // std::cout << "diff1 = " << diff1 << std::endl;
      // std::cout << "diff2 = " << diff2 << std::endl;

      if( isModdable(diff1, M_PI) && isModdable(diff2, 2.0*M_PI) ){
	phi = [=](const double s){ return phi0-diff1 + std::atan( sign_phi * absk * std::tan(s-s0) ); };
	Dphi = [=](const double s){
	  const double r = sign_phi * absk;
	  return r / (std::cos(s-s0)*std::cos(s-s0)) / ( 1.0 + r*r * std::tan(s-s0)*std::tan(s-s0) );
	};

	const double diff3 = phi(ell) - x2.phi();
	if( isModdable(diff3) ){
	  is_ok = true;
	  break;
	}
      }
    } // for br
  } // for phi0

  assert( is_ok );

  return Sol( VF{theta, phi, Dtheta, Dphi}, ell, x1, x2 );
}





Sol SolveGeodesicsAltering( const Pt& x1, const Pt& x2 ){
  I2 sign1, sign2;
  getSign( sign1, sign2, x1, x2 );
  assert( sign1(0) == -sign2(0) );
  assert( sign1(1) == sign2(1) );

  const int sign_theta1 = sign1(0);
  const int sign_theta2 = sign2(0);
  const int sign_phi = sign1(1);

  const double ell = geodesicLength( x1, x2 );

  std::vector<double> phi0s;
  for(const int sign1 : {-1,1}){
    for(const int sign2 : {-1,1}){
      double phi0;
      bool is_sol;
      getPhi0WithRatioCheck( phi0, is_sol, x1, x2, sign1, sign2 );
      if(is_sol) phi0s.push_back(phi0);
    }}

  // for(double phi0 : phi0s) std::cout << "debug. phi0 = " << phi0 << std::endl;

  F theta1, phi1, theta2, phi2, Dtheta1, Dphi1, Dtheta2, Dphi2;
  double sE;
  bool is_ok = false;
  for(const double phi0 : phi0s){
    // { double phi0 = phi0s[1];
    const double tmp = std::tan(x1.theta())*std::tan(x1.theta()) * std::sin(x1.phi()-phi0)*std::sin(x1.phi()-phi0);
    const double absk = std::sqrt( 1.0/( 1.0 + 1.0/tmp ) );

    for(int br1=-4; br1<=4; br1++){
      for(int br2=-4; br2<=4; br2++){
	// {{ int br1 = 2; int br2 = -2;
	const double s01 = -std::atan( sign_phi * std::tan(x1.phi()-phi0)/absk ) + br1*M_PI;
	const double s02 = ell-std::atan( sign_phi * std::tan(x2.phi()-phi0)/absk ) + br2*M_PI;

	sE = Mod( s01 + 0.5*M_PI + M_PI ) - M_PI;
	// std::cout << "debug. sE = " << sE << std::endl;
	assert( !isModdable(sE) );

	theta1 = [=](const double s){ return std::acos( -sign_theta1 * std::sqrt(1.0-absk*absk) * std::sin(s-s01) ); };
	theta2 = [=](const double s){ return std::acos( -sign_theta2 * std::sqrt(1.0-absk*absk) * std::sin(s-s02) ); };
	Dtheta1 = [=](const double s){
	  const double r = sign_theta1 * std::sqrt(1.0-absk*absk);
	  return r * std::cos(s-s01) / std::sqrt( 1.0 - r*r*std::sin(s-s01)*std::sin(s-s01) );
	};
	Dtheta2 = [=](const double s){
	  const double r = sign_theta2 * std::sqrt(1.0-absk*absk);
	  return r * std::cos(s-s02) / std::sqrt( 1.0 - r*r*std::sin(s-s02)*std::sin(s-s02) );
	};

	F phi_tmp1 = [=](const double s){ return phi0 + std::atan( sign_phi * absk * std::tan(s-s01) ); };
	F phi_tmp2 = [=](const double s){ return phi0 + std::atan( sign_phi * absk * std::tan(s-s02) ); };

	const double diff1 = phi_tmp1(0.0) - x1.phi();
	const double diff2 = phi_tmp2(ell) - x2.phi();

	const double diff3 = theta1(0.0) - x1.theta();
	const double diff4 = theta2(ell) - x2.theta();

	// std::cout << "debug. diff1 = " << diff1 << std::endl;
	// std::cout << "debug. diff2 = " << diff2 << std::endl;
	// std::cout << "debug. diff3 = " << diff3 << std::endl;
	// std::cout << "debug. diff4 = " << diff4 << std::endl;

	if( isModdable(diff1, M_PI) && isModdable(diff2, M_PI)
	    && isModdable(diff3, 2.0*M_PI) && isModdable(diff4, 2.0*M_PI) ){
	  phi1 = [=](const double s){ return phi0-diff1 + std::atan( sign_phi * absk * std::tan(s-s01) ); };
	  phi2 = [=](const double s){ return phi0-diff2 + std::atan( sign_phi * absk * std::tan(s-s02) ); };
	  Dphi1 = [=](const double s){
	    const double r = sign_phi * absk;
	    return r / (std::cos(s-s01)*std::cos(s-s01)) / ( 1.0 + r*r * std::tan(s-s01)*std::tan(s-s01) );
	  };
	  Dphi2 = [=](const double s){
	    const double r = sign_phi * absk;
	    return r / (std::cos(s-s02)*std::cos(s-s02)) / ( 1.0 + r*r * std::tan(s-s02)*std::tan(s-s02) );
	  };

	  const double diff5 = phi1(0) - x1.phi();
	  const double diff6 = phi2(ell) - x2.phi();
	  const double diff7 = phi1(sE) - phi2(sE);

	  const double diff8 = theta1(0) - x1.theta();
	  const double diff9 = theta2(ell) - x2.theta();
	  const double diff10 = theta1(sE) - theta2(sE);

	  // std::cout << "debug. diffs = "
	  //           << diff5 << ", "
	  //           << diff6 << ", "
	  //           << diff7 << ", "
	  //           << diff8 << ", "
	  //           << diff9 << ", "
	  //           << diff10 << std::endl;

	  if( isModdable(diff5)&&isModdable(diff6)&&isModdable(diff7)
	      &&isModdable(diff8)&&isModdable(diff9)&&isModdable(diff10)){
	    is_ok = true;
	    break;
	  }
	}
      }} // for br
  } // for phi0

  assert( is_ok );

  return Sol( VF{theta1, phi1, Dtheta1, Dphi1}, VF{theta2, phi2, Dtheta2, Dphi2}, ell, sE, x1, x2 );
}


Sol SolveGeodesics( const Pt& x1, const Pt& x2, const bool is_verbose=true ){
  if( isModdable( x2.phi()-x1.phi() ) ){
    if(is_verbose) std::cout << "delta phi = 0" << std::endl;
    return SolveGeodesicsConstPhi( x1, x2 );
  }

  if( isModdable( x2.phi()-x1.phi()+M_PI ) ){
    if(is_verbose) std::cout << "delta phi = pi" << std::endl;
    return SolveGeodesicsDeltaPhiEqPi( x1, x2 );
  }

  if( isModdable( x1.theta()-0.5*M_PI ) || isModdable( x2.theta()-0.5*M_PI ) ){
    if(is_verbose) std::cout << "theta end = pi/2" << std::endl;
    return SolveGeodesicsEndPtPiHalf( x1, x2 );
  }

  I2 sign1, sign2;
  getSign( sign1, sign2, x1, x2 );
  if(sign1(0)==sign2(0)){
    if(is_verbose) std::cout << "monotonic" << std::endl;
    return SolveGeodesicsMonotonic( x1, x2 );
  }
  else{
    if(is_verbose) std::cout << "altering" << std::endl;
    return SolveGeodesicsAltering( x1, x2 );
  }

  assert(false);
}

