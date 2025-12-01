#include <iostream>
// #include <cmath>
#include <math.h>
#include <vector>
#include <stdfloat>
#include <numbers>

#include <Eigen/Dense>

namespace Geodesic{


  constexpr Double TOLMOD=1.0e-4;
  constexpr Double TOL3 = 1.0e-6;

  // constexpr Double EPSNUMDER=1.0e-5;
  constexpr Double EPSNUMDER=1.0e-6;

  constexpr int BRMAX = 8;

  constexpr Double _M_PI = std::numbers::pi_v<Double>;


using I2=Eigen::Vector2i;

constexpr int DIM = 2;
// using V2=Eigen::Vector2d;
using V2=Eigen::Matrix<Double, 1, 2>;

constexpr int EDIM = 3;
using V3=Eigen::Matrix<Double, 1, 3>;
// using V3=Eigen::Vector3d;



Double Mod(Double a, Double b=2.0*_M_PI){
  int p = int(std::floor(a / b));
  // if(a<0) p += 1;
  Double r = a - p*b;
  return r;
}

Double Mod2(Double a){
  Double tmp = Mod(a);
  if(tmp>_M_PI) tmp -= 2.0*_M_PI;
  return tmp;
}

template<typename V>
V Mod2(const V& a){
  Double b=2.0*_M_PI;

  V res = a;
  for(auto& elem : res){
    // int p = int(std::floor((elem+2.0*b) / b));
    // int p = int(std::floor((elem+2.0*b) / b));
    // // if(a<0) p += 1;
    // Double r = elem - (p+2)*b;
    // // if(r>_M_PI) r -= 2.0*_M_PI;
    elem = Mod(elem);
  }
  return res;
}


int sgn(const Double a){
  int res = 1;
  if(a<0) res *= -1;
  return res;
}


bool isModdable( const Double v, const Double q=2.0*_M_PI, const Double TOL=TOLMOD ){
  const Double tmp1 = Mod(std::abs(v),q);
  const Double tmp2 = Mod(-std::abs(v),q);
  return tmp1<TOL || tmp2<TOL;
}


  int decide_branch( Double a, Double b=_M_PI ){
    int res;
    assert( isModdable(a,b) );
    for(res=-BRMAX; res<=BRMAX+1; res++){
      if( std::abs(a-res*b)<TOLMOD ) break;
    }
    if(res==BRMAX+1) assert(false);
    return res;
  }




V3 embedding3D( const V2& xi ){
  return V3(std::sin(xi[0])*std::cos(xi[1]),
	    std::sin(xi[0])*std::sin(xi[1]),
	    std::cos(xi[0]));
}


Double _acos(const Double arg){
  double res;
  if(arg<-1.0) res=_M_PI;
  else if(arg>1.0) res=0.0;
  else {
    res = std::acos( arg );
  }
  return res;
}


// Double _sin(const Double arg){
//   double argd = arg;
//   return std::sin(arg);
// }



V2 projectionS2( const V3& x ){
  const Double r = x.norm();

  // const Double theta = std::acos(x(2)/r);
  // std::cout << "r = " << r << std::endl;
  const Double theta = _acos(x(2)/r);
  // std::cout << "theta = " << theta << std::endl;
  // std::cout << "debug. ismoddable = " << isModdable(theta, _M_PI) << std::endl;
  Double phi = 0.0;
  if(!isModdable(theta, _M_PI)){
    // std::cout << "debug. arg = " << x(0)/(r*std::sin(theta)) << std::endl;
    const Double arg = x(0)/(r*std::sin(theta));
    // Double arg = x(0);
    // std::cout << "sine = " << std::sin(theta) << std::endl;
    // arg /= r;
    // std::cout << "arg3 = " << arg << std::endl;
    // std::cout << "theta = " << theta << std::endl;
    // double thetad = theta;
    // // Double sine = std::sin( thetad );
    // Double sine = sinl( theta );
    // std::cout << "sine = " << sine << std::endl;
    // arg /= sine;
    // std::cout << "arg3 = " << arg << std::endl;
    // if(arg<=-1.0) phi=_M_PI;
    // else if(arg>=1.0) phi=0.0;
    // else {
    //   phi = std::acos( x(0)/(r*std::sin(theta)) );
    //   if(x(1)<0.0) phi *= -1.0;
    // }
    // std::cout << "x0 = " << x(0) << std::endl;
    // std::cout << "sin = " << r*std::sin(theta) << std::endl;
    // std::cout << "sin = " << 1.0/r/std::sin(theta) << std::endl;
    // std::cout << "arg1 = " << x(0)/(r*std::sin(theta)) << std::endl;
    // std::cout << "arg2 = " << 1.0/(r*_sin(theta)) << std::endl;
    // std::cout << "arg3 = " << arg << std::endl;
    phi = _acos(arg);
    if(x(1)<0.0) phi *= -1.0;
  }
  // std::cout << "phi = " << phi << std::endl;
  return V2(theta, phi);
}


struct Pt{
  const V3 x;
  const V2 xi;
  const bool is_singular;

  Pt() = delete;

  explicit Pt( const V3& x )
    : x(x)
    , xi(projectionS2(x))
    , is_singular( isModdable(xi(0), _M_PI) )
  {}

  // explicit Pt( const V2& xi )
  //   : x(embedding3D(xi))
  //   , xi(xi)
  //   , is_singular( isModdable(xi(0), _M_PI) )
  // {}

  explicit Pt(const Double x1,
              const Double x2,
              const Double x3)
    : x(V3(x1,x2,x3))
    , xi(projectionS2(x))
    , is_singular( isModdable(xi(0), _M_PI) )
  {}

  explicit Pt( const Double xi1,
               const Double xi2 )
    : x(embedding3D(V2(xi1,xi2)))
    , xi(V2(xi1,xi2))
    , is_singular( isModdable(xi(0), _M_PI) )
  {}


  Double theta() const { return xi(0); }
  Double phi() const {
    // assert( !is_singular );
    return xi(1);
  }

  // V3 e0() const {
  //   return V3(std::cos(xi[0])*std::cos(xi[1]),
  //             std::cos(xi[0])*std::sin(xi[1]),
  //             -std::sin(xi[0]));
  // }

  // V3 e1() const {
  //   return V3(-std::sin(xi[0])*std::sin(xi[1]),
  //             std::sin(xi[0])*std::cos(xi[1]),
  //             0.0);
  // }
  V3 Delta0() const {
    // return V3(std::vector<Double>({std::cos(xi[0])*std::cos(xi[1]),
    //       std::cos(xi[0])*std::sin(xi[1]),
    //       -std::sin(xi[0])})
    //   );
    return V3(std::cos(xi[0])*std::cos(xi[1]),
              std::cos(xi[0])*std::sin(xi[1]),
              -std::sin(xi[0]));
  }

  V3 Delta1() const {
    Double zero = 0.0;
    return V3(-std::sin(xi[0])*std::sin(xi[1]),
          std::sin(xi[0])*std::cos(xi[1]),
          zero);
    // return V3(std::vector<Double>({-std::sin(xi[0])*std::sin(xi[1]),
    //       std::sin(xi[0])*std::cos(xi[1]),
    //       zero}));
  }

  V3 e0() const {
    return Delta0();
  }

  V3 e1() const {
    return Delta1()/std::sin(xi[0]);
  }


  Double operator[](const int mu) const { return x(mu); }
};


Double geodesicLength(const Pt& x1, const Pt& x2){
  const Double inner = x1.x.dot( x2.x );
  return _acos( inner );
}

Double sphericalarea( const Pt& x1, const Pt& x2, const Pt& x3 ){
  Double a = geodesicLength( x2, x3 );
  Double b = geodesicLength( x3, x1 );
  Double c = geodesicLength( x1, x2 );
  Double s = 0.5*(a+b+c);
  Double tantan = std::tan(0.5*s) * std::tan(0.5*(s-a)) * std::tan(0.5*(s-b)) * std::tan(0.5*(s-c));
  return 4.0*std::atan(std::sqrt(tantan));
}


using F=std::function<Double(const Double)>;
using VF=std::array<F, 2*DIM>;

struct Sol{
  const bool is_split;
  const Double ell;
  const Double sE;
  const Pt p1;
  const Pt p2;

  std::vector<VF> solutions;

  Sol() = delete;

  Sol(const VF sol,
      const Double ell,
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
      const Double ell,
      const Double sE,
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

  Double theta(const Double s, const int level=0){
    assert(level<solutions.size());
    return solutions[level][0](s);
  }
  Double phi(const Double s, const int level=0){
    assert(level<solutions.size());
    return solutions[level][1](s);
  }
  Double Dtheta(const Double s, const int level=0){
    assert(level<solutions.size());
    return solutions[level][2](s);
  }
  Double Dphi(const Double s, const int level=0){
    assert(level<solutions.size());
    return solutions[level][3](s);
  }
};



void getSign(I2& sign1, I2& sign2,
             const Pt& x1, const Pt& x2, const Double eps=EPSNUMDER){
  const V3 p = x1.x;
  const V3 q = x2.x;

  V2 deriv1, deriv2;
  // const V3 p2 = embedding3D(projectionS2(p + eps*(q-p)));
  // const V3 diff1 = p2-p;
  // std::cout << diff1.transpose() << std::endl;
  // deriv1 << diff1.dot( x1.e0() ), diff1.dot( x1.e1() );
  // std::cout << deriv1.transpose() << std::endl;
  // const V3 q2 = embedding3D(projectionS2(q - eps*(q-p)));
  // const V3 diff2 = q-q2;
  // deriv2 << diff2.dot( x2.e0() ), diff2.dot( x2.e1() );
  // std::cout << deriv2.transpose() << std::endl;
  const V3 p2 = embedding3D(projectionS2(p + eps*(q-p)/(q-p).norm()));
  const V3 diff1 = p2-p;
  deriv1 << diff1.dot( x1.Delta0() ), diff1.dot( x1.Delta1() );
  // const V2 deriv1( {diff1.dot( x1.Delta0() ), diff1.dot( x1.Delta1() )});
  const V3 q2 = embedding3D(projectionS2(q - eps*(q-p)));
  const V3 diff2 = q-q2;
  // const V2 deriv2( {diff2.dot( x2.Delta0() ), diff2.dot( x2.Delta1() )});
  deriv2 << diff2.dot( x2.Delta0() ), diff2.dot( x2.Delta1() );
  std::cout << deriv1 << std::endl;
  std::cout << deriv2 << std::endl;

  sign1 = deriv1.array().sign().matrix().cast<int>();
  sign2 = deriv2.array().sign().matrix().cast<int>();
}


Double getCosPhi0(const Pt& x1, const Pt& x2, const int sign1, const int sign2){
  const Double t12 = std::tan( x1.theta() )/std::tan( x2.theta() );
  const Double numer = sign2 * ( t12*std::cos( x1.phi() ) - sign1*std::cos( x2.phi() ) );
  const Double denom = std::sqrt( t12*t12 + 1.0 - 2.0*sign1*t12*std::cos( x1.phi() - x2.phi() ) );
  return numer/denom;
}


Double getPhi0(const Pt& x1, const Pt& x2, const int sign1, const int sign2){
  // return std::acos( getCosPhi0(x1, x2, sign1, sign2) );
  return _acos( getCosPhi0(x1, x2, sign1, sign2) );
}

  void getPhi0WithRatioCheck(Double& phi0, bool& is_sol,
                             const Pt& x1, const Pt& x2,
                             const int sign1, const int sign2,
                             const Double TOL=TOL3){
    phi0 = getPhi0(x1, x2, sign1, sign2);
    const Double diff = sign1*std::tan(x1.theta())*std::sin(x1.phi()-phi0) - ( std::tan(x2.theta())*std::sin(x2.phi()-phi0) );
  std::cout << std::setprecision(16) << "phi. diff = " << diff << std::endl;
  is_sol = std::abs(diff) < TOL;
  //   const Double ratio = sign1*std::tan(x1.theta())*std::sin(x1.phi()-phi0) / ( std::tan(x2.theta())*std::sin(x2.phi()-phi0) );
  // std::cout << std::setprecision(16) << "ratio = " << ratio << std::endl;
  // is_sol = std::abs(std::abs(ratio)-1.0) < TOL;
}


Sol SolveGeodesicsConstPhi( const Pt& x1, const Pt& x2 ){
  assert( isModdable( x1.phi()-x2.phi() ) || isModdable( x1.theta(), _M_PI ) || isModdable( x2.theta(), _M_PI ) );

  const Double ell = geodesicLength( x1, x2 ); // std::abs( x2.theta()-x1.theta() );

  Double cphi = x1.phi();
  if( isModdable( x1.theta(), _M_PI ) ) {
    cphi = x2.phi();
    assert( !isModdable( x2.theta(), _M_PI ) );
  }

  F theta = [=](const Double s){ return x1.theta() + (x2.theta()-x1.theta()) * s/ell; };
  F phi = [=](const Double s){ return cphi; };
  F Dtheta = [=](const Double s){ return (x2.theta()-x1.theta()) * 1.0/ell; };
  F Dphi = [=](const Double s){ return 0.0; };

  return Sol( VF{theta, phi, Dtheta, Dphi}, ell, x1, x2 );
}


Sol SolveGeodesicsDeltaPhiEqPi( const Pt& x1, const Pt& x2 ){
  assert( (!isModdable( x1.phi()-x2.phi(), 2.0*_M_PI )) && isModdable( x1.phi()-x2.phi(), _M_PI ) );

  const Double ell = geodesicLength( x1, x2 );

  const Double sN = geodesicLength( x1, Pt(0., 0., 1.) );
  const Double sS = geodesicLength( x2, Pt(0., 0., -1.) );

  Double sE, thetaE;
  if(sN<sS){
    sE = sN;
    thetaE = 0.0;
  }
  else{
    sE = sS;
    thetaE = _M_PI;
  }

  F theta1 = [=](const Double s){ return x1.theta() + (thetaE-x1.theta()) * s/sE; };
  F theta2 = [=](const Double s){ return thetaE + (x2.theta()-thetaE) * (s-sE)/(ell-sE); };
  F phi1 = [=](const Double s){ return x1.phi(); };
  F phi2 = [=](const Double s){ return x2.phi(); };

  F Dtheta1 = [=](const Double s){ return (thetaE-x1.theta()) * 1.0/sE; };
  F Dtheta2 = [=](const Double s){ return (x2.theta()-thetaE) * 1.0/(ell-sE); };
  F Dphi1 = [=](const Double s){ return 0.0; };
  F Dphi2 = [=](const Double s){ return 0.0; };

  return Sol( VF{theta1, phi1, Dtheta1, Dphi1},
	      VF{theta2, phi2, Dtheta2, Dphi2},
	      ell, sE, x1, x2 );
}


Sol SolveGeodesicsEndPtPiHalf( const Pt& x1, const Pt& x2 ){
  assert( isModdable( x1.theta()-0.5*_M_PI ) || isModdable( x2.theta()-0.5*_M_PI ) );

  const Double ell = geodesicLength( x1, x2 );

  Double phi0, phip, thetap, s0;
  if( isModdable( x1.theta()-0.5*_M_PI ) ){
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

  const int sign_phi = sgn( Mod( x2.phi()-x1.phi()+_M_PI ) - _M_PI );
  const int sign_theta = sgn( x2.theta()-x1.theta() );

  const Double tmp = std::tan(thetap)*std::tan(thetap) * std::sin(phip-phi0)*std::sin(phip-phi0);
  const Double absk = std::sqrt( 1.0/( 1.0 + 1.0/tmp ) );
  
  F theta = [=](const Double s){ return std::acos( -sign_theta * std::sqrt(1.0-absk*absk) * std::sin(s-s0) ); };
  F phi = [=](const Double s){ return phi0 + std::atan( sign_phi * absk * std::tan(s-s0) ); };

  F Dtheta = [=](const Double s){
    const Double r = sign_theta * std::sqrt(1.0-absk*absk);
    return r * std::cos(s-s0) / std::sqrt( 1.0 - r*r*std::sin(s-s0)*std::sin(s-s0) );
  };
  F Dphi = [=](const Double s){
    const Double r = sign_phi * absk;
    return r / (std::cos(s-s0)*std::cos(s-s0)) / ( 1.0 + r*r * std::tan(s-s0)*std::tan(s-s0) );
  };

  return Sol( VF{theta, phi, Dtheta, Dphi}, ell, x1, x2 );
}


Sol SolveGeodesicsMonotonic( const Pt& x1, const Pt& x2 ){
  I2 sign1, sign2;
  getSign( sign1, sign2, x1, x2 );
  assert( sign1(0) == sign2(0) );
  // assert( sign1(1) == sign2(1) );

  const int sign_theta = sign1(0);
  const int sign_phi = sign1(1);

  const Double ell = geodesicLength( x1, x2 );

  std::vector<Double> phi0s;
  for(const int sign1 : {-1,1}){
    for(const int sign2 : {-1,1}){
      Double phi0;
      bool is_sol;
      getPhi0WithRatioCheck( phi0, is_sol, x1, x2, sign1, sign2 );
      if(is_sol) {
        phi0s.push_back(phi0);
        std::cout << "yes" << std::endl;
      }
    }}

  F theta, phi, Dtheta, Dphi;
  bool is_ok = false;
  for(const Double phi0 : phi0s){
    const Double tmp = std::tan(x1.theta())*std::tan(x1.theta()) * std::sin(x1.phi()-phi0)*std::sin(x1.phi()-phi0);
    const Double absk = std::sqrt( 1.0/( 1.0 + 1.0/tmp ) );

    for(int br=-8; br<=8; br++){
      const Double s0 = -std::atan( sign_phi * std::tan(x1.phi()-phi0)/absk ) + br*_M_PI;

      theta = [=](const Double s){ return std::acos( -sign_theta * std::sqrt(1.0-absk*absk) * std::sin(s-s0) ); };
      Dtheta = [=](const Double s){
	const Double r = sign_theta * std::sqrt(1.0-absk*absk);
	return r * std::cos(s-s0) / std::sqrt( 1.0 - r*r*std::sin(s-s0)*std::sin(s-s0) );
      };
      F phi_tmp = [=](const Double s){ return phi0 + std::atan( sign_phi * absk * std::tan(s-s0) ); };

      const Double diff1 = phi_tmp(0.0) - x1.phi();
      const Double diff2 = theta(ell) - x2.theta();
      std::cout << "diff1 = " << diff1 << std::endl;
      std::cout << "diff2 = " << diff2 << std::endl;

      if( isModdable(diff1, _M_PI) && isModdable(diff2, 2.0*_M_PI) ){
	phi = [=](const Double s){ return phi0-diff1 + std::atan( sign_phi * absk * std::tan(s-s0) ); };
	Dphi = [=](const Double s){
	  const Double r = sign_phi * absk;
	  return r / (std::cos(s-s0)*std::cos(s-s0)) / ( 1.0 + r*r * std::tan(s-s0)*std::tan(s-s0) );
	};

	const Double diff3 = phi(ell) - x2.phi();
        // @@@
        std::cout << "diff3 = " << diff3 << std::endl;
	if( isModdable(diff3, _M_PI) ){
          phi = [=](const Double s){ return phi0-diff1+diff3 + std::atan( sign_phi * absk * std::tan(s-s0) ); };
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
  // assert( sign1(1) == sign2(1) );

  const int sign_theta1 = sign1(0);
  const int sign_theta2 = sign2(0);
  const int sign_phi = sign1(1);

  const Double ell = geodesicLength( x1, x2 );

  std::vector<Double> phi0s;
  for(const int sign1 : {-1,1}){
    for(const int sign2 : {-1,1}){
      Double phi0;
      bool is_sol;
      getPhi0WithRatioCheck( phi0, is_sol, x1, x2, sign1, sign2 );
      if(is_sol) phi0s.push_back(phi0);
    }}

  // for(Double phi0 : phi0s) std::cout << "debug. phi0 = " << phi0 << std::endl;

  F theta1, phi1, theta2, phi2, Dtheta1, Dphi1, Dtheta2, Dphi2;
  Double sE;
  bool is_ok = false;
  for(const Double phi0 : phi0s){
    // { Double phi0 = phi0s[1];
    const Double tmp = std::tan(x1.theta())*std::tan(x1.theta()) * std::sin(x1.phi()-phi0)*std::sin(x1.phi()-phi0);
    const Double absk = std::sqrt( 1.0/( 1.0 + 1.0/tmp ) );

    for(int br1=-8; br1<=8; br1++){
      for(int br2=-8; br2<=8; br2++){
	// {{ int br1 = 2; int br2 = -2;
	const Double s01 = -std::atan( sign_phi * std::tan(x1.phi()-phi0)/absk ) + br1*_M_PI;
	const Double s02 = ell-std::atan( sign_phi * std::tan(x2.phi()-phi0)/absk ) + br2*_M_PI;

	sE = Mod( s01 + 0.5*_M_PI + _M_PI ) - _M_PI;
        // std::cout << "debug. sE = " << sE << std::endl;
        // assert( !isModdable(sE) ); // @@@

	theta1 = [=](const Double s){ return std::acos( -sign_theta1 * std::sqrt(1.0-absk*absk) * std::sin(s-s01) ); };
	theta2 = [=](const Double s){ return std::acos( -sign_theta2 * std::sqrt(1.0-absk*absk) * std::sin(s-s02) ); };
	Dtheta1 = [=](const Double s){
	  const Double r = sign_theta1 * std::sqrt(1.0-absk*absk);
	  return r * std::cos(s-s01) / std::sqrt( 1.0 - r*r*std::sin(s-s01)*std::sin(s-s01) );
	};
	Dtheta2 = [=](const Double s){
	  const Double r = sign_theta2 * std::sqrt(1.0-absk*absk);
	  return r * std::cos(s-s02) / std::sqrt( 1.0 - r*r*std::sin(s-s02)*std::sin(s-s02) );
	};

	F phi_tmp1 = [=](const Double s){ return phi0 + std::atan( sign_phi * absk * std::tan(s-s01) ); };
	F phi_tmp2 = [=](const Double s){ return phi0 + std::atan( sign_phi * absk * std::tan(s-s02) ); };

	const Double diff1 = phi_tmp1(0.0) - x1.phi();
	const Double diff2 = phi_tmp2(ell) - x2.phi();

	const Double diff3 = theta1(0.0) - x1.theta();
	const Double diff4 = theta2(ell) - x2.theta();

	// std::cout << "debug. diff1 = " << diff1 << std::endl;
	// std::cout << "debug. diff2 = " << diff2 << std::endl;
	// std::cout << "debug. diff3 = " << diff3 << std::endl;
	// std::cout << "debug. diff4 = " << diff4 << std::endl;

	if( isModdable(diff1, _M_PI) && isModdable(diff2, _M_PI)
	    && isModdable(diff3, 2.0*_M_PI) && isModdable(diff4, 2.0*_M_PI) ){
	  phi1 = [=](const Double s){ return phi0-diff1 + std::atan( sign_phi * absk * std::tan(s-s01) ); };
	  phi2 = [=](const Double s){ return phi0-diff2 + std::atan( sign_phi * absk * std::tan(s-s02) ); };
	  Dphi1 = [=](const Double s){
	    const Double r = sign_phi * absk;
	    return r / (std::cos(s-s01)*std::cos(s-s01)) / ( 1.0 + r*r * std::tan(s-s01)*std::tan(s-s01) );
	  };
	  Dphi2 = [=](const Double s){
	    const Double r = sign_phi * absk;
	    return r / (std::cos(s-s02)*std::cos(s-s02)) / ( 1.0 + r*r * std::tan(s-s02)*std::tan(s-s02) );
	  };

	  const Double diff5 = phi1(0) - x1.phi();
	  const Double diff6 = phi2(ell) - x2.phi();
	  const Double diff7 = phi1(sE) - phi2(sE);

	  const Double diff8 = theta1(0) - x1.theta();
	  const Double diff9 = theta2(ell) - x2.theta();
	  const Double diff10 = theta1(sE) - theta2(sE);


	  std::cout << std::setprecision(16)
                    << "debug. diffs = "
	            << diff5 << ", "
	            << diff6 << ", "
	            << diff7 << ", "
	            << diff8 << ", "
	            << diff9 << ", "
	            << diff10 << std::endl;

	  if( isModdable(diff5, _M_PI)&&isModdable(diff6)&&isModdable(diff7, _M_PI)
	      &&isModdable(diff8)&&isModdable(diff9)&&isModdable(diff10)){
            phi1 = [=](const Double s){ return phi0-diff1-diff5 + std::atan( sign_phi * absk * std::tan(s-s01) ); };
            phi2 = [=](const Double s){ return phi0-diff2-diff7 + std::atan( sign_phi * absk * std::tan(s-s02) ); };

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
  if( isModdable( x2.phi()-x1.phi() ) || isModdable( x1.theta(), _M_PI ) || isModdable( x2.theta(), _M_PI ) ){
    if(is_verbose) std::cout << "delta phi = 0" << std::endl;
    return SolveGeodesicsConstPhi( x1, x2 );
  }

  if( isModdable( x2.phi()-x1.phi()+_M_PI ) ){
    if(is_verbose) std::cout << "delta phi = pi" << std::endl;
    return SolveGeodesicsDeltaPhiEqPi( x1, x2 );
  }

  if( isModdable( x1.theta()-0.5*_M_PI ) || isModdable( x2.theta()-0.5*_M_PI ) ){
    if(is_verbose) std::cout << "theta end = pi/2" << std::endl;
    return SolveGeodesicsEndPtPiHalf( x1, x2 );
  }

  I2 sign1, sign2;
  getSign( sign1, sign2, x1, x2 );
  std::cout << sign1.transpose() << " "
            << sign2.transpose() << std::endl;
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

};
