#include <iostream>
#include <cmath>
// #include <math.h>
#include <vector>
#include <limits>
#include <cstdlib>
#include <cassert>
#include <functional>
#include <numbers>

#include <cstdint>
// #include <Eigen/Dense>


namespace Geodesic{


  // using Double = long double;
  using Double = double;
  // using Idx = std::size_t; // std::int32_t;
  using Idx = std::int32_t;
// using Complex = std::complex<Double>;

// constexpr Double TOLLOOSE=1.0e-10;
// constexpr Double TOLLOOSE=1.0e-6;

  // constexpr Double TOLMOD=1.0e-7;
  constexpr Double TOLMOD=1.0e-6;
  constexpr Double TOLPHI = 1.0e-8;
  constexpr Double EPSNUMDER=1.0e-5;

  constexpr Double _M_PI = std::numbers::pi;

  constexpr int BRMAX = 8;



  constexpr Idx DIM = 2;
  constexpr Idx EDIM = 3;

  using F=std::function<Double(const Double)>;
  using VF=std::array<F, 2*DIM>;

  constexpr Double ZERO = 0.0;
  constexpr Double ONE = 1.0;
  constexpr Double TWO = 2.0;
  constexpr Double HALF = 0.5;



// using I2=Eigen::Vector2i;

// constexpr int DIM = 2;
// // using V2=Eigen::Vector2d;
// using V2=Eigen::Matrix<Double, 1, 2>;

// constexpr int EDIM = 3;
// using V3=Eigen::Matrix<Double, 1, 3>;
// // using V3=Eigen::Vector3d;

  template<typename T, int N>
  struct V {
    std::vector<T> v;

    V()
      : v(N)
    {}

    V(const std::vector<T>& other)
      : v(other)
    {
      // for(T elem : {args...}) v.push_back(elem);
      // assert(v.size()==N);
    }

    V(const std::initializer_list<T>& other)
      : v(std::vector<T>(other))
    {
      // for(T elem : {args...}) v.push_back(elem);
      // assert(v.size()==N);
    }

    template<typename... Args>
    V(Args... args)
    {
      for(T elem : {args...}) v.push_back(elem);
      assert(v.size()==N);
    }

    auto begin(){ return v.begin(); }
    auto end(){ return v.end(); }
    auto begin() const { return v.begin(); }
    auto end() const { return v.end(); }


    template<typename... Args>
    void operator << (Args... args)
    {
      v.clear();
      for(const T& elem : {args...}) v.push_back(elem);
      assert(v.size()==N);
    }

    Double norm() const {
      Double res = ZERO;
      for(const T& elem : v) res += std::abs(elem)*std::abs(elem);
      res = std::sqrt(res);
      return res;
    }

    static V<T,N> Zero(){
      V<T,N> res;
      for(T& elem : res.v) elem = ZERO;
      return res;
    }

    inline void push_back( const T elem ) { v.push_back(elem); }

    // V<int, N> sign() const {
    //   V<int, N> res;
    //   for(const T& elem : v) {
    //     int tmp=0;
    //     if(elem>=0) tmp = 1;
    //     else tmp = -1;
    //     res.push_back(tmp);
    //   }
    //   return res;
    // }

    V<Double, N> sign() const {
      V<Double, N> res;
      for(Idx i=0; i<N; i++) {
        Double tmp=ZERO;
        if(v[i]>=0) tmp = 1.;
        else tmp = -1.;
        res[i] = tmp;
      }
      return res;
    }

    Double dot(const V& rhs) const {
      Double res = ZERO;
      for(Idx i=0; i<N; i++) res += v[i]*rhs[i];
      return res;
    }

    Double mean() const {
      Double res = ZERO;
      for(Idx i=0; i<N; i++) res += v[i];
      res /= N;
      return res;
    }

    V<T,N> cross(const V& rhs) const {
      static_assert( N==3 );
      V<T, N> res;
      res[0] = v[1]*rhs[2] - v[2]*rhs[1];
      res[1] = v[2]*rhs[0] - v[0]*rhs[2];
      res[2] = v[0]*rhs[1] - v[1]*rhs[0];
      return res;
    }

    T operator[](const Idx i) const { return v[i]; }
    T& operator[](const Idx i) { return v[i]; }
    T operator()(const Idx i) const { return v[i]; }
    T& operator()(const Idx i) { return v[i]; }

    V& operator+=(const V& rhs){
      for(Idx i=0; i<v.size(); i++) v[i] += rhs.v[i];
      return *this;
    }
    friend V operator+(V v, const V& w) { v += w; return v; }

    V operator-(){
      V<T,N> res;
      for(Idx i=0; i<v.size(); i++) res.v[i] = -v[i];
      return res;
    }


    V& operator-=(const V& rhs){
      for(Idx i=0; i<v.size(); i++) v[i] -= rhs.v[i];
      return *this;
    }
    friend V operator-(V v, const V& w) { v -= w; return v; }

    template<typename R>
    V& operator*=(const R rhs){
      for(Idx i=0; i<v.size(); i++) v[i] *= rhs;
      return *this;
    }
    template<typename R>
    friend V operator*(V v, const R a) { v *= a; return v; }
    template<typename R>
    friend V operator*(const R a, V v) { v *= a; return v; }

    template<typename R>
    V& operator/=(const R rhs){
      for(Idx i=0; i<v.size(); i++) v[i] /= rhs;
      return *this;
    }
    template<typename R>
    friend V operator/(V v, const R a) { v /= a; return v; }
    template<typename R>
    friend V operator/(const R a, V v) { v /= a; return v; }

    friend std::ostream& operator <<(std::ostream& s, const V& v)
    {
      for(Idx i=0; i<N; i++) s << v[i] << " ";
      return s;
    }

  };


  using V2 = V<Double, 2>;
  using V3 = V<Double, 3>;


  // using I2=Eigen::Vector2i;
  // using I2 = V<Idx, 2>;

  // // using V2=Eigen::Vector2d;
  // using V2=Eigen::Matrix<Double, 1, 2>;

  // using V3=Eigen::Matrix<Double, 1, 3>;
  // // using V3=Eigen::Vector3d;




Double Mod(Double a, Double b=TWO*_M_PI){
  int p = int(std::floor(a / b));
  // if(a<0) p += 1;
  Double r = a - p*b;
  return r;
  // constexpr Double ten = 10.0;
  // // return std::fmodl(a+ten*b, b);
  // return fmodl(a+ten*b, b);
}

Double Mod2(Double a){
  Double tmp = Mod(a);
  if(tmp>_M_PI) tmp -= TWO*_M_PI;
  return tmp;
}

// template<typename V>
// V Mod2(const V& a){
//   Double b=TWO*M_PI;

//   V res = a;
//   for(auto& elem : res){
//     // int p = int(std::floor((elem+TWO*b) / b));
//     // int p = int(std::floor((elem+TWO*b) / b));
//     // // if(a<0) p += 1;
//     // Double r = elem - (p+2)*b;
//     // // if(r>M_PI) r -= TWO*M_PI;
//     elem = Mod(elem);
//   }
//   return res;
// }


Double sgn(const Double a){
  Double res = ONE;
  if(a<0) res *= -1;
  return res;
}


bool isModdable( const Double v, const Double q=TWO*_M_PI, const Double TOL=TOLMOD ){
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
  // return V3(std::sin(xi[0])*std::cos(xi[1]),
  //           std::sin(xi[0])*std::sin(xi[1]),
  //           std::cos(xi[0]));
  return V3(std::vector<Double>({std::sin(xi[0])*std::cos(xi[1]),
        std::sin(xi[0])*std::sin(xi[1]),
        std::cos(xi[0])})
    );
}


Double _acos(const Double arg){
  Double res;
  if(arg<-ONE) {
    // std::clog << "# cutoff found" << std::endl;
    res=_M_PI;
  }
  else if(arg>ONE) {
    // std::clog << "# cutoff found" << std::endl;
    res=ZERO;
  }
  else {
    res = std::acos( arg );
  }
  return res;
}



V2 projectionS2( const V3& x ){
  const Double r = x.norm();
  const Double theta = _acos(x(2)/r);
  Double phi = ZERO;
  if(!isModdable(theta, _M_PI)){
    const Double arg = x(0)/(r*std::sin(theta));
    phi = _acos(arg);
    if(x(1)<ZERO) phi *= -ONE;
  }
  return V2(std::vector<Double>({theta, phi}));
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
    : x(V3(std::vector<Double>({x1,x2,x3})))
    , xi(projectionS2(x))
    , is_singular( isModdable(xi(0), _M_PI) )
  {}

  explicit Pt( const Double xi1,
               const Double xi2 )
    : x(embedding3D(V2(std::vector<Double>({xi1,xi2}))))
    , xi(V2(std::vector<Double>({xi1,xi2})))
    , is_singular( isModdable(xi(0), _M_PI) )
  {}


  Double theta() const { return xi(0); }
  Double phi() const {
    assert( !is_singular );
    return xi(1);
  }

  V3 Delta0() const {
    return V3(std::vector<Double>({std::cos(xi[0])*std::cos(xi[1]),
          std::cos(xi[0])*std::sin(xi[1]),
          -std::sin(xi[0])})
      );
  }

  V3 Delta1() const {
    Double zero = ZERO;
    return V3(std::vector<Double>({-std::sin(xi[0])*std::sin(xi[1]),
          std::sin(xi[0])*std::cos(xi[1]),
          zero}));
  }

  V3 e0() const {
    return Delta0();
  }

  V3 e1() const {
    return Delta1()/std::sin(xi[0]);
  }

  Double operator[](const Idx mu) const { return x(mu); }
};


Double geodesicLength(const Pt& x1, const Pt& x2){
  const Double inner = x1.x.dot( x2.x );
  return _acos( inner );
}

Double sphericalarea( const Pt& x1, const Pt& x2, const Pt& x3 ){
  Double a = geodesicLength( x2, x3 );
  Double b = geodesicLength( x3, x1 );
  Double c = geodesicLength( x1, x2 );
  Double s = HALF*(a+b+c);
  Double tantan = std::tan(HALF*s) * std::tan(HALF*(s-a)) * std::tan(HALF*(s-b)) * std::tan(HALF*(s-c));
  return 4.0*std::atan(std::sqrt(tantan));
}



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



void getSign(V2& sign1, V2& sign2,
             const Pt& x1, const Pt& x2, const Double eps=EPSNUMDER){
  const V3 p = x1.x;
  const V3 q = x2.x;

  const V3 p2 = embedding3D(projectionS2(p + eps*(q-p)/(q-p).norm()));
  const V3 diff1 = p2-p;
  const V2 deriv1( {diff1.dot( x1.Delta0() ), diff1.dot( x1.Delta1() )});
  const V3 q2 = embedding3D(projectionS2(q - eps*(q-p)));
  const V3 diff2 = q-q2;
  const V2 deriv2( {diff2.dot( x2.Delta0() ), diff2.dot( x2.Delta1() )});

  sign1 = deriv1.sign();
  sign2 = deriv2.sign();
}


Double getCosPhi0(const Pt& x1, const Pt& x2, const Double sign1, const Double sign2){
  const Double t12 = std::tan( x1.theta() )/std::tan( x2.theta() );
  const Double numer = sign2 * ( t12*std::cos( x1.phi() ) - sign1*std::cos( x2.phi() ) );
  const Double denom = std::sqrt( t12*t12 + ONE - TWO*sign1*t12*std::cos( x1.phi() - x2.phi() ) );
  return numer/denom;
}


Double getPhi0(const Pt& x1, const Pt& x2, const Double sign1, const Double sign2){
  return _acos( getCosPhi0(x1, x2, sign1, sign2) );
}

void getPhi0WithRatioCheck(Double& phi0, bool& is_sol,
			   const Pt& x1, const Pt& x2,
			   const Double sign1, const Double sign2,
			   const Double TOL=TOLPHI){
  phi0 = getPhi0(x1, x2, sign1, sign2);
  const Double ratio = sign1*std::tan(x1.theta())*std::sin(x1.phi()-phi0) / ( std::tan(x2.theta())*std::sin(x2.phi()-phi0) );
  is_sol = std::abs(std::abs(ratio)-ONE) < TOL;
}


Sol SolveGeodesicsConstPhi( const Pt& x1, const Pt& x2 ){
  assert( isModdable( x1.phi()-x2.phi() ) );

  const Double ell = std::abs( x2.theta()-x1.theta() );

  F theta = [=](const Double s){ return x1.theta() + (x2.theta()-x1.theta()) * s/ell; };
  F phi = [=](const Double s){ return x1.phi(); };
  F Dtheta = [=](const Double s){ return (x2.theta()-x1.theta()) * ONE/ell; };
  F Dphi = [=](const Double s){ return ZERO; };

  return Sol( VF{theta, phi, Dtheta, Dphi}, ell, x1, x2 );
}


Sol SolveGeodesicsDeltaPhiEqPi( const Pt& x1, const Pt& x2 ){
  assert( (!isModdable( x1.phi()-x2.phi(), TWO*_M_PI )) && isModdable( x1.phi()-x2.phi(), _M_PI ) );

  const Double ell = geodesicLength( x1, x2 );

  const Double sN = geodesicLength( x1, Pt(ZERO, ZERO, ONE) );
  const Double sS = geodesicLength( x2, Pt(ZERO, ZERO, -ONE) );

  Double sE, thetaE;
  if(sN<sS){
    sE = sN;
    thetaE = ZERO;
  }
  else{
    sE = sS;
    thetaE = _M_PI;
  }

  F theta1 = [=](const Double s){ return x1.theta() + (thetaE-x1.theta()) * s/sE; };
  F theta2 = [=](const Double s){ return thetaE + (x2.theta()-thetaE) * (s-sE)/(ell-sE); };
  F phi1 = [=](const Double s){ return x1.phi(); };
  F phi2 = [=](const Double s){ return x2.phi(); };

  F Dtheta1 = [=](const Double s){ return (thetaE-x1.theta()) * ONE/sE; };
  F Dtheta2 = [=](const Double s){ return (x2.theta()-thetaE) * ONE/(ell-sE); };
  F Dphi1 = [=](const Double s){ return ZERO; };
  F Dphi2 = [=](const Double s){ return ZERO; };

  return Sol( VF{theta1, phi1, Dtheta1, Dphi1},
	      VF{theta2, phi2, Dtheta2, Dphi2},
	      ell, sE, x1, x2 );
}


Sol SolveGeodesicsEndPtPiHalf( const Pt& x1, const Pt& x2 ){
  assert( isModdable( x1.theta()-HALF*_M_PI ) || isModdable( x2.theta()-HALF*_M_PI ) );

  const Double ell = geodesicLength( x1, x2 );

  Double phi0, phip, thetap, s0;
  if( isModdable( x1.theta()-HALF*_M_PI ) ){
    phi0 = x1.phi();
    phip = x2.phi();
    thetap = x2.theta();
    s0 = ZERO;
  }
  else{
    phi0 = x2.phi();
    phip = x1.phi();
    thetap = x1.theta();
    s0 = ell;
  }

  const Double sign_phi = sgn( Mod( x2.phi()-x1.phi()+_M_PI ) - _M_PI );
  const Double sign_theta = sgn( x2.theta()-x1.theta() );

  const Double tmp = std::tan(thetap)*std::tan(thetap) * std::sin(phip-phi0)*std::sin(phip-phi0);
  const Double absk = std::sqrt( ONE/( ONE + ONE/tmp ) );

  F theta = [=](const Double s){ return _acos( -sign_theta * std::sqrt(ONE-absk*absk) * std::sin(s-s0) ); };
  F phi = [=](const Double s){ return phi0 + std::atan( sign_phi * absk * std::tan(s-s0) ); };

  F Dtheta = [=](const Double s){
    const Double r = sign_theta * std::sqrt(ONE-absk*absk);
    return r * std::cos(s-s0) / std::sqrt( ONE - r*r*std::sin(s-s0)*std::sin(s-s0) );
  };
  F Dphi = [=](const Double s){
    const Double r = sign_phi * absk;
    return r / (std::cos(s-s0)*std::cos(s-s0)) / ( ONE + r*r * std::tan(s-s0)*std::tan(s-s0) );
  };

  return Sol( VF{theta, phi, Dtheta, Dphi}, ell, x1, x2 );
}


Sol SolveGeodesicsMonotonic( const Pt& x1, const Pt& x2 ){
  V2 sign1, sign2;
  getSign( sign1, sign2, x1, x2 );
  assert( std::abs(sign1(0)-sign2(0))<1.0e-10 );
  assert( std::abs(sign1(1)-sign2(1))<1.0e-10 );
  // assert( sign1(1) == sign2(1) );

  const Double sign_theta = sign1(0);
  const Double sign_phi = sign1(1);

  const Double ell = geodesicLength( x1, x2 );

  std::vector<Double> phi0s;
  for(const Double sign1 : {-ONE,ONE}){
    for(const Double sign2 : {-ONE,ONE}){
      Double phi0;
      bool is_sol;
      getPhi0WithRatioCheck( phi0, is_sol, x1, x2, sign1, sign2 );
      if(is_sol) phi0s.push_back(phi0);
    }}

  F theta, phi, Dtheta, Dphi;
  bool is_ok = false;
  for(const Double phi0 : phi0s){
    const Double tmp = std::tan(x1.theta())*std::tan(x1.theta()) * std::sin(x1.phi()-phi0)*std::sin(x1.phi()-phi0);
    const Double absk = std::sqrt( ONE/( ONE + ONE/tmp ) );

    for(int br=-BRMAX; br<=BRMAX; br++){
      const Double s0 = -std::atan( sign_phi * std::tan(x1.phi()-phi0)/absk ) + br*_M_PI;

      theta = [=](const Double s){ return _acos( -sign_theta * std::sqrt(ONE-absk*absk) * std::sin(s-s0) ); };
      Dtheta = [=](const Double s){
	const Double r = sign_theta * std::sqrt(ONE-absk*absk);
	return r * std::cos(s-s0) / std::sqrt( ONE - r*r*std::sin(s-s0)*std::sin(s-s0) );
      };
      F phi_tmp = [=](const Double s){ return phi0 + std::atan( sign_phi * absk * std::tan(s-s0) ); };

      const Double diff1 = phi_tmp(ZERO) - x1.phi();
      const Double diff2 = theta(ell) - x2.theta();

      // std::cout << "diff1 = " << diff1 << std::endl;
      // std::cout << "diff2 = " << diff2 << std::endl;
      // std::cout << "diff1 = " << isModdable(diff1, _M_PI) << std::endl;
      // std::cout << "diff2 = " << isModdable(diff2, TWO*_M_PI) << std::endl;

      if( isModdable(diff1, _M_PI) && isModdable(diff2, TWO*_M_PI) ){
	phi = [=](const Double s){ return phi0-diff1 + std::atan( sign_phi * absk * std::tan(s-s0) ); };
	Dphi = [=](const Double s){
	  const Double r = sign_phi * absk;
	  return r / (std::cos(s-s0)*std::cos(s-s0)) / ( ONE + r*r * std::tan(s-s0)*std::tan(s-s0) );
	};

	const Double diff3 = phi(ell) - x2.phi();
        // std::cout << "diff3 = " << diff3 << std::endl;
        // std::cout << "diff3 = " << isModdable(diff3, TWO*_M_PI) << std::endl;
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
  V2 sign1, sign2;
  getSign( sign1, sign2, x1, x2 );
  assert( std::abs(sign1(0)+sign2(0))<1.0e-10 );
  assert( std::abs(sign1(1)-sign2(1))<1.0e-10 );
  // assert( sign1(0) == -sign2(0) );
  // assert( sign1(1) == sign2(1) );

  const Double sign_theta1 = sign1(0);
  const Double sign_theta2 = sign2(0);
  const Double sign_phi = sign1(1);

  const Double ell = geodesicLength( x1, x2 );

  std::vector<Double> phi0s;
  for(const Double sign1 : {-ONE,ONE}){
    for(const Double sign2 : {-ONE,ONE}){
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
    const Double absk = std::sqrt( ONE/( ONE + ONE/tmp ) );

    for(int br1=-BRMAX; br1<=BRMAX; br1++){
      for(int br2=-BRMAX; br2<=BRMAX; br2++){
	// {{ int br1 = 2; int br2 = -2;
	const Double s01 = -std::atan( sign_phi * std::tan(x1.phi()-phi0)/absk ) + br1*_M_PI;
	const Double s02 = ell-std::atan( sign_phi * std::tan(x2.phi()-phi0)/absk ) + br2*_M_PI;

	sE = Mod( s01 + HALF*_M_PI + _M_PI ) - _M_PI;
	// std::cout << "debug. sE = " << sE << std::endl;
	assert( !isModdable(sE) );

	theta1 = [=](const Double s){ return _acos( -sign_theta1 * std::sqrt(ONE-absk*absk) * std::sin(s-s01) ); };
	theta2 = [=](const Double s){ return _acos( -sign_theta2 * std::sqrt(ONE-absk*absk) * std::sin(s-s02) ); };
	Dtheta1 = [=](const Double s){
	  const Double r = sign_theta1 * std::sqrt(ONE-absk*absk);
	  return r * std::cos(s-s01) / std::sqrt( ONE - r*r*std::sin(s-s01)*std::sin(s-s01) );
	};
	Dtheta2 = [=](const Double s){
	  const Double r = sign_theta2 * std::sqrt(ONE-absk*absk);
	  return r * std::cos(s-s02) / std::sqrt( ONE - r*r*std::sin(s-s02)*std::sin(s-s02) );
	};

	F phi_tmp1 = [=](const Double s){ return phi0 + std::atan( sign_phi * absk * std::tan(s-s01) ); };
	F phi_tmp2 = [=](const Double s){ return phi0 + std::atan( sign_phi * absk * std::tan(s-s02) ); };

	const Double diff1 = phi_tmp1(ZERO) - x1.phi();
	const Double diff2 = phi_tmp2(ell) - x2.phi();

	const Double diff3 = theta1(ZERO) - x1.theta();
	const Double diff4 = theta2(ell) - x2.theta();

	// std::cout << "debug. diff1 = " << diff1 << std::endl;
	// std::cout << "debug. diff2 = " << diff2 << std::endl;
	// std::cout << "debug. diff3 = " << diff3 << std::endl;
	// std::cout << "debug. diff4 = " << diff4 << std::endl;

	if( isModdable(diff1, _M_PI) && isModdable(diff2, _M_PI)
	    && isModdable(diff3, TWO*_M_PI) && isModdable(diff4, TWO*_M_PI) ){
	  phi1 = [=](const Double s){ return phi0-diff1 + std::atan( sign_phi * absk * std::tan(s-s01) ); };
	  phi2 = [=](const Double s){ return phi0-diff2 + std::atan( sign_phi * absk * std::tan(s-s02) ); };
	  Dphi1 = [=](const Double s){
	    const Double r = sign_phi * absk;
	    return r / (std::cos(s-s01)*std::cos(s-s01)) / ( ONE + r*r * std::tan(s-s01)*std::tan(s-s01) );
	  };
	  Dphi2 = [=](const Double s){
	    const Double r = sign_phi * absk;
	    return r / (std::cos(s-s02)*std::cos(s-s02)) / ( ONE + r*r * std::tan(s-s02)*std::tan(s-s02) );
	  };

	  const Double diff5 = phi1(ZERO) - x1.phi();
	  const Double diff6 = phi2(ell) - x2.phi();
	  const Double diff7 = phi1(sE) - phi2(sE);

	  const Double diff8 = theta1(ZERO) - x1.theta();
	  const Double diff9 = theta2(ell) - x2.theta();
	  const Double diff10 = theta1(sE) - theta2(sE);

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

  if( isModdable( x2.phi()-x1.phi()+_M_PI ) ){
    if(is_verbose) std::cout << "delta phi = pi" << std::endl;
    return SolveGeodesicsDeltaPhiEqPi( x1, x2 );
  }

  if( isModdable( x1.theta()-HALF*_M_PI ) || isModdable( x2.theta()-HALF*_M_PI ) ){
    if(is_verbose) std::cout << "theta end = pi/2" << std::endl;
    return SolveGeodesicsEndPtPiHalf( x1, x2 );
  }

  V2 sign1, sign2;
  getSign( sign1, sign2, x1, x2 );
  // std::cout << sign1 << " "
  //           << sign2 << std::endl;

  // std::cout << sign1.transpose() << " "
  //           << sign2.transpose() << std::endl;
  // assert(  );
  // assert( std::abs(sign1(1)-sign2(1))<ONEe-10 );

  if(std::abs(sign1(0)-sign2(0))<1.0e-10){
    if(is_verbose) std::cout << "monotonic" << std::endl;
    return SolveGeodesicsMonotonic( x1, x2 );
  }
  else{
    if(is_verbose) std::cout << "altering" << std::endl;
    return SolveGeodesicsAltering( x1, x2 );
  }

  assert(false);
}

}
