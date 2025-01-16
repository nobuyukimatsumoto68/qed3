#pragma once

#include <cmath>
#include <vector>

// #include "cg_cuda.h"

struct Zolotarev{
  using FLOAT = double;
  static constexpr FLOAT ONE = 1.0;
  static constexpr FLOAT TWO = 2.0;
  static constexpr FLOAT HALF = ONE/TWO;

  const int n;
  const double k;
  const double kp;
  const int size;

  std::vector<double> c;
  std::vector<double> cp;
  double M;
  double lambda_inv;

  double E;
  std::vector<double> A;


  Zolotarev(const double k_=0.01,
	    const int n_=21 )
    : n(n_)
    , k(k_)
    , kp(std::sqrt(1.0-k*k))
    , size( int(n/2)+1 )
    , c(size, 0.0)
    , cp(size, 0.0)
    , A(size, 0.0)
  {
    get_coeffs();
    partial_fraction();
  }

  // A.D.Kennedy 2004; https://arxiv.org/abs/hep-lat/0402038
  static void sncndnK(const FLOAT z, const FLOAT k,
		      FLOAT& sn, FLOAT& cn, FLOAT& dn,
		      FLOAT& K) {
    FLOAT agm;
    int sgn;
    sn = arithgeom(z, ONE, std::sqrt(ONE - k*k), agm);
    K = M_PI / (TWO * agm);
    sgn = ((int) (std::abs(z) / K)) % 4; /* sgn = 0, 1, 2, 3 */
    sgn ^= sgn >> 1; /* (sgn & 1) = 0, 1, 1, 0 */
    sgn = 1 - ((sgn & 1) << 1); /* sgn = 1, -1, -1, 1 */
    cn = ((FLOAT) sgn) * std::sqrt(ONE - sn * sn);
    dn = std::sqrt(ONE - k*k* sn * sn);
  }

  static FLOAT arithgeom(const FLOAT z, FLOAT a, FLOAT b, FLOAT& agm) {
    static FLOAT pb = -ONE;
    FLOAT xi;
    if (b <= pb) {
      pb = -ONE;
      agm = a;
      return std::sin(z * a);
    }
    pb = b;
    xi = arithgeom(z, HALF*(a+b), std::sqrt(a*b), agm);
    return 2*a*xi / ((a+b) + (a-b)*xi*xi);
  }

  void get_coeffs() {
    double Kp = 1.0, xibar;
    for(int m=0; m<size; m++){
      double sn, cn, dn;
      double z = 2.0 * Kp * m / n;
      sncndnK( z, kp, sn, cn, dn, Kp );
      if(m==0) continue;
      c[m] = - std::pow( cn / sn, 2 );
      z = 2.0 * Kp * (m-0.5) / n;
      sncndnK( z, kp, sn, cn, dn, Kp );
      cp[m] = - std::pow( cn / sn, 2 );
      if(m==1) xibar = 1.0/dn;
    }

    M = 1.0;
    for(int m=1; m<size; m++) M *= (1.0-c[m]) / (1.0-cp[m]);

    lambda_inv = xibar / M;
    for(int m=1; m<size; m++) lambda_inv *= (1.0-c[m]*xibar*xibar) / (1.0-cp[m]*xibar*xibar);
  }

  void partial_fraction() {
    E = 1.0;
    for(int m=1; m<size; m++) {
      double numer = 1.0, denom = 1.0;
      for(int ell=1; ell<size; ell++) {
	numer *= k*k/cp[m] - k*k/c[ell];
	if(m!=ell) denom *= k*k/cp[m] - k*k/cp[ell];
      }
      A[m] = numer/denom;
      E*=c[m]/cp[m];
    }
  }

  double operator[]( const double x ) const {
    double res = 2.0 / (1.0+lambda_inv) * x / (k*M);
    for(int m=1; m<size; m++) res *= (k*k - c[m]*x*x) / (k*k - cp[m]*x*x);
    return res;
  }

  double operator()( const double x ) const {
    double res = 1.0;
    for(int m=1; m<size; m++) res += A[m] / (x*x - k*k/cp[m]);
    res *= E * 2.0 / (1.0+lambda_inv) * x / (k*M);
    return res;
  }

  inline double Delta() const {
    return (lambda_inv - 1.0) / (lambda_inv + 1.0);
  }
};



struct Overlap : private Zolotarev {
  using Gauge=U1onS2;
  using WilsonDirac=Dirac1fonS2;

  static constexpr Idx N = CompilationConst::N;

  SparseDW SDW; // actual data used in M_DW, M_DWH
  CSR M_DW;
  CSR M_DWH;
  double lambda_max;

  Overlap( const WilsonDirac& DW,
	   // const double lambda_max_=12.0,
	   const double k_=0.01,
	   const int n_=21,
	   const bool locate_on_gpu=true)
    : Zolotarev(k_, n_)
    , SDW(DW)
      // , lambda_max(lambda_max_)
  {
    SDW.associate( M_DW, false );
    SDW.associate( M_DWH, true );
  }

  void compute( const Gauge& U ) {
    SDW.update( U );
    compute_lambda_max();
  }

  void compute_lambda_max( const double TOL=1.0e-4, const int MAXITER=500 ) {
    std::vector<Complex> q(N, 0.0);
    std::vector<Complex> x(N, 0.0);

    for(int i=0; i<N; i++) q[i] = (1.0*i+1.0)/N;

    MatPoly Op;
    Op.push_back ( cplx(1.0), {&M_DW, &M_DWH} );

    CuC *d_x, *d_q, *d_scalar; // , *d_tmp, *d_tmp2;
    CUDA_CHECK(cudaMalloc(&d_x, N*CD));
    CUDA_CHECK(cudaMalloc(&d_q, N*CD));
    CUDA_CHECK(cudaMalloc(&d_scalar, CD));
    
    CUDA_CHECK(cudaMemset(d_x, 0, N*CD));
    CUDA_CHECK(cudaMemset(d_q, 0, N*CD));
    CUDA_CHECK(cudaMemset(d_scalar, 0, CD));

    Complex dot;
    double norm=1.0, mu_0=1.0, mu_m1=1.0, mu_m2=1.0;

    CUDA_CHECK(cudaMemcpy(d_q, reinterpret_cast<const CuC*>(q.data()), N*CD, H2D));
    Op.dot2self<N>(norm, d_scalar, d_q);
    for(int i=0; i<N; i++) q[i] = q[i]/std::sqrt(norm);

    double lambda=100.0, lambda_old=1000.0;

    for(int i=0; i<MAXITER; i++){
      Op.from_cpu<N>( x, q );
      CUDA_CHECK(cudaMemcpy(d_x, reinterpret_cast<const CuC*>(x.data()), N*CD, H2D));
      //
      Op.dot2self<N>(norm, d_scalar, d_x);
      for(int i=0; i<N; i++) q[i] = x[i]/std::sqrt(norm);
      CUDA_CHECK(cudaMemcpy(d_q, reinterpret_cast<const CuC*>(q.data()), N*CD, H2D));
      
      Op.dot<N>(reinterpret_cast<CuC&>(dot), d_scalar, d_x, d_q);
      mu_m2=mu_m1;
      mu_m1=mu_0;
      mu_0=dot.real();

      const double r = (mu_0-mu_m1)/(mu_m1-mu_m2);
      const double a = (mu_0-mu_m1)/std::pow(r,i-1)/(r-1);
      lambda_old = lambda;
      lambda = mu_0 - a*std::pow(r,i);

      if(std::abs(lambda_old-lambda)/lambda<TOL) {
#ifdef IsVerbose
	std::clog << "# lambda_max estimate escaped in i = " << i << std::endl;
#endif
	break;
      }
    }

    CUDA_CHECK(cudaFree(d_x));
    CUDA_CHECK(cudaFree(d_q));
    CUDA_CHECK(cudaFree(d_scalar));

    lambda_max = std::sqrt( (1.0+100.0*TOL)*lambda );
  }


  void operator()(std::vector<Complex>& res, const std::vector<Complex>& xi) const {
    res = xi;
    std::vector<Complex> tmp(xi.size());

    for(int m=1; m<size; m++) {
      MatPoly Op;
      Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
      const CuC a = cplx(-k*k/cp[m]);
      Op.push_back ( a, {} );
      Op.solve<N>( tmp, xi );
      for(Idx i=0; i<res.size(); i++) res[i] += A[m] * tmp[i];
    }

    for(Idx i=0; i<res.size(); i++) tmp[i] = E * 2.0 / (1.0+lambda_inv) / (k*M) * res[i];
    MatPoly Op;
    Op.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
    Op.from_cpu<N>( res, tmp );

    for(Idx i=0; i<res.size(); i++) res[i] += xi[i];
  }


};




// struct OverlapPseudoFermion {
//   using Complex = std::complex<double>;
//   using Link = std::array<int,2>; // <int,int>;
//   using Face = std::vector<int>;
//   using Gauge = U1onS2;
//   using Force=U1onS2;
//   using Rng = ParallelRng;
//   using Idx = long int;

//   static constexpr Complex I = Complex(0.0, 1.0);
//   const int NS=2;

//   const Dirac1fonS2& D;
//   const CGCUDA cg;
//   const Sparse& sparse;
//   Zolotarev sgn;

//   std::vector<Complex> phi;
//   std::vector<Complex> eta;

//   OverlapPseudoFermion()=delete;

//   OverlapPseudoFermion(const Dirac1fonS2& D_, const double k_=0.01, const int n_=21)
//     : D(D_)
//     , cg(D)
//     , sparse(cg.sparse)
//     , sgn(k_, n_)
//     , phi(D.lattice.n_sites*NS, 0.0)
//     , eta(D.lattice.n_sites*NS, 0.0)
//   {}

//   Complex operator[](const int i) const { return phi[i]; }
//   Complex& operator[](const int i) { return phi[i]; }

//   Complex operator()(const int ix, const int a) const { return phi[NS*ix+a]; }
//   Complex& operator()(const int ix, const int a) { return phi[NS*ix+a]; }

//   void H( Complex* res, const Complex* v, const U1onS2& U,
// 	  const double lambda_max = 1.0) const {
//     for(Idx i=0; i<sparse.N; i++) res[i] = v[i];

//     for(int m=1; m<sgn.size; m++) {
//       std::vector<Complex> tmp(sparse.N);
//       // cg.general_op_inv( tmp.data(), v, U, 1.0/(lambda_max*lambda_max), sgn.k*sgn.k/sgn.cp[m] );
//       cg.generalH_op_inv( tmp.data(), v, U, 1.0, sgn.k*sgn.k/sgn.cp[m], lambda_max );
//       for(Idx i=0; i<sparse.N; i++) res[i] += sgn.A[m]*tmp[i];
//     }

//     std::vector<Complex> tmp(sparse.N);
//     std::vector<Complex> tmp2(sparse.N);
//     for(Idx i=0; i<sparse.N; i++) tmp[i] = res[i];
//     multHW( tmp2, tmp, U, lambda_max );
//     const double aa = sgn.E * 2.0 / (1.0+sgn.lambda_inv) / (sgn.k*sgn.M);
//     for(Idx i=0; i<sparse.N; i++) res[i] = aa*tmp2[i];
//   }


//   void multD( std::vector<Complex>& Dxi, const std::vector<Complex>& xi, const Gauge& U ) const {
//     assert( Dxi.size()==D.lattice.n_sites*NS );
//     assert( xi.size()==D.lattice.n_sites*NS );

//     Complex D_coo[sparse.len], D_csr[sparse.len];
//     D.coo_format(D_coo, U);
//     sparse.coo2csr( D_csr, D_coo );
//     sparse.mult<Complex>( Dxi.data(), xi.data(), D_csr );
//   }

//   void multHW( std::vector<Complex>& Dxi, const std::vector<Complex>& xi, const Gauge& U,
// 	       const double lambda_max ) const {
//     assert( Dxi.size()==D.lattice.n_sites*NS );
//     assert( xi.size()==D.lattice.n_sites*NS );

//     Complex D_coo[sparse.len], D_csr[sparse.len];
//     D.H_coo_format(D_coo, U, lambda_max);
//     sparse.coo2csr( D_csr, D_coo );
//     sparse.mult<Complex>( Dxi.data(), xi.data(), D_csr );
//   }


//   void multDH( std::vector<Complex>& DHxi, const std::vector<Complex>& xi, const Gauge& U ) const {
//     assert( DHxi.size()==D.lattice.n_sites*NS );
//     assert( xi.size()==D.lattice.n_sites*NS );

//     Complex D_coo[sparse.len], D_csrH[sparse.len];
//     D.coo_format(D_coo, U);
//     sparse.coo2csrH( D_csrH, D_coo );
//     sparse.multT<Complex>( DHxi.data(), xi.data(), D_csrH );    
//   }


//   void multDHD( std::vector<Complex>& DHDxi, const std::vector<Complex>& xi, const Gauge& U ) const {
//     assert( DHDxi.size()==D.lattice.n_sites*NS );
//     assert( xi.size()==D.lattice.n_sites*NS );

//     Complex D_coo[sparse.len], D_csr[sparse.len], D_csrH[sparse.len];
//     D.coo_format(D_coo, U);
//     sparse.coo2csr_csrH( D_csr, D_csrH, D_coo );

//     std::vector<Complex> tmp( D.lattice.n_sites*NS );
//     sparse.mult<Complex>( tmp.data(), xi.data(), D_csr );
//     sparse.multT<Complex>( DHDxi.data(), tmp.data(), D_csrH );    
//   }


// //   void gen( const Gauge& U, Rng& rng ) {
// //     const int N = D.lattice.n_sites*NS;
// //     std::vector<Complex> xi(N, 0.0);

// //     for(int ix=0; ix<D.lattice.n_sites; ix++) for(int a=0; a<NS; a++) xi[NS*ix+a] = ( rng.gaussian_site(ix) + I*rng.gaussian_site(ix) ) / std::sqrt(2.0);

// //     multDH( phi, xi, U );

// //     update_eta(U);
// //   }


// //   void update_eta( const Gauge& U ) { cg( eta.data(), phi.data(), U ); }


// //   auto begin(){ return phi.begin(); }
// //   auto end(){ return phi.end(); }
// //   auto begin() const { return phi.begin(); }
// //   auto end() const { return phi.end(); }


// //   Complex dot( const std::vector<Complex>& eta1, const std::vector<Complex>& xi) const {
// //     assert( eta1.size()==xi.size() );
// //     Complex res = 0.0;
// //     for(int i=0; i<eta1.size(); i++) res += std::conj(eta1[i]) * xi[i];
// //     return res;
// //   }

// //   Complex dot( const std::vector<Complex>& eta1 ) const {
// //     return dot(this->phi, eta1);
// //   }

// //   double S() const { return dot( eta ).real(); }


// //   double get_force( const Gauge& U, const Link& ell ) const {
// //     const int N = D.lattice.n_sites*NS;

// //     std::vector<Complex> dD;
// //     std::vector<int> is;
// //     std::vector<int> js;
// //     D.d_coo_format( dD, is, js, U, ell );

// //     std::vector<Complex> dD_eta(N);
// //     cg.sparse.multcoo( dD_eta, eta, dD, is, js );

// //     std::vector<Complex> DH_dD_eta(N);
// //     multDH( DH_dD_eta, dD_eta, U );

// //     return -2.0 * dot( eta, DH_dD_eta ).real();
// //   }


// //   Force dS( const Gauge& U ) const {
// //     Force pi( U.lattice ); // 0 initialized
// // #ifdef _OPENMP
// // #pragma omp parallel for num_threads(nparallel)
// // #endif
// //     for(int ell=0; ell<U.lattice.n_links; ell++) pi[ell] = get_force( U, U.lattice.links[ell] );
// //     return pi;
// //   }


// };
