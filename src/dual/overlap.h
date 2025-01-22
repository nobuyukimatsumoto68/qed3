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
  double C;
  std::vector<double> A;

  // double *d_c, *d_A, *d_C;

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
    C = E * 2.0 / (1.0+lambda_inv) / (k*M);

    // CUDA_CHECK(cudaMalloc(&d_c, size*DB));
    // CUDA_CHECK(cudaMalloc(&d_A, size*DB));
    // CUDA_CHECK(cudaMalloc(&d_C, DB));

    // A[0] = 1.0;
    // CUDA_CHECK(cudaMemcpy(d_c, c.data(), size*DB, H2D));
    // CUDA_CHECK(cudaMemcpy(d_A, A.data(), size*DB, H2D));
    // CUDA_CHECK(cudaMemcpy(d_C, &C, DB, H2D));
  }

  ~Zolotarev(){
    // CUDA_CHECK(cudaFree(d_c));
    // CUDA_CHECK(cudaFree(d_A));
    // CUDA_CHECK(cudaFree(d_C));
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
    // res *= E * 2.0 / (1.0+lambda_inv) * x / (k*M);
    res *= C * x;
    return res;
  }

  inline double Delta() const {
    return (lambda_inv - 1.0) / (lambda_inv + 1.0);
  }
};



struct Overlap : public Zolotarev {
  using Gauge=U1onS2;
  using WilsonDirac=Dirac1fonS2;
  using Link = std::array<Idx,2>; // <int,int>;

  static constexpr Idx N = CompilationConst::N;

  const WilsonDirac& DW;
  DWDevice d_DW; // actual data used in M_DW, M_DWH
  CSR M_DW;
  CSR M_DWH;
  double lambda_max;

  Overlap( const WilsonDirac& DW_,
	   // const double lambda_max_=12.0,
	   const double k_=0.01,
	   const int n_=21,
	   const bool locate_on_gpu=true)
    : Zolotarev(k_, n_)
    , DW(DW_)
    , d_DW(DW)
      // , lambda_max(lambda_max_)
  {
    d_DW.associateCSR( M_DW, false );
    d_DW.associateCSR( M_DWH, true );
  }

  void compute( const Gauge& U ) {
    d_DW.update( U );
    compute_lambda_max();
  }

  void compute_lambda_max( const double TOL=1.0e-8, const int MAXITER=500 ) {
    std::vector<Complex> q(N, 0.0);

    q[0] = std::sqrt(1.0/2.0);
    q[1] = std::sqrt(1.0/2.0);

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

    double lambda=100.0, lambda_old=1000.0;

    for(int i=0; i<MAXITER; i++){
      Op.on_gpu<N>( d_x, d_q );
      //
      Op.dot2self<N>(&norm, d_x);
      CUDA_CHECK(cudaMemcpy(d_q, d_x, N*CD, H2D));
      Op.Zdscal<N>( 1.0/std::sqrt(norm), d_q );

      Op.dot<N>(reinterpret_cast<CuC*>(&dot), d_x, d_q);
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


  void mult_device(CuC* d_res, const CuC* d_xi) const {
    CUDA_CHECK(cudaMemcpy(d_res, d_xi, N*CD, D2D));

    CuC* d_tmp;
    CUDA_CHECK(cudaMalloc(&d_tmp, N*CD));

    for(int m=1; m<size; m++) {
      MatPoly Op;
      Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
      const CuC a = cplx(-k*k/cp[m]);
      Op.push_back ( a, {} );
      Op.solve<N>( d_tmp, d_xi );
      Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, A[m], d_tmp, d_res);
    }

    MatPoly Op;
    CUDA_CHECK(cudaMemcpy(d_tmp, d_res, N*CD, D2D));
    Op.Zdscal<N>(C, d_tmp);

    Op.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
    Op.on_gpu<N>( d_res, d_tmp );
    Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_xi, d_res); // A[0]=1.0

    CUDA_CHECK(cudaFree(d_tmp));
  }


  void adj_device(CuC* d_res, const CuC* d_xi) const {
    CuC* d_DHxi;
    CUDA_CHECK(cudaMalloc(&d_DHxi, N*CD));

    CUDA_CHECK(cudaMemcpy(d_res, d_xi, N*CD, D2D));
    CUDA_CHECK(cudaMemcpy(d_DHxi, d_xi, N*CD, D2D));
    {
      MatPoly Op;
      Op.push_back ( cplx(1.0/(lambda_max)), {&M_DWH} );
      Op.on_gpu<N>( d_DHxi, d_xi );
      Op.Zdscal<N>( C, d_DHxi );
    }
    CUDA_CHECK(cudaMemcpy(d_res, d_DHxi, N*CD, D2D));

    CuC* d_tmp;
    CUDA_CHECK(cudaMalloc(&d_tmp, N*CD));

    for(int m=1; m<size; m++) {
      MatPoly Op;
      Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
      const CuC a = cplx(-k*k/cp[m]);
      Op.push_back ( a, {} );
      Op.solve<N>( d_tmp, d_DHxi );
      Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, A[m], d_tmp, d_res);
    }

    Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_xi, d_res); // A[0]=1.0

    CUDA_CHECK(cudaFree(d_tmp));
    CUDA_CHECK(cudaFree(d_DHxi));
  }



  void sq_device( CuC* d_res, const CuC* d_xi) const {
    CuC *d_tmp1, *d_tmp2;
    CUDA_CHECK(cudaMalloc(&d_tmp1, N*CD));
    CUDA_CHECK(cudaMalloc(&d_tmp2, N*CD));

    this->mult_device(d_tmp1, d_xi);
    this->adj_device(d_tmp2, d_xi);

    CUDA_CHECK(cudaMemcpy(d_res, d_tmp1, N*CD, D2D));
    Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_tmp2, d_res); // A[0]=1.0

    CUDA_CHECK(cudaFree(d_tmp1));
    CUDA_CHECK(cudaFree(d_tmp2));
  }



  void mult_device3(CuC* d_res, const CuC* d_xi) const {
    const int mmax = size; // size;
    std::vector<CuC*> d_Zs(mmax);
    for(int m=1; m<mmax; m++) CUDA_CHECK(cudaMalloc(&d_Zs[m], N*CD));

    for(int m=1; m<mmax; m++) {
      MatPoly Op;
      Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );

      const CuC a = cplx(-k*k/cp[m]);
      Op.push_back ( a, {} );
      Op.solve<N>( d_Zs[m], d_xi );
      Op.Zdscal<N>( A[m], d_Zs[m] );
    }

    CUDA_CHECK(cudaMemcpy(d_res, d_xi, N*CD, D2D));
    for(int m=1; m<mmax; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_Zs[m], d_res);

    CuC* d_tmp; CUDA_CHECK(cudaMalloc(&d_tmp, N*CD));
    MatPoly Op; Op.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
    CUDA_CHECK(cudaMemcpy(d_tmp, d_res, N*CD, D2D));
    Op.on_gpu<N>( d_res, d_tmp );

    Op.Zdscal<N>(C, d_res);
    Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_xi, d_res); // 1+V

    for(int m=1; m<mmax; m++) CUDA_CHECK(cudaFree(d_Zs[m]));
    CUDA_CHECK(cudaFree(d_tmp));
  }


  double grad_device( const Link& link, const Gauge& U, CuC* d_eta ) const {
    CuC inner;
    double res = 0.0;

    std::vector<CuC*> d_Ys(size), d_Zs(size);
    for(int m=1; m<size; m++) {
      CUDA_CHECK(cudaMalloc(&d_Zs[m], N*CD));
      CUDA_CHECK(cudaMalloc(&d_Ys[m], N*CD));
    }

    MatPoly X;     X.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );

    CuC* d_Xdag_eta; CUDA_CHECK(cudaMalloc(&d_Xdag_eta, N*CD));
    {
      MatPoly XH;   XH.push_back ( cplx(1.0/(lambda_max)), {&M_DWH} );
      CUDA_CHECK(cudaMemcpy(d_Xdag_eta, d_eta, N*CD, D2D)); // @@@@@@
      XH.on_gpu<N>(d_Xdag_eta, d_eta);
    }

    for(int m=1; m<size; m++) {
      MatPoly Op;
      Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
      const CuC a = cplx(-k*k/cp[m]); Op.push_back ( a, {} );

      Op.solve<N>( d_Zs[m], d_eta );
      Op.solve<N>( d_Ys[m], d_Xdag_eta );
    }

    COO coo;
    DW.d_coo_format(coo.en, U, link);
    coo.do_it();

    {
      CuC *d_XZm, *d_dD_Ym;
      CUDA_CHECK(cudaMalloc(&d_XZm, N*CD));
      CUDA_CHECK(cudaMalloc(&d_dD_Ym, N*CD));

      for(int m=1; m<size; m++) {
        X.on_gpu<N>(d_XZm, d_Zs[m]);
        coo( d_dD_Ym, d_Ys[m] );

        X.Zdscal<N>( A[m], d_XZm );
        X.dot<N>( &inner, d_XZm, d_dD_Ym );
        res -= real(inner);
      }

      CUDA_CHECK(cudaFree(d_XZm));
      CUDA_CHECK(cudaFree(d_dD_Ym));
    }
    {
      CuC *d_XYm, *d_dD_Zm;
      CUDA_CHECK(cudaMalloc(&d_XYm, N*CD));
      CUDA_CHECK(cudaMalloc(&d_dD_Zm, N*CD));

      for(int m=1; m<size; m++) {
        X.on_gpu<N>(d_XYm, d_Ys[m]);
        coo( d_dD_Zm, d_Zs[m] );

        X.Zdscal<N>( A[m], d_XYm );
        X.dot<N>( &inner, d_XYm, d_dD_Zm );
        res -= real(inner);
      }

      CUDA_CHECK(cudaFree(d_XYm));
      CUDA_CHECK(cudaFree(d_dD_Zm));
    }
    {
      CuC *d_sum, *d_dD_sum;
      CUDA_CHECK(cudaMalloc(&d_sum, N*CD));
      CUDA_CHECK(cudaMalloc(&d_dD_sum, N*CD));

      CUDA_CHECK(cudaMemcpy(d_sum, d_eta, N*CD, D2D));
      for(int m=1; m<size; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_sum, A[m], d_Zs[m], d_sum);

      coo( d_dD_sum, d_sum );

      X.dot<N>( &inner, d_eta, d_dD_sum );
      res += real(inner);

      CUDA_CHECK(cudaFree(d_sum));
      CUDA_CHECK(cudaFree(d_dD_sum));
    }

    CUDA_CHECK(cudaFree(d_Xdag_eta));
    for(int m=1; m<size; m++) {
      CUDA_CHECK(cudaFree(d_Zs[m]));
      CUDA_CHECK(cudaFree(d_Ys[m]));
    }

    res *= -2.0*C/lambda_max;

    return res;
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
















  // void adj(std::vector<Complex>& res, const std::vector<Complex>& xi) const {
  //   std::vector<Complex> DHxi(xi.size());

  //   {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max)), {&M_DWH} );
  //     Op.from_cpu<N>( DHxi, xi );
  //     for(Idx i=0; i<res.size(); i++) DHxi[i] *= C;
  //   }

  //   res = DHxi;

  //   std::vector<Complex> tmp(xi.size());
  //   for(int m=1; m<size; m++) {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
  //     const CuC a = cplx(-k*k/cp[m]);
  //     Op.push_back ( a, {} );
  //     Op.solve<N>( tmp, DHxi );
  //     for(Idx i=0; i<res.size(); i++) res[i] += A[m] * tmp[i];
  //   }

  //   for(Idx i=0; i<res.size(); i++) res[i] += xi[i];
  // }


  // void sq(std::vector<Complex>& res, const std::vector<Complex>& xi) const {
  //   std::vector<Complex> tmp1(xi.size()), tmp2(xi.size());
  //   this->mult(tmp1, xi);
  //   this->adj(tmp2, xi);
  //   for(Idx i=0;  i<res.size(); i++) res[i] = tmp1[i] + tmp2[i];
  // }



  // double grad( const Link& link, const Gauge& U, const std::vector<Complex>& eta ) const {
  //   double res = 0.0;

  //   std::vector<Complex> Xdag_eta(eta.size());
  //   {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max)), {&M_DWH} );
  //     Op.from_cpu<N>(Xdag_eta, eta);
  //   }

  //   std::vector<std::vector<Complex>> Zs(size, std::vector<Complex>(eta.size()) );
  //   std::vector<std::vector<Complex>> Ys(size, std::vector<Complex>(eta.size()) );

  //   for(int m=1; m<size; m++) {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );

  //     const CuC a = cplx(-k*k/cp[m]);
  //     Op.push_back ( a, {} );
  //     Op.solve<N>( Zs[m], eta );
  //     Op.solve<N>( Ys[m], Xdag_eta );
  //   }

  //   COO coo;
  //   DW.d_coo_format(coo.en, U, link);
  //   coo.do_it();

  //   {
  //     std::vector<Complex> sum(eta.size(), 0.0);
  //     for(int m=1; m<size; m++) {
  //       for(Idx i=0; i<eta.size(); i++) sum[i] += Zs[m][i];
  //     }
  //     std::vector<Complex> dD_sum(eta.size(), 0.0);
  //     matmulcoo<N>( reinterpret_cast<CuC*>(dD_sum.data()),
  //                   reinterpret_cast<const CuC*>(sum.data()),
  //                    coo.en );
  //     for(Idx i=0; i<eta.size(); i++) res += std::real( std::conj(eta[i]) * dD_sum[i] );
  //   }
  //   {
  //     for(int m=1; m<size; m++) {
  //       std::vector<Complex> XZm(eta.size(), 0.0);
  //       {
  //         MatPoly Op;
  //         Op.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
  //         Op.from_cpu<N>(XZm, Zs[m]);
  //       }
  //       std::vector<Complex> dD_Ym(eta.size(), 0.0);
  //       matmulcoo<N>( reinterpret_cast<CuC*>(dD_Ym.data()),
  //                     reinterpret_cast<const CuC*>(Ys[m].data()),
  //                     coo.en );
  //       for(Idx i=0; i<eta.size(); i++) res -= std::real( std::conj(XZm[i]) * dD_Ym[i] );
  //     }
  //   }
  //   {
  //     for(int m=1; m<size; m++) {
  //       std::vector<Complex> XYm(eta.size(), 0.0);
  //       {
  //         MatPoly Op;
  //         Op.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
  //         Op.from_cpu<N>(XYm, Ys[m]);
  //       }
  //       std::vector<Complex> dD_Zm(eta.size(), 0.0);
  //       matmulcoo<N>( reinterpret_cast<CuC*>(dD_Zm.data()),
  //                     reinterpret_cast<const CuC*>(Zs[m].data()),
  //                     coo.en );
  //       for(Idx i=0; i<eta.size(); i++) res -= std::real( std::conj(XYm[i]) * dD_Zm[i] );
  //     }
  //   }

  //   res *= -2.0 * E * 2.0 / (1.0+lambda_inv) / (k*M);

  //   return res;
  // }


  // void mult(std::vector<Complex>& res, const std::vector<Complex>& xi) const {
  //   res = xi;
  //   std::vector<Complex> tmp(xi.size());

  //   for(int m=1; m<size; m++) {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
  //     const CuC a = cplx(-k*k/cp[m]);
  //     Op.push_back ( a, {} );
  //     Op.solve<N>( tmp, xi );
  //     for(Idx i=0; i<res.size(); i++) res[i] += A[m] * tmp[i];
  //   }

  //   for(Idx i=0; i<res.size(); i++) tmp[i] = E * 2.0 / (1.0+lambda_inv) / (k*M) * res[i];
  //   MatPoly Op;
  //   Op.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
  //   Op.from_cpu<N>( res, tmp );

  //   for(Idx i=0; i<res.size(); i++) res[i] += xi[i];
  // }

  // MatPoly Dummy;
  // double mu;
  // Dummy.dot2self<N>(&mu, d_xi);
  // std::cout << "|d_xi| = " << mu << std::endl;



  // double grad( const Link& link, const Gauge& U, CuC* d_eta ) const {
  //   double res = 0.0;

  //   double mu;
  //   CuC inner;
  //   MatPoly Dummy;

  //   CuC* d_Xdag_eta;
  //   CUDA_CHECK(cudaMalloc(&d_Xdag_eta, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_Xdag_eta, d_eta, N*CD, D2D)); // @@@@@@
  //   // {
  //   //   MatPoly Op;
  //   //   Op.push_back ( cplx(1.0/(lambda_max)), {&M_DWH} );
  //   //   Op.on_gpu<N>(d_Xdag_eta, d_eta);
  //   // }

  //   std::vector<CuC*> d_Zs(size);
  //   std::vector<CuC*> d_Ys(size);
  //   for(int m=1; m<size; m++) {
  //     CUDA_CHECK(cudaMalloc(&d_Zs[m], N*CD));
  //     CUDA_CHECK(cudaMalloc(&d_Ys[m], N*CD));
  //   }

  //   for(int m=1; m<size; m++) {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );

  //     const CuC a = cplx(-k*k/cp[m]);
  //     Op.push_back ( a, {} );
  //     Op.solve<N>( d_Zs[m], d_eta );
  //     Op.solve<N>( d_Ys[m], d_Xdag_eta );

  //     Op.Zdscal<N>( A[m], d_Zs[m] );
  //     Op.Zdscal<N>( A[m], d_Ys[m] );
  //   }

  //   COO coo;
  //   DW.d_coo_format(coo.en, U, link);
  //   coo.do_it();

  //   // {
  //   //   CuC *d_sum, *d_dD_sum;
  //   //   CUDA_CHECK(cudaMalloc(&d_sum, N*CD));
  //   //   CUDA_CHECK(cudaMalloc(&d_dD_sum, N*CD));
  //   //   CUDA_CHECK(cudaMemcpy(d_sum, d_eta, N*CD, D2D));
  //   //   for(int m=1; m<size; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_sum, 1.0, d_Zs[m], d_sum);

  //   //   coo( d_dD_sum, d_sum );

  //   //   Dummy.dot<N>( &inner, d_eta, d_dD_sum );
  //   //   res += real(inner);

  //   //   CUDA_CHECK(cudaFree(d_sum));
  //   //   CUDA_CHECK(cudaFree(d_dD_sum));
  //   // }
  //   // {
  //   //   CuC *d_XZm, *d_dD_Ym;
  //   //   CUDA_CHECK(cudaMalloc(&d_XZm, N*CD));
  //   //   // CUDA_CHECK(cudaMalloc(&d_dD_Ym, N*CD));

  //   //   for(int m=1; m<size; m++) {
  //   //     {
  //   //       MatPoly Op;
  //   //       Op.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
  //   //       Op.on_gpu<N>(d_XZm, d_Zs[m]);
  //   //     }
  //   //     coo( d_dD_Ym, d_Ys[m] );

  //   //     Dummy.dot<N>( &inner, d_XZm, d_dD_Ym );
  //   //     res -= real(inner);
  //   //   }

  //   //   CUDA_CHECK(cudaFree(d_XZm));
  //   //   CUDA_CHECK(cudaFree(d_dD_Ym));
  //   // }
  //   {
  //     CuC *d_XYm, *d_dD_Zm;
  //     CUDA_CHECK(cudaMalloc(&d_XYm, N*CD));
  //     CUDA_CHECK(cudaMalloc(&d_dD_Zm, N*CD));

  //     for(int m=1; m<size; m++) {
  //       {
  //         MatPoly Op;
  //         Op.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
  //         Op.on_gpu<N>(d_XYm, d_Ys[m]);
  //       }
  //       coo( d_dD_Zm, d_Zs[m] );

  //       Dummy.dot<N>( &inner, d_XYm, d_dD_Zm );
  //       res -= real(inner);
  //     }
  //     CUDA_CHECK(cudaFree(d_XYm));
  //     CUDA_CHECK(cudaFree(d_dD_Zm));
  //   }

  //   res *= 2.0*C/lambda_max;

  //   CUDA_CHECK(cudaFree(d_Xdag_eta));
  //   for(int m=1; m<size; m++) {
  //     CUDA_CHECK(cudaFree(d_Zs[m]));
  //     CUDA_CHECK(cudaFree(d_Ys[m]));
  //   }

  //   return res;
  // }


  // void mult_device2(CuC* d_res, const CuC* d_xi) const {
  //   const int mmax = size; // size;
  //   std::vector<CuC*> d_Zs(mmax);
  //   for(int m=1; m<mmax; m++) CUDA_CHECK(cudaMalloc(&d_Zs[m], N*CD));

  //   for(int m=1; m<mmax; m++) {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );

  //     const CuC a = cplx(-k*k/cp[m]);
  //     Op.push_back ( a, {} );
  //     Op.solve<N>( d_Zs[m], d_xi );
  //     Op.Zdscal<N>( A[m], d_Zs[m] );
  //   }

  //   // CUDA_CHECK(cudaMemcpy(d_res, d_xi, N*CD, D2D));
  //   CUDA_CHECK(cudaMemset(d_res, 0, N*CD));
  //   for(int m=1; m<mmax; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_Zs[m], d_res);

  //   for(int m=1; m<mmax; m++) CUDA_CHECK(cudaFree(d_Zs[m]));
  // }



  // double grad2( const Link& link, const Gauge& U, CuC* d_xi ) const {
  //   const int mmax = size; // size
  //   std::vector<CuC*> d_Zs(mmax);
  //   for(int m=1; m<mmax; m++) CUDA_CHECK(cudaMalloc(&d_Zs[m], N*CD));

  //   for(int m=1; m<mmax; m++) {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );

  //     const CuC a = cplx(-k*k/cp[m]);
  //     Op.push_back ( a, {} );
  //     Op.solve<N>( d_Zs[m], d_xi );
  //   }

  //   COO coo;
  //   DW.d_coo_format(coo.en, U, link);
  //   coo.do_it();

  //   CuC *d_tmp1, *d_tmp2;
  //   CUDA_CHECK(cudaMalloc(&d_tmp1, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_tmp2, N*CD));

  //   double res = 0.0;
  //   MatPoly Dummy;

  //   for(int m=1; m<mmax; m++){
  //     coo( d_tmp1, d_Zs[m] ); // DH
  //     M_DWH( d_tmp2, d_tmp1 );

  //     CuC inner;
  //     Dummy.Zdscal<N>( A[m], d_Zs[m] );
  //     Dummy.dot<N>( &inner, d_Zs[m], d_tmp2 );
  //     res += -2.0 * real(inner) / (lambda_max * lambda_max);
  //   }

  //   CUDA_CHECK(cudaFree(d_tmp1));
  //   CUDA_CHECK(cudaFree(d_tmp2));
  //   for(int m=1; m<mmax; m++) CUDA_CHECK(cudaFree(d_Zs[m]));

  //   return res;
  // }

