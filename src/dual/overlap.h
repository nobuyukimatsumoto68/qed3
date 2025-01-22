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
  double k;
  double kp;
  const int size;

  std::vector<double> c;
  std::vector<double> cp;
  double M;
  double lambda_inv;

  double E;
  double C;
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
    C = E * 2.0 / (1.0+lambda_inv) / (k*M);
  }

  ~Zolotarev(){
    // CUDA_CHECK(cudaFree(d_c));
    // CUDA_CHECK(cudaFree(d_A));
    // CUDA_CHECK(cudaFree(d_C));
  }

  void update(const double k_)
  {
    k = k_;
    kp = std::sqrt(1.0-k*k);
    get_coeffs();
    partial_fraction();
    C = E * 2.0 / (1.0+lambda_inv) / (k*M);
  }


  // A.D.Kennedy 2004; https://arxiv.org/abs/hep-lat/0402038
  static void sncndnK(const FLOAT z, const FLOAT k,
		      FLOAT& sn, FLOAT& cn, FLOAT& dn,
		      FLOAT& K) {
    FLOAT agm;
    int sgn;
    sn = arithgeom(z, ONE, std::sqrt(ONE - k*k), &agm);
    K = M_PI / (TWO * agm);
    sgn = ((int) (std::abs(z) / K)) % 4; /* sgn = 0, 1, 2, 3 */
    sgn ^= sgn >> 1; /* (sgn & 1) = 0, 1, 1, 0 */
    sgn = 1 - ((sgn & 1) << 1); /* sgn = 1, -1, -1, 1 */
    cn = ((FLOAT) sgn) * std::sqrt(ONE - sn * sn);
    dn = std::sqrt(ONE - k*k* sn * sn);
  }

  static FLOAT arithgeom(const FLOAT z, FLOAT a, FLOAT b, FLOAT* agm) {
    static FLOAT pb = -ONE;
    FLOAT xi;
    if (b <= pb) {
      pb = -ONE;
      *agm = a;
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
  double lambda_max, lambda_min;

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
    Zolotarev::update(lambda_min/lambda_max);
  }

  void compute_lambda_max( const double TOL=1.0e-8, const int MAXITER=500 ) {
    std::vector<Complex> q(N, 0.0);

    q[0] = std::sqrt(1.0/2.0);
    q[1] = std::sqrt(1.0/2.0);

    MatPoly Op;
    Op.push_back ( cplx(1.0), {&M_DW, &M_DWH} );

    CuC *d_x, *d_q;
    CUDA_CHECK(cudaMalloc(&d_x, N*CD));
    CUDA_CHECK(cudaMalloc(&d_q, N*CD));
    CUDA_CHECK(cudaMemcpy(d_q, reinterpret_cast<const CuC*>(q.data()), N*CD, H2D));

    Complex dot;
    double norm=1.0, mu_0=1.0, mu_m1=1.0, mu_m2=1.0;

    double lambda=100.0, lambda_old=1000.0;

    for(int i=0; i<MAXITER; i++){
      Op.on_gpu<N>( d_x, d_q );
      //
      Op.dot2self<N>(&norm, d_x);
      CUDA_CHECK(cudaMemcpy(d_q, d_x, N*CD, D2D));
      Op.Zdscal<N>( 1.0/std::sqrt(norm), d_q );

      Op.dot<N>(reinterpret_cast<CuC*>(&dot), d_x, d_q);
      mu_m2=mu_m1;
      mu_m1=mu_0;
      mu_0=dot.real();

      const double r = (mu_0-mu_m1)/(mu_m1-mu_m2);
      const double a = (mu_0-mu_m1)/std::pow(r,i-1)/(r-1);
      lambda_old = lambda;
      lambda = mu_0 - a*std::pow(r,i);

      if(std::abs(lambda_old-lambda)/std::abs(lambda)<TOL) {
#ifdef IsVerbose
	std::clog << "# lambda_max estimate escaped in i = " << i << std::endl;
#endif
	break;
      }
    }

    CUDA_CHECK(cudaMemcpy(d_q, reinterpret_cast<const CuC*>(q.data()), N*CD, H2D));
    double lambda2=100.0, lambda2_old=1000.0;

    for(int i=0; i<MAXITER; i++){
      Op.solve<N>( d_x, d_q );
      //
      Op.dot2self<N>(&norm, d_x);
      CUDA_CHECK(cudaMemcpy(d_q, d_x, N*CD, D2D));
      Op.Zdscal<N>( 1.0/std::sqrt(norm), d_q );

      Op.dot<N>(reinterpret_cast<CuC*>(&dot), d_x, d_q);
      mu_m2=mu_m1;
      mu_m1=mu_0;
      mu_0=dot.real();

      const double r = (mu_0-mu_m1)/(mu_m1-mu_m2);
      const double a = (mu_0-mu_m1)/std::pow(r,i-1)/(r-1);
      lambda2_old = lambda2;
      lambda2 = mu_0 - a*std::pow(r,i);

      if(std::abs(lambda2_old-lambda2)/std::abs(lambda2)<TOL) {
#ifdef IsVerbose
	std::clog << "# lambda_min estimate escaped in i = " << i << std::endl;
#endif
	break;
      }
    }

    CUDA_CHECK(cudaFree(d_x));
    CUDA_CHECK(cudaFree(d_q));

    lambda_min = std::sqrt( (1.0-100*TOL)/lambda2 );
    lambda_max = std::sqrt( (1.0+100*TOL)*lambda );
  }

  void mult(std::vector<Complex>& x, const std::vector<Complex>& b) const {
    CuC *d_x, *d_r; // , *d_tmp, *d_tmp2;
    CUDA_CHECK(cudaMalloc(&d_x, N*CD));
    CUDA_CHECK(cudaMalloc(&d_r, N*CD));
    CUDA_CHECK(cudaMemcpy(d_r, reinterpret_cast<const CuC*>(b.data()), N*CD, H2D));

    mult_device(d_x, d_r);

    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(x.data()), d_x, N*CD, D2H));
    CUDA_CHECK(cudaFree(d_x));
    CUDA_CHECK(cudaFree(d_r));
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


  double grad_device( const Link& link, const Gauge& U, CuC* d_eta ) const {
    CuC inner;
    double res = 0.0;

    std::vector<CuC*> d_Ys(size), d_Zs(size);
    for(int m=1; m<size; m++) {
      CUDA_CHECK(cudaMalloc(&d_Zs[m], N*CD));
      CUDA_CHECK(cudaMalloc(&d_Ys[m], N*CD));
    }

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

    MatPoly X; X.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
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







// potentially factor 2
  // void sq_device2( CuC* d_res, const CuC* d_xi) const {
  //   CuC *d_res2;
  //   CUDA_CHECK(cudaMalloc(&d_res2, N*CD));

  //   cublasHandle_t handle, handle2;
  //   cudaStream_t stream, stream2;
  //   CUBLAS_CHECK(cublasCreate(&handle));
  //   CUBLAS_CHECK(cublasCreate(&handle2));
  //   CUDA_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
  //   CUDA_CHECK(cudaStreamCreateWithFlags(&stream2, cudaStreamNonBlocking));
  //   CUBLAS_CHECK(cublasSetStream(handle, stream));
  //   CUBLAS_CHECK(cublasSetStream(handle2, stream2));

  //   CuC *d_DHxi, *d_tmp, *d_tmp2;
  //   CUDA_CHECK(cudaMalloc(&d_DHxi, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_tmp, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_tmp2, N*CD));

  //   // this->mult_device(d_tmp1, d_xi);
  //   {
  //     CUDA_CHECK(cudaMemcpy(d_res, d_xi, N*CD, D2D));

  //     for(int m=1; m<size; m++) {
  //       // MatPoly Op(handle);
  //       MatPoly Op;
  //       Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
  //       const CuC a = cplx(-k*k/cp[m]);
  //       Op.push_back ( a, {} );
  //       Op.solve<N>( d_tmp, d_xi );
  //       Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, A[m], d_tmp, d_res);
  //     }

  //     CUDA_CHECK(cudaMemcpy(d_tmp, d_res, N*CD, D2D));
  //     MatPoly Op; Op.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
  //     Op.on_gpu<N>( d_res, d_tmp );

  //     Op.Zdscal<N>(C, d_res);
  //     Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_xi, d_res); // A[0]=1.0
  //   }
  //   // this->adj_device(d_tmp2, d_xi);
  //   {
  //     CUDA_CHECK(cudaMemcpy(d_DHxi, d_xi, N*CD, D2D));
  //     {
  //       MatPoly Op(handle2); Op.push_back ( cplx(1.0/(lambda_max)), {&M_DWH} );
  //       Op.on_gpu<N>( d_DHxi, d_xi );
  //       Op.Zdscal<N>( C, d_DHxi );
  //     }
  //     CUDA_CHECK(cudaMemcpy(d_res2, d_DHxi, N*CD, D2D));

  //     for(int m=1; m<size; m++) {
  //       MatPoly Op;
  //       Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
  //       const CuC a = cplx(-k*k/cp[m]);
  //       Op.push_back ( a, {} );
  //       Op.solve<N>( d_tmp2, d_DHxi );
  //       Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res2, A[m], d_tmp2, d_res2);
  //     }

  //     Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res2, 1.0, d_xi, d_res2); // A[0]=1.0
  //   }
  //   CUDA_CHECK(cudaFree(d_tmp));
  //   CUDA_CHECK(cudaFree(d_tmp2));
  //   CUDA_CHECK(cudaFree(d_DHxi));
  //   CUBLAS_CHECK(cublasDestroy(handle));
  //   CUBLAS_CHECK(cublasDestroy(handle2));

  //   // sync();
  //   CUDA_CHECK(cudaStreamSynchronize(stream));
  //   CUDA_CHECK(cudaStreamSynchronize(stream2));
  //   CUDA_CHECK(cudaStreamDestroy(stream));
  //   CUDA_CHECK(cudaStreamDestroy(stream2));

  //   // CUDA_CHECK(cudaMemcpy(d_res, d_res1, N*CD, D2D));
  //   Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_res2, d_res); // A[0]=1.0

  //   // CUDA_CHECK(cudaFree(d_tmp1));
  //   CUDA_CHECK(cudaFree(d_res2));
  // }
