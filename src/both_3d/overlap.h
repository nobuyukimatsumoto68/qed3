#pragma once

#include <cmath>
#include <vector>


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


template<typename WilsonDirac>
struct Overlap : public Zolotarev {
  static constexpr Idx N = Comp::N;
  static constexpr int nstreams = Comp::NSTREAMS;

  const WilsonDirac& DW;
  DWDevice<WilsonDirac> d_DW; // actual data used in M_DW, M_DWH
  CSR M_DW;
  CSR M_DWH;
  double lambda_max, lambda_min;

  std::vector<cudaStream_t> stream;
  std::vector<cublasHandle_t> handle;

  std::vector<CuC*> d_Ys, d_Zs, d_XYs, d_XZs;

  bool is_precalc;

  explicit Overlap( const WilsonDirac& DW_,
                    const int n_=21,
                    const double k_=0.01,
                    const bool locate_on_gpu=true)
    : Zolotarev(k_, n_)
    , DW(DW_)
    , d_DW(DW)
    , stream(nstreams)
    , handle(nstreams)
    , d_Ys(size)
    , d_Zs(size)
    , d_XYs(size)
    , d_XZs(size)
    , is_precalc(false)
  {
    d_DW.associateCSR( M_DW, false );
    d_DW.associateCSR( M_DWH, true );

    for(int istream=0; istream<nstreams; istream++) {
      // CUDA_CHECK(cudaStreamCreateWithFlags(&stream[istream], cudaStreamNonBlocking));
      CUDA_CHECK(cudaStreamCreate( &stream[istream] ));
      CUBLAS_CHECK(cublasCreate(&handle[istream]));
      CUBLAS_CHECK(cublasSetStream(handle[istream], stream[istream]));
    }

    for(int m=0; m<size; m++) {
      const int istream = m%nstreams;
      CUDA_CHECK(cudaMallocAsync(&d_Zs[m], N*CD, stream[istream]));
      CUDA_CHECK(cudaMallocAsync(&d_Ys[m], N*CD, stream[istream]));
      CUDA_CHECK(cudaMallocAsync(&d_XZs[m], N*CD, stream[istream]));
      CUDA_CHECK(cudaMallocAsync(&d_XYs[m], N*CD, stream[istream]));
    }
    CUDA_CHECK(cudaDeviceSynchronize());
  }

  ~Overlap()
  {
    for(int m=0; m<size; m++) {
      const int istream = m%nstreams;
      CUDA_CHECK(cudaFreeAsync(d_Zs[m], stream[istream]));
      CUDA_CHECK(cudaFreeAsync(d_Ys[m], stream[istream]));
      CUDA_CHECK(cudaFreeAsync(d_XZs[m], stream[istream]));
      CUDA_CHECK(cudaFreeAsync(d_XYs[m], stream[istream]));
    }
    for(int istream=0; istream<nstreams; istream++) {
      CUDA_CHECK(cudaStreamSynchronize(stream[istream]));
      CUDA_CHECK(cudaStreamDestroy(stream[istream]));
      CUBLAS_CHECK(cublasDestroy(handle[istream]));
    }
    CUDA_CHECK(cudaDeviceSynchronize());
  }

  template<typename Gauge>
  void update( const Gauge& U ) {
    d_DW.update( U );
    compute_lambda_max();
    Zolotarev::update(lambda_min/lambda_max);
#ifdef InfoDelta
    std::clog << "# Delta : " << Delta() << std::endl;
#endif
    is_precalc = false;
  }

  void compute_lambda_max( const double TOL=1.0e-8, const int MAXITER=500 ) {
    std::vector<Complex> q(N, 0.0);

    q[0] = std::sqrt(1.0/2.0);
    q[1] = std::sqrt(1.0/2.0);

    MatPoly Op( handle[0], stream[0] );
    Op.push_back ( cplx(1.0), {&M_DW, &M_DWH} );

    CuC *d_x, *d_q;
    CUDA_CHECK(cudaMallocAsync(&d_x, N*CD, stream[0]));
    CUDA_CHECK(cudaMallocAsync(&d_q, N*CD, stream[nstreams-1]));
    CUDA_CHECK(cudaMemcpyAsync(d_q, reinterpret_cast<const CuC*>(q.data()), N*CD, H2D, stream[nstreams-1]));
    CUDA_CHECK(cudaStreamSynchronize(stream[nstreams-1]));

    Complex dot;
    double norm=1.0, mu_0=1.0, mu_m1=1.0, mu_m2=1.0;
    double lambda=100.0, lambda_old=1000.0;

    for(int i=0; i<MAXITER; i++){
      Op.on_gpuAsync<N>( d_x, d_q );
      //
      CUDA_CHECK(cudaMemcpyAsync(d_q, d_x, N*CD, D2D, stream[nstreams-1])); // stream 2
      Op.dot2selfAsync<N>(&norm, d_x); // stream 1
      CUDA_CHECK(cudaStreamSynchronize(stream[nstreams-1]));
      //
      Op.ZdscalAsync<N>( 1.0/std::sqrt(norm), d_q );

      Op.dotAsync<N>(reinterpret_cast<CuC*>(&dot), d_x, d_q);
      CUDA_CHECK(cudaStreamSynchronize(stream[0]));
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
      Op.solveAsync<N>( d_x, d_q, Comp::TOL_OUTER );
      //
      CUDA_CHECK(cudaMemcpyAsync(d_q, d_x, N*CD, D2D, stream[nstreams-1])); // stream 2
      Op.dot2selfAsync<N>(&norm, d_x);
      CUDA_CHECK(cudaStreamSynchronize(stream[nstreams-1]));
      //
      Op.ZdscalAsync<N>( 1.0/std::sqrt(norm), d_q );

      Op.dotAsync<N>(reinterpret_cast<CuC*>(&dot), d_x, d_q);
      CUDA_CHECK(cudaStreamSynchronize(stream[0]));
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

    CUDA_CHECK(cudaFreeAsync(d_x, stream[0]));
    CUDA_CHECK(cudaFreeAsync(d_q, stream[nstreams-1]));

    // lambda_min = 0.5*std::sqrt( 1.0/lambda2 );
    // lambda_max = 2.0*std::sqrt( lambda );
    lambda_min = std::sqrt( (1.0-100*TOL)/lambda2 );
    lambda_max = std::sqrt( (1.0+100*TOL)*lambda );
    // lambda_min = 0.01; // std::sqrt( (1.0-100*TOL)/lambda2 );
    // lambda_max = 16; // std::sqrt( (1.0+100*TOL)*lambda );

    CUDA_CHECK(cudaDeviceSynchronize());
  }


  void mult_deviceAsyncLaunch(CuC* d_res, const CuC* d_xi) const {
    // can parallelize
#ifdef _OPENMP
#pragma omp parallel for num_threads(nstreams) schedule(static)
#endif
    for(int m=1; m<size; m++) {
      const int istream = m%nstreams;
      MatPoly Op(handle[istream], stream[istream]);
      Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
      const CuC a = cplx(-k*k/cp[m]);
      Op.push_back ( a, {} );
      Op.solveAsync<N>( d_Zs[m], d_xi, Comp::TOL_INNER );
    }

    // reduction
    CUDA_CHECK(cudaMemcpy(d_Zs[0], d_xi, N*CD, D2D)); // E(1+Z)
    for(int m=1; m<size; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_Zs[0], A[m], d_Zs[m], d_Zs[0]);

    MatPoly Op( handle[0], stream[0] );
    Op.push_back( cplx(1.0/(lambda_max)), {&M_DW} );
    Op.on_gpuAsync<N>( d_res, d_Zs[0] );

    Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, C, d_res, d_xi); // 1+V
  }

  void adj_deviceAsyncLaunch(CuC* d_res, const CuC* d_xi) const {
    CUDA_CHECK(cudaMemcpy(d_Ys[0], d_xi, N*CD, D2D)); // E(1+Y)
    {
      MatPoly OpGlob( handle[0], stream[0] );
      OpGlob.push_back ( cplx(1.0/(lambda_max)), {&M_DWH} );
      OpGlob.on_gpuAsync<N>( d_Ys[0], d_xi );
    }

    CUDA_CHECK(cudaMemcpy(d_res, d_Ys[0], N*CD, D2D));

    // can parallelize
#ifdef _OPENMP
#pragma omp parallel for num_threads(nstreams) schedule(static)
#endif
    for(int m=1; m<size; m++) {
      const int istream = m%nstreams;
      MatPoly Op(handle[istream], stream[istream]);
      Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
      const CuC a = cplx(-k*k/cp[m]);
      Op.push_back ( a, {} );
      Op.solveAsync<N>( d_Ys[m], d_Ys[0], Comp::TOL_INNER );
    }

    // reduction
    for(int m=1; m<size; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, A[m], d_Ys[m], d_res);
    Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, C, d_res, d_xi); // A[0]=1.0
  }


  void sq_deviceAsyncLaunch( CuC* d_res, const CuC* d_xi ) const {
    CuC *d_tmp1, *d_tmp2;
    CUDA_CHECK(cudaMalloc(&d_tmp1, N*CD));
    CUDA_CHECK(cudaMalloc(&d_tmp2, N*CD));

    this->mult_deviceAsyncLaunch(d_tmp1, d_xi);
    this->adj_deviceAsyncLaunch(d_tmp2, d_xi);

    CUDA_CHECK(cudaMemcpy(d_res, d_tmp1, N*CD, D2D));
    Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_tmp2, d_res); // A[0]=1.0
    CUDA_CHECK(cudaFree(d_tmp1));
    CUDA_CHECK(cudaFree(d_tmp2));
  }


  template<typename Gauge>
  void precalc_grad_deviceAsyncLaunch( const Gauge& U, const CuC* d_eta ) {
    {
      MatPoly XH(handle[0], stream[0]);   XH.push_back ( cplx(1.0/(lambda_max)), {&M_DWH} );
      XH.on_gpuAsync<N>(d_Ys[0], d_eta);
    }

    // can parallelize omp parallel nparallel=nstreams static assert(nparallel>=nstreams)
#ifdef _OPENMP
#pragma omp parallel for num_threads(nstreams) schedule(static)
#endif
    for(int m=1; m<size; m++) {
      const int istream = m%nstreams;
      MatPoly Op(handle[istream], stream[istream]);
      Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
      const CuC a = cplx(-k*k/cp[m]); Op.push_back ( a, {} );

      Op.solveAsync<N>( d_Zs[m], d_eta, Comp::TOL_INNER );
      Op.solveAsync<N>( d_Ys[m], d_Ys[0], Comp::TOL_INNER );
    }

    is_precalc = true;
    CUDA_CHECK(cudaDeviceSynchronize());
  }

  template<typename Link, typename Gauge>
  double grad_deviceAsyncLaunch( const Link& link, const Gauge& U, const CuC* d_eta ) const {
    assert( is_precalc );

    COO coo;
    DW.d_coo_format(coo.en, U, link);
    coo.do_it();

    std::vector<double> tmp2reduce(size, 0.0);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nstreams) schedule(static)
#endif
    for(int m=1; m<size; m++) {
      const int istream = m%nstreams;
      MatPoly X(handle[istream], stream[istream]);
      X.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );

      X.on_gpuAsync<N>(d_XZs[m], d_Zs[m]);
      coo.Async( d_XYs[m], d_Ys[m], stream[istream] );

      CuC inner;
      X.dotAsync<N>( &inner, d_XZs[m], d_XYs[m] );
      tmp2reduce[m] -= real(inner);

      X.on_gpuAsync<N>(d_XYs[m], d_Ys[m]);
      coo.Async( d_XZs[m], d_Zs[m], stream[istream] );

      X.dotAsync<N>( &inner, d_XYs[m], d_XZs[m] );
      tmp2reduce[m] -= real(inner);

      tmp2reduce[m] *= A[m];
    }

    double res = 0.0;
    // reductions
    for(int m=1; m<size; m++) res+=tmp2reduce[m];

    CUDA_CHECK(cudaMemcpy(d_Zs[0], d_eta, N*CD, D2D));
    for(int m=1; m<size; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_Zs[0], A[m], d_Zs[m], d_Zs[0]);
    coo( d_Ys[0], d_Zs[0] );
    CuC inner;
    {
      MatPoly XH(handle[0], stream[0]);
      XH.dotAsync<N>( &inner, d_eta, d_Ys[0] );
    }
    res += real(inner);
    res *= -2.0*C/lambda_max;
    return res;
  }



    // void mult(std::vector<Complex>& x, const std::vector<Complex>& b) const {
  //   CuC *d_x, *d_r; // , *d_tmp, *d_tmp2;
  //   CUDA_CHECK(cudaMalloc(&d_x, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_r, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_r, reinterpret_cast<const CuC*>(b.data()), N*CD, H2D));

  //   mult_device(d_x, d_r);

  //   CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(x.data()), d_x, N*CD, D2H));
  //   CUDA_CHECK(cudaFree(d_x));
  //   CUDA_CHECK(cudaFree(d_r));
  // }


  // void mult_device(CuC* d_res, const CuC* d_xi) const {
  //   std::vector<CuC*> d_Xs(size);
  //   for(int m=0; m<size; m++) CUDA_CHECK(cudaMalloc(&d_Xs[m], N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_Xs[0], d_xi, N*CD, D2D)); // E(1+Z)

  //   // can parallelize
  //   for(int m=1; m<size; m++) {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
  //     const CuC a = cplx(-k*k/cp[m]);
  //     Op.push_back ( a, {} );
  //     Op.solve<N>( d_Xs[m], d_xi );
  //   }
  //   for(int m=1; m<size; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_Xs[0], A[m], d_Xs[m], d_Xs[0]);

  //   MatPoly Op;
  //   // Op.Zdscal<N>(C, d_Xs[0]);

  //   Op.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
  //   Op.on_gpu<N>( d_res, d_Xs[0] );
  //   // Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_xi, d_res); // 1+V
  //   Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, C, d_res, d_xi); // 1+V

  //   for(int m=0; m<size; m++) CUDA_CHECK(cudaFree(d_Xs[m]));
  // }


  // void adj_device(CuC* d_res, const CuC* d_xi) const {
  //   std::vector<CuC*> d_Ys(size);
  //   for(int m=0; m<size; m++) CUDA_CHECK(cudaMalloc(&d_Ys[m], N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_Ys[0], d_xi, N*CD, D2D)); // E(1+Y)

  //   MatPoly OpGlob;
  //   OpGlob.push_back ( cplx(1.0/(lambda_max)), {&M_DWH} );
  //   OpGlob.on_gpu<N>( d_Ys[0], d_xi );
  //   CUDA_CHECK(cudaMemcpy(d_res, d_Ys[0], N*CD, D2D));

  //   // can parallelize
  //   for(int m=1; m<size; m++) {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
  //     const CuC a = cplx(-k*k/cp[m]);
  //     Op.push_back ( a, {} );
  //     Op.solve<N>( d_Ys[m], d_Ys[0] );
  //   }

  //   for(int m=1; m<size; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, A[m], d_Ys[m], d_res);
  //   OpGlob.Zdscal<N>( C, d_res );
  //   Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_xi, d_res); // A[0]=1.0

  //   for(int m=0; m<size; m++) CUDA_CHECK(cudaFree(d_Ys[m]));
  // }


    // void sq_device( CuC* d_res, const CuC* d_xi) const {
  //   CuC *d_tmp1, *d_tmp2;
  //   CUDA_CHECK(cudaMalloc(&d_tmp1, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_tmp2, N*CD));

  //   this->mult_device(d_tmp1, d_xi);
  //   this->adj_device(d_tmp2, d_xi);

  //   CUDA_CHECK(cudaMemcpy(d_res, d_tmp1, N*CD, D2D));
  //   Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_tmp2, d_res); // A[0]=1.0

  //   CUDA_CHECK(cudaFree(d_tmp1));
  //   CUDA_CHECK(cudaFree(d_tmp2));
  // }



    // double grad_device( const Link& link, const Gauge& U, CuC* d_eta ) const {
  //   std::vector<CuC*> d_Ys(size), d_Zs(size);
  //   for(int m=0; m<size; m++) {
  //     CUDA_CHECK(cudaMalloc(&d_Zs[m], N*CD));
  //     CUDA_CHECK(cudaMalloc(&d_Ys[m], N*CD));
  //   }
  //   {
  //     MatPoly XH;   XH.push_back ( cplx(1.0/(lambda_max)), {&M_DWH} );
  //     XH.on_gpu<N>(d_Ys[0], d_eta);
  //   }

  //   // can parallelize omp parallel nparallel=nstreams static assert(nparallel>=nstreams)
  //   for(int m=1; m<size; m++) {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
  //     const CuC a = cplx(-k*k/cp[m]); Op.push_back ( a, {} );

  //     Op.solve<N>( d_Zs[m], d_eta );
  //     Op.solve<N>( d_Ys[m], d_Ys[0] );
  //   }

  //   COO coo;
  //   DW.d_coo_format(coo.en, U, link);
  //   coo.do_it();

  //   CuC inner;
  //   double res = 0.0;

  //   MatPoly X; X.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
  //   for(int m=1; m<size; m++) {
  //     X.on_gpu<N>(d_Zs[0], d_Zs[m]);
  //     coo( d_Ys[0], d_Ys[m] );

  //     X.Zdscal<N>( A[m], d_Zs[0] );
  //     X.dot<N>( &inner, d_Zs[0], d_Ys[0] );
  //     res -= real(inner);
  //   }
  //   for(int m=1; m<size; m++) {
  //     X.on_gpu<N>(d_Ys[0], d_Ys[m]);
  //     coo( d_Zs[0], d_Zs[m] );

  //     X.Zdscal<N>( A[m], d_Ys[0] );
  //     X.dot<N>( &inner, d_Ys[0], d_Zs[0] );
  //     res -= real(inner);
  //   }

  //   CUDA_CHECK(cudaMemcpy(d_Zs[0], d_eta, N*CD, D2D));
  //   for(int m=1; m<size; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_Zs[0], A[m], d_Zs[m], d_Zs[0]);
  //   coo( d_Ys[0], d_Zs[0] );
  //   X.dot<N>( &inner, d_eta, d_Ys[0] );
  //   res += real(inner);

  //   for(int m=0; m<size; m++) {
  //     CUDA_CHECK(cudaFree(d_Zs[m]));
  //     CUDA_CHECK(cudaFree(d_Ys[m]));
  //   }
  //   res *= -2.0*C/lambda_max;
  //   return res;
  // }


};


