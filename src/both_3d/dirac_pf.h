#pragma once

template<typename WilsonDirac>
struct DiracPf {

  static constexpr Idx N = Comp::N;
  static constexpr int nstreams = Comp::NSTREAMS;

  const WilsonDirac& DW;
  DWDevice<WilsonDirac> d_DW; // actual data used in M_DW, M_DWH
  CSR M_DW;
  CSR M_DWH;

  std::vector<cudaStream_t> stream;
  std::vector<cublasHandle_t> handle;

  std::vector<CuC*> d_dDeta;
  CuC* d_Deta;

  bool is_precalc;

  explicit DiracPf( const WilsonDirac& DW_,
                    const bool locate_on_gpu=true)
    : DW(DW_)
    , d_DW(DW)
    , stream(nstreams)
    , handle(nstreams)
    , d_dDeta(nstreams)
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

    CUDA_CHECK(cudaMallocAsync(&d_Deta, N*CD, stream[0]));
    for(int m=0; m<nstreams; m++) CUDA_CHECK(cudaMallocAsync(&d_dDeta[m], N*CD, stream[m]));

    CUDA_CHECK(cudaDeviceSynchronize());
  }

  ~DiracPf()
  {
    for(int m=0; m<nstreams; m++) CUDA_CHECK(cudaFreeAsync(d_dDeta[m], stream[m]));
    CUDA_CHECK(cudaFreeAsync(d_Deta, stream[0]));

    for(int istream=0; istream<nstreams; istream++) {
      CUDA_CHECK(cudaStreamSynchronize(stream[istream]));
      CUDA_CHECK(cudaStreamDestroy(stream[istream]));
      CUBLAS_CHECK(cublasDestroy(handle[istream]));
    }
    CUDA_CHECK(cudaDeviceSynchronize());
  }

  template<typename Gauge>
  void update( const Gauge& U ){
    d_DW.update( U );
    is_precalc = false;
  }

  void sq_deviceAsyncLaunch( CuC* d_res, const CuC* d_xi ) const {
    MatPoly X(handle[0], stream[0]);
    X.push_back ( cplx(1.0), {&M_DW, &M_DWH} );
    X.on_gpuAsync<N>(d_res, d_xi);
    CUDA_CHECK(cudaDeviceSynchronize());
  }

  void mult_deviceAsyncLaunch(CuC* d_res, const CuC* d_xi) const {
    MatPoly OpGlob( handle[0], stream[0] );
    OpGlob.push_back ( cplx(1.0), {&M_DW} );
    OpGlob.on_gpuAsync<N>( d_res, d_xi );
    CUDA_CHECK(cudaDeviceSynchronize());
  }

  void adj_deviceAsyncLaunch(CuC* d_res, const CuC* d_xi) const {
    MatPoly OpGlob( handle[0], stream[0] );
    OpGlob.push_back ( cplx(1.0), {&M_DWH} );
    OpGlob.on_gpuAsync<N>( d_res, d_xi );
    CUDA_CHECK(cudaDeviceSynchronize());
  }

  template<typename Gauge>
  void precalc_grad_deviceAsyncLaunch( const Gauge& U, const CuC* d_eta ) {
    is_precalc = true;
    MatPoly X(handle[0], stream[0]);

    X.push_back ( cplx(1.0), {&M_DW} );
    X.on_gpuAsync<N>(d_Deta, d_eta);

    CUDA_CHECK(cudaDeviceSynchronize());
  }


  template<typename Link, typename Gauge>
  double grad_deviceAsyncLaunch( const Link& link, const Gauge& U, const CuC* d_eta ) const {
    assert( is_precalc );
    const int m = omp_get_thread_num();

    COO coo;
    DW.d_coo_format(coo.en, U, link);
    coo.do_it();
    coo.Async( d_dDeta[m], d_eta, stream[m] );

    CuC inner;
    MatPoly XH(handle[m], stream[m]);
    XH.dotAsync<N>( &inner, d_Deta, d_dDeta[m] );

    double res = -2.0 * real(inner);
    return res;
  }


  // for looped wrapper of grad returning Gauge


};
