#pragma once


template<class WilsonDirac>
struct DWDevice{
  using T = CuC;

  const WilsonDirac& D;
  const Idx N;

  static constexpr int NS = 2;

  Idx len;

  // @@@
  // std::vector<Idx> ell2em;
  // std::vector<Idx> ell2emT;
  // std::vector<Idx> cols_csr;
  // std::vector<Idx> rows_csr;
  // std::vector<Idx> cols_csrT;
  // std::vector<Idx> rows_csrT;
  // std::vector<Complex> v_coo, v_csr, v_csrH;
  Idx *ell2em, *ell2emT, *cols_csr, *rows_csr, *cols_csrT, *rows_csrT;
  Complex *v_coo, *v_csr, *v_csrH;

  Idx *d_cols, *d_rows, *d_colsT, *d_rowsT;
  CuC *d_val, *d_valH;

  std::vector<Idx> is;
  std::vector<Idx> js;


  DWDevice(const WilsonDirac& D_)
    : D(D_)
    , N(Comp::N)
  {
    initialize();
  } // end of constructor


  ~DWDevice()
  {
    CUDA_CHECK(cudaFree(d_cols));
    CUDA_CHECK(cudaFree(d_rows));
    CUDA_CHECK(cudaFree(d_colsT));
    CUDA_CHECK(cudaFree(d_rowsT));

    CUDA_CHECK(cudaFree(d_val));
    CUDA_CHECK(cudaFree(d_valH));

    CUDA_CHECK(cudaFreeHost(ell2em));
    CUDA_CHECK(cudaFreeHost(ell2emT));
    CUDA_CHECK(cudaFreeHost(cols_csr));
    CUDA_CHECK(cudaFreeHost(rows_csr));
    CUDA_CHECK(cudaFreeHost(cols_csrT));
    CUDA_CHECK(cudaFreeHost(rows_csrT));
    CUDA_CHECK(cudaFreeHost(v_coo));
    CUDA_CHECK(cudaFreeHost(v_csr));
    CUDA_CHECK(cudaFreeHost(v_csrH));
  }

  // @@@
  void initialize(){
    // ========= COO ========= //

    D.coo_structure(is, js);
    len = js.size();
    assert( is.size()==len );

    // ========= CSR ========= //

    Idx em=0;
    Idx emT=0;

    // ell2em.resize(len);
    // cols_csr.resize(len);
    CUDA_CHECK( cudaMallocHost( &ell2em, len*sizeof(Idx) ) );
    CUDA_CHECK( cudaMallocHost( &cols_csr, len*sizeof(Idx) ) );
    memset(ell2em, 0, len*sizeof(Idx));
    memset(cols_csr, 0, len*sizeof(Idx));

    // rows_csr.clear();
    CUDA_CHECK( cudaMallocHost( &rows_csr, (N+1)*sizeof(Idx) ) );
    memset(rows_csr, 0, (N+1)*sizeof(Idx));

    // ell2emT.resize(len);
    // cols_csrT.resize(len);
    CUDA_CHECK( cudaMallocHost( &ell2emT, len*sizeof(Idx) ) );
    CUDA_CHECK( cudaMallocHost( &cols_csrT, len*sizeof(Idx) ) );
    memset(ell2emT, 0, len*sizeof(Idx));
    memset(cols_csrT, 0, len*sizeof(Idx));

    // rows_csrT.clear();
    CUDA_CHECK( cudaMallocHost( &rows_csrT, (N+1)*sizeof(Idx) ) );
    memset(rows_csrT, 0, (N+1)*sizeof(Idx));

    Idx counter=0;
    // rows_csr.push_back( em );
    // rows_csrT.push_back( em );
    rows_csr[counter] = em;
    rows_csrT[counter] = em;
    counter++;


    std::vector<std::vector<Idx>> ell_record(N);
    std::vector<std::vector<Idx>> ell_recordT(N);

#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL_SORT)
#endif
    for(Idx i=0; i<N; i++){
      for(Idx ell=0; ell<len; ell++){
	if( is[ell]==i ) ell_record[i].push_back(ell);
	if( js[ell]==i ) ell_recordT[i].push_back(ell);
      }
    }

    for(Idx i=0; i<N; i++){
      for(const Idx ell : ell_record[i]){
        ell2em[ell] = em;
        cols_csr[em] = js[ell];
        ++em;
      }
      for(const Idx ell : ell_recordT[i]){
        ell2emT[ell] = emT;
        cols_csrT[emT] = is[ell];
        ++emT;
      }
      // rows_csr.push_back( em );
      // rows_csrT.push_back( emT );
      rows_csr[counter] = em;
      rows_csrT[counter] = emT;
      counter++;
    }

    assert(counter==N+1);
    // assert( rows_csr.size()==N+1 );
    // assert( rows_csrT.size()==N+1 );

    // const std::vector<Idx>& cols = cols_csr;
    // const std::vector<Idx>& rows = rows_csr;
    // const std::vector<Idx>& colsT = cols_csrT;
    // const std::vector<Idx>& rowsT = rows_csrT;

    CUDA_CHECK(cudaMalloc(&d_cols, len*sizeof(Idx)));
    CUDA_CHECK(cudaMalloc(&d_rows, (N+1)*sizeof(Idx)));
    // CUDA_CHECK(cudaMemcpy(d_cols, cols.data(), len*sizeof(Idx), H2D));
    // CUDA_CHECK(cudaMemcpy(d_rows, rows.data(), (N+1)*sizeof(Idx), H2D));
    CUDA_CHECK(cudaMemcpy(d_cols, cols_csr, len*sizeof(Idx), H2D));
    CUDA_CHECK(cudaMemcpy(d_rows, rows_csr, (N+1)*sizeof(Idx), H2D));
    //
    CUDA_CHECK(cudaMalloc(&d_colsT, len*sizeof(Idx)));
    CUDA_CHECK(cudaMalloc(&d_rowsT, (N+1)*sizeof(Idx)));
    // CUDA_CHECK(cudaMemcpy(d_colsT, colsT.data(), len*sizeof(Idx), H2D));
    // CUDA_CHECK(cudaMemcpy(d_rowsT, rowsT.data(), (N+1)*sizeof(Idx), H2D));
    CUDA_CHECK(cudaMemcpy(d_colsT, cols_csrT, len*sizeof(Idx), H2D));
    CUDA_CHECK(cudaMemcpy(d_rowsT, rows_csrT, (N+1)*sizeof(Idx), H2D));
    //
    CUDA_CHECK(cudaMalloc(&d_val, len*CD));
    CUDA_CHECK(cudaMalloc(&d_valH, len*CD));

    // v_coo.resize(len);
    // v_csr.resize(len);
    // v_csrH.resize(len);
    CUDA_CHECK( cudaMallocHost( &v_coo, len*CD ) );
    CUDA_CHECK( cudaMallocHost( &v_csr, len*CD ) );
    CUDA_CHECK( cudaMallocHost( &v_csrH, len*CD ) );
    memset(v_coo, 0, len*CD);
    memset(v_csr, 0, len*CD);
    memset(v_csrH, 0, len*CD);
  }

  // @@@
  // void coo2csr_csrH( std::vector<Complex>& v_csr,
  //       	     std::vector<Complex>& v_csrH,
  //       	     const std::vector<Complex>& v_coo) const {
  void coo2csr_csrH() const {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL)
#endif
    for(Idx ell=0; ell<len; ell++) {
      v_csr[ ell2em[ell] ] = v_coo[ell];
      v_csrH[ ell2emT[ell] ] = std::conj( v_coo[ell] );
    }
  }

  template<typename Gauge>
  void update( const Gauge& U ){
    D.coo_format( v_coo, U );
    // coo2csr_csrH( v_csr, v_csrH, v_coo );
    coo2csr_csrH();

    // CUDA_CHECK(cudaMemcpy(d_val, reinterpret_cast<const CuC*>(v_csr.data()), len*CD, H2D));
    // CUDA_CHECK(cudaMemcpy(d_valH, reinterpret_cast<const CuC*>(v_csrH.data()), len*CD, H2D));
    CUDA_CHECK(cudaMemcpy(d_val, v_csr, len*CD, H2D));
    CUDA_CHECK(cudaMemcpy(d_valH, v_csrH, len*CD, H2D));
  }


  void associateCSR( CSR<Comp::N>& M, const bool is_transpose=false ){
    if(!is_transpose){
      M.cols = this->d_cols;
      M.rows = this->d_rows;
      M.val = this->d_val;
    }
    else{
      M.cols = this->d_colsT;
      M.rows = this->d_rowsT;
      M.val = this->d_valH;
    }
  }
};
