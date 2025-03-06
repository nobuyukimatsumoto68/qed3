#pragma once

// template<class WilsonDirac, class Lattice>
template<class WilsonDirac>
struct DWDevice{
  using T = CuC;
  // using WilsonDirac=Dirac1fonS2;
  // using Lattice=S2Trivalent;

  const WilsonDirac& D;
  // const Lattice& lattice;
  const Idx N;

  static constexpr int NS = 2;

  Idx len;

  std::vector<Idx> ell2em;
  std::vector<Idx> ell2emT;

  std::vector<Idx> cols_csr;
  std::vector<Idx> rows_csr;

  std::vector<Idx> cols_csrT;
  std::vector<Idx> rows_csrT;

  std::vector<Complex> v_coo, v_csr, v_csrH;

  Idx *d_cols, *d_rows, *d_colsT, *d_rowsT;
  CuC *d_val, *d_valH;

  DWDevice(const WilsonDirac& D_)
    : D(D_)
    // , lattice(D.lattice)
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
  }

  void initialize(){
    std::vector<Idx> is;
    std::vector<Idx> js;

    // ========= COO ========= //

    D.coo_structure(is, js);
    len = js.size();
    assert( is.size()==len );

    // ========= CSR ========= //

    Idx em=0;
    Idx emT=0;

    ell2em.resize(len);
    cols_csr.resize(len);
    rows_csr.clear();

    ell2emT.resize(len);
    cols_csrT.resize(len);
    rows_csrT.clear();

    rows_csr.push_back( em );
    rows_csrT.push_back( em );
    for(Idx i=0; i<N; i++){
      for(Idx ell=0; ell<len; ell++){
	if( is[ell]==i ){
	  ell2em[ell] = em;
	  cols_csr[em] = js[ell];
	  ++em;
	}
	if( js[ell]==i ){
	  ell2emT[ell] = emT;
	  cols_csrT[emT] = is[ell];
	  ++emT;
	}
      }
      rows_csr.push_back( em );
      rows_csrT.push_back( emT );
    }

    assert( rows_csr.size()==N+1 );
    assert( rows_csrT.size()==N+1 );

    const std::vector<Idx>& cols = cols_csr;
    const std::vector<Idx>& rows = rows_csr;
    const std::vector<Idx>& colsT = cols_csrT;
    const std::vector<Idx>& rowsT = rows_csrT;

    CUDA_CHECK(cudaMalloc(&d_cols, len*sizeof(Idx)));
    CUDA_CHECK(cudaMalloc(&d_rows, (N+1)*sizeof(Idx)));
    CUDA_CHECK(cudaMemcpy(d_cols, cols.data(), len*sizeof(Idx), H2D));
    CUDA_CHECK(cudaMemcpy(d_rows, rows.data(), (N+1)*sizeof(Idx), H2D));
    //
    CUDA_CHECK(cudaMalloc(&d_colsT, len*sizeof(Idx)));
    CUDA_CHECK(cudaMalloc(&d_rowsT, (N+1)*sizeof(Idx)));
    CUDA_CHECK(cudaMemcpy(d_colsT, colsT.data(), len*sizeof(Idx), H2D));
    CUDA_CHECK(cudaMemcpy(d_rowsT, rowsT.data(), (N+1)*sizeof(Idx), H2D));
    //
    CUDA_CHECK(cudaMalloc(&d_val, len*CD));
    CUDA_CHECK(cudaMalloc(&d_valH, len*CD));

    v_coo.resize(len);
    v_csr.resize(len);
    v_csrH.resize(len);

    // is_set = true;
  }

  // // template<typename T>
  // void coo2csr( T* v_csr,
  //               const T* v_coo) const {
  //   for(Idx ell=0; ell<len; ell++) v_csr[ ell2em[ell] ] = v_coo[ell];
  // }

  // // template<typename T>
  // void coo2csrT( T* v_csrT,
  //       	 const T* v_coo) const {
  //   for(Idx ell=0; ell<len; ell++) v_csrT[ ell2emT[ell] ] = v_coo[ell];
  // }

  // // template<typename T>
  // void coo2csrH( T* v_csrH,
  //       	 const T* v_coo) const {
  //   for(Idx ell=0; ell<len; ell++) v_csrH[ ell2emT[ell] ] = conj( v_coo[ell] );
  // }

  // template<typename T>
  void coo2csr_csrH( std::vector<Complex>& v_csr,
		     std::vector<Complex>& v_csrH,
		     const std::vector<Complex>& v_coo) const {
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
    coo2csr_csrH( v_csr, v_csrH, v_coo );

    CUDA_CHECK(cudaMemcpy(d_val, reinterpret_cast<const CuC*>(v_csr.data()), len*CD, H2D));
    CUDA_CHECK(cudaMemcpy(d_valH, reinterpret_cast<const CuC*>(v_csrH.data()), len*CD, H2D));
    //    }
  }

  // // ugliness
  // void associate( SparseMatrix& M, const bool is_transpose=false ){
  //   assert(!locate_on_gpu);
  //   M.on_gpu = false;
  //   if(!is_transpose){
  //     M.cols = this->cols_csr.data();
  //     M.rows = this->rows_csr.data();
  //     M.val = this->v_csr.data();
  //   }
  //   else{
  //     M.cols = this->cols_csrT.data();
  //     M.rows = this->rows_csrT.data();
  //     M.val = this->v_csrH.data();
  //   }
  // }

  void associateCSR( CSR& M, const bool is_transpose=false ){
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
