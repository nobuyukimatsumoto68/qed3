#pragma once

struct SparseHelper{
  const Lattice& lattice;
  const Idx N;

  const bool locate_on_gpu;
  bool is_set;

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

  SparseHelper(const Lattice& lattice_, const bool locate_on_gpu_=true)
    : lattice(lattice_)
    , N(NS*lattice.n_sites)
    , locate_on_gpu(locate_on_gpu_)
    , is_set(false)
  {
  } // end of constructor


  ~SparseHelper()
  {
    if(locate_on_gpu&&is_set){
      cudacheck(cudaFree(d_cols));
      cudacheck(cudaFree(d_rows));
      cudacheck(cudaFree(d_colsT));
      cudacheck(cudaFree(d_rowsT));

      cudacheck(cudaFree(d_val));
      cudacheck(cudaFree(d_valH));
    }
  }

  void set_dirac(){
    std::vector<Idx> is;
    std::vector<Idx> js;

    // ========= COO ========= //
    for(Idx ix=0; ix<lattice.n_sites; ix++){
      for(int jj=0; jj<3; jj++){
	Idx iy = lattice.nns[ix][jj];

	is.push_back( NS*ix ); js.push_back( NS*iy );
	is.push_back( NS*ix ); js.push_back( NS*iy+1 );
	is.push_back( NS*ix ); js.push_back( NS*ix );
	is.push_back( NS*ix ); js.push_back( NS*ix+1 );
	is.push_back( NS*ix+1 ); js.push_back( NS*iy );
	is.push_back( NS*ix+1 ); js.push_back( NS*iy+1 );
	is.push_back( NS*ix+1 ); js.push_back( NS*ix );
	is.push_back( NS*ix+1 ); js.push_back( NS*ix+1 );
      }
    }
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

    if(locate_on_gpu){
      const std::vector<Idx>& cols = cols_csr;
      const std::vector<Idx>& rows = rows_csr;
      const std::vector<Idx>& colsT = cols_csrT;
      const std::vector<Idx>& rowsT = rows_csrT;

      cudacheck(cudaMalloc(&d_cols, len*sizeof(Idx)));
      cudacheck(cudaMalloc(&d_rows, (N+1)*sizeof(Idx)));
      cudacheck(cudaMemcpy(d_cols, cols.data(), len*sizeof(Idx), H2D));
      cudacheck(cudaMemcpy(d_rows, rows.data(), (N+1)*sizeof(Idx), H2D));
      //
      cudacheck(cudaMalloc(&d_colsT, len*sizeof(Idx)));
      cudacheck(cudaMalloc(&d_rowsT, (N+1)*sizeof(Idx)));
      cudacheck(cudaMemcpy(d_colsT, colsT.data(), len*sizeof(Idx), H2D));
      cudacheck(cudaMemcpy(d_rowsT, rowsT.data(), (N+1)*sizeof(Idx), H2D));
      //
      cudacheck(cudaMalloc(&d_val, len*CD));
      cudacheck(cudaMalloc(&d_valH, len*CD));
    }

    v_coo.resize(len);
    v_csr.resize(len);
    v_csrH.resize(len);

    is_set = true;
  }


  void set_diag(){
    len = N;

    // ========= CSR ========= //

    cols_csr.resize(N);
    for(Idx i=0; i<N; i++){
      cols_csr[i] = i;
      rows_csr.push_back( i );
    }
    rows_csr.push_back( N );

    cols_csrT = cols_csr;
    rows_csrT = rows_csr;

    assert( rows_csr.size()==N+1 );
    assert( rows_csrT.size()==N+1 );

    if(locate_on_gpu){
      const std::vector<Idx>& cols = cols_csr;
      const std::vector<Idx>& rows = rows_csr;
      const std::vector<Idx>& colsT = cols_csrT;
      const std::vector<Idx>& rowsT = rows_csrT;

      cudacheck(cudaMalloc(&d_cols, len*sizeof(Idx)));
      cudacheck(cudaMalloc(&d_rows, (N+1)*sizeof(Idx)));
      cudacheck(cudaMemcpy(d_cols, cols.data(), len*sizeof(Idx), H2D));
      cudacheck(cudaMemcpy(d_rows, rows.data(), (N+1)*sizeof(Idx), H2D));
      //
      cudacheck(cudaMalloc(&d_colsT, len*sizeof(Idx)));
      cudacheck(cudaMalloc(&d_rowsT, (N+1)*sizeof(Idx)));
      cudacheck(cudaMemcpy(d_colsT, colsT.data(), len*sizeof(Idx), H2D));
      cudacheck(cudaMemcpy(d_rowsT, rowsT.data(), (N+1)*sizeof(Idx), H2D));
    }

    v_coo.resize(len);
    v_csr.resize(len);
    v_csrH.resize(len);

    is_set = true;
  }

  template<typename T>
  void coo2csr( T* v_csr,
	        const T* v_coo) const {
    for(Idx ell=0; ell<len; ell++) v_csr[ ell2em[ell] ] = v_coo[ell];
  }

  template<typename T>
  void coo2csrT( T* v_csrT,
		 const T* v_coo) const {
    for(Idx ell=0; ell<len; ell++) v_csrT[ ell2emT[ell] ] = v_coo[ell];
  }

  template<typename T>
  void coo2csrH( T* v_csrH,
		 const T* v_coo) const {
    for(Idx ell=0; ell<len; ell++) v_csrH[ ell2emT[ell] ] = std::conj( v_coo[ell] );
  }

  template<typename T>
  void coo2csr_csrH( T* v_csr,
		     T* v_csrH,
		     const T* v_coo) const {
    for(Idx ell=0; ell<len; ell++) {
      v_csr[ ell2em[ell] ] = v_coo[ell];
      v_csrH[ ell2emT[ell] ] = std::conj( v_coo[ell] );
    }
  }

  void reset_DU( const Dirac1fonS2& D, const U1onS2& U ){
    D.coo_format( v_coo.data(), U );
    coo2csr_csrH<Complex>( v_csr.data(), v_csrH.data(), v_coo.data() );

    if(locate_on_gpu){
      cudacheck(cudaMemcpy(d_val, reinterpret_cast<const CuC*>(v_csr.data()), len*CD, H2D));
      cudacheck(cudaMemcpy(d_valH, reinterpret_cast<const CuC*>(v_csrH.data()), len*CD, H2D));
    }
  }
  
};
