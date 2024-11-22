#pragma once

struct Sparse{
  const Lattice& lattice;
  const int N;

  static constexpr int NS = 2;


  int len;

  std::vector<int> is;
  std::vector<int> js;

  std::vector<int> ell2em;
  std::vector<int> cols_csr;
  std::vector<int> rows_csr;

  std::vector<int> ell2emT;
  std::vector<int> cols_csrT;
  std::vector<int> rows_csrT;

  Sparse(const Lattice& lattice_)
    : lattice(lattice_)
    , N(NS*lattice.n_sites)
  {
    // ========= COO ========= //

    for(int ix=0; ix<lattice.n_sites; ix++){
      for(int jj=0; jj<3; jj++){
	int iy = lattice.nns[ix][jj];

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

    int em=0;
    int emT=0;

    ell2em.resize(len);
    cols_csr.resize(len);
    rows_csr.clear();

    ell2emT.resize(len);
    cols_csrT.resize(len);
    rows_csrT.clear();

    rows_csr.push_back( em );
    rows_csrT.push_back( em );
    for(int i=0; i<N; i++){
      for(int ell=0; ell<len; ell++){
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
  } // end of constructor


  template<typename T>
  void coo2csr( T* v_csr,
	        const T* v_coo) const {
    for(int ell=0; ell<len; ell++) v_csr[ ell2em[ell] ] = v_coo[ell];
  }

  template<typename T>
  void coo2csrT( T* v_csrT,
		 const T* v_coo) const {
    for(int ell=0; ell<len; ell++) v_csrT[ ell2emT[ell] ] = v_coo[ell];
  }

  template<typename T>
  void coo2csrH( T* v_csrH,
		 const T* v_coo) const {
    for(int ell=0; ell<len; ell++) v_csrH[ ell2emT[ell] ] = std::conj( v_coo[ell] );
  }

  template<typename T>
  void coo2csr_csrH( T* v_csr,
		     T* v_csrH,
		     const T* v_coo) const {
    for(int ell=0; ell<len; ell++) {
      v_csr[ ell2em[ell] ] = v_coo[ell];
      v_csrH[ ell2emT[ell] ] = std::conj( v_coo[ell] );
    }
  }

  template<typename T>
  void mult( T* res, const T* v,
	     const T* v_csr ) const {
    assert(res!=v);
    for(int i=0; i<N; i++){
      res[i] = 0.0;
      const int row_start = rows_csr[i];
      const int row_end = rows_csr[i+1];
      for(int jj=row_start; jj<row_end; jj++){
	res[i] += v_csr[jj] * v[ cols_csr[jj] ];
      }
    }
  }

  template<typename T>
  void multT( T* res, const T* v,
	      const T* v_csrT ) const {
    assert(res!=v);
    for(int i=0; i<N; i++){
      res[i] = 0.0;
      const int row_start = rows_csrT[i];
      const int row_end = rows_csrT[i+1];
      for(int jj=row_start; jj<row_end; jj++){
	res[i] += v_csrT[jj] * v[ cols_csrT[jj] ];
      }
    }
  }

  // for customized sparse structure such as for the derivatives
  template <typename T>
  void multcoo( std::vector<T>& res, const std::vector<T>& v,
		const std::vector<T>& val, const std::vector<int>& isC, const std::vector<int>& jsC ) const {
    res.resize( v.size() );
    for(int i=0; i<v.size(); i++) res[i] = 0.0;
    for(int i=0; i<isC.size(); i++) res[isC[i]] += val[i] * v[jsC[i]];
  }




};
