#pragma once


struct CGCUDA{ // wrapper
  using Complex=std::complex<double>;

  const Sparse sparse;
  const Dirac1fonS2& D;

  int *d_cols, *d_rows, *d_colsT, *d_rowsT;

  constexpr Idx N;
  constexpr Idx len;

  CGCUDA( const Dirac1fonS2& D )
    : sparse(D.lattice)
    , D(D)
    , N(sparse.N)
    , len(sparse.len)
  {
    const std::vector<int>& cols = sparse.cols_csr;
    const std::vector<int>& rows = sparse.rows_csr;
    const std::vector<int>& colsT = sparse.cols_csrT;
    const std::vector<int>& rowsT = sparse.rows_csrT;

    cudacheck(cudaMalloc(&d_cols, len*sizeof(int)));
    cudacheck(cudaMalloc(&d_rows, (N+1)*sizeof(int)));
    cudacheck(cudaMemcpy(d_cols, cols.data(), len*sizeof(int), H2D));
    cudacheck(cudaMemcpy(d_rows, rows.data(), (N+1)*sizeof(int), H2D));
    //
    cudacheck(cudaMalloc(&d_colsT, len*sizeof(int)));
    cudacheck(cudaMalloc(&d_rowsT, (N+1)*sizeof(int)));
    cudacheck(cudaMemcpy(d_colsT, colsT.data(), len*sizeof(int), H2D));
    cudacheck(cudaMemcpy(d_rowsT, rowsT.data(), (N+1)*sizeof(int), H2D));
  }

  ~CGCUDA()
  {
    cudacheck(cudaFree(d_cols));
    cudacheck(cudaFree(d_rows));
    cudacheck(cudaFree(d_colsT));
    cudacheck(cudaFree(d_rowsT));
  }

  void operator()( Complex* res, const Complex* v, const U1onS2& U ) const {

    Complex v_coo[ len ], v_csr[ len ], v_csrH[ len ];
    D.coo_format( v_coo, U );
    sparse.coo2csr_csrH( v_csr, v_csrH, v_coo );

    solve<N,len>(reinterpret_cast<CuC*>(res),
		 reinterpret_cast<const CuC*>(v),
		 reinterpret_cast<const CuC*>(v_csr),
		 // sparse.cols_csr, sparse.rows_csr,
		 reinterpret_cast<const CuC*>(v_csrH),
		 // sparse.cols_csrT, sparse.rows_csrT,
		 d_cols, d_rows, d_colsT, d_rowsT
		 );
  }


  void general_op_inv( Complex* res, const Complex* v, const U1onS2& U,
		       const double aa, const double bb ) const {

    Complex v_coo[ len ], v_csr[ len ], v_csrH[ len ];
    D.coo_format( v_coo, U );
    sparse.coo2csr_csrH( v_csr, v_csrH, v_coo );

    solve<N,len>(reinterpret_cast<CuC*>(res),
		 reinterpret_cast<const CuC*>(v),
		 reinterpret_cast<const CuC*>(v_csr),
		 // sparse.cols_csr, sparse.rows_csr,
		 reinterpret_cast<const CuC*>(v_csrH),
		 // sparse.cols_csrT, sparse.rows_csrT,
		 d_cols, d_rows, d_colsT, d_rowsT,
		 aa, bb
		 );
  }

  void generalH_op_inv( Complex* res, const Complex* v, const U1onS2& U,
			const double aa, const double bb,
			const double lambda_max ) const {

    Complex v_coo[ len ], v_csr[ len ], v_csrH[ len ];
    D.H_coo_format( v_coo, U, lambda_max);
    sparse.coo2csr_csrH( v_csr, v_csrH, v_coo );

    solve<N,len>(reinterpret_cast<CuC*>(res),
		 reinterpret_cast<const CuC*>(v),
		 reinterpret_cast<const CuC*>(v_csr),
		 // sparse.cols_csr, sparse.rows_csr,
		 reinterpret_cast<const CuC*>(v_csrH),
		 // sparse.cols_csrT, sparse.rows_csrT,
		 d_cols, d_rows, d_colsT, d_rowsT,
		 aa, bb
		 );
  }


  void general_op( Complex* res, const Complex* v, const U1onS2& U,
		   const double aa, const double bb ) const {

    Complex v_coo[ len ], v_csr[ len ], v_csrH[ len ];
    D.coo_format( v_coo, U );
    sparse.coo2csr_csrH( v_csr, v_csrH, v_coo );

    multA_wrapper<N,len>(reinterpret_cast<CuC*>(res),
			 reinterpret_cast<const CuC*>(v),
			 reinterpret_cast<const CuC*>(v_csr),
			 reinterpret_cast<const CuC*>(v_csrH),
			 d_cols, d_rows, d_colsT, d_rowsT,
			 aa, bb
			 );
  }

};





