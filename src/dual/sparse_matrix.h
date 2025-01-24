#pragma once


struct LinOp{
  using T = CuC;
  virtual void operator()( T* d_res, const T* d_v ) const = 0;
};


struct CSR : public LinOp {
  using T = CuC;

  T* val;
  Idx* cols;
  Idx* rows;

  void operator()( T* d_res, const T* d_v ) const {
    constexpr Idx N = Comp::N;
    mult<T,N><<<NBlocks, NThreadsPerBlock>>>(d_res,
					     d_v,
					     this->val,
					     this->cols,
					     this->rows);
  }

};


template<typename F>
struct LinOpWrapper : public LinOp {
  using T = CuC;
  const F& f;

  LinOpWrapper(const F& f_ )
    : f(f_)
  {}

  void operator()( T* d_res, const T* d_v ) const {
    f( d_res, d_v );
  }
};


struct COOEntry {
  using T = CuC;
  T v;
  Idx i;
  Idx j;

  // COOEntry() = delete;
  COOEntry( const T v, const Idx i, const Idx j )
    : v(v)
    , i(i)
    , j(j)
  {}

  COOEntry( const Complex v_, const Idx i, const Idx j )
    : v(cplx(v_))
    , i(i)
    , j(j)
  {}

  bool operator<(const COOEntry& other) const {
    bool res = false;

    if(this->i < other.i) {
      res = true;
      return res;
    }
    else if(this->j == other.j){
      if(this->j < other.j) {
        res = true;
        return res;
      }
    }
    return res;
  }

  friend std::ostream& operator<<(std::ostream& stream, const COOEntry& obj) {
    stream << real(obj.v) << " "<< imag(obj.v) << " " << obj.i << " " << obj.j << std::flush; // << std::endl
    return stream;
  }

};


// for gradients
struct COO : public LinOp {
  using T = CuC;

  std::vector<COOEntry> en;

  //

  T* d_val;
  Idx* d_cols;
  Idx* d_rows;

  bool is_set;

  COO()
    : is_set(false)
  {}

  ~COO()
  {
    if(is_set){
      CUDA_CHECK(cudaFree(d_val));
      CUDA_CHECK(cudaFree(d_cols));
      CUDA_CHECK(cudaFree(d_rows));
    }
  }

  void do_it(){
    constexpr Idx N = Comp::N;
    std::sort( en.begin(), en.end() );

    const Idx len=en.size();
    std::vector<CuC> v(len);
    std::vector<Idx> cols(len);
    std::vector<Idx> rows;

    for(Idx k=0; k<len; k++){
      v[k] = en[k].v;
      cols[k] = en[k].j;
    }

    Idx k=0;
    rows.push_back(k);
    for(Idx i=0; i<N; i++){
      while(en[k].i == i) k++;
      rows.push_back(k);
    }

    CUDA_CHECK(cudaMalloc(&d_cols, len*sizeof(Idx)));
    CUDA_CHECK(cudaMalloc(&d_rows, (N+1)*sizeof(Idx)));
    CUDA_CHECK(cudaMalloc(&d_val, len*CD));

    CUDA_CHECK(cudaMemcpy(d_cols, cols.data(), len*sizeof(Idx), H2D));
    CUDA_CHECK(cudaMemcpy(d_rows, rows.data(), (N+1)*sizeof(Idx), H2D));
    CUDA_CHECK(cudaMemcpy(d_val, reinterpret_cast<const CuC*>(v.data()), len*CD, H2D));

    is_set = true;
  }

  void do_conjugate(){
    constexpr Idx N = Comp::N;
    const Idx len=en.size();

    std::vector<COOEntry> enH;
    for(Idx k=0; k<len; k++) enH.push_back( COOEntry( conj(en[k].v), en[k].j, en[k].i ) );

    std::sort( enH.begin(), enH.end() );

    std::vector<CuC> v(len);
    std::vector<Idx> cols(len);
    std::vector<Idx> rows;

    for(Idx k=0; k<len; k++){
      v[k] = enH[k].v;
      cols[k] = enH[k].j;
    }

    Idx k=0;
    rows.push_back(k);
    for(Idx i=0; i<N; i++){
      while(enH[k].i == i) k++;
      rows.push_back(k);
    }

    CUDA_CHECK(cudaMalloc(&d_cols, len*sizeof(Idx)));
    CUDA_CHECK(cudaMalloc(&d_rows, (N+1)*sizeof(Idx)));
    CUDA_CHECK(cudaMalloc(&d_val, len*CD));

    CUDA_CHECK(cudaMemcpy(d_cols, cols.data(), len*sizeof(Idx), H2D));
    CUDA_CHECK(cudaMemcpy(d_rows, rows.data(), (N+1)*sizeof(Idx), H2D));
    CUDA_CHECK(cudaMemcpy(d_val, reinterpret_cast<const CuC*>(v.data()), len*CD, H2D));

    is_set = true;
  }


  void operator()( T* d_res, const T* d_v ) const {
    constexpr Idx N = Comp::N;
    mult<T,N><<<NBlocks, NThreadsPerBlock>>>(d_res,
					     d_v,
					     this->d_val,
					     this->d_cols,
					     this->d_rows);
  }

};



template <Idx N>
void matmulcoo( CuC* res, const CuC* v,
		const std::vector<COOEntry>& coo) {
  for(int i=0; i<N; i++) res[i] = cplx(0.0);
  for(int k=0; k<coo.size(); k++) res[coo[k].i] = res[coo[k].i] + coo[k].v * v[coo[k].j];
}


template <Idx N>
void matmul( CuC* res, const CuC* v,
	     const std::vector<CuC>& val,
	     const std::vector<Idx>& cols,
	     const std::vector<Idx>& rows ) {
  for(int i=0; i<N; i++){
    res[i] = cplx(0.0);
    const int row_start = rows[i];
    const int row_end = rows[i+1];
    for(int jj=row_start; jj<row_end; jj++){
      res[i] = res[i] + val[jj] * v[cols[jj]];
    }
  }
}

