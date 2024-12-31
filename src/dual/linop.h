#pragma once

template<typename T>
struct LinOp{
  std::vector<std::vector<const SparseMatrix<T>*>> vec_mats;
  std::vector<T> coeffs;

  LinOp( const std::initializer_list<int> struc={1} )
    : vec_mats(struc.size())
    , coeffs(struc.size())
  {
    auto itr = struc.begin();
    for(int i=0; i<struc.size(); i++) {
      vec_mats[i].resize( *itr );
      ++itr;
    }
  }

  int size() const { return vec_mats.size(); }
  std::vector<const SparseMatrix<T>*> operator[](const int i) const { return vec_mats[i]; }

  void set_coeff( const int i, const T coeff ){ coeffs[i] = coeff; }

  void set_matrix( const int i, const int j,
		   const SparseMatrix<T>* mat_ptr ){
    vec_mats[i][j] = mat_ptr;
  }


  template<Idx N>
  void on_cpu( std::vector<Complex>& v, const std::vector<Complex>& v0 ) const {
    std::vector<Complex> tmp(N), Mv0(N);

    for(int i=0; i<vec_mats.size(); i++){
      tmp = v0;
      Mv0 = v0;
      for(int j=0; j<vec_mats[i].size(); j++){
	vec_mats[i][j]->act_cpu( Mv0, tmp );
	tmp = Mv0;
      }
      for(Idx j=0; j<N; j++) v[j] += coeffs[i] * Mv0[j];
    }
  }


  template<Idx N>
  void from_cpu( std::vector<Complex>& v, const std::vector<Complex>& v0 ) const {
    CuC *d_v, *d_v0;
    cudacheck(cudaMalloc(&d_v, N*CD));
    cudacheck(cudaMalloc(&d_v0, N*CD));

    cudacheck(cudaMemcpy(d_v0, reinterpret_cast<const CuC*>(v0.data()), N*CD, H2D));
    cudacheck(cudaMemset(d_v, 0, N*CD));
    
    on_gpu<N>( d_v, d_v0 );

    cudacheck(cudaMemcpy(reinterpret_cast<CuC*>(v.data()), d_v, N*CD, D2H));

    cudacheck(cudaFree(d_v));
    cudacheck(cudaFree(d_v0));
  }

  template<Idx N> __host__
  void on_gpu(CuC* d_v, const CuC* d_v0) const {
    CuC *d_tmp, *d_Mv0, *d_coeffs;
    cudacheck(cudaMalloc(&d_tmp, N*CD));
    cudacheck(cudaMalloc(&d_Mv0, N*CD));
    cudacheck(cudaMalloc(&d_coeffs, coeffs.size()*CD));

    cudacheck(cudaMemset(d_tmp, 0, N*CD));
    cudacheck(cudaMemset(d_Mv0, 0, N*CD));
    cudacheck(cudaMemcpy(d_coeffs, coeffs.data(), coeffs.size()*CD, H2D));

    for(int i=0; i<vec_mats.size(); i++){
      cudacheck(cudaMemcpy(d_tmp, d_v0, N*CD, D2D));
      cudacheck(cudaMemcpy(d_Mv0, d_v0, N*CD, D2D));
      for(int j=0; j<vec_mats[i].size(); j++){
	mult<T, N><<<NBlocks, NThreadsPerBlock>>>(d_Mv0,
						  d_tmp,
						  vec_mats[i][j]->val,
						  vec_mats[i][j]->cols,
						  vec_mats[i][j]->rows);
	cudacheck(cudaMemcpy(d_tmp, d_Mv0, N*CD, D2D));
      }
      ScalarMult<CuC, N><<<NBlocks, NThreadsPerBlock>>>(d_v,
							d_coeffs+i,
							d_Mv0);
    }

    cudacheck(cudaFree(d_tmp));
    cudacheck(cudaFree(d_Mv0));
    cudacheck(cudaFree(d_coeffs));
  }


};


