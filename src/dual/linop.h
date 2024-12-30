#pragma once

template<typename T>
struct LinOp{
  std::vector<std::vector<SparseMatrix<T>>> vec_mats;
  std::vector<T> coeffs;

  LinOp( const std::initializer_list<int> struc )
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
  std::vector<SparseMatrix<T>> operator[](const int i) const { return vec_mats[i]; }

  void set_coeff( const int i, const T coeff ){ coeffs[i] = coeff; }

  void set_matrix( const int i, const int j,
		   const SparseMatrix<T>& mat ){
    vec_mats[i][j] = mat;
  }

};


