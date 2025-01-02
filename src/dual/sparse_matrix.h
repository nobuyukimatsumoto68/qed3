#pragma once

template<typename T>
struct SparseMatrix{
  bool on_gpu = false;
  T* val;
  Idx* cols;
  Idx* rows;

  void act_cpu( std::vector<T>& res, const std::vector<T>& v ) const {
    assert( !on_gpu );

    constexpr Idx N = CompilationConst::N;
    for(Idx i=0; i<N; i++) {
      res[i] = 0.0;
      const int row_start = rows[i];
      const int row_end = rows[i+1];
      for(int jj=row_start; jj<row_end; jj++) res[i] = res[i] + val[jj] * v[ cols[jj] ];
    }
  }


  void act_gpu( T* d_res, const T* d_v ) const {
    assert( on_gpu );

    constexpr Idx N = CompilationConst::N;
    mult<T,N><<<NBlocks, NThreadsPerBlock>>>(d_res,
					     d_v,
					     this->val,
					     this->cols,
					     this->rows);
  }



};

