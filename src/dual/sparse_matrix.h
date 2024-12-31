#pragma once

template<typename T>
struct SparseMatrix{
  bool on_gpu = false;
  T* val;
  Idx* cols;
  Idx* rows;

  void assign_from_helper_DW_gpu( const SparseHelper& H_DW ){
    on_gpu = true;
    cols = H_DW.d_cols;
    rows = H_DW.d_rows;
    val = H_DW.d_val;
  }

  void assign_from_helper_DWH_gpu( const SparseHelper& H_DW ){
    on_gpu = true;
    cols = H_DW.d_colsT;
    rows = H_DW.d_rowsT;
    val = H_DW.d_valH;
  }

  void assign_from_helper_DW_cpu( SparseHelper& H_DW ){
    on_gpu = false;
    cols = H_DW.cols_csr.data();
    rows = H_DW.rows_csr.data();
    val = H_DW.v_csr.data();
  }

  void assign_from_helper_DWH_cpu( SparseHelper& H_DW ){
    on_gpu = false;
    cols = H_DW.cols_csrT.data();
    rows = H_DW.rows_csrT.data();
    val = H_DW.v_csrH.data();
  }

  void act_cpu( std::vector<Complex>& res, const std::vector<Complex>& v ) const {
    assert( !on_gpu );

    constexpr Idx N = CompilationConst::N;
    for(Idx i=0; i<N; i++) {
      res[i] = 0.0;
      const int row_start = rows[i];
      const int row_end = rows[i+1];
      for(int jj=row_start; jj<row_end; jj++) res[i] = res[i] + val[jj] * v[ cols[jj] ];
    }
  }

};

