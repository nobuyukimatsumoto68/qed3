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

  // void act_gpu2( T* d_res, const T* d_v ) const {
  //   assert( on_gpu );

  //   cusparseHandle_t handle = NULL;

  //   const cusparseDirection_t dirA = CUSPARSE_DIRECTION_COLUMN;
  //   constexpr Idx blockDim = ComilationConst::N;
  //   constexpr int mb = 1;
  //   constexpr int nb = 1;
  //   constexpr int nnzb = 1;

  //   int* bsrRowPtr;
  //   cudaMalloc((void**)&bsrRowPtrC, sizeof(int) *(mb+1));

  //   CuC* descrA
  //   cusparseXcsr2bsrNnz(handle, dirA, m, n,

  // 			descrA, csrRowPtrA, csrColIndA, blockDim,
  // 			descrC, bsrRowPtrC, &nnzb);

  //   cudaMalloc((void**)&bsrColIndC, sizeof(int)*nnzb);
  //   cudaMalloc((void**)&bsrValC, sizeof(float)*(blockDim*blockDim)*nnzb);
  //   cusparseScsr2bsr(handle, dirA, m, n,
  // 		     descrA, csrValA, csrRowPtrA, csrColIndA, blockDim,
  // 		     descrC, bsrValC, bsrRowPtrC, bsrColIndC);
  //   // step 2: allocate vector x and vector y large enough for bsrmv
  //   cudaMalloc((void**)&x, sizeof(float)*(nb*blockDim));
  //   cudaMalloc((void**)&y, sizeof(float)*(mb*blockDim));
  //   cudaMemcpy(x, hx, sizeof(float)*n, cudaMemcpyHostToDevice);
  //   cudaMemcpy(y, hy, sizeof(float)*m, cudaMemcpyHostToDevice);
  //   // step 3: perform bsrmv
  //   cusparseSbsrmv(handle, dirA, transA, mb, nb, nnzb, &alpha,
  // 		   descrC, bsrValC, bsrRowPtrC, bsrColIndC, blockDim, x, &beta, y);

  //     cusparseZbsrmv(cusparseHandle_t         handle,
  //              cusparseDirection_t      dir,
  //              cusparseOperation_t      trans,
  //              int                      mb,
  //              int                      nb,
  //              int                      nnzb,
  //              const cuDoubleComplex*   alpha,
  //              const cusparseMatDescr_t descr,
  //              const cuDoubleComplex*   bsrVal,
  //              const int*               bsrRowPtr,
  //              const int*               bsrColInd,
  //              int                      blockDim,
  //              const cuDoubleComplex*   x,
  //              const cuDoubleComplex*   beta,
  //              cuDoubleComplex*         y)


  //   constexpr Idx N = CompilationConst::N;
  //   mult<T,N><<<NBlocks, NThreadsPerBlock>>>(d_res,
  // 					     d_v,
  // 					     this->val,
  // 					     this->cols,
  // 					     this->rows);
  // }



};

