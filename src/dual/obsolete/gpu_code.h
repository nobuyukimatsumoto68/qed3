#pragma once


__device__ __host__ cuDoubleComplex cplx(const double c) { return make_cuDoubleComplex(c, 0.0); }
__device__ __host__ cuDoubleComplex cplx(const std::complex<double> c) { return make_cuDoubleComplex(c.real(), c.imag()); }
__device__ __host__ cuDoubleComplex  operator+(cuDoubleComplex a, cuDoubleComplex b) { return cuCadd(a,b); }
__device__ __host__ cuDoubleComplex  operator+(cuDoubleComplex a, double b) { return cuCadd(a,cplx(b)); }
__device__ __host__ cuDoubleComplex  operator+(double a, cuDoubleComplex b) { return cuCadd(cplx(a),b); }
__device__ __host__ cuDoubleComplex  operator-(cuDoubleComplex a, cuDoubleComplex b) { return cuCsub(a,b); }
__device__ __host__ cuDoubleComplex  operator-(cuDoubleComplex a, double b) { return cuCsub(a,cplx(b)); }
__device__ __host__ cuDoubleComplex  operator-(double a, cuDoubleComplex b) { return cuCsub(cplx(a),b); }
__device__ __host__ cuDoubleComplex  operator-(cuDoubleComplex b) { return cplx(0.0)-b; }
__device__ __host__ cuDoubleComplex  operator*(cuDoubleComplex a, cuDoubleComplex b) { return cuCmul(a,b); }
__device__ __host__ cuDoubleComplex  operator*(cuDoubleComplex a, double b) { return cuCmul(a,cplx(b)); }
__device__ __host__ cuDoubleComplex  operator*(double a, cuDoubleComplex b) { return cuCmul(cplx(a),b); }
__device__ __host__ cuDoubleComplex  operator/(cuDoubleComplex a, cuDoubleComplex b) { return cuCdiv(a,b); }
__device__ __host__ cuDoubleComplex  operator/(cuDoubleComplex a, double b) { return cuCdiv(a,cplx(b)); }
__device__ __host__ cuDoubleComplex  operator/(double a, cuDoubleComplex b) { return cuCdiv(cplx(a),b); }
__device__ __host__ inline double real(const cuDoubleComplex c ){ return cuCreal(c); }
__device__ __host__ inline double imag(const cuDoubleComplex c ){ return cuCimag(c); }
__device__ __host__ inline cuDoubleComplex conj(cuDoubleComplex c ){ return cuConj(c); }



__host__
void cudacheck( cudaError status ){
  if(status!=0) std::cout << status << std::endl;
  assert(cudaSuccess == status);
}

// __host__
// void set2zero( double* v, const Idx size ){ for(Idx i=0; i<size; i++) v[i] = 0.0; }
// __host__
// void set2zero( CuC* v, const Idx size ){ for(Idx i=0; i<size; i++) v[i] = cplx(0.0); }
// __host__
// void set2zero( unsigned long* v, const Idx size ){ for(Idx i=0; i<size; i++) v[i] = 0; }
// __host__ __device__
// int mod(const int a, const int b){ return (b +(a%b))%b; }

template<Idx N> __global__
void daxpy(CuC* d_res, const CuC* d_a, const CuC* d_x, const CuC* d_y){
  Idx i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<N) d_res[i] = *d_a * d_x[i] + d_y[i];
}


// __global__
// void daxpby(CuC* d_res, double a, CuC* d_x, double b, CuC* d_y, const int N){
//   Idx i = blockIdx.x*blockDim.x + threadIdx.x;
//   if(i<N) d_res[i] = a * d_x[i] + b * d_y[i];
// }

// __global__
// template<Idx N>
// void add(CuC* d_res, const CuC* d_x, const CuC* d_y){
//   Idx i = blockIdx.x*blockDim.x + threadIdx.x;
//   if(i<N) d_res[i] = d_x[i] + d_y[i];
// }



template<Idx N> __global__
void mult ( CuC* res,
	    const CuC* v,
	    const CuC* v_csr,
	    const Idx* cols,
	    const Idx* rows
	    ){
  Idx i = blockIdx.x*blockDim.x + threadIdx.x;

  if(i<N) {
    res[i] = cplx(0.0);
    const int row_start = rows[i];
    const int row_end = rows[i+1];
    for(int jj=row_start; jj<row_end; jj++) res[i] = res[i] + v_csr[jj] * v[ cols[jj] ];
  }
}



// __host__
// template<Idx N>
// void multA(CuC* d_v, CuC* d_tmp, CuC* d_v0,
// 	   CuC* d_val, int* d_cols, int* d_rows,
// 	   CuC* d_valH, int* d_colsT, int* d_rowsT
// 	   ){
//   cudacheck(cudaMemset(d_tmp, 0, N*CD));
//   mult<N><<<NBlocks, NThreadsPerBlock>>>(d_tmp, d_v0, d_val, d_cols, d_rows);

//   cudacheck(cudaMemset(d_v, 0, N*CD));
//   mult<N><<<NBlocks, NThreadsPerBlock>>>(d_v, d_tmp, d_valH, d_colsT, d_rowsT);
// }


// __host__
// template<Idx N>
// void multD(CuC* d_v, CuC* d_v0,
// 	   CuC* d_val, int* d_cols, int* d_rows
// 	   ){
//   cudacheck(cudaMemset(d_v, 0, N*CD));
//   mult<N><<<NBlocks, NThreadsPerBlock>>>(d_v, d_v0, d_val, d_cols, d_rows);
// }


// __host__
// template<Idx N>
// void multA(CuC* d_v, CuC* d_tmp, CuC* d_v0,
// 	   CuC* d_val, int* d_cols, int* d_rows,
// 	   CuC* d_valH, int* d_colsT, int* d_rowsT,
// 	   const double a, const double b
// 	   ){
//   cudacheck(cudaMemset(d_tmp, 0, N*CD));
//   mult<N><<<NBlocks, NThreadsPerBlock>>>(d_tmp, d_v0, d_val, d_cols, d_rows);

//   cudacheck(cudaMemset(d_v, 0, N*CD));
//   mult<N><<<NBlocks, NThreadsPerBlock>>>(d_v, d_tmp, d_valH, d_colsT, d_rowsT);

//   cudacheck(cudaMemset(d_tmp, 0, N*CD));
//   daxpby<N><<<NBlocks, NThreadsPerBlock>>>(d_v, a, d_v, b, d_v0);
// }


template<Idx N> __host__
void linop_wrapper(CuC* d_v, CuC* d_tmp, CuC* d_Mv0, CuC* d_v0,
		   // const std::vector<std::vector<SparseMatrix>>& d_mats
		   const LinOp<CuC>& d_op
		   ){
  cudacheck(cudaMemset(d_tmp, 0, N*CD));
  cudacheck(cudaMemset(d_Mv0, 0, N*CD));
  cudacheck(cudaMemset(d_v, 0, N*CD));

  for(int i=0; i<d_op.size(); i++){
    cudacheck(cudaMemcpy(d_tmp, d_v0, CD, D2D));
    for(int j=0; j<d_op[i].size(); j++){
      mult<N><<<NBlocks, NThreadsPerBlock>>>(d_Mv0, d_tmp,
					     d_op[i][j].val,
					     d_op[i][j].cols,
					     d_op[i][j].rows);
      cudacheck(cudaMemcpy(d_tmp, d_Mv0, CD, D2D));
    }
    daxpy<N><<<NBlocks, NThreadsPerBlock>>>(d_v, d_op.coeffs[i], d_Mv0, d_v);
  }

  // cudacheck(cudaMemset(d_v, 0, N*CD));
  // mult<N><<<NBlocks, NThreadsPerBlock>>>(d_v, d_tmp, d_valH, d_colsT, d_rowsT);

  // cudacheck(cudaMemset(d_tmp, 0, N*CD));
  // daxpby<N><<<NBlocks, NThreadsPerBlock>>>(d_v, a, d_v, b, d_v0);
}



// https://forums.developer.nvidia.com/t/atomic-add-for-complex-numbers/39757
__device__
void atomicAddCuC(CuC* a, CuC b){
  //transform the addresses of real and imag. parts to double pointers
  double *x = (double*)a;
  double *y = x+1;
  //use atomicAdd for double variables
  atomicAdd(x, real(b));
  atomicAdd(y, imag(b));
}


template<Idx N> __global__
void dot_normalized(CuC* d_res, CuC* d_p, CuC* d_q){
  __shared__ CuC tmp[NThreadsPerBlock];
  if (threadIdx.x == 0) for(int j=0; j<NThreadsPerBlock; j++) tmp[j] = cplx(0.0);
  __syncthreads();

  Idx i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<N){
    tmp[threadIdx.x] = conj(d_p[i])*d_q[i]/N;
    __syncthreads();

    if(threadIdx.x == 0){
      CuC sum = cplx(0.0);
      for(int j=0; j<NThreadsPerBlock; j++) sum = sum+tmp[j];
      atomicAddCuC(d_res, sum);
    }
  }
}


template<Idx N> __host__
void dot_normalized_wrapper(CuC& scalar, CuC* d_scalar, CuC* d_p, CuC* d_q){
  scalar = cplx(0.0);
  cudacheck(cudaMemcpy(d_scalar, &scalar, CD, H2D));
  dot_normalized<N><<<NBlocks, NThreadsPerBlock>>>(d_scalar, d_p, d_q);
  cudacheck(cudaMemcpy(&scalar, d_scalar, CD, D2H));
}


template<Idx N> __host__
void dot2self_normalized_wrapper(double& scalar, CuC* d_scalar, CuC* d_p){
  scalar = 0.0;
  CuC dummy = cplx(scalar);
  cudacheck(cudaMemcpy(d_scalar, &dummy, CD, H2D));
  dot_normalized<N><<<NBlocks, NThreadsPerBlock>>>(d_scalar, d_p, d_p);
  cudacheck(cudaMemcpy(&dummy, d_scalar, CD, D2H));
  assert( abs( imag(dummy) )<1.0e-13*std::sqrt(N) );
  scalar = real(dummy);
}



template<Idx N, Idx len> __host__
void solve(CuC* x, const CuC* b,
	   // const CuC* val,
	   // const CuC* valH,
	   // int *d_cols, int *d_rows, int *d_colsT, int *d_rowsT,
	   const LinOp<CuC>& d_op,
	   const double tol=1.0e-13, const int maxiter=1e8){

  // sparse matrix
  // CuC *d_val, *d_valH;
  // cudacheck(cudaMalloc(&d_val, len*CD));
  // cudacheck(cudaMemcpy(d_val, val, len*CD, H2D));
  // //
  // cudacheck(cudaMalloc(&d_valH, len*CD));
  // cudacheck(cudaMemcpy(d_valH, valH, len*CD, H2D));
  //
  // int *d_cols, *d_rows, *d_colsT, *d_rowsT;
  // cudacheck(cudaMalloc(&d_cols, len*sizeof(int)));
  // cudacheck(cudaMalloc(&d_rows, (N+1)*sizeof(int)));
  // cudacheck(cudaMemcpy(d_cols, cols.data(), len*sizeof(int), H2D));
  // cudacheck(cudaMemcpy(d_rows, rows.data(), (N+1)*sizeof(int), H2D));
  // //
  // cudacheck(cudaMalloc(&d_colsT, len*sizeof(int)));
  // cudacheck(cudaMalloc(&d_rowsT, (N+1)*sizeof(int)));
  // cudacheck(cudaMemcpy(d_colsT, colsT.data(), len*sizeof(int), H2D));
  // cudacheck(cudaMemcpy(d_rowsT, rowsT.data(), (N+1)*sizeof(int), H2D));

  // CG
  CuC *d_x, *d_r, *d_p, *d_q, *d_tmp, *d_tmp2;
  cudacheck(cudaMalloc(&d_x, N*CD));
  cudacheck(cudaMalloc(&d_r, N*CD));
  cudacheck(cudaMalloc(&d_p, N*CD));
  cudacheck(cudaMalloc(&d_q, N*CD));
  cudacheck(cudaMalloc(&d_tmp, N*CD));
  cudacheck(cudaMalloc(&d_tmp2, N*CD));
  //
  cudacheck(cudaMemset(d_x, 0, N*CD)); // added @@
  // cudacheck(cudaMemset(d_r, 0, N*CD)); // added @@
  // cudacheck(cudaMemset(d_p, 0, N*CD)); // added @@
  cudacheck(cudaMemset(d_q, 0, N*CD)); // added @@
  cudacheck(cudaMemset(d_tmp, 0, N*CD)); // added @@
  cudacheck(cudaMemset(d_tmp2, 0, N*CD)); // added @@

  CuC *d_scalar;
  cudacheck(cudaMalloc(&d_scalar, CD));
  cudacheck(cudaMemset(d_scalar, 0, CD)); // added @@

  cudacheck(cudaMemcpy(d_r, b, N*CD, H2D));
  cudacheck(cudaMemcpy(d_p, d_r, N*CD, D2D));

  double mu; dot2self_normalized_wrapper<N>(mu, d_scalar, d_r);
  assert(mu>=0.0);
  double mu_old = mu;

  double b_norm_sq; dot2self_normalized_wrapper<N>(b_norm_sq, d_scalar, d_r);
  assert(b_norm_sq>=0.0);
  double mu_crit = tol*tol*b_norm_sq;

  if(mu<mu_crit) {
#ifdef IsVerbose
    std::clog << "NO SOLVE" << std::endl;
#endif
  }
  else{
    int k=0;
    CuC gam;

    for(; k<maxiter; ++k){
      linop_wrapper<N>(d_q, d_tmp, d_tmp2, d_p,
		       d_op );

      dot_normalized_wrapper<N>(gam, d_scalar, d_p, d_q);

      CuC al = mu/gam;
      cudacheck(cudaMemcpy(d_scalar, &al, CD, H2D));
      daxpy<N><<<NBlocks, NThreadsPerBlock>>>(d_x, d_scalar, d_p, d_x);

      al = -al;
      cudacheck(cudaMemcpy(d_scalar, &al, CD, H2D));
      daxpy<N><<<NBlocks, NThreadsPerBlock>>>(d_r, d_scalar, d_q, d_r);

      dot2self_normalized_wrapper<N>(mu, d_scalar, d_r);
      assert(mu>=0.0);

      if(mu<mu_crit || std::isnan(mu)) break;
      CuC bet = cplx(mu/mu_old);
      mu_old = mu;

      cudacheck(cudaMemcpy(d_scalar, &bet, CD, H2D));
      daxpy<N><<<NBlocks, NThreadsPerBlock>>>(d_p, d_scalar, d_p, d_r);

      if(k%100==0) {
#ifdef IsVerbose
	std::clog << "SOLVER:       #iterations: " << k << ", mu =         " << mu << std::endl;
#endif
      }
    }
#ifdef IsVerbose
    std::clog << "SOLVER:       #iterations: " << k << std::endl;
    std::clog << "SOLVER:       mu =         " << mu << std::endl;
#endif
  }

  cudacheck(cudaMemcpy(x, d_x, N*CD, D2H));

  cudacheck(cudaFree(d_x));
  cudacheck(cudaFree(d_r));
  cudacheck(cudaFree(d_p));
  cudacheck(cudaFree(d_q));
  cudacheck(cudaFree(d_tmp));
  cudacheck(cudaFree(d_tmp2));
  //
  cudacheck(cudaFree(d_scalar));

  // cudacheck(cudaFree(d_val));
  // cudacheck(cudaFree(d_valH));

  // cudacheck(cudaFree(d_cols));
  // cudacheck(cudaFree(d_rows));
  // cudacheck(cudaFree(d_colsT));
  // cudacheck(cudaFree(d_rowsT));
}




// __host__
// template<Idx N, Idx len>
// void solve(CuC* x, const CuC* b,
// 	   const CuC* val,
// 	   const CuC* valH,
// 	   int *d_cols, int *d_rows, int *d_colsT, int *d_rowsT,
// 	   const double aa, const double bb,
// 	   const double tol=1.0e-13, const int maxiter=1e8){

//   // sparse matrix
//   CuC *d_val, *d_valH;
//   cudacheck(cudaMalloc(&d_val, len*CD));
//   cudacheck(cudaMemcpy(d_val, val, len*CD, H2D));
//   //
//   cudacheck(cudaMalloc(&d_valH, len*CD));
//   cudacheck(cudaMemcpy(d_valH, valH, len*CD, H2D));
//   //
//   // int *d_cols, *d_rows, *d_colsT, *d_rowsT;
//   // cudacheck(cudaMalloc(&d_cols, len*sizeof(int)));
//   // cudacheck(cudaMalloc(&d_rows, (N+1)*sizeof(int)));
//   // cudacheck(cudaMemcpy(d_cols, cols.data(), len*sizeof(int), H2D));
//   // cudacheck(cudaMemcpy(d_rows, rows.data(), (N+1)*sizeof(int), H2D));
//   // //
//   // cudacheck(cudaMalloc(&d_colsT, len*sizeof(int)));
//   // cudacheck(cudaMalloc(&d_rowsT, (N+1)*sizeof(int)));
//   // cudacheck(cudaMemcpy(d_colsT, colsT.data(), len*sizeof(int), H2D));
//   // cudacheck(cudaMemcpy(d_rowsT, rowsT.data(), (N+1)*sizeof(int), H2D));

//   // CG
//   CuC *d_x, *d_r, *d_p, *d_q, *d_tmp;
//   cudacheck(cudaMalloc(&d_x, N*CD));
//   cudacheck(cudaMalloc(&d_r, N*CD));
//   cudacheck(cudaMalloc(&d_p, N*CD));
//   cudacheck(cudaMalloc(&d_q, N*CD));
//   cudacheck(cudaMalloc(&d_tmp, N*CD));
//   cudacheck(cudaMemset(d_x, 0, N*CD)); // added @@
//   // cudacheck(cudaMemset(d_r, 0, N*CD)); // added @@
//   // cudacheck(cudaMemset(d_p, 0, N*CD)); // added @@
//   cudacheck(cudaMemset(d_q, 0, N*CD)); // added @@
//   cudacheck(cudaMemset(d_tmp, 0, N*CD)); // added @@

//   CuC *d_scalar;
//   cudacheck(cudaMalloc(&d_scalar, CD));
//   cudacheck(cudaMemset(d_scalar, 0, CD)); // added @@

//   cudacheck(cudaMemcpy(d_r, b, N*CD, H2D));
//   cudacheck(cudaMemcpy(d_p, d_r, N*CD, D2D));

//   double mu; dot2self_normalized_wrapper(mu, d_scalar, d_r, N);
//   assert(mu>=0.0);
//   double mu_old = mu;

//   double b_norm_sq; dot2self_normalized_wrapper(b_norm_sq, d_scalar, d_r, N);
//   assert(b_norm_sq>=0.0);
//   double mu_crit = tol*tol*b_norm_sq;

//   if(mu<mu_crit) {
// #ifdef IsVerbose
//     std::clog << "NO SOLVE" << std::endl;
// #endif
//   }
//   else{
//     int k=0;
//     CuC gam;

//     for(; k<maxiter; ++k){
//       multA(d_q, d_tmp, d_p,
// 	    d_val, d_cols, d_rows,
// 	    d_valH, d_colsT, d_rowsT,
// 	    aa, bb,
// 	    N
// 	    );

//       dot_normalized_wrapper(gam, d_scalar, d_p, d_q, N);

//       CuC al = mu/gam;
//       cudacheck(cudaMemcpy(d_scalar, &al, CD, H2D));
//       daxpy<<<NBlocks, NThreadsPerBlock>>>(d_x, d_scalar, d_p, d_x, N);

//       al = -al;
//       cudacheck(cudaMemcpy(d_scalar, &al, CD, H2D));
//       daxpy<<<NBlocks, NThreadsPerBlock>>>(d_r, d_scalar, d_q, d_r, N);

//       dot2self_normalized_wrapper(mu, d_scalar, d_r, N);
//       assert(mu>=0.0);

//       if(mu<mu_crit || std::isnan(mu)) break;
//       CuC bet = cplx(mu/mu_old);
//       mu_old = mu;

//       cudacheck(cudaMemcpy(d_scalar, &bet, CD, H2D));
//       daxpy<<<NBlocks, NThreadsPerBlock>>>(d_p, d_scalar, d_p, d_r, N);

//       if(k%100==0) {
// #ifdef IsVerbose
// 	std::clog << "SOLVER:       #iterations: " << k << ", mu =         " << mu << std::endl;
// #endif
//       }
//     }
// #ifdef IsVerbose
//     std::clog << "SOLVER:       #iterations: " << k << std::endl;
//     std::clog << "SOLVER:       mu =         " << mu << std::endl;
// #endif
//   }

//   cudacheck(cudaMemcpy(x, d_x, N*CD, D2H));

//   cudacheck(cudaFree(d_x));
//   cudacheck(cudaFree(d_r));
//   cudacheck(cudaFree(d_p));
//   cudacheck(cudaFree(d_q));
//   cudacheck(cudaFree(d_tmp));
//   cudacheck(cudaFree(d_scalar));

//   cudacheck(cudaFree(d_val));
//   cudacheck(cudaFree(d_valH));

//   // cudacheck(cudaFree(d_cols));
//   // cudacheck(cudaFree(d_rows));
//   // cudacheck(cudaFree(d_colsT));
//   // cudacheck(cudaFree(d_rowsT));
// }




// __host__
// template<Idx N, Idx len>
// void multD_wrapper(CuC* x, const CuC* b,
// 		   const CuC* val,
// 		   int *d_cols, int *d_rows,
// 		   const double aa, const double bb
// 		   ){
//   // sparse matrix
//   CuC *d_val;
//   cudacheck(cudaMalloc(&d_val, len*CD));
//   cudacheck(cudaMemcpy(d_val, val, len*CD, H2D));

//   CuC *d_x, *d_b;
//   cudacheck(cudaMalloc(&d_x, N*CD));
//   cudacheck(cudaMalloc(&d_b, N*CD));
//   cudacheck(cudaMemset(d_x, 0, N*CD)); // added @@
//   cudacheck(cudaMemset(d_b, 0, N*CD)); // added @@

//   cudacheck(cudaMemcpy(d_b, b, N*CD, H2D));

//   multD_wrapper(d_x, d_b,
// 		d_val, d_cols, d_rows,
// 		aa, bb,
// 		N, len
// 		);

//   cudacheck(cudaMemcpy(x, d_x, N*CD, D2H));

//   cudacheck(cudaFree(d_x));
//   cudacheck(cudaFree(d_b));
//   // cudacheck(cudaFree(d_tmp));

//   cudacheck(cudaFree(d_val));
//   // cudacheck(cudaFree(d_valH));
// }



// __host__
// template<Idx N, Idx len>
// void multA_wrapper(CuC* x, const CuC* b,
// 		   const CuC* val, const CuC* valH,
// 		   int *d_cols, int *d_rows,
// 		   int *d_colsT, int *d_rowsT,
// 		   const double aa, const double bb
// 		   ){
//   // sparse matrix
//   CuC *d_val, *d_valH;
//   cudacheck(cudaMalloc(&d_val, len*CD));
//   cudacheck(cudaMemcpy(d_val, val, len*CD, H2D));
//   //
//   cudacheck(cudaMalloc(&d_valH, len*CD));
//   cudacheck(cudaMemcpy(d_valH, valH, len*CD, H2D));

//   CuC *d_x, *d_b, *d_tmp;
//   cudacheck(cudaMalloc(&d_x, N*CD));
//   cudacheck(cudaMalloc(&d_b, N*CD));
//   cudacheck(cudaMalloc(&d_tmp, N*CD));
//   cudacheck(cudaMemset(d_x, 0, N*CD)); // added @@
//   cudacheck(cudaMemset(d_b, 0, N*CD)); // added @@
//   cudacheck(cudaMemset(d_tmp, 0, N*CD)); // added @@

//   cudacheck(cudaMemcpy(d_b, b, N*CD, H2D));

//   multA(d_x, d_tmp, d_b,
// 	d_val, d_cols, d_rows,
// 	d_valH, d_colsT, d_rowsT,
// 	aa, bb,
// 	N
// 	);

//   cudacheck(cudaMemcpy(x, d_x, N*CD, D2H));

//   cudacheck(cudaFree(d_x));
//   cudacheck(cudaFree(d_b));
//   cudacheck(cudaFree(d_tmp));

//   cudacheck(cudaFree(d_val));
//   cudacheck(cudaFree(d_valH));
// }







