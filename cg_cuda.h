#pragma once

#include <cuComplex.h>
#include <cuda_runtime.h>


// ======================================


#define NThreadsPerBlock (128) // 256, 1024
#define NBlocks (N+NThreadsPerBlock)/NThreadsPerBlock
#define H2D (cudaMemcpyHostToDevice)
#define D2H (cudaMemcpyDeviceToHost)
#define D2D (cudaMemcpyDeviceToDevice)
#define DB (sizeof(double))
#define CD (sizeof(cuDoubleComplex))
#define cuI ( make_cuDoubleComplex(0.0,1.0) )


__device__ __host__ cuDoubleComplex cplx(const double c) { return make_cuDoubleComplex(c, 0.0); }
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

using CuC = cuDoubleComplex;
using Idx = long int;

__host__
void cudacheck( cudaError status ){
  if(status!=0) std::cout << status << std::endl;
  assert(cudaSuccess == status);
}
__host__
void set2zero( double* v, const Idx size ){ for(Idx i=0; i<size; i++) v[i] = 0.0; }
__host__
void set2zero( CuC* v, const Idx size ){ for(Idx i=0; i<size; i++) v[i] = cplx(0.0); }
__host__
void set2zero( unsigned long* v, const Idx size ){ for(Idx i=0; i<size; i++) v[i] = 0; }
__host__ __device__
int mod(const int a, const int b){ return (b +(a%b))%b; }
__global__
void daxpy(CuC* d_res, CuC* d_a, CuC* d_x, CuC* d_y, const int N){
  Idx i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<N) d_res[i] = *d_a * d_x[i] + d_y[i];
}



__global__
void mult ( CuC* res,
	    const CuC* v,
	    const CuC* v_csr,
	    const int* cols,
	    const int* rows,
	    const int N
	    ){
  Idx i = blockIdx.x*blockDim.x + threadIdx.x;

  if(i<N) {
    res[i] = cplx(0.0);
    const int row_start = rows[i];
    const int row_end = rows[i+1];
    for(int jj=row_start; jj<row_end; jj++) res[i] = res[i] + v_csr[jj] * v[ cols[jj] ];
  }
}


// __global__
// void mult_gam5_M( CuC* res,
// 		  const CuC* v,
// 		  const CuC* v_csr,
// 		  const int* cols,
// 		  const int* rows,
// 		  const int N
// 		  ){
//   Idx i = blockIdx.x*blockDim.x + threadIdx.x;

//   if(i<N) {
//     res[i] = cplx(0.0);
//     const int row_start = rows[i];
//     const int row_end = rows[i+1];
//     for(int jj=row_start; jj<row_end; jj++) res[i] = res[i] + v_csr[jj] * v[ cols[jj] ];
//     res[i] = ( -2*(i%2) + 1 ) * res[i];
//   }
// }


// __global__
// void multgam5_self ( CuC* v,
// 		     const int N
// 		){
//   Idx i = blockIdx.x*blockDim.x + threadIdx.x;
//   if(i<N) v[i] = ( -2*(i%2) + 1 ) * v[i];
// }


__host__
void multA(CuC* d_v, CuC* d_tmp, CuC* d_v0,
	   CuC* d_val, int* d_cols, int* d_rows,
	   CuC* d_valH, int* d_colsT, int* d_rowsT,
	   const int N
	   ){
  cudacheck(cudaMemset(d_tmp, 0, N*CD));
  mult<<<NBlocks, NThreadsPerBlock>>>(d_tmp, d_v0, d_val, d_cols, d_rows, N);

  cudacheck(cudaMemset(d_v, 0, N*CD));
  mult<<<NBlocks, NThreadsPerBlock>>>(d_v, d_tmp, d_valH, d_colsT, d_rowsT, N);
}


// __host__
// void multA(CuC* d_v, CuC* d_tmp, CuC* d_v0,
// 	   CuC* d_val, int* d_cols, int* d_rows,
// 	   // CuC* d_valH, int* d_colsT, int* d_rowsT,
// 	   const int N
// 	   ){
//   cudacheck(cudaMemset(d_tmp, 0, N*CD));
//   mult<<<NBlocks, NThreadsPerBlock>>>(d_tmp, d_v0, d_val, d_cols, d_rows, N);

//   cudacheck(cudaMemset(d_v, 0, N*CD));
//   // mult_gam5_M<<<NBlocks, NThreadsPerBlock>>>(d_tmp, d_v0, d_val, d_cols, d_rows, N);
//   mult<<<NBlocks, NThreadsPerBlock>>>(d_v, d_tmp, d_valH, d_colsT, d_rowsT, N);
// }




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
__global__
void dot_normalized(CuC* d_res, CuC* d_p, CuC* d_q, const int N){
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
__host__
void dot_normalized_wrapper(CuC& scalar, CuC* d_scalar, CuC* d_p, CuC* d_q, const int N){
  scalar = cplx(0.0);
  cudacheck(cudaMemcpy(d_scalar, &scalar, CD, H2D));
  dot_normalized<<<NBlocks, NThreadsPerBlock>>>(d_scalar, d_p, d_q, N);
  cudacheck(cudaMemcpy(&scalar, d_scalar, CD, D2H));
}
__host__
void dot2self_normalized_wrapper(double& scalar, CuC* d_scalar, CuC* d_p, const int N){
  scalar = 0.0;
  CuC dummy = cplx(scalar);
  cudacheck(cudaMemcpy(d_scalar, &dummy, CD, H2D));
  dot_normalized<<<NBlocks, NThreadsPerBlock>>>(d_scalar, d_p, d_p, N);
  cudacheck(cudaMemcpy(&dummy, d_scalar, CD, D2H));
  std::cout << abs( imag(dummy) ) << std::endl;
  assert( abs( imag(dummy) )<1.0e-13*std::sqrt(N) );
  scalar = real(dummy);
}



__host__
void solve(CuC* x, const CuC* b,
	   const CuC* val, const std::vector<int>& cols, const std::vector<int>& rows,
	   const CuC* valH, const std::vector<int>& colsT, const std::vector<int>& rowsT,
	   const int N, const int len,
	   const double tol=1.0e-13, const int maxiter=1e8){

  // sparse matrix
  CuC *d_val, *d_valH;
  cudacheck(cudaMalloc(&d_val, len*CD));
  cudacheck(cudaMemcpy(d_val, val, len*CD, H2D));
  //
  cudacheck(cudaMalloc(&d_valH, len*CD));
  cudacheck(cudaMemcpy(d_valH, valH, len*CD, H2D));
  //
  int *d_cols, *d_rows, *d_colsT, *d_rowsT;
  cudacheck(cudaMalloc(&d_cols, len*sizeof(int)));
  cudacheck(cudaMalloc(&d_rows, (N+1)*sizeof(int)));
  cudacheck(cudaMemcpy(d_cols, cols.data(), len*sizeof(int), H2D));
  cudacheck(cudaMemcpy(d_rows, rows.data(), (N+1)*sizeof(int), H2D));
  //
  cudacheck(cudaMalloc(&d_colsT, len*sizeof(int)));
  cudacheck(cudaMalloc(&d_rowsT, (N+1)*sizeof(int)));
  cudacheck(cudaMemcpy(d_colsT, colsT.data(), len*sizeof(int), H2D));
  cudacheck(cudaMemcpy(d_rowsT, rowsT.data(), (N+1)*sizeof(int), H2D));

  // CG
  CuC *d_x, *d_r, *d_p, *d_q, *d_tmp;
  cudacheck(cudaMalloc(&d_x, N*CD));
  cudacheck(cudaMalloc(&d_r, N*CD));
  cudacheck(cudaMalloc(&d_p, N*CD));
  cudacheck(cudaMalloc(&d_q, N*CD));
  cudacheck(cudaMalloc(&d_tmp, N*CD));
  cudacheck(cudaMemset(d_x, 0, N*CD)); // added @@
  cudacheck(cudaMemset(d_r, 0, N*CD)); // added @@
  cudacheck(cudaMemset(d_p, 0, N*CD)); // added @@
  cudacheck(cudaMemset(d_q, 0, N*CD)); // added @@
  cudacheck(cudaMemset(d_tmp, 0, N*CD)); // added @@

  CuC *d_scalar;
  cudacheck(cudaMalloc(&d_scalar, CD));
  cudacheck(cudaMemset(d_scalar, 0, CD)); // added @@

  cudacheck(cudaMemcpy(d_r, b, N*CD, H2D));
  cudacheck(cudaMemcpy(d_p, d_r, N*CD, D2D));

  double mu; dot2self_normalized_wrapper(mu, d_scalar, d_r, N);
  assert(mu>=0.0);
  double mu_old = mu;

  double b_norm_sq; dot2self_normalized_wrapper(b_norm_sq, d_scalar, d_r, N);
  assert(b_norm_sq>=0.0);
  double mu_crit = tol*tol*b_norm_sq;

  if(mu<mu_crit) std::clog << "NO SOLVE" << std::endl;
  else{
    int k=0;
    CuC gam;

    for(; k<maxiter; ++k){
      // multA(d_q, d_tmp, d_p, nu);
      multA(d_q, d_tmp, d_p,
	    d_val, d_cols, d_rows,
	    d_valH, d_colsT, d_rowsT,
	    N
	    );

      dot_normalized_wrapper(gam, d_scalar, d_p, d_q, N);

      CuC al = mu/gam;
      cudacheck(cudaMemcpy(d_scalar, &al, CD, H2D));
      daxpy<<<NBlocks, NThreadsPerBlock>>>(d_x, d_scalar, d_p, d_x, N);

      al = -al;
      cudacheck(cudaMemcpy(d_scalar, &al, CD, H2D));
      daxpy<<<NBlocks, NThreadsPerBlock>>>(d_r, d_scalar, d_q, d_r, N);

      dot2self_normalized_wrapper(mu, d_scalar, d_r, N);
      assert(mu>=0.0);

      if(mu<mu_crit || std::isnan(mu)) break;
      CuC bet = cplx(mu/mu_old);
      mu_old = mu;

      cudacheck(cudaMemcpy(d_scalar, &bet, CD, H2D));
      daxpy<<<NBlocks, NThreadsPerBlock>>>(d_p, d_scalar, d_p, d_r, N);

      if(k%100==0) {
	std::clog << "SOLVER:       #iterations: " << k << ", mu =         " << mu << std::endl;
      }
    }
    std::clog << "SOLVER:       #iterations: " << k << std::endl;
    std::clog << "SOLVER:       mu =         " << mu << std::endl;
  }

  cudacheck(cudaMemcpy(x, d_x, N*CD, D2H));

  cudacheck(cudaFree(d_x));
  cudacheck(cudaFree(d_r));
  cudacheck(cudaFree(d_p));
  cudacheck(cudaFree(d_q));
  cudacheck(cudaFree(d_tmp));
  cudacheck(cudaFree(d_scalar));

  cudacheck(cudaFree(d_val));
  cudacheck(cudaFree(d_cols));
  cudacheck(cudaFree(d_rows));
  cudacheck(cudaFree(d_valH));
  cudacheck(cudaFree(d_colsT));
  cudacheck(cudaFree(d_rowsT));
}


// __host__
// void solve(CuC* x, const CuC* b,
// 	   const CuC* val, const std::vector<int>& cols, const std::vector<int>& rows,
// 	   // const CuC* valH, const std::vector<int>& colsT, const std::vector<int>& rowsT,
// 	   const int N, const int len,
// 	   const double tol=1.0e-13, const int maxiter=1e8){

//   // sparse matrix
//   CuC *d_val, *d_valH;
//   cudacheck(cudaMalloc(&d_val, len*CD));
//   cudacheck(cudaMemcpy(d_val, val, len*CD, H2D));
//   //
//   // cudacheck(cudaMalloc(&d_valH, len*CD));
//   // cudacheck(cudaMemcpy(d_valH, valH, len*CD, H2D));
//   //
//   int *d_cols, *d_rows; // , *d_colsT, *d_rowsT;
//   cudacheck(cudaMalloc(&d_cols, len*sizeof(int)));
//   cudacheck(cudaMalloc(&d_rows, (N+1)*sizeof(int)));
//   cudacheck(cudaMemcpy(d_cols, cols.data(), len*sizeof(int), H2D));
//   cudacheck(cudaMemcpy(d_rows, rows.data(), (N+1)*sizeof(int), H2D));
//   //
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
//   cudacheck(cudaMemset(d_r, 0, N*CD)); // added @@
//   cudacheck(cudaMemset(d_p, 0, N*CD)); // added @@
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

//   if(mu<mu_crit) std::clog << "NO SOLVE" << std::endl;
//   else{
//     int k=0;
//     CuC gam;

//     for(; k<maxiter; ++k){
//       // multA(d_q, d_tmp, d_p, nu);
//       multA(d_q, d_tmp, d_p,
// 	    d_val, d_cols, d_rows,
// 	    // d_valH, d_colsT, d_rowsT,
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
// 	std::clog << "SOLVER:       #iterations: " << k << ", mu =         " << mu << std::endl;
//       }
//     }
//     std::clog << "SOLVER:       #iterations: " << k << std::endl;
//     std::clog << "SOLVER:       mu =         " << mu << std::endl;
//   }

//   cudacheck(cudaMemcpy(x, d_x, N*CD, D2H));

//   cudacheck(cudaFree(d_x));
//   cudacheck(cudaFree(d_r));
//   cudacheck(cudaFree(d_p));
//   cudacheck(cudaFree(d_q));
//   cudacheck(cudaFree(d_tmp));
//   cudacheck(cudaFree(d_scalar));

//   cudacheck(cudaFree(d_val));
//   cudacheck(cudaFree(d_cols));
//   cudacheck(cudaFree(d_rows));
//   // cudacheck(cudaFree(d_valH));
//   // cudacheck(cudaFree(d_colsT));
//   // cudacheck(cudaFree(d_rowsT));
// }




struct CGCUDA{ // wrapper
  using Complex=std::complex<double>;

  const Sparse sparse;
  const Dirac1fonS2& D;

  CGCUDA( const Lattice& lattice,
	  const Dirac1fonS2& D )
    : sparse(lattice)
    , D(D)
  {}

  void operator()( Complex* res, const Complex* v, const U1onS2& U ) const {

    Complex v_coo[ sparse.len ], v_csr[ sparse.len ], v_csrH[ sparse.len ];
    D.coo_format( v_coo, sparse.N, U );
    sparse.coo2csr_csrH( v_csr, v_csrH, v_coo );

    solve(reinterpret_cast<CuC*>(res),
	  reinterpret_cast<const CuC*>(v),
	  reinterpret_cast<const CuC*>(v_csr), sparse.cols_csr, sparse.rows_csr,
	  reinterpret_cast<const CuC*>(v_csrH), sparse.cols_csrT, sparse.rows_csrT,
	  sparse.N, sparse.len
	  );
  }

};




// struct CGCUDA{
//   const int N;
//   const int len;

//   CuC *d_val, *d_valH;
//   int *d_cols, *d_rows, *d_colsT, *d_rowsT;
//   CuC *d_x, *d_r, *d_p, *d_q, *d_tmp;
//   CuC *d_scalar;


//   CGCUDA( const Sparse& sparse )
//     : N(sparse.N)
//     , len(sparse.len)
//   {
//     int device;
//     cudacheck(cudaGetDeviceCount(&device));
//     cudaDeviceProp device_prop[device];
//     cudaGetDeviceProperties(&device_prop[0], 0);
//     cudacheck(cudaSetDevice(0));// "TITAN V"

//     // sparse matrix
//     cudacheck(cudaMalloc(&d_val, len*CD));
//     cudacheck(cudaMalloc(&d_valH, len*CD));
//     //
//     cudacheck(cudaMalloc(&d_cols, len*sizeof(int)));
//     cudacheck(cudaMalloc(&d_rows, (N+1)*sizeof(int)));
//     cudacheck(cudaMalloc(&d_colsT, len*sizeof(int)));
//     cudacheck(cudaMalloc(&d_rowsT, (N+1)*sizeof(int)));

//     cudacheck(cudaMemcpy(d_cols, sparse.cols_csr.data(), len*sizeof(int), H2D));
//     cudacheck(cudaMemcpy(d_rows, sparse.rows_csr.data(), (N+1)*sizeof(int), H2D));
//     cudacheck(cudaMemcpy(d_colsT, sparse.cols_csrT.data(), len*sizeof(int), H2D));
//     cudacheck(cudaMemcpy(d_rowsT, sparse.rows_csrT.data(), (N+1)*sizeof(int), H2D));

//     // ------------------

//     // CG
//     cudacheck(cudaMalloc(&d_x, N*CD));
//     cudacheck(cudaMalloc(&d_r, N*CD));
//     cudacheck(cudaMalloc(&d_p, N*CD));
//     cudacheck(cudaMalloc(&d_q, N*CD));
//     cudacheck(cudaMalloc(&d_tmp, N*CD));

//     cudacheck(cudaMalloc(&d_scalar, CD));

//   }

//   ~CGCUDA(){
//     cudacheck(cudaFree(d_x));
//     cudacheck(cudaFree(d_r));
//     cudacheck(cudaFree(d_p));
//     cudacheck(cudaFree(d_q));
//     cudacheck(cudaFree(d_tmp));
//     cudacheck(cudaFree(d_scalar));

//     cudacheck(cudaFree(d_val));
//     cudacheck(cudaFree(d_valH));
//     cudacheck(cudaFree(d_cols));
//     cudacheck(cudaFree(d_rows));
//     cudacheck(cudaFree(d_colsT));
//     cudacheck(cudaFree(d_rowsT));

//     cudacheck(cudaDeviceReset());
//   }
  


//   __host__
//   void solve(CuC* x, CuC* b,
// 	     CuC* val, CuC* valH,
// 	     const double tol=1.0e-13, const int maxiter=1e8){
//     // CG
//     cudacheck(cudaMemset(d_x, 0, N*CD)); // added @@
//     cudacheck(cudaMemset(d_r, 0, N*CD)); // added @@
//     cudacheck(cudaMemset(d_p, 0, N*CD)); // added @@
//     cudacheck(cudaMemset(d_q, 0, N*CD)); // added @@
//     cudacheck(cudaMemset(d_tmp, 0, N*CD)); // added @@
//     cudacheck(cudaMemset(d_scalar, 0, CD)); // added @@

//     cudacheck(cudaMemcpy(d_r, b, N*CD, H2D));
//     cudacheck(cudaMemcpy(d_p, d_r, N*CD, D2D));

//     double mu; dot2self_normalized_wrapper(mu, d_scalar, d_r, N);
//     assert(mu>=0.0);
//     double mu_old = mu;

//     double b_norm_sq; dot2self_normalized_wrapper(b_norm_sq, d_scalar, d_r, N);
//     assert(b_norm_sq>=0.0);
//     double mu_crit = tol*tol*b_norm_sq;

//     if(mu<mu_crit) std::clog << "NO SOLVE" << std::endl;
//     else{
//       int k=0;
//       CuC gam;

//       for(; k<maxiter; ++k){
// 	multA(d_q, d_tmp, d_p,
// 	      d_val, d_cols, d_rows,
// 	      d_valH, d_colsT, d_rowsT,
// 	      N
// 	      );

// 	dot_normalized_wrapper(gam, d_scalar, d_p, d_q, N);

// 	CuC al = mu/gam;
// 	cudacheck(cudaMemcpy(d_scalar, &al, CD, H2D));
// 	daxpy<<<NBlocks, NThreadsPerBlock>>>(d_x, d_scalar, d_p, d_x, N);

// 	al = -al;
// 	cudacheck(cudaMemcpy(d_scalar, &al, CD, H2D));
// 	daxpy<<<NBlocks, NThreadsPerBlock>>>(d_r, d_scalar, d_q, d_r, N);

// 	dot2self_normalized_wrapper(mu, d_scalar, d_r, N);
// 	assert(mu>=0.0);

// 	if(mu<mu_crit || std::isnan(mu)) break;
// 	CuC bet = cplx(mu/mu_old);
// 	mu_old = mu;

// 	cudacheck(cudaMemcpy(d_scalar, &bet, CD, H2D));
// 	daxpy<<<NBlocks, NThreadsPerBlock>>>(d_p, d_scalar, d_p, d_r, N);

// 	if(k%100==0) {
// 	  std::clog << "SOLVER:       #iterations: " << k << ", mu =         " << mu << std::endl;
// 	}
//       }
//       std::clog << "SOLVER:       #iterations: " << k << std::endl;
//       std::clog << "SOLVER:       mu =         " << mu << std::endl;
//     }

//     cudacheck(cudaMemcpy(x, d_x, N*CD, D2H));
//   }

// };
