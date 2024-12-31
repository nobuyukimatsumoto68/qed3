#pragma once


#include <cuComplex.h>
#include <cuda_runtime.h>



#define NThreadsPerBlock (256) // 256, 1024
#define NBlocks (N+NThreadsPerBlock)/NThreadsPerBlock
#define H2D (cudaMemcpyHostToDevice)
#define D2H (cudaMemcpyDeviceToHost)
#define D2D (cudaMemcpyDeviceToDevice)
#define DB (sizeof(double))
#define CD (sizeof(cuDoubleComplex))
#define cuI ( make_cuDoubleComplex(0.0,1.0) )




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



template<typename T, Idx N> __global__
void mult ( T* res,
	    const T* v,
	    const T* v_csr,
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


template<typename T, Idx N> __global__
void Taxpy(T* d_res, const T* d_a, const T* d_x, const T* d_y){
  Idx i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<N) d_res[i] = *d_a * d_x[i] + d_y[i];
}


template<typename T, Idx N> __global__
void ScalarMult(T* d_res, const T* d_a, const T* d_x){
  Idx i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<N) d_res[i] = *d_a * d_x[i];
}
