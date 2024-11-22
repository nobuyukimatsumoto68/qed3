#pragma once

#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>

#define H2D (cudaMemcpyHostToDevice)
#define D2H (cudaMemcpyDeviceToHost)
#define D2D (cudaMemcpyDeviceToDevice)
#define DB (sizeof(double))
#define CD (sizeof(Complex))
#define cuI ( make_cuDoubleComplex(0.0,1.0) )

using complex = cuDoubleComplex;

__device__ __host__
complex cplx(const double c) { return make_cuDoubleComplex(c, 0.0); }
complex cplx(const std::complex<double> c) { return make_cuDoubleComplex(c.real(), c.imag()); }


__device__ __host__ complex  operator+(complex a, complex b) { return cuCadd(a,b); }
__device__ __host__ complex  operator+(complex a, double b) { return cuCadd(a,cplx(b)); }
__device__ __host__ complex  operator+(double a, complex b) { return cuCadd(cplx(a),b); }
__device__ __host__ complex  operator-(complex a, complex b) { return cuCsub(a,b); }
__device__ __host__ complex  operator-(complex a, double b) { return cuCsub(a,cplx(b)); }
__device__ __host__ complex  operator-(double a, complex b) { return cuCsub(cplx(a),b); }
__device__ __host__ complex  operator-(complex b) { return cplx(0.0)-b; }
__device__ __host__ complex  operator*(complex a, complex b) { return cuCmul(a,b); }
__device__ __host__ complex  operator*(complex a, double b) { return cuCmul(a,cplx(b)); }
__device__ __host__ complex  operator*(double a, complex b) { return cuCmul(cplx(a),b); }
__device__ __host__ complex  operator/(complex a, complex b) { return cuCdiv(a,b); }
__device__ __host__ complex  operator/(complex a, double b) { return cuCdiv(a,cplx(b)); }
__device__ __host__ complex  operator/(double a, complex b) { return cuCdiv(cplx(a),b); }


__device__ __host__ inline double real(const complex c ){ return cuCreal(c); }
__device__ __host__ inline double imag(const complex c ){ return cuCimag(c); }
__device__ __host__ inline complex conj(complex c ){ return cuConj(c); }


__host__
void cudacheck( cudaError status ){
  if(status!=0) std::cout << status << std::endl;
  assert(cudaSuccess == status);
}

__host__
void cudacheck( cusolverStatus_t status ){
  if(status!=CUSOLVER_STATUS_SUCCESS) std::cout << status << std::endl;
  assert(CUSOLVER_STATUS_SUCCESS == status);
}
