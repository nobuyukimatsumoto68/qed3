#pragma once



#define NThreadsPerBlock (512) // 256, 1024
#define NBlocks (N+NThreadsPerBlock)/NThreadsPerBlock
#define H2D (cudaMemcpyHostToDevice)
#define D2H (cudaMemcpyDeviceToHost)
#define D2D (cudaMemcpyDeviceToDevice)
#define DB (sizeof(double))
#define CD (sizeof(cuDoubleComplex))
#define cuI ( make_cuDoubleComplex(0.0,1.0) )




__device__ __host__ cuDoubleComplex cplx(const double c) { return make_cuDoubleComplex(c, 0.0); }
// __device__ __host__ cuDoubleComplex cplx(const std::complex<double> c) { return make_cuDoubleComplex(c.real(), c.imag()); }
__host__ cuDoubleComplex cplx(const std::complex<double> c) { return make_cuDoubleComplex(c.real(), c.imag()); }
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
__device__ __host__ inline double abs(const cuDoubleComplex c ){ return cuCabs(c); }
__device__ __host__ inline cuDoubleComplex conj(cuDoubleComplex c ){ return cuConj(c); }



// __host__
// void CUDA_CHECK( cudaError status ){
//   if(status!=0) std::cout << status << std::endl;
//   assert(cudaSuccess == status);
// }


// __host__
// void cublascheck( cublasStatus_t status ){
//   if(status!=0) std::cout << status << std::endl;
//   assert(cudaSuccess == status);
// }
#define CUDA_CHECK(err)                                                                            \
    do {                                                                                           \
        cudaError_t err_ = (err);                                                                  \
        if (err_ != cudaSuccess) {                                                                 \
            std::printf("CUDA error %d at %s:%d\n", err_, __FILE__, __LINE__);                     \
            throw std::runtime_error("CUDA error");                                                \
        }                                                                                          \
    } while (0)

// cublas API error checking
#define CUBLAS_CHECK(err)                                                                          \
    do {                                                                                           \
        cublasStatus_t err_ = (err);                                                               \
        if (err_ != CUBLAS_STATUS_SUCCESS) {                                                       \
            std::printf("cublas error %d at %s:%d\n", err_, __FILE__, __LINE__);                   \
            throw std::runtime_error("cublas error");                                              \
        }                                                                                          \
    } while (0)

// cusolver API error checking
#define CUSOLVER_CHECK(err)                                                                        \
    do {                                                                                           \
        cusolverStatus_t err_ = (err);                                                             \
        if (err_ != CUSOLVER_STATUS_SUCCESS) {                                                     \
            printf("cusolver error %d at %s:%d\n", err_, __FILE__, __LINE__);                      \
            throw std::runtime_error("cusolver error");                                            \
        }                                                                                          \
    } while (0)

#define CUSPARSE_CHECK(err)                                                                        \
    do {                                                                                           \
        cusparseStatus_t err_ = (err);                                                             \
        if (err_ != CUSPARSE_STATUS_SUCCESS) {                                                     \
            printf("cusparse error %d at %s:%d\n", err_, __FILE__, __LINE__);                      \
            throw std::runtime_error("cusparse error");                                            \
        }                                                                                          \
    } while (0)


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
void Taxpy(T* d_res, const T d_a, const T* d_x, const T* d_y){
  Idx i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<N) d_res[i] = d_a * d_x[i] + d_y[i];
}


template<typename T, typename T2, Idx N> __global__
void Taxpy_gen(T* d_res, const T2* d_a, const T* d_x, const T* d_y){
  Idx i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<N) d_res[i] = *d_a * d_x[i] + d_y[i];
}

template<typename T, typename T2, Idx N> __global__
void Taxpy_gen(T* d_res, const T2 d_a, const T* d_x, const T* d_y){
  Idx i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<N) d_res[i] = d_a * d_x[i] + d_y[i];
}


// template<typename T, Idx N> __global__
// void ScalarMult(T* d_res, const T* d_a, const T* d_x){
//   Idx i = blockIdx.x*blockDim.x + threadIdx.x;
//   if(i<N) d_res[i] = *d_a * d_x[i];
// }







// https://forums.developer.nvidia.com/t/atomic-add-for-complex-numbers/39757
__device__
void atomicAddCuC(CuC* a, const CuC b){
  //transform the addresses of real and imag. parts to double pointers
  double *x = reinterpret_cast<double*>(a);
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

  int i = blockIdx.x*blockDim.x + threadIdx.x;
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

// template<Idx N> __global__
// void dot_normalized(CuC* a, CuC* b, CuC* c) {
//   __shared__ CuC cache[threadsPerBlock];
//   int tid = threadIdx.x + blockIdx.x * blockDim.x;
//   int cacheIndex = threadIdx.x;

//   CuC temp = cplx(0.);
//   while (tid < N){
//     temp += conj(a[tid]) * b[tid] / N;
//     tid += blockDim.x * gridDim.x;
//   }
	
//   // set the cache values
//   cache[cacheIndex] = temp;
	
//   // synchronize threads in this block
//   __syncthreads();
	
//   // for reductions, threadsPerBlock must be a power of 2
//   // because of the following code
//   int i = blockDim.x/2;
//   while (i != 0){
//     if (cacheIndex < i) cache[cacheIndex] += cache[cacheIndex + i];
//     __syncthreads();
//     i /= 2;
//   }
	
//   if (cacheIndex == 0) c[blockIdx.x] = cache[0];
// }


template<Idx N> __host__
void dot_normalized_wrapper(CuC& scalar, CuC* d_scalar, CuC* d_p, CuC* d_q){
  scalar = cplx(0.0);
  CUDA_CHECK(cudaMemcpy(d_scalar, &scalar, CD, H2D));
  dot_normalized<N><<<NBlocks, NThreadsPerBlock>>>(d_scalar, d_p, d_q);
  CUDA_CHECK(cudaMemcpy(&scalar, d_scalar, CD, D2H));
}


template<Idx N> __host__
void dot2self_normalized_wrapper(double& scalar, CuC* d_scalar, CuC* d_p, const double TOL=1.0e-12){
  scalar = 0.0;
  CuC dummy = cplx(scalar);
  CUDA_CHECK(cudaMemcpy(d_scalar, &dummy, CD, H2D));
  dot_normalized<N><<<NBlocks, NThreadsPerBlock>>>(d_scalar, d_p, d_p);
  CUDA_CHECK(cudaMemcpy(&dummy, d_scalar, CD, D2H));
  assert( abs( imag(dummy) )<TOL*std::sqrt(N) );
  scalar = real(dummy);
}


