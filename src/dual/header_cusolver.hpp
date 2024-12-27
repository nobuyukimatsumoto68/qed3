#pragma once

// #include "header_cuda.hpp"
#include <cusolverDn.h>


__host__
void cudacheck( cusolverStatus_t status ){
  if(status!=CUSOLVER_STATUS_SUCCESS) std::cout << status << std::endl;
  assert(CUSOLVER_STATUS_SUCCESS == status);
}
