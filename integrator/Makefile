CXX=g++
CXXFLAGS=-Wall -Wextra -g
INCLUDES=-I"/opt/eigen/" -I"/opt/gsl-2.8/"
LDFLAGS=-L"/usr/local/lib" -lgsl -lgslcblas -lm

NVCC=/usr/local/cuda-12.6/bin/nvcc # nvcc
NVCCFLAGS = -arch=sm_70 -O3 -lcusolver -std=c++17 # -diag-suppress<1650-D>
INCLUDES_CUDA =
LDFLAGS_CUDA = -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/

INCLUDES_CUDA += $(INCLUDES)


GRP_DIR=$(shell pwd)/newQFE/grp/
CXXFLAGS+=-DGRP_DIR="\"$(GRP_DIR)\""
NVCCFLAGS+=-DGRP_DIR="\"$(GRP_DIR)\""


all:
	$(CXX) test_integrator.cc $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)
