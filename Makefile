
CXX=g++
CXXFLAGS=-std=c++17 -g -O3
INCLUDES=-I"./newQFE/include"

NVCC=/usr/local/cuda-12.6/bin/nvcc # nvcc
NVCCFLAGS = -arch=sm_70 -O3 -lcusolver -std=c++17 # -diag-suppress<1650-D>
INCLUDES_CUDA =
LDFLAGS_CUDA = -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/

INCLUDES_CUDA += $(INCLUDES)


GRP_DIR=$(shell pwd)/newQFE/grp/
CXXFLAGS+=-DGRP_DIR="\"$(GRP_DIR)\""
NVCCFLAGS+=-DGRP_DIR="\"$(GRP_DIR)\""


all:
	# $(CXX) get_geometry.cc $(INCLUDES) $(CXXFLAGS)
	# $(CXX) test_geometry.cc $(INCLUDES) $(CXXFLAGS)
	# $(CXX) test_geometry_flat.cc $(INCLUDES) $(CXXFLAGS)
	# $(CXX) test_u1.cc $(INCLUDES) $(CXXFLAGS)
	$(CXX) test_dirac.cc $(INCLUDES) $(CXXFLAGS)
	# $(CXX) test_dirac_flat.cc $(INCLUDES) $(CXXFLAGS)

eig: test_dirac.cu header_cusolver.hpp
	$(NVCC) $< $(NVCCFLAGS) $(INCLUDES_CUDA) -o b.out

geometry:
	$(CXX) get_geometry.cc $(INCLUDES) $(CXXFLAGS)
