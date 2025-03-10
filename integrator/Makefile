CXX=g++
# CXX=clang++
CXXFLAGS=-Wall -Wextra -std=c++23 -fuse-ld=lld -march=native -O3 -g
INCLUDES=-I"/usr/include/c++/13" -I"/usr/include/x86_64-linux-gnu/c++/13/" -I"/opt/eigen-3.4.0/" -I"/opt/gsl-2.8/"
# INCLUDES=-I"/opt/eigen-3.4.0/" -I"/opt/gsl-2.8/"

QFEINCLUDES=-I"../newQFE/include"

GRP_DIR="/mnt/hdd_barracuda/qed3/newQFE/grp/"
QFECXXFLAGS=-DGRP_DIR="\"$(GRP_DIR)\""

LDFLAGS=-L"/lib/x86_64-linux-gnu/" -L"/lib/x86_64-linux-gnu/" -L"/lib64/" -L"/usr/local/lib" -lgsl -lgslcblas -lm

NVCC=/usr/local/cuda-12.6/bin/nvcc # nvcc
NVCCFLAGS = -arch=sm_70 -O3 -lcusolver -std=c++17 # -diag-suppress<1650-D>
INCLUDES_CUDA =
LDFLAGS_CUDA = -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/

INCLUDES_CUDA += $(INCLUDES)



all:
	# $(CXX) test_integrator.cc $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)
	# $(CXX) test_sorting.cc $(INCLUDES) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)
	$(CXX) test_sorting.cc $(INCLUDES) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)
	$(CXX) test_sorting_simp.cc -o b.out $(CXXFLAGS) $(QFECXXFLAGS) $(QFEINCLUDES) $(INCLUDES) $(LDFLAGS)
	$(CXX) get_geometry.cc -o geom $(CXXFLAGS) $(QFECXXFLAGS) $(QFEINCLUDES) $(INCLUDES) $(LDFLAGS)
	# $(CXX) get_coords.cc -o coord $(CXXFLAGS) $(QFECXXFLAGS) $(QFEINCLUDES) $(INCLUDES) $(LDFLAGS)
