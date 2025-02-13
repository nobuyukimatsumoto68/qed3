# CXX=g++
# CXXFLAGS=-std=c++17 -g -O3
INCLUDES=-I"../../newQFE/include"

NVCC=/usr/local/cuda-12.6/bin/nvcc # nvcc
# NVCCFLAGS = -arch=sm_70 -O3 -lcusolver -std=c++17 # -diag-suppress<1650-D>
# NVCCFLAGS=-w -arch=sm_70 -O3 -g -std=c++17 -lcublas -lcusolver -lgomp -Xcompiler -fopenmp,-mavx2,-ftree-vectorize # -diag-suppress<1650-D>
NVCCFLAGS=-arch=sm_70 -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp,-march=native # -diag-suppress<1650-D>

INCLUDES_CUDA=-I/usr/local/cuda-12.6/include/
# LDFLAGS_CUDA=-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/
LDFLAGS_CUDA=-L/usr/local/cuda-12.6/lib64/

INCLUDES_CUDA += $(INCLUDES)


GRP_DIR="/mnt/hdd_barracuda/qed3/newQFE/grp/"
CXXFLAGS+=-DGRP_DIR="\"$(GRP_DIR)\""
NVCCFLAGS+=-DGRP_DIR="\"$(GRP_DIR)\""

simp:
	# $(NVCC) test_dirac.cu $(NVCCFLAGS) $(INCLUDES_CUDA) $(LDFLAGS_CUDA) # $(INCLUDES) $(CXXFLAGS)
	$(NVCC) prop.cu $(NVCCFLAGS) $(INCLUDES_CUDA) $(LDFLAGS_CUDA) # $(INCLUDES) $(CXXFLAGS)

SRCDIR="./"
