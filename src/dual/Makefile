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


GRP_DIR=$(shell pwd)/newQFE/grp/
CXXFLAGS+=-DGRP_DIR="\"$(GRP_DIR)\""
NVCCFLAGS+=-DGRP_DIR="\"$(GRP_DIR)\""


# simp:
	#$(CXX) get_geometry.cc $(INCLUDES) $(CXXFLAGS)
	# $(CXX) test_geometry.cc $(INCLUDES) $(CXXFLAGS)
	# $(CXX) test_geometry_flat.cc $(INCLUDES) $(CXXFLAGS)
	# $(CXX) test_u1.cc $(INCLUDES) $(CXXFLAGS) -o c.out
	# $(CXX) test_dirac.cc $(INCLUDES) $(CXXFLAGS)
	# $(CXX) test_dirac_flat.cc $(INCLUDES) $(CXXFLAGS)
	# $(CXX) test_dirac_dual.cc $(INCLUDES) # $(CXXFLAGS)
	# $(CXX) test_dirac_tim.cc $(INCLUDES) $(CXXFLAGS)

SRCDIR="./"

dual:
	# $(CXX) get_geometry.cc $(INCLUDES) $(CXXFLAGS)
	# $(CXX) test_geometry.cc $(INCLUDES) $(CXXFLAGS)
	# $(CXX) test_geometry_flat.cc $(INCLUDES) $(CXXFLAGS)
	# $(CXX) test_u1_dual.cc $(INCLUDES) $(CXXFLAGS)
	# $(CXX) test_dirac.cc $(INCLUDES) $(CXXFLAGS)
	# $(CXX) test_dirac_flat.cc $(INCLUDES) $(CXXFLAGS)
	# $(CXX) test_dirac_dual.cc $(INCLUDES) # $(CXXFLAGS)
	# $(CXX) test_sparse.cc $(INCLUDES) $(CXXFLAGS)
	# $(NVCC) test_cg.cu $(NVCCFLAGS) $(INCLUDES_CUDA) $(LDFLAGS_CUDA)
	# $(NVCC) $(SRCDIR)main.cu $(NVCCFLAGS) $(INCLUDES_CUDA) $(LDFLAGS_CUDA)
	# $(NVCC) $(SRCDIR)test.cu $(NVCCFLAGS) $(INCLUDES_CUDA) $(LDFLAGS_CUDA)
	# $(NVCC) eig.cu -o eig $(NVCCFLAGS) $(INCLUDES_CUDA) $(LDFLAGS_CUDA)
	# $(NVCC) grad.cu $(NVCCFLAGS) $(INCLUDES_CUDA) $(LDFLAGS_CUDA)
	$(NVCC) prop.cu $(NVCCFLAGS) $(INCLUDES_CUDA) $(LDFLAGS_CUDA)
	# $(NVCC) prop_dirac.cu -o c.out $(NVCCFLAGS) $(INCLUDES_CUDA) $(LDFLAGS_CUDA)
	# $(NVCC) test.cu $(NVCCFLAGS) $(INCLUDES_CUDA) $(LDFLAGS_CUDA)
	# $(CXX) test_hmc_scalar.cc $(INCLUDES) # $(CXXFLAGS)
	# $(CXX) hmc_dirac_dual.cc $(INCLUDES) $(CXXFLAGS)
	# $(CXX) test_dirac_tim.cc $(INCLUDES) $(CXXFLAGS)

# eig: test_dirac.cu header_cusolver.hpp
# 	$(NVCC) $< $(NVCCFLAGS) $(INCLUDES_CUDA) -o b.out
# # eig: test_dirac_dual.cu header_cusolver.hpp
# # 	$(NVCC) $< $(NVCCFLAGS) $(INCLUDES_CUDA) -o b.out

# geometry:
# 	$(CXX) get_geometry.cc $(INCLUDES) $(CXXFLAGS)
