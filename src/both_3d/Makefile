CXX=g++
# INCLUDES=-I"/mnt/hdd_barracuda/qfe_mod//include"
# GRP_DIR="/mnt/hdd_barracuda/qfe_mod/grp/"
# INCLUDES+=-DGRP_DIR="\"$(GRP_DIR)\""

# NVCC=/usr/local/cuda-12.6/bin/nvcc # nvcc
# NVCCFLAGS=-arch=sm_70 -pg -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp,-march=native # -diag-suppress<1650-D>
# =======
INCLUDES=-I"/projectnb/qfe/nmatsum/qed3/opt/eigen" # -I"/mnt/hdd_barracuda/qed3/newQFE/include"
# GRP_DIR="/mnt/hdd_barracuda/qed3/newQFE/grp/"
# INCLUDES+=-DGRP_DIR="\"$(GRP_DIR)\""

NVCC=nvcc # /usr/local/cuda-12.6/bin/nvcc # nvcc
NVCCFLAGS=-arch=sm_70 -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp # -diag-suppress<1650-D>
# CXXFLAGS=-O3 -std=c++20 omp -Xcompiler -fopenmp,-march=native # -diag-suppress<1650-D>
# >>>>>>> d7b97ca92d97e875a0371c9b0a724cd0e6e17e9d

INCLUDES_CUDA= #-I/usr/local/cuda-12.6/include/
LDFLAGS_CUDA= # -L/usr/local/cuda-12.6/lib64/
INCLUDES_CUDA += $(INCLUDES)

NVCCFLAGS+= # -DGRP_DIR="\"$(GRP_DIR)\""

SRCS := $(wildcard *.cu)
# SRCS := $(wildcard *.h)
OBJS := $(SRCS:%.cu=%.o)
# OBJS += $(SRCS:%.h=%.o)

all: $(OBJS)

%.d: %.cu
	$(NVCC) $(INCLUDES_CUDA) -M $< -o $@
	# $(CXX) $(INCLUDES) -M $< -o $@

# %.d: %.h
# 	$(NVCC) $(INCLUDES_CUDA) -M $< -o $@

include $(SRCS:.cu=.d)

# %.o: %.h
# 	$(NVCC) $(NVCCFLAGS) $(INCLUDES_CUDA) $(LDFLAGS_CUDA) -c $< -o $@ # $(<:.d=.cu)

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) $(INCLUDES_CUDA) $(LDFLAGS_CUDA) $< -o $@ # $(<:.d=.cu)
	# $(CXX) $(NVCCFLAGS) $(INCLUDES_CUDA) $(LDFLAGS_CUDA) $< -o $@ # $(<:.d=.cu)

