#!/bin/bash
# Build bench_nthreads_claude.cu -> bench_256_claude.o / bench_512_claude.o / bench_1024_claude.o
# Run from benchmark/ directory (or anywhere; uses absolute paths).
set -e

SRCDIR="$(cd "$(dirname "$0")/.." && pwd)"
BENCHDIR="$(cd "$(dirname "$0")" && pwd)"

source /lustre2/qed3/env.sh

INCLUDES="-I\"${SRCDIR}/includes/\" \
          -I\"/project/qed3//qed3/qfe_mod/include/\" \
          -I\"/project/qed3//highfive/include/\" \
          -I\"/srv/software/el8/x86_64/eb/HDF5/1.14.2-GCC-12.3.0-serial/include/\" \
          -I\"/project/qed3//gsl/include/\""

NVCCFLAGS="-arch=sm_80 -g -O3 -std=c++20 \
  -L${CUDA_HOME}/lib64 -L${CUDA_MATH}/lib64 -L${CUDA_COMP}/lib -L${CUDA_HOME}/lib64/stubs \
  -lcudart -lcuda -lnvidia-ml -lcublas -lcufft -ldl -lgomp -Xcompiler -fopenmp"

LDFLAGS="-L\"/srv/software/el8/x86_64/eb/HDF5/1.14.2-GCC-12.3.0-serial/lib/\" \
         -L\"/project/qed3//gsl/lib/\" -lhdf5 -lgsl -lgslcblas -lm"

SRC="${BENCHDIR}/bench_nthreads_claude.cu"

for NT in 256 512 1024; do
  OUT="${BENCHDIR}/bench_${NT}_claude.o"
  echo "==> Compiling NThreadsPerBlock=${NT} -> ${OUT}"
  eval nvcc ${NVCCFLAGS} -DNThreadsPerBlock=${NT} ${INCLUDES} ${LDFLAGS} "${SRC}" -o "${OUT}"
done

echo "==> Done. Binaries:"
ls -lh "${BENCHDIR}"/bench_*_claude.o
