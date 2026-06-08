#!/bin/bash
# Build bench_nthreads_claude.cu (lives in both_3d/) -> bench_256/512/1024_claude.o in benchmark/.
# Compiles from the both_3d/ source dir so bare includes resolve via -Iincludes/ (same as production).
set -e

SRCDIR="$(cd "$(dirname "$0")/.." && pwd)"      # both_3d/
BENCHDIR="$(cd "$(dirname "$0")" && pwd)"       # both_3d/benchmark/

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

SRC="${SRCDIR}/bench_nthreads_claude.cu"

cd "${SRCDIR}"
for NT in 256 512 1024; do
  OUT="${BENCHDIR}/bench_${NT}_claude.o"
  echo "==> Compiling NThreadsPerBlock=${NT} -> ${OUT}"
  eval nvcc ${NVCCFLAGS} -DNThreadsPerBlock=${NT} ${INCLUDES} ${LDFLAGS} "${SRC}" -o "${OUT}"
done

echo "==> Done. Binaries:"
ls -lh "${BENCHDIR}"/bench_*_claude.o
