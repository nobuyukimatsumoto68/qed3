#!/usr/bin/env bash
# _claude: Chunk-1 validation of the single-shift mixed-precision reliable-update CG
# (matpoly_mixed_claude.h / mixedprec_impl_plan_claude.md). Builds + RUNS test_solve_mixed_claude.cu
# at L1/L2/L4: for the seed/mid/largest Zolotarev pole it checks mixed-vs-fp64 solution agreement,
# the mixed solution's TRUE fp64 residual (reliable update reaches tol), and the wall-time speedup.
# Reads back via the .log. No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=tmp_solve_mixed_build_claude.log
export OMP_NUM_THREADS=4
export CUDA_VISIBLE_DEVICES=0

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_solve_mixed_claude.cu

LVALS=( 1 2 4 )

{
  echo "######## single-shift mixed-prec CG validation  $(date) ########"
  for LV in "${LVALS[@]}"; do
    bin="test_solve_mixed_L${LV}.x"
    echo ""
    echo "==================== BUILD LREF=${LV} ===================="
    if $NVCC $NVCCFLAGS $INCLUDES -DLREF=${LV} $SRC $LDFLAGS -o "$bin" ; then
      echo ">>> BUILD OK: LREF=${LV}"
      echo "-------------------- RUN LREF=${LV} --------------------"
      ./"$bin"
    else
      echo "!!! BUILD FAILED: LREF=${LV}"
    fi
  done
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
