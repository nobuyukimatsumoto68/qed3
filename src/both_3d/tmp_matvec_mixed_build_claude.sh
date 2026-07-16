#!/usr/bin/env bash
# _claude: Chunk-0 GATE of the mixed-precision inner-solver plan (mixedprec_impl_plan_claude.md).
# Builds + RUNS test_matvec_mixed_claude.cu at L1/L2/L4 and reports, per L, the fp32-vs-fp64
# D_W^dag D_W matvec rel error + speedup. Purpose: confirm the fp32 bandwidth win is worthwhile
# before building the reliable-update CG. Reads back via the .log. No rm.
#
# Run on ONE GPU (set CUDA_VISIBLE_DEVICES below if needed). USER runs this; Claude reads the log.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=tmp_matvec_mixed_build_claude.log
export OMP_NUM_THREADS=4
export CUDA_VISIBLE_DEVICES=0

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_matvec_mixed_claude.cu

LVALS=( 1 2 4 )

{
  echo "######## fp32-vs-fp64 D_W^dag D_W matvec GATE  $(date) ########"
  for LV in "${LVALS[@]}"; do
    bin="test_matvec_mixed_L${LV}.x"
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
