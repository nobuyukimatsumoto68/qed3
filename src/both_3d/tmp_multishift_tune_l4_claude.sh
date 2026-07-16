#!/usr/bin/env bash
# _claude: L4 tuning of the mixed-prec MULTISHIFT -- joint (tol_f x n_nested) scan on REAL configs.
# n_nested = # nested fp32 defect-correction stages inside EACH per-shift cleanup (nn=0 = pure fp64
# cleanup = the ~1.35x baseline). Goal: push L4 net past 1.35x toward the ~1.5x floor by doing the
# cleanup in fp32 too. Nf6 gsq8 massless, Nt=128. Reads back via the .log. No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=tmp_multishift_tune_l4_v2_claude.log   # fp32-cleanup stage now targets abs crit_d (over-solve fixed)
export OMP_NUM_THREADS=4
export CUDA_VISIBLE_DEVICES=0

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_multishift_mixed_claude.cu

L4DIR="Nf6_gsq8.000000at0.200000nu01.000000nt128L4"
CFGS=( "$L4DIR/ckpoint_lat.40" "$L4DIR/ckpoint_lat.60" "$L4DIR/ckpoint_lat.80" )
BIN=test_multishift_tune_l4.x

{
  echo "######## L4 mixed-multishift (tol_f x n_nested) tuning  $(date) ########"
  echo "==================== BUILD LREF=4 Nt=128 ===================="
  if $NVCC $NVCCFLAGS $INCLUDES -DLREF=4 -DNTIME=128 $SRC $LDFLAGS -o "$BIN" ; then
    echo ">>> BUILD OK"
    for c in "${CFGS[@]}"; do
      echo "-------------------- RUN config=$c --------------------"
      if [ -f "$c" ]; then ./"$BIN" "$c"; else echo "!!! MISSING config $c"; fi
    done
  else
    echo "!!! BUILD FAILED"
  fi
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
