#!/usr/bin/env bash
# Compile + run jj_exact_diag_deter_free_claude.cu: FREE single-insertion exact <J(a,0)J(a,t)>, tp+sp,
# for L=1,2 with the LATTICE overlap propagator (reads EXISTING prop_deter_L<L>/Dinv.0.h5).  ins=0.
# Runs on GPU 1 (binary calls cudaSetDevice(0) on the only visible device).
#
# CACHING: the dense K(a,0) (tp+sp) is built ONCE (N op_K applies each) and saved to
#   data_free_Kcache_L<L>/K_ins0.h5.  On any later run it is READ and the op_K applies are SKIPPED
#   (only the cheap dense matmul + time-shift trace run).  To force a rebuild, delete K_ins0.h5.
# Output: data_free_vmRe0.000000vmIm0.000000/corr_deter_exact1_L<L>/corr.0.h5, keys h0/t0_b/{tp,sp}/Vpp
#   = single-insertion <J(0,0)J(0,t)> (NO sum over insertions).
#
# COST (first run, builds K): 2*N op_K applies.  L=1 (N=3072) ~ minutes ; L=2 (N=10752) ~ 10-40 min.
#   Cached reruns: seconds.   (Continuum prop: add --prop-file cont_prop_L<L>/Dinv.0.h5 --out-tag cont.)
# RERUN of the correlator: corr_deter_exact1_L*/corr.0.h5 is "complete"-gated -> rm it to redo
#   (the K cache is separate and is reused).
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES=1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

for L in 1 2; do
  BIN=jj_exact_diag_deter_free_L${L}.o
  LOG=jj_exact_diag_deter_free_L${L}_claude.log
  echo "### compile L=${L} ###" | tee "$LOG"
  $NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=${L} jj_exact_diag_deter_free_claude.cu -o "$BIN" \
    2>&1 | tee -a "$LOG" || { echo "BUILD FAILED L=${L}"; exit 1; }
  echo "### run L=${L}: ins=0, LATTICE prop, n-t0=2, OMP=${OMP_NUM_THREADS} ###" | tee -a "$LOG"
  ./"$BIN" --ins 0 --n-t0 2 \
    2>&1 | tee -a "$LOG" || { echo "RUN FAILED L=${L}"; exit 1; }
  echo "### done L=${L} ###" | tee -a "$LOG"
done
echo "### all L done ###"
