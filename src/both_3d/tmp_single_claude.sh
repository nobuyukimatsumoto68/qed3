#!/usr/bin/env bash
# Compile + run the SINGLE-INSERTION counterparts of loc and disp (--ins 0) for L=1,2,4 with the LATTICE
# overlap propagator (reads EXISTING prop_deter_L<L>/Dinv.0.h5; no LU rebuild).  Free field.
# L=4 loads the dense P (D_ov^-1, ~41 GB RAM via the lean load_mat); the trace itself is cheap.
# (exact1 has no L=4 counterpart -- the dense K would need ~41k op_K applies -- so the threesome ratio
#  at L=4 shows loc1 and disp1 only, both -> -2.)
#
# These mirror jj_exact_diag_deter_free (corr_deter_exact1, also ins=0): NO sum over sites/links and
# NO summation weight (RAW trace tr[W(t0) P W(t) P]), so loc1/disp1/exact1 are apples-to-apples.
#   loc1  (jj_local_deter --ins 0): single SITE 0, channels s1,s2,s3 -> corr_deter_local1_L<L>
#   disp1 (jj_disp_deter  --ins 0): single LINK 0, channel sp        -> corr_deter_disp1_L<L>
# Compare with corr_deter_exact1_L<L> (tp: site 0 ; sp: link 0) in comp_threesome_jj_claude.ipynb.
#
# Runs on GPU 0 (CUDA_VISIBLE_DEVICES=0; the binary calls cudaSetDevice(0) on the only visible device).
# RERUN NOTE: corr_deter_{local1,disp1}_L*/corr.0.h5 are "complete"-gated -> the program SKIPS them if
#   present.  To regenerate, delete them yourself (no rm in this script by policy):
#     rm data_free_vmRe0.000000vmIm0.000000/corr_deter_local1_L{1,2,4}/corr.0.h5
#     rm data_free_vmRe0.000000vmIm0.000000/corr_deter_disp1_L{1,2,4}/corr.0.h5
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES=0
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

for L in 1 2 4; do
  LOG=jj_single_L${L}_claude.log
  echo "### compile loc1 + disp1 L=${L} ###" | tee "$LOG"
  $NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=${L} jj_local_deter_claude.cu -o jj_local_deter_L${L}.o \
    2>&1 | tee -a "$LOG" || { echo "BUILD FAILED loc1 L=${L}"; exit 1; }
  $NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=${L} jj_disp_deter_claude.cu  -o jj_disp_deter_L${L}.o \
    2>&1 | tee -a "$LOG" || { echo "BUILD FAILED disp1 L=${L}"; exit 1; }
  echo "### run loc1 L=${L}: --ins 0, LATTICE prop, n-t0=2 ###" | tee -a "$LOG"
  ./jj_local_deter_L${L}.o --ins 0 --n-t0 2 \
    2>&1 | tee -a "$LOG" || { echo "RUN FAILED loc1 L=${L}"; exit 1; }
  echo "### run disp1 L=${L}: --ins 0, LATTICE prop, n-t0=2 ###" | tee -a "$LOG"
  ./jj_disp_deter_L${L}.o --ins 0 --n-t0 2 \
    2>&1 | tee -a "$LOG" || { echo "RUN FAILED disp1 L=${L}"; exit 1; }
  echo "### done L=${L} ###" | tee -a "$LOG"
done
echo "### all L done ###"
