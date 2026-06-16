#!/usr/bin/env bash
# Compile + run jj_exact_diag_deter_claude.cu (INSERTION-DIAGONAL exact conserved current, tp+sp) for
# L=1,2 with the LATTICE overlap propagator (reads the EXISTING prop_deter_L<L>/Dinv.0.h5).  Free field.
# Writes data_free_vmRe0.000000vmIm0.000000/corr_deter_exactdiag_L<L>/corr.0.h5, keys h0/t0_b/{tp,sp}/Vpp.
# Compare with loc (corr_deter_local_L*: s3=tp, s1+s2=sp) and disp (corr_deter_disp_L*: sp) in
# comp_trio_jj_claude.ipynb.
#
# COST: op_K applies = N * n_ins (== jj_kbuild_exact).  L=1 (N=3072, n_ins=42) ~ 1.3e5 solves: minutes.
#   L=2 (N=10752, n_ins=162) ~ 1.7e6 overlap solves: HOURS -- run L=1 first; do L=2 only if patient /
#   on a fast GPU.  (Edit the `for L in` line to "1" to do L=1 only.)
# OMP_NUM_THREADS speeds the conn_shift time-shift (CPU); set to the cores you want.
# (Continuum prop instead: add  --prop-file cont_prop_L<L>/Dinv.0.h5 --out-tag cont .)
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-16}

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

for L in 1 2; do
  BIN=jj_exact_diag_deter_L${L}.o
  LOG=jj_exact_diag_deter_lat_L${L}_claude.log
  echo "### compile L=${L} ###" | tee "$LOG"
  $NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=${L} jj_exact_diag_deter_claude.cu -o "$BIN" \
    2>&1 | tee -a "$LOG" || { echo "BUILD FAILED L=${L}"; exit 1; }
  echo "### run L=${L}: LATTICE prop (reads prop_deter_L${L}/Dinv.0.h5), n-t0=2, OMP=${OMP_NUM_THREADS} ###" | tee -a "$LOG"
  ./"$BIN" --n-t0 2 \
    2>&1 | tee -a "$LOG" || { echo "RUN FAILED L=${L}"; exit 1; }
  echo "### done L=${L} ###" | tee -a "$LOG"
done
echo "### all L done ###"
