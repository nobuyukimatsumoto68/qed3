#!/usr/bin/env bash
# FAST rerun of the diagonally-n-summed DISP current (sp+tp) for L=1,2, lattice prop.  This regenerates
# corr_deter_disp_L{1,2} with the NEW temporal channel (tp) and the W_ov_kappa kernel (sp numerics
# unchanged).  "Fast" = no overlap solves: disp builds the bare kernel + traces against the dense P.
# (loc summed corr_deter_local_L{1,2} is already current -- s1,s2,s3,ylm present, summed path unchanged
#  -- so it is NOT rerun here.  exact summed is the SLOW one: separate tmp_exactsum_free script.)
#
# PRE-REQ: delete the stale (tp-less) disp outputs first (run yourself; no rm in this script):
#   rm data_free_vmRe0.000000vmIm0.000000/corr_deter_disp_L1/corr.0.h5
#   rm data_free_vmRe0.000000vmIm0.000000/corr_deter_disp_L2/corr.0.h5
# (Otherwise the "complete"-gate makes the run SKIP and you keep the old sp-only file.)
# Reads the EXISTING data_free_*/prop_deter_L<L>/Dinv.0.h5 (read-only).
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES=0
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

for L in 1 2; do
  BIN=jj_disp_deter_L${L}.o
  LOG=jj_disp_summed_L${L}_claude.log
  echo "### compile disp L=${L} ###" | tee "$LOG"
  $NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=${L} jj_disp_deter_claude.cu -o "$BIN" \
    2>&1 | tee -a "$LOG" || { echo "BUILD FAILED L=${L}"; exit 1; }
  echo "### run disp SUMMED L=${L} (no --ins => sum over links/sites; sp+tp), lattice prop ###" | tee -a "$LOG"
  ./"$BIN" --n-t0 2 \
    2>&1 | tee -a "$LOG" || { echo "RUN FAILED L=${L}"; exit 1; }
  echo "### done L=${L} ###" | tee -a "$LOG"
done
echo "### all L done -- open comp_threesome_jj_claude.ipynb (summed) for the Gs/Gt -> -2 check ###"
