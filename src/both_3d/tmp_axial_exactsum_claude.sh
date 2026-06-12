#!/usr/bin/env bash
# SLOW: diagonally-n-summed EXACT AXIAL current (jj_exact_axial_deter_free --sum), Eq.(3.55)+(4.29), free,
# lattice prop.  Output corr_deter_exactsum_axial_L<L> (keys h0/t0_b/{tp,sp}/Apm + ylm/l{l}/Apm + disc).
# Reads prop_deter_L<L>/Dinv.0.h5.  Builds the G=(1-D_ov) cache once (data_free_Gcache_L<L>/G.h5).
#
# COST: like the vector exactsum -- build-use-discard dense K per insertion => (n_sites+n_links)*N op_K
#   solves, plus a one-time G build (N op_oneMinusD) and 4 dense N x N matmuls per insertion.
#   L=1 (N=3072): ~minutes.  L=2 (N=10752): ~HOURS.  Edit the `for L` line to add 2 when ready to wait.
# GPU1, OMP=4.
#
# RERUN NOTE (no rm in this script): corr_deter_exactsum_axial_L<L>/corr.0.h5 is written atomically but is
# NOT complete-gated for skip; to redo, delete it yourself first:
#   rm data_free_vmRe0.000000vmIm0.000000/corr_deter_exactsum_axial_L1/corr.0.h5
# (The G cache is reusable across runs/insertions -- leave it.)
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES=1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

for L in 1; do   # add 2 here when ready for the multi-hour L=2 run
  BIN=jj_exact_axial_deter_free_L${L}.o
  LOG=jj_exactsum_axial_free_L${L}_claude.log
  echo "### compile L=${L} ###" | tee "$LOG"
  $NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=${L} jj_exact_axial_deter_free_claude.cu -o "$BIN" \
    2>&1 | tee -a "$LOG" || { echo "BUILD FAILED L=${L}"; exit 1; }
  echo "### run --sum L=${L} (lattice prop, n-t0=2) ###" | tee -a "$LOG"
  ./"$BIN" --sum --n-t0 2 \
    2>&1 | tee -a "$LOG" || { echo "RUN FAILED L=${L}"; exit 1; }
  echo "### done L=${L} ###" | tee -a "$LOG"
done
echo "### all done -- corr_deter_exactsum_axial_L*; check comp_trio_jj_axial_claude.ipynb ###"
