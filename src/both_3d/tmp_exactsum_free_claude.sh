#!/usr/bin/env bash
# SLOW: diagonally-n-summed EXACT current (jj_exact_diag_deter_free --sum), Eq.(4.29), free, lattice prop.
# Output corr_deter_exactsum_L<L> (keys h0/t0_b/{tp,sp}/{Vpp,Vmm} + disc).  Reads prop_deter_L<L>/Dinv.0.h5.
#
# COST: build-use-discard dense K per insertion => (n_sites+n_links)*N op_K (multishift) solves.
#   L=1 (N=3072, 12 sites + 30 links): ~minutes.
#   L=2 (N=10752, 42 sites + 120 links): ~HOURS.  Edit the `for L` line to add 2 when ready to wait.
# No caching (all-n K cache would be ~300 GB at L=2).  GPU0, OMP=4.
#
# RERUN NOTE: corr_deter_exactsum_L<L>/corr.0.h5 is written atomically (tmp->rename) but is NOT
# complete-gated for skip here; to redo, delete it yourself first (no rm in this script):
#   rm data_free_vmRe0.000000vmIm0.000000/corr_deter_exactsum_L1/corr.0.h5
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES=1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

for L in 2; do   # add 2 here when ready for the multi-hour L=2 run
  BIN=jj_exact_diag_deter_free_L${L}.o
  LOG=jj_exactsum_free_L${L}_claude.log
  echo "### compile L=${L} ###" | tee "$LOG"
  $NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=${L} jj_exact_diag_deter_free_claude.cu -o "$BIN" \
    2>&1 | tee -a "$LOG" || { echo "BUILD FAILED L=${L}"; exit 1; }
  echo "### run --sum L=${L} (lattice prop, n-t0=2) ###" | tee -a "$LOG"
  ./"$BIN" --sum --n-t0 2 \
    2>&1 | tee -a "$LOG" || { echo "RUN FAILED L=${L}"; exit 1; }
  echo "### done L=${L} ###" | tee -a "$LOG"
done
echo "### all done -- corr_deter_exactsum_L*; check ratio in comp_threesome_jj_claude.ipynb ###"
