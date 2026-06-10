#!/usr/bin/env bash
# D(L2): local op + CONTINUUM G at L=2 -> corr_deter_local_cont_L2
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
module load cuda/12.8 2>/dev/null; module load gcc/13.2.0 2>/dev/null
export CUDA_VISIBLE_DEVICES=0
[ -f jj_local_deter_L2.o ] || { echo "missing jj_local_deter_L2.o"; exit 1; }
[ -f cont_prop_L2/Dinv.0.h5 ] || { echo "missing cont_prop_L2"; exit 1; }
rm -f data_free_vmRe0.000000vmIm0.000000/corr_deter_local_cont_L2/corr.0.h5*
./jj_local_deter_L2.o --mass-re 0.0 --mass-im 0.0 --n-t0 2 \
  --prop-file cont_prop_L2/Dinv.0.h5 --out-tag cont || { echo "D(L2) FAILED"; exit 1; }
echo "DONE D(L2)"
