#!/usr/bin/env bash
# compile + run jj_local_deter_claude.cu for L=1,2,4 with the analytic continuum propagator.
# CFT check: Gt (s3, Eq.4.31), Gs (s1+s2, Eq.4.28), Ylm tower (Eq.4.35).
# Writes data_free_vmRe0.000000vmIm0.000000/corr_deter_local_cont_L<L>/corr.0.h5.
#
# NOTE 1 (rerun): if corr_deter_local_cont_L<L>/corr.0.h5 already exists with its "complete" flag,
#   the program SKIPS it (no regenerate).  Delete the stale L=1 and L=2 files first (run yourself):
#     rm data_free_vmRe0.000000vmIm0.000000/corr_deter_local_cont_L1/corr.0.h5
#     rm data_free_vmRe0.000000vmIm0.000000/corr_deter_local_cont_L2/corr.0.h5
#   (L=4 has no prior output.)
# NOTE 2 (RAM): the deterministic trace loads the full dense P.  L=4 (cont_prop_L4 = 26 GB) needs
#   ~41 GB host RAM at load peak -- run it when the machine is otherwise idle; if it OOMs, run L=4
#   on a larger-RAM node.
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

for L in 1 2 4; do
  BIN=jj_local_deter_L${L}.o
  LOG=jj_local_deter_cont_L${L}_claude.log
  echo "### compile L=${L} ###" | tee "$LOG"
  $NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=${L} jj_local_deter_claude.cu -o "$BIN" \
    2>&1 | tee -a "$LOG" || { echo "BUILD FAILED L=${L}"; exit 1; }
  echo "### run L=${L}: cont prop, out-tag=cont, n-t0=2 ###" | tee -a "$LOG"
  ./"$BIN" --prop-file cont_prop_L${L}/Dinv.0.h5 --out-tag cont --n-t0 2 \
    2>&1 | tee -a "$LOG" || { echo "RUN FAILED L=${L}"; exit 1; }
  echo "### done L=${L} ###" | tee -a "$LOG"
done
echo "### all L done ###"
