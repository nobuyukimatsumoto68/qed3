#!/usr/bin/env bash
# Run the LOCAL deterministic jj with the LATTICE overlap propagator for L=1,2 (free field).
# Reads the EXISTING data_free_.../prop_deter_L<L>/Dinv.0.h5 (no --prop-file) and writes
# data_free_.../corr_deter_local_L<L>/corr.0.h5 (no --out-tag).  Does NOT rebuild the propagator
# (the dense LU is expensive; prop_deter_L1=302M, L2=3.7G, L4=55G already exist -- the L4 file is big
#  because the free propagator also dumps the dense D_ov for debug; load_mat reads only Dm_inv).
# Companion of tmp_claude.sh (which feeds the CONTINUUM propagator -> corr_deter_local_cont_L<L>).
#
# NOTE 1 (rerun): the existing corr_deter_local_L{1,2,4}/corr.0.h5 carry the OLD tp/sp keys (pre
#   s1/s2/s3 redesign) with their "complete" flag, so the program SKIPS them.  Delete first (run
#   yourself):
#     rm data_free_vmRe0.000000vmIm0.000000/corr_deter_local_L1/corr.0.h5
#     rm data_free_vmRe0.000000vmIm0.000000/corr_deter_local_L2/corr.0.h5
#     rm data_free_vmRe0.000000vmIm0.000000/corr_deter_local_L4/corr.0.h5
# NOTE 2 (RAM): L=4 loads the dense Dm_inv (N=41472) -> ~41 GB host RAM at the load peak (one N^2
#   double buffer + the N^2 complex matrix).  Run when the machine is idle; if it OOMs, run L=4 on a
#   larger-RAM node.  L=1,2 are comfortable (<3 GB).
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

for L in 1 2 4; do
  BIN=jj_local_deter_L${L}.o
  LOG=jj_local_deter_lat_L${L}_claude.log
  echo "### compile L=${L} ###" | tee "$LOG"
  $NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=${L} jj_local_deter_claude.cu -o "$BIN" \
    2>&1 | tee -a "$LOG" || { echo "BUILD FAILED L=${L}"; exit 1; }
  echo "### run L=${L}: LATTICE prop (reads prop_deter_L${L}/Dinv.0.h5), n-t0=2 ###" | tee -a "$LOG"
  ./"$BIN" --n-t0 2 \
    2>&1 | tee -a "$LOG" || { echo "RUN FAILED L=${L}"; exit 1; }
  echo "### done L=${L} ###" | tee -a "$LOG"
done
echo "### all L done ###"
