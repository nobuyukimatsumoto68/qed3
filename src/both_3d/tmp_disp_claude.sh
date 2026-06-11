#!/usr/bin/env bash
# Compile + run jj_disp_deter_claude.cu (DISPLACED link sp current) for L=1,2 with the LATTICE overlap
# propagator (reads the EXISTING prop_deter_L<L>/Dinv.0.h5; no LU rebuild).  Free field.
# Writes data_free_vmRe0.000000vmIm0.000000/corr_deter_disp_L<L>/corr.0.h5, key h0/t0_b/sp/{Vpp,Vmm}.
# Compare against the loc result corr_deter_local_L<L> (s1+s2 = G_s) in comp_loc_disp_jj_claude.ipynb.
#
# (Continuum prop instead: add  --prop-file cont_prop_L<L>/Dinv.0.h5 --out-tag cont  -> corr_deter_disp_cont_L<L>.)
# Runs on GPU 1 (CUDA_VISIBLE_DEVICES=1; the binary calls cudaSetDevice(0) on the only visible device).
#
# RERUN NOTE: corr_deter_disp_L{1,2}/corr.0.h5 already exist from the first run and are "complete"-gated,
#   so the program will SKIP them.  Delete both first (run yourself) to regenerate with the i-fix
#   (W_d = C/kappa, no Dirac-i => sign now matches loc):
#     rm data_free_vmRe0.000000vmIm0.000000/corr_deter_disp_L1/corr.0.h5
#     rm data_free_vmRe0.000000vmIm0.000000/corr_deter_disp_L2/corr.0.h5
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES=1

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

for L in 1 2; do
  BIN=jj_disp_deter_L${L}.o
  LOG=jj_disp_deter_lat_L${L}_claude.log
  echo "### compile L=${L} ###" | tee "$LOG"
  $NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=${L} jj_disp_deter_claude.cu -o "$BIN" \
    2>&1 | tee -a "$LOG" || { echo "BUILD FAILED L=${L}"; exit 1; }
  echo "### run L=${L}: LATTICE prop (reads prop_deter_L${L}/Dinv.0.h5), n-t0=2 ###" | tee -a "$LOG"
  ./"$BIN" --n-t0 2 \
    2>&1 | tee -a "$LOG" || { echo "RUN FAILED L=${L}"; exit 1; }
  echo "### done L=${L} ###" | tee -a "$LOG"
done
echo "### all L done ###"
