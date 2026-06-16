#!/usr/bin/env bash
# FAST: diagonally-n-summed AXIAL loc + disp currents (no overlap solves; bare kernels traced against the
# dense lattice P) for L=1,2.  Writes corr_deter_local_axial_L<L> (s1,s2,s3,ylm) and
# corr_deter_disp_axial_L<L> (sp,tp,ylm), single complex channel "Apm".  Reads the EXISTING
# data_free_*/prop_deter_L<L>/Dinv.0.h5 (read-only).  GPU1, OMP=4.
#
# PRE-REQ (run yourself; no rm in this script): delete any stale outputs so the complete-gate does not skip:
#   rm data_free_vmRe0.000000vmIm0.000000/corr_deter_local_axial_L1/corr.0.h5
#   rm data_free_vmRe0.000000vmIm0.000000/corr_deter_local_axial_L2/corr.0.h5
#   rm data_free_vmRe0.000000vmIm0.000000/corr_deter_disp_axial_L1/corr.0.h5
#   rm data_free_vmRe0.000000vmIm0.000000/corr_deter_disp_axial_L2/corr.0.h5
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES=1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

for L in 1 2; do
  for SRC in jj_local_axial_deter jj_disp_axial_deter; do
    BIN=${SRC}_L${L}.o
    LOG=${SRC}_L${L}_claude.log
    echo "### compile ${SRC} L=${L} ###" | tee "$LOG"
    $NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=${L} ${SRC}_claude.cu -o "$BIN" \
      2>&1 | tee -a "$LOG" || { echo "BUILD FAILED ${SRC} L=${L}"; exit 1; }
    echo "### run ${SRC} SUMMED L=${L} (no --ins; lattice prop, n-t0=2) ###" | tee -a "$LOG"
    ./"$BIN" --n-t0 2 \
      2>&1 | tee -a "$LOG" || { echo "RUN FAILED ${SRC} L=${L}"; exit 1; }
    echo "### done ${SRC} L=${L} ###" | tee -a "$LOG"
  done
done
echo "### all done -- corr_deter_{local,disp}_axial_L{1,2}; open comp_trio_jj_axial_claude.ipynb ###"
