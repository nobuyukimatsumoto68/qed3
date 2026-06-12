#!/usr/bin/env bash
# L=4 loc + disp AXIAL (trace-only, NO solves) -- extend the disp-vs-loc L-trend to L=4 to test whether
# the disp sp ratio keeps converging toward loc's -1 (discretization) or stalls.  Reads the EXISTING
# data_free_*/prop_deter_L4/Dinv.0.h5 ("Dm_inv", 27.5 GB dense; lean load -> ~41 GB PEAK host RAM).
# Output corr_deter_{local,disp}_axial_L4 (s1,s2,s3 / sp,tp + ylm; single complex channel Apm).  GPU1, OMP=4.
#
# *** RUN ONLY AFTER the L=2 exact1 run (on GPU1) finishes: that frees GPU1 AND releases its ~15 GB host RAM.
#     L=4 load peaks ~41 GB on a 62 GB box -- needs that RAM back to avoid OOM/swap.  Reads 27.5 GB from
#     disk per binary -> several min I/O each. ***
#
# RERUN NOTE (no rm in this script): delete stale outputs yourself first if re-running:
#   rm data_free_vmRe0.000000vmIm0.000000/corr_deter_local_axial_L4/corr.0.h5
#   rm data_free_vmRe0.000000vmIm0.000000/corr_deter_disp_axial_L4/corr.0.h5
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES=1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

for SRC in jj_local_axial_deter jj_disp_axial_deter; do
  BIN=${SRC}_L4.o
  LOG=${SRC}_L4_claude.log
  echo "### compile ${SRC} L=4  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
  $NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=4 ${SRC}_claude.cu -o "$BIN" \
    2>&1 | tee -a "$LOG" || { echo "BUILD FAILED ${SRC} L=4"; exit 1; }
  echo "### run ${SRC} SUMMED L=4 (no --ins; lattice prop ~27.5 GB, n-t0=2)  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
  ./"$BIN" --n-t0 2 2>&1 | tee -a "$LOG" || { echo "RUN FAILED ${SRC} L=4"; exit 1; }
  echo "### done ${SRC} L=4  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
done
echo "### all done -- corr_deter_{local,disp}_axial_L4; re-enable L=4 in comp_trio_jj_axial_claude.ipynb ###"
