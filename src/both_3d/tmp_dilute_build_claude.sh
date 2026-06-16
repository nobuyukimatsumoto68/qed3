#!/usr/bin/env bash
# CHECKPOINT for jj_corr_dilute_claude.cu after Chunk 1 (header + seed_from_string + dir_out tag) and
# Chunk 1b (ylm stripped via N_ELL=0).  Confirms the renamed/ylm-stripped base BUILDS, RUNS, and writes
# exact tp/sp/disc/axial with NO ylm keys to data_free_.../corr_dil_nt02_nhits1/.
# It is NOT yet diluted / superposed / local -- those are the next chunks (3 superposition, 3L local).
#
# Usage:  bash tmp_dilute_build_claude.sh [GPU]     (default GPU 0)
# Read back: jj_corr_dilute_build_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
GPU="${1:-0}"
export CUDA_VISIBLE_DEVICES=$GPU
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
LOG=jj_corr_dilute_build_claude.log

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

echo "### compile jj_corr_dilute_claude.cu (Chunk 1+1b) [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_corr_dilute_claude.cu -o jj_corr_dilute.o 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}
if [ "$st" -ne 0 ]; then echo "### BUILD FAILED (status $st) ###" | tee -a "$LOG"; exit 1; fi

echo "### build OK; FREE L=1 run (nhits=1, n-t0=2) [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
./jj_corr_dilute.o --nhits 1 --n-t0 2 2>&1 | tee -a "$LOG"
echo "### run status ${PIPESTATUS[0]} [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
echo "### output -> data_free_vmRe0.000000vmIm0.000000/corr_dil_nt02_nhits1/corr.0.h0.h5 ###" | tee -a "$LOG"
