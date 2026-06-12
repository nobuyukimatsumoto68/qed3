#!/usr/bin/env bash
# Build + run the same-site spin-diagonality probe (test_samesite_spindiag_claude.cu): does the LOCAL
# s1==s2 property (same-site propagator block spin-diagonal) survive a Gaussian gauge field?
# Cheap: 4 gauge cases (U=1 + Gaussian w=0.1,0.3,1.0) x 2 point-source solves at L=1.
# Usage:  bash tmp_spindiag_claude.sh [GPU]    (default GPU 0)   ; read back jj_spindiag_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
GPU="${1:-0}"
export CUDA_VISIBLE_DEVICES=$GPU
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
LOG=jj_spindiag_claude.log

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

echo "### compile test_samesite_spindiag_claude.cu [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS test_samesite_spindiag_claude.cu -o test_samesite_spindiag.o 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}
if [ "$st" -ne 0 ]; then echo "### BUILD FAILED (status $st) ###" | tee -a "$LOG"; exit 1; fi
echo "### run [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
./test_samesite_spindiag.o 2>&1 | tee -a "$LOG"
echo "### run status ${PIPESTATUS[0]} [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
