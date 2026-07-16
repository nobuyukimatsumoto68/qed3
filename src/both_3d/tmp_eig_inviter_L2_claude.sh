#!/usr/bin/env bash
# tmp_eig_inviter_L2_claude.sh
# _claude: INVERSE-ITERATION ground-truth check (independent of the Chebyshev/IRL).  Applies A^{-1} (multishift
# solve on the SAME frozen-window Dov the IRL uses) repeatedly -> converges to the SMALLEST eigenvalue of
# A = D_m^dag D_m.  If lambda -> ~0.05 the small mode is REAL (IRL is buggy); if it stays >1 the operator has
# no small modes.  Runs on thermalized ckpoint_lat.740 (Wilson lambda_min=0.0643, INSIDE frozen window 0.06).
# Single process on GPU 1.  Reads back: eig_inviter_L2_claude.log .
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

GPU=1
LREF=2
NT=128
CONFIG="Nf2_gsq8.000000at0.200000nu01.000000nt128L2/ckpoint_lat.740"

export CONFIG_LAT="$CONFIG"
export INV_ITER=1                 # driver does inverse iteration + returns (no IRL)

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

LOG=eig_inviter_L2_claude.log

echo "### INVERSE ITERATION L2  Nt=${NT} (N=$((2*(10*LREF*LREF+2)*NT)))  config=${CONFIG}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"

BIN="eig_lanczos_L${LREF}_nt${NT}.o"
echo "### BUILD -> $BIN ###" | tee -a "$LOG"
$NVCC $NVCCBASE -DLREF=${LREF} -DNT_CLI=${NT} $INCLUDES $LDFLAGS eig_lanczos_claude.cu -o "$BIN" 2>>"$LOG" \
  || { echo "### BUILD FAILED (see $LOG) ###" | tee -a "$LOG"; exit 1; }

echo "### RUN inverse iteration (mass=0)  GPU$GPU  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
CUDA_VISIBLE_DEVICES=$GPU ./"$BIN" 0.0 0.0 2>&1 | tee -a "$LOG" \
  || { echo "### RUN FAILED ###" | tee -a "$LOG"; exit 1; }

echo "### inverse iteration done  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
