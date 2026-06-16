#!/usr/bin/env bash
# Production runner for jj_corr_dilute_claude.cu (exact-K + local current; spin x time DILUTION + master-field
# source-origin SUPERPOSITION).  Compiles then runs; FREE (U=1) by default, or interacting with ENSDIR.
# Output: data_<ESNID>/corr_dil_nt0<N>_nhits<H>_s<0|1>td<TD>/corr.<k>.h<h>.h5 (per hit; atomic + resume-skip).
#
# Usage:
#   bash run_jj_dilute_claude.sh [GPU] [NHITS] [NT0] [ENSDIR] [SPIN_DIL] [TD]
#     GPU      : CUDA device (default 0)
#     NHITS    : stochastic hits (default 16)
#     NT0      : superposed origins (default 2 -> t0={0,Nt/2})
#     ENSDIR   : sea config dir (ends with '/'); OMIT or "" => FREE field (U=1) validation
#     SPIN_DIL : 1 => --spin-dilution ON (NS=2 spin classes); 0 => OFF (default 0)
#     TD       : --time-dilution (number of interval time classes; require (Nt/2)%TD==0; default 2)
#   => dilution patterns per hit = (SPIN_DIL? 2:1) x TD.  L=1 cost ~ (#patterns) x ~196 s/pattern.
#   Read back: jj_dilute_run_claude.log
#
# SPIN-ONLY (recommended for the connected JJ): SPIN_DIL=1 TD=1 -> 2 patterns -> dir ..._s1td1.
# Per-hit files are atomic + complete-gated, so the run is resumable (rerun the same command to continue).
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
GPU="${1:-0}"
NHITS="${2:-16}"
NT0="${3:-2}"
ENSDIR="${4:-}"
SPIN_DIL="${5:-0}"
TD="${6:-2}"
export CUDA_VISIBLE_DEVICES=$GPU
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
LOG=jj_dilute_run_claude.log

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

echo "### compile jj_corr_dilute_claude.cu  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_corr_dilute_claude.cu -o jj_corr_dilute.o 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}
if [ "$st" -ne 0 ]; then echo "### BUILD FAILED (status $st) ###" | tee -a "$LOG"; exit 1; fi

ENSARG=""
if [ -n "$ENSDIR" ]; then ENSARG="--ens-dir $ENSDIR"; fi
SPINARG=""
if [ "$SPIN_DIL" -eq 1 ]; then SPINARG="--spin-dilution"; fi
echo "### run  GPU=$GPU NHITS=$NHITS NT0=$NT0 SPIN_DIL=$SPIN_DIL TD=$TD  ${ENSARG:-[FREE U=1]}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
./jj_corr_dilute.o --nhits "$NHITS" --n-t0 "$NT0" --time-dilution "$TD" $SPINARG $ENSARG 2>&1 | tee -a "$LOG"
echo "### run status ${PIPESTATUS[0]}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
echo "### output -> data_<ESNID>/corr_dil_nt0${NT0}_nhits${NHITS}_s${SPIN_DIL}td${TD}/  (validate notebooks) ###" | tee -a "$LOG"
