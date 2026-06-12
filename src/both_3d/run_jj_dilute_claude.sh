#!/usr/bin/env bash
# Production runner for jj_corr_dilute_claude.cu (exact-K + local current; 4-pattern spin x e/o DILUTION;
# master-field source-origin SUPERPOSITION).  Compiles then runs; FREE (U=1) by default, or interacting
# with --ens-dir.  Output: data_<ESNID>/corr_dil_nt0<N>_nhits<H>/corr.<k>.h<h>.h5 (per hit; t0s + rng_seed).
#
# Usage:
#   bash run_jj_dilute_claude.sh [GPU] [NHITS] [NT0] [ENSDIR]
#     GPU    : CUDA device (default 0)
#     NHITS  : stochastic hits (default 16).  Each hit = 4 dilution patterns; at L=1 ~13 min/hit.
#     NT0    : superposed origins (default 2 -> t0={0,Nt/2}).
#     ENSDIR : sea config dir (ends with '/'); OMIT or "" => FREE field (U=1) validation.
#   Read back: jj_dilute_run_claude.log
#
# FAIR-COST NOTE: 1 diluted hit (4 patterns) ~ 4 plain (volume) hits.  To compare variance vs the volume
# run corr_nt02_nhits64 at matched cost, use NHITS ~ 16 (= 64 volume-hit equivalent).
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
GPU="${1:-0}"
NHITS="${2:-16}"
NT0="${3:-2}"
ENSDIR="${4:-}"
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
echo "### run  GPU=$GPU NHITS=$NHITS NT0=$NT0  ${ENSARG:-[FREE U=1]}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
./jj_corr_dilute.o --nhits "$NHITS" --n-t0 "$NT0" $ENSARG 2>&1 | tee -a "$LOG"
echo "### run status ${PIPESTATUS[0]}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
echo "### output -> data_<ESNID>/corr_dil_nt0${NT0}_nhits${NHITS}/  (validate: jj_dilute_validate_claude.ipynb) ###" | tee -a "$LOG"
