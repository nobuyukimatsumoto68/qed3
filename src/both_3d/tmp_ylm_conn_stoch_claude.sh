#!/usr/bin/env bash
# Chunk B1: stochastic CONNECTED per-m Y_lm tower (jj_local_ylm_conn_stoch_claude.cu).
# FREE field, L=1, single origin t0=0, spin dilution.  Validates vs the Chunk A deterministic
# ground truth corr_deter_local_ylmA_L1 (h0/t0_0/ylm/s{a}/l{l}/m{m}): stoch (hit-averaged) -> determ.
#
# Usage: bash tmp_ylm_conn_stoch_claude.sh [GPU] [NHITS]   (default GPU 0, NHITS 16).
#   A first NHITS=16 check is enough to see stoch ~ determ; the full validation uses NHITS=140.
# Per-hit atomic + complete-gated (resume-safe).  Output:
#   data_free_vmRe0.000000vmIm0.000000/corr_ylm_conn_t00_nhits<NHITS>_s1/corr.0.h<h>.h5
# Read back: jj_ylm_conn_stoch_claude.log .  Then jj_local_ylm_validate_claude.ipynb (Chunk C).
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES="${1:-0}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
NHITS="${2:-16}"
BIN=jj_local_ylm_conn_stoch.o
LOG=jj_ylm_conn_stoch_claude.log

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

echo "### compile -> $BIN  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_local_ylm_conn_stoch_claude.cu -o "$BIN" 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}; if [ "$st" -ne 0 ]; then echo "### BUILD FAILED ($st) ###" | tee -a "$LOG"; exit 1; fi

echo "### run FREE L=1, t0=0, spin dilution, nhits=$NHITS  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
./"$BIN" --nhits "$NHITS" --t0 0 --spin-dilution 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}; if [ "$st" -ne 0 ]; then echo "### RUN FAILED ($st) ###" | tee -a "$LOG"; exit 1; fi

echo "### done -> data_free_vmRe0.000000vmIm0.000000/corr_ylm_conn_t00_nhits${NHITS}_s1/  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
echo "### CHECK in jj_local_ylm_validate: hit-avg stoch g^a_{lm} -> determ corr_deter_local_ylmA_L1 (t0_0) ###" | tee -a "$LOG"
