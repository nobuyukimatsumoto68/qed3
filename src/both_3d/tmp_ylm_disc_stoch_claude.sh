#!/usr/bin/env bash
# Chunk B2: stochastic DISCONNECTED per-m Y_lm one-point loops + sigma_PS (jj_local_ylm_disc_stoch_claude.cu).
# FREE field, L=1, time+spin dilution (disc_tblock=8 => 16 time x 2 spin classes/hit).
# Validation (Chunk C notebook): hit-averaged J^a_{lm} -> 0 for l>=1 (=> C_disc -> 0); sigma_PS density
#   etadag_xi/(Nt*sum A_n) -> -1, sigma_PS -> -2 (cf. dense condensate_deter_free).
#
# Usage: bash tmp_ylm_disc_stoch_claude.sh [GPU] [NHITS]   (default GPU 0, NHITS 16).
# Per-hit atomic + complete-gated.  Output:
#   data_free_vmRe0.000000vmIm0.000000/corr_ylm_disc_tb8_nhits<NHITS>/corr.0.h<h>.h5
# Read back: jj_ylm_disc_stoch_claude.log .
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES="${1:-0}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
NHITS="${2:-16}"
BIN=jj_local_ylm_disc_stoch.o
LOG=jj_ylm_disc_stoch_claude.log

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

echo "### compile -> $BIN  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_local_ylm_disc_stoch_claude.cu -o "$BIN" 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}; if [ "$st" -ne 0 ]; then echo "### BUILD FAILED ($st) ###" | tee -a "$LOG"; exit 1; fi

echo "### run FREE L=1, disc_tblock=8, nhits=$NHITS  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
./"$BIN" --nhits "$NHITS" --disc-tblock 8 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}; if [ "$st" -ne 0 ]; then echo "### RUN FAILED ($st) ###" | tee -a "$LOG"; exit 1; fi

echo "### done -> data_free_vmRe0.000000vmIm0.000000/corr_ylm_disc_tb8_nhits${NHITS}/  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
