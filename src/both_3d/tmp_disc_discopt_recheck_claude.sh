#!/usr/bin/env bash
# Re-validate the mP disc-- ride-along optimization (cyclicity: disc-- now rides the connected V-- sink,
# standalone [disc --] sweep commented out -> ~20% faster mP).  The disc-- estimator CHANGED ordering
# (eta^dag K^dag phimm vs the old K(n,t) tilphi), so it is NOT byte-identical to the old run -- the check is
# STATISTICAL: over a few hits the disc-- two-point should still -> 0 (free), and the physical sum /
# V-- mean should match the determ reference, exactly as before the change.
#
# Builds a DISTINCT binary (jj_corr_dilute_discopt.o) so no ETXTBSY with the running validation jobs, and
# writes to the standard nhits16 dir (corr_dil_nt02_nhits16_s1td2) -- DISTINCT from the live nhits140 run.
# Usage: bash tmp_disc_discopt_recheck_claude.sh [GPU] [NHITS]   (default GPU 0, NHITS 16).
#   GPUs are busy with the 140-hit runs; under MPS this co-runs (slower) or run it once one frees up.
# Read back: disc_discopt_recheck_claude.log .  Then in jj_validate_mP_claude.ipynb set
#   MASS='mP', DIL='corr_dil_nt02_nhits16_s1td2' and confirm: disc tp/sp -> ~0, physical sum & V-- vs determ.
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES="${1:-0}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
NHITS="${2:-16}"
LOG=disc_discopt_recheck_claude.log

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

echo "### compile -> jj_corr_dilute_discopt.o (optimized disc--)  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_corr_dilute_claude.cu -o jj_corr_dilute_discopt.o 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}; if [ "$st" -ne 0 ]; then echo "### BUILD FAILED ($st) ###" | tee -a "$LOG"; exit 1; fi

echo "### run FREE m_P=0.1i  nhits=$NHITS  n-t0 2 spin+td2  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
echo "### (per-hit log should NO LONGER show a '[disc --] tilde trace' line; mP/hit ~1330s vs ~1678s) ###" | tee -a "$LOG"
./jj_corr_dilute_discopt.o --nhits "$NHITS" --n-t0 2 --spin-dilution --time-dilution 2 \
  --mass-re 0.0 --mass-im 0.1 2>&1 | tee -a "$LOG"
echo "### run status ${PIPESTATUS[0]}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
echo "### output -> data_free_vmRe0.000000vmIm0.100000/corr_dil_nt02_nhits${NHITS}_s1td2/ ###" | tee -a "$LOG"
echo "### CHECK in jj_validate_mP: DIL='corr_dil_nt02_nhits${NHITS}_s1td2' -> disc tp/sp ~0, sum & V-- match determ ###" | tee -a "$LOG"
