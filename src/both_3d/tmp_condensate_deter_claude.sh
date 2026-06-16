#!/usr/bin/env bash
# Free-field cross-check for the condensate (Eqs. 1.23 PS, 1.55 FS).  Two phases:
#   (1) DENSE EXACT reference: condensate_deter_free_claude.cu (complete unit-source basis) ->
#       data_free_vmRe0.000000vmIm0.000000/condensate_deter_L1/cond.h5
#   (2) STOCHASTIC: jj_corr_dilute_claude.cu in FREE mode (omit --ens-dir), built to a DISTINCT binary
#       (jj_corr_dilute_cond.o) so it does NOT clash with the live production jj_corr_dilute.o (ETXTBSY) ->
#       data_free_vmRe0.000000vmIm0.000000/corr_dil_nt02_nhits${NHITS}_s1td2/  (h0/condensate/* per hit)
# Compare in condensate_validate_claude.ipynb (jackknife stoch over hits, /Vst average, vs dense exact).
# Massless run (mass=0).  For m_F / m_P: re-run BOTH phases with --mass-re / --mass-im set to the sea mass.
# GPU0, OMP=4.  Read back: condensate_deter_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES=0
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
LOG=condensate_deter_claude.log

NHITS=${1:-16}                 # stochastic free hits (default 16); bump for tighter errors
MRE=${2:-0.0}                  # sea mass (real = m_F)
MIM=${3:-0.0}                  # sea mass (imag = m_P)

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

echo "### (1) compile dense reference -> condensate_deter.o  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS condensate_deter_free_claude.cu -o condensate_deter.o 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}; if [ "$st" -ne 0 ]; then echo "### BUILD FAILED ($st) ###" | tee -a "$LOG"; exit 1; fi

echo "### (1) run dense exact reference (mass=($MRE,$MIM))  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
./condensate_deter.o --mass-re "$MRE" --mass-im "$MIM" 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}; if [ "$st" -ne 0 ]; then echo "### DENSE RUN FAILED ($st) ###" | tee -a "$LOG"; exit 1; fi

echo "### (2) compile dilute -> jj_corr_dilute_cond.o (DISTINCT; avoid prod ETXTBSY)  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_corr_dilute_claude.cu -o jj_corr_dilute_cond.o 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}; if [ "$st" -ne 0 ]; then echo "### BUILD FAILED ($st) ###" | tee -a "$LOG"; exit 1; fi

echo "### (2) run STOCHASTIC free dilute nhits=$NHITS (mass=($MRE,$MIM))  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
./jj_corr_dilute_cond.o --nhits "$NHITS" --n-t0 2 --spin-dilution --time-dilution 2 \
  --mass-re "$MRE" --mass-im "$MIM" 2>&1 | tee -a "$LOG"
echo "### run status ${PIPESTATUS[0]}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
echo "### dense -> condensate_deter_L1/cond.h5 ; stoch -> corr_dil_nt02_nhits${NHITS}_s1td2/ ; compare in condensate_validate_claude.ipynb ###" | tee -a "$LOG"
