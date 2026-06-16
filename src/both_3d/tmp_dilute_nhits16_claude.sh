#!/usr/bin/env bash
# Quick FREE spin-only run (nhits=16) of jj_corr_dilute WITH the new local-axial code, to validate local axial.
# Builds to a DISTINCT binary (jj_corr_dilute_la.o) so it does NOT clash with the interacting run's
# jj_corr_dilute.o (which is executing on GPU1 -- recompiling that file would fail with "text file busy").
# GPU0, OMP=4.  Output: data_free_*/corr_dil_nt02_nhits16_s1td1/ (fresh; has h0/axial/s{1,2,3}/Apm now).
# ~16 hits x 2 patterns x ~196 s ~ 1.7 h.  Read back: jj_dilute_nhits16_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES=0
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
LOG=jj_dilute_nhits16_claude.log

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

echo "### compile -> jj_corr_dilute_la.o (distinct binary)  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_corr_dilute_claude.cu -o jj_corr_dilute_la.o 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}
if [ "$st" -ne 0 ]; then echo "### BUILD FAILED (status $st) ###" | tee -a "$LOG"; exit 1; fi

echo "### run FREE spin-only nhits=16 (GPU0)  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
./jj_corr_dilute_la.o --nhits 16 --n-t0 2 --spin-dilution --time-dilution 1 2>&1 | tee -a "$LOG"
echo "### run status ${PIPESTATUS[0]}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
echo "### output -> data_free_*/corr_dil_nt02_nhits16_s1td1/  (set DIL there in jj_dilute_validate_axial) ###" | tee -a "$LOG"
