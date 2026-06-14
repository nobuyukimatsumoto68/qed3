#!/usr/bin/env bash
# FREE-FIELD full-validation run, MASSLESS (m=0).  GPU 0 (free; the m_F/m_P 140-hit jobs are packed on GPU 1).
# Full obs: exact-K tp/sp/disc + axial tp/sp + local s1,s2,s3 (vector) + local axial + condensates (PS/FS).
# Master-field (n_t0=2) + SPIN x e/o-TIME dilution (--spin-dilution --time-dilution 2 => 4 patterns/hit).
# nhits=140.  Free mode (no --ens-dir => U=1; no --mass args => massless, non-parity path).
# Compare vs the existing MASSLESS deterministic ground truth (L=1) in the validation notebooks
# (corr_deter_exact1[_axial], corr_deter_local[_axial], condensate_deter_L1 -- all already present).
#
# Builds to a DISTINCT binary (jj_corr_dilute_msl.o) so it does NOT clash (ETXTBSY) with the production
# jj_corr_dilute.o or the other validation binaries (mF/mP/cond/msltest).
# Output: data_free_vmRe0.000000vmIm0.000000/corr_dil_nt02_nhits140_s1td2/corr.0.h<h>.h5  (per hit, resumable)
# Runtime: ~13 min/hit (4 patterns) -> 140 hits ~30 h.  Read back: jj_dilute_valid_massless_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES=0
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
LOG=jj_dilute_valid_massless_claude.log

NHITS=140
NT0=2

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

echo "### compile -> jj_corr_dilute_msl.o (distinct binary)  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_corr_dilute_claude.cu -o jj_corr_dilute_msl.o 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}; if [ "$st" -ne 0 ]; then echo "### BUILD FAILED ($st) ###" | tee -a "$LOG"; exit 1; fi

echo "### run FREE massless nhits=$NHITS nt0=$NT0 spin+td2 (GPU0)  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
./jj_corr_dilute_msl.o --nhits "$NHITS" --n-t0 "$NT0" --spin-dilution --time-dilution 2 2>&1 | tee -a "$LOG"
echo "### run status ${PIPESTATUS[0]}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
echo "### output -> data_free_vmRe0.000000vmIm0.000000/corr_dil_nt0${NT0}_nhits${NHITS}_s1td2/ ###" | tee -a "$LOG"
