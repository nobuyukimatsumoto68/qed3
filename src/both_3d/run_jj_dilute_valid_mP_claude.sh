#!/usr/bin/env bash
# FREE-FIELD full-validation run, PARITY-breaking mass m_P = 0.1 i (imaginary).  GPU 1 (PACKED with m_F).
# Full obs: exact-K tp/sp/disc + axial tp/sp + local s1,s2,s3 (vector) + local axial + condensates (PS/FS).
# Master-field (n_t0=2) + SPIN x e/o-TIME dilution (--spin-dilution --time-dilution 2 => 4 patterns/hit).
# nhits=140.  Free mode (no --ens-dir => U=1).  m_P uses the tilde-D dagger leg (Eq. 3.61) -> the dilute's
# parity path; the condensate xidag_eta / xidag_1mDdag_eta take their m_P forms.  Compare vs the
# DETERMINISTIC ground truth at m=0.1i (L=1) -- NOTE: that determ must be generated at this mass first.
#
# Builds to a DISTINCT binary (jj_corr_dilute_mP.o) so it does NOT clash (ETXTBSY) with the production
# jj_corr_dilute.o or the condensate-check jj_corr_dilute_cond.o.
# Output: data_free_vmRe0.000000vmIm0.100000/corr_dil_nt02_nhits140_s1td2/corr.0.h<h>.h5  (per hit, resumable)
# PACKING (jj_dilute_multijob_local_handoff_claude.md): shares GPU 1 with the m_F job (disjoint output dirs
#   + distinct binaries -> safe; one process only fills ~50-60% of the GPU).  Launch both in the background.
# Runtime: ~free hit (4 patterns; parity does extra tilde legs) ~20 min -> 140 hits ~45 h (slower when packed).
# Read back: jj_dilute_valid_mP_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES=1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
LOG=jj_dilute_valid_mP_claude.log

NHITS=140
NT0=2
MRE=0.0
MIM=0.1                            # parity-breaking imaginary mass m_P = 0.1 i

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

echo "### compile -> jj_corr_dilute_mP.o (distinct binary)  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_corr_dilute_claude.cu -o jj_corr_dilute_mP.o 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}; if [ "$st" -ne 0 ]; then echo "### BUILD FAILED ($st) ###" | tee -a "$LOG"; exit 1; fi

echo "### run FREE m_P=$MIM i  nhits=$NHITS nt0=$NT0  spin+td2 (GPU1)  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
./jj_corr_dilute_mP.o --nhits "$NHITS" --n-t0 "$NT0" --spin-dilution --time-dilution 2 \
  --mass-re "$MRE" --mass-im "$MIM" 2>&1 | tee -a "$LOG"
echo "### run status ${PIPESTATUS[0]}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
echo "### output -> data_free_vmRe0.000000vmIm0.100000/corr_dil_nt0${NT0}_nhits${NHITS}_s1td2/ ###" | tee -a "$LOG"
