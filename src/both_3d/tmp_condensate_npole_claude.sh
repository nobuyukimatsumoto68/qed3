#!/usr/bin/env bash
# Contact-term probe: dense free condensate sigma_PS/sigma_FS at HIGH Zolotarev pole count (n=51) vs the
# current n=11, for all three mass cases (m0, m_F=0.1, m_P=0.1i), CONTACT-SUBTRACTED (md Sec 10).
# Goal: the symmetry-protected residuals should be ~0; whatever residual does NOT shrink from n=11 -> n=51
# is a genuine lattice O(m) (term-2) effect, not a finite-Zolotarev sign-function artifact.
#
# Builds a DISTINCT binary (condensate_deter_npole.o) and writes to pole-tagged dirs
# (condensate_deter_n{11,51}_L1/) so it never clobbers the existing n=11 reference (condensate_deter_L1/).
# The contact-subtracted numbers are PRINTED to the log -- that is the deliverable.
# Usage: bash tmp_condensate_npole_claude.sh [GPU]   (default 0).  Read back: condensate_npole_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES="${1:-0}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
LOG=condensate_npole_claude.log

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

echo "### compile -> condensate_deter_npole.o  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS condensate_deter_free_npole_claude.cu -o condensate_deter_npole.o 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}; if [ "$st" -ne 0 ]; then echo "### BUILD FAILED ($st) ###" | tee -a "$LOG"; exit 1; fi

# mass cases: m0 (0,0), mF (0.1,0), mP (0,0.1).  pole counts: 11 (current) and 51 (high-precision).
for NP in 11 51; do
  for MC in "m0 0.0 0.0" "mF 0.1 0.0" "mP 0.0 0.1"; do
    set -- $MC; NAME=$1; MRE=$2; MIM=$3
    echo "### run $NAME (m=($MRE,$MIM))  npole=$NP  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
    ./condensate_deter_npole.o --mass-re "$MRE" --mass-im "$MIM" --npole "$NP" --out-tag "_n${NP}" 2>&1 | tee -a "$LOG"
    st=${PIPESTATUS[0]}; if [ "$st" -ne 0 ]; then echo "### RUN FAILED ($NAME n=$NP, $st) ###" | tee -a "$LOG"; exit 1; fi
  done
done

echo "### done.  Compare the 'CONTACT-SUBTRACTED' lines: n=11 vs n=51 per mass.  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
echo "### (residual that is unchanged n11->n51 = genuine lattice O(m); residual that shrinks = Zolotarev) ###" | tee -a "$LOG"
