#!/usr/bin/env bash
# PRODUCTION driver: jj_corr_dilute_claude.cu on the INTERACTING Nf=2 sea ensemble (gsq=2, L=1), with
# the developed dilution code.  Spin x e/o-time dilution (--spin-dilution --time-dilution 2 => 4 patterns/hit),
# master-field superposition, exact-K (tp/sp/disc/axial) + local currents + condensates (PS/FS).
#
# Sweep: configs ckpoint_lat.k for k = KMIN, KMIN+STRIDE, ... (NCONF configs).  Here KMIN=40, STRIDE=8,
#        NCONF=140 => k=40,48,...,1152 ; KMAX=KMIN+STRIDE*NCONF=1160 (exclusive).  nhits=1 per config.
# Output: data_Nf2_gsq2...nt128L1_vmRe0.000000vmIm0.000000/corr_dil_nt0<NT0>_nhits1_s1td2/corr.<k>.h0.h5
#   (per config, atomic + complete-gated -> RESUMABLE: rerun the same command to continue where it stopped).
# Read back: jj_dilute_prod_claude.log
#
# Usage:  bash run_jj_dilute_prod_claude.sh [GPU]      (default GPU 1)
#
# NOTE on dilution: spin-only (TD=1) is the recommended per-cost config for the connected JJ.  For the
#   4-pattern (spin x e/o) instead, set TD=2 below.
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1

# ----- ensemble + sweep -----
ENS_DIR="Nf2_gsq2.000000at0.200000nu01.000000nt128L1/"   # bare dir (no data_ prefix); holds ckpoint_lat.k
KMIN=40
STRIDE=8
NCONF=140
KMAX=$(( KMIN + STRIDE*NCONF ))    # exclusive -> 1160 ; k=40,48,...,1152 (140 configs)
NHITS=1                            # interacting: 1 hit/config (statistics come from the configs)
NT0=2                              # master-field superposed origins t0={0,Nt/2}
SPIN_DIL=1                         # --spin-dilution ON
# TD=1                             # was: --time-dilution 1 => spin-only (2 patterns)
TD=2                               # --time-dilution 2 => spin x e/o time (4 patterns); also lowers condensate variance

GPU="${1:-1}"
export CUDA_VISIBLE_DEVICES=$GPU
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
LOG=jj_dilute_prod_claude.log

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

echo "### compile jj_corr_dilute_claude.cu  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_corr_dilute_claude.cu -o jj_corr_dilute.o 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}
if [ "$st" -ne 0 ]; then echo "### BUILD FAILED (status $st) ###" | tee -a "$LOG"; exit 1; fi

SPINARG=""; if [ "$SPIN_DIL" -eq 1 ]; then SPINARG="--spin-dilution"; fi
echo "### run  GPU=$GPU  ens=$ENS_DIR  k=$KMIN..$((KMAX-STRIDE)) step $STRIDE ($NCONF configs)" \
     "nhits=$NHITS nt0=$NT0 spin=$SPIN_DIL td=$TD  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
./jj_corr_dilute.o --ens-dir "$ENS_DIR" --kmin "$KMIN" --stride "$STRIDE" --kmax "$KMAX" \
  --nhits "$NHITS" --n-t0 "$NT0" --time-dilution "$TD" $SPINARG 2>&1 | tee -a "$LOG"
echo "### run status ${PIPESTATUS[0]}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
echo "### output -> data_Nf2_gsq2...nt128L1_vmRe0.000000vmIm0.000000/corr_dil_nt0${NT0}_nhits${NHITS}_s${SPIN_DIL}td${TD}/ ###" | tee -a "$LOG"
