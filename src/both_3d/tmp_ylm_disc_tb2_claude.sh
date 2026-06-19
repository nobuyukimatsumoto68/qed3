#!/usr/bin/env bash
# DISC-ONLY rerun with --disc-tblock 2 (interval = Nt/tblock = 64).  Each dilution class fills only
# 2 timeslices, spaced by 64, with its own independent Z2 source + its own solve.  The same-source
# self-contraction (see disc_dilution_bias_claude.md) is then pushed to dt = multiples of the period 64,
# i.e. ONLY dt=64 (lattice midpoint) inside 1..127 -- the small-dt disc signal is fully clean.
# interval = 64 classes x 2 spin = 128 solves/config (~4x the old tb8 disc; ~72 s/config).
#
# CONN is already done -- this does NOT touch it.  Output -> data_Nf<nf>_..._/corr_ylm_disc_tb2_nhits1/
# (distinct dir, no overwrite).  Per-config atomic .h5 + resume-skip; safe to re-launch.
# MPS must be up (nvidia-cuda-mps-control -d).  Nf2 -> GPU 0 ; Nf4, Nf6 -> GPU 1 (co-resident under MPS).
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

KMIN=40; STRIDE=8; NCONF=140
KMAX=$(( KMIN + STRIDE*NCONF ))     # exclusive -> 1160 ; k=40,48,...,1152 (140 configs)
NHITS=1
DISC_TB=2                            # <<< interval = Nt/2 = 64 ; bias pushed to dt=64 only

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

DISC=jj_local_ylm_disc_tb2.o         # distinct binary name (no ETXTBSY vs the tb8/tb1 binaries)

echo "### compile -> $DISC  [$(date +%F_%H:%M:%S)] ###"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_local_ylm_disc_stoch_claude.cu -o "$DISC" || { echo "### DISC BUILD FAILED ###"; exit 1; }

ENS () { echo "Nf${1}_gsq8.000000at0.200000nu01.000000nt128L1/"; }

run_disc () {   # $1=Nf  $2=GPU
  local NF=$1 GPU=$2 LOG="ylm_disc_tb2_nf${1}_claude.log"
  echo "### DISC tb2 Nf=$NF -> GPU $GPU  ($NCONF cfg, interval=Nt/2=64)  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
  CUDA_VISIBLE_DEVICES=$GPU ./"$DISC" --ens-dir "$(ENS $NF)" --kmin "$KMIN" --stride "$STRIDE" --kmax "$KMAX" \
    --nhits "$NHITS" --disc-tblock "$DISC_TB" 2>&1 | tee -a "$LOG"
  echo "### DISC tb2 Nf=$NF done (status ${PIPESTATUS[0]})  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
}

echo "### DISC tb2: Nf2 GPU0 ; Nf4, Nf6 GPU1 (MPS)  [$(date +%F_%H:%M:%S)] ###"
run_disc 2 0 &
run_disc 4 1 &
run_disc 6 1 &
wait
echo "### ALL disc tb2 done  [$(date +%F_%H:%M:%S)] ###"
