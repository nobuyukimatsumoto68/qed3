#!/usr/bin/env bash
# Production jj_corr_dilute on the MASSLESS INTERACTING ensembles, gsq=8.0, nt128 L1, Nf=2,4,6.
# Distribution: Nf2 -> GPU 0 ; Nf4, Nf6 -> GPU 1 (packed via the running MPS daemon).
# Spin x e/o-time dilution (--spin-dilution --time-dilution 2 => 4 patterns/hit), 1 hit/config, n-t0=2.
# Configs: ckpoint_lat.k for k = KMIN, KMIN+STRIDE, ... (NCONF configs); per-config atomic .h5 + resume-skip,
# so re-running continues where it stopped.  Massless -> no --mass args; disc-- optimization is inert here.
# Output: data_Nf<nf>_gsq8...nt128L1_vmRe0.000000vmIm0.000000/corr_dil_nt02_nhits1_s1td2/corr.<k>.h0.h5
# Runtime: ~950 s/config (interacting) x NCONF -> ~37 h per ensemble at NCONF=140 (MPS keeps GPU1's two ~solo).
#
# ONE shared binary (jj_corr_dilute_massint.o) built once, run as 3 instances (read-only exec is fine).
# Do NOT re-run this script's BUILD step while these are live (ETXTBSY); to resume, the binary already exists.
# Read back: jj_massint_nf{2,4,6}_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

KMIN=40; STRIDE=8; NCONF=140
KMAX=$(( KMIN + STRIDE*NCONF ))     # exclusive -> 1160 ; k=40,48,...,1152 (140 configs)
NT0=2

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
BIN=jj_corr_dilute_massint.o

# ---- build once (skip if you are resuming and the binary already exists + is newer than the source) ----
echo "### compile -> $BIN  [$(date +%F_%H:%M:%S)] ###"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_corr_dilute_claude.cu -o "$BIN"
if [ $? -ne 0 ]; then echo "### BUILD FAILED ###"; exit 1; fi

# ---- launch: Nf2 on GPU0 ; Nf4, Nf6 on GPU1 (background) ----
run_one () {   # $1=Nf  $2=GPU
  local NF=$1 GPU=$2
  local ENS="Nf${NF}_gsq8.000000at0.200000nu01.000000nt128L1/"   # bare config dir (ckpoint_lat.k)
  local LOG="jj_massint_nf${NF}_claude.log"
  echo "### Nf=$NF -> GPU $GPU  ens=$ENS  k=$KMIN..$((KMAX-STRIDE)) step $STRIDE ($NCONF configs)  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
  CUDA_VISIBLE_DEVICES=$GPU ./"$BIN" --ens-dir "$ENS" --kmin "$KMIN" --stride "$STRIDE" --kmax "$KMAX" \
    --nhits 1 --n-t0 "$NT0" --spin-dilution --time-dilution 2 2>&1 | tee -a "$LOG"
  echo "### Nf=$NF done (status ${PIPESTATUS[0]})  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
}

run_one 2 0 &
run_one 4 1 &
run_one 6 1 &
wait
echo "### all three massless-interacting jj runs finished  [$(date +%F_%H:%M:%S)] ###"
