#!/usr/bin/env bash
# PRODUCTION local-current Y_lm tower on the L=2 MASSLESS INTERACTING ensembles (gsq=8.0, nt128 L2,
# Nf=2,4,6).  L2 variant of run_ylm_prod_claude.sh.  Differences vs L1:
#   - compile with -DN_REFINE_CLI=2  (sites 10*4+2=42, links 120)
#   - ensembles  Nf<nf>_..._nt128L2/  (massless sea: Nf2 k=1..1601, Nf4 k=1..518, Nf6 k=1..289)
#   - CONFIG SAMPLING INTERVAL changed: KMIN=40, STRIDE=4, NCONF per chain (Nf2=140, Nf4=120, Nf6=63)
#   - DISC --disc-tblock 2 unchanged (interval Nt/2=64; bias fix is temporal, L-independent)
# Output AUTO-TAGS from --ens-dir -> data_Nf<nf>_..._nt128L2_vmRe0.000000vmIm0.000000/
#   {corr_ylm_conn_t00_nhits1_s1, corr_ylm_disc_tb2_nhits1}/  (per-config atomic .h5 + resume-skip).
# Two programs (conn = heavy vector+axial tower, t0=0 wall + spin dilution; disc = one-point + sigma_PS,
# time+spin dilution).  Physical tower (analysis) = -C_conn + C_disc.
# ALL ON GPU 0 (dev=1 not used).  Nf list = positional args (default "2 4"); Nf6 run later, e.g.
#   bash run_ylm_prod_L2_claude.sh 6 .   Listed Nf run SIMULTANEOUSLY on GPU 0, co-resident under MPS.
# MPS must be up:  nvidia-cuda-mps-control -d .   Read back: ylm_{conn,disc}_L2_nf{2,4,6}_claude.log .
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

GPU=0
NF_LIST="${*:-2 4}"

KMIN=40
STRIDE=4
NHITS=1
DISC_TB=2

# per-Nf number of configs (chain lengths: Nf2=1601, Nf4=518, Nf6=289)
nconf () {   # $1=Nf -> echoes NCONF
  case "$1" in
    2) echo 140 ;;
    4) echo 120 ;;
    6) echo 63  ;;
    *) echo 0   ;;
  esac
}

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -DN_REFINE_CLI=2 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

CONN=jj_local_ylm_conn_stoch_L2.o
DISC=jj_local_ylm_disc_stoch_L2.o

# ---- build once (distinct L2 binaries; no clobber of the L1 .o) ----
echo "### compile L2 (-DN_REFINE_CLI=2) -> $CONN, $DISC  [$(date +%F_%H:%M:%S)] ###"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_local_ylm_conn_stoch_claude.cu -o "$CONN" || { echo "### CONN BUILD FAILED ###"; exit 1; }
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_local_ylm_disc_stoch_claude.cu -o "$DISC" || { echo "### DISC BUILD FAILED ###"; exit 1; }

ENS () { echo "Nf${1}_gsq8.000000at0.200000nu01.000000nt128L2/"; }

run_conn () {   # $1=Nf  $2=GPU
  local NF=$1 GPU=$2 LOG="ylm_conn_L2_nf${1}_claude.log"
  local NC KMAX
  NC=$(nconf "$NF")
  KMAX=$(( KMIN + STRIDE*NC ))     # exclusive
  echo "### CONN L2 Nf=$NF -> GPU $GPU  ($NC cfg, k=$KMIN..$((KMAX-STRIDE)) step $STRIDE)  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
  CUDA_VISIBLE_DEVICES=$GPU ./"$CONN" --ens-dir "$(ENS $NF)" --kmin "$KMIN" --stride "$STRIDE" --kmax "$KMAX" \
    --nhits "$NHITS" --t0 0 --spin-dilution 2>&1 | tee -a "$LOG"
  echo "### CONN L2 Nf=$NF done (status ${PIPESTATUS[0]})  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
}

run_disc () {   # $1=Nf  $2=GPU
  local NF=$1 GPU=$2 LOG="ylm_disc_L2_nf${1}_claude.log"
  local NC KMAX
  NC=$(nconf "$NF")
  KMAX=$(( KMIN + STRIDE*NC ))
  echo "### DISC L2 tb2 Nf=$NF -> GPU $GPU  ($NC cfg, interval=Nt/2=64)  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
  CUDA_VISIBLE_DEVICES=$GPU ./"$DISC" --ens-dir "$(ENS $NF)" --kmin "$KMIN" --stride "$STRIDE" --kmax "$KMAX" \
    --nhits "$NHITS" --disc-tblock "$DISC_TB" 2>&1 | tee -a "$LOG"
  echo "### DISC L2 Nf=$NF done (status ${PIPESTATUS[0]})  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
}

# ---- PHASE 1: CONN (the heavy tower).  All listed Nf simultaneously on GPU 0 (<=2 procs under MPS) ----
echo "### PHASE 1: CONN L2  (Nf: $NF_LIST -> GPU $GPU)  [$(date +%F_%H:%M:%S)] ###"
for nf in $NF_LIST; do
  run_conn "$nf" "$GPU" &
done
wait
echo "### PHASE 1 (conn L2) done  [$(date +%F_%H:%M:%S)] ###"

# ---- PHASE 2: DISC (lighter).  Same distribution ----
echo "### PHASE 2: DISC L2  (Nf: $NF_LIST -> GPU $GPU)  [$(date +%F_%H:%M:%S)] ###"
for nf in $NF_LIST; do
  run_disc "$nf" "$GPU" &
done
wait
echo "### ALL ylm L2 production done  [$(date +%F_%H:%M:%S)] ###"
