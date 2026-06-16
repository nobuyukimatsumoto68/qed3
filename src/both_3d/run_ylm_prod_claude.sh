#!/usr/bin/env bash
# PRODUCTION local-current Y_lm tower on the MASSLESS INTERACTING ensembles (gsq=8.0, nt128 L1, Nf=2,4,6).
# Two programs:
#   CONN (jj_local_ylm_conn_stoch.o): per-m vector+axial tower, single t0=0 wall source, spin dilution.
#   DISC (jj_local_ylm_disc_stoch.o): per-m disc one-point J^a_{lm} + sigma_PS, time+spin dilution (tb=8).
# Massless -> no --mass args (vector/axial, no parity).  Physical tower (analysis) = -C_conn + C_disc.
#
# Distribution (mirror run_jj_massint): Nf2 -> GPU 0 ; Nf4, Nf6 -> GPU 1, co-resident under the running
# CUDA MPS daemon.  TWO PHASES (conn, then disc) so each GPU has at most 2 procs at once (12 GB TITAN V).
# Configs: ckpoint_lat.k, k = KMIN, KMIN+STRIDE, ... (NCONF configs); per-config atomic .h5 + resume-skip.
# Output: data_Nf<nf>_gsq8...nt128L1_vmRe0.000000vmIm0.000000/{corr_ylm_conn_t00_nhits1_s1,corr_ylm_disc_tb8_nhits1}/
# Read back: ylm_{conn,disc}_nf{2,4,6}_claude.log .   (MPS must be up: `nvidia-cuda-mps-control -d`.)
#
# nhits=1/config (stats from the 140 configs).  Disc is intrinsically noisy -- bump --nhits if needed.
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

KMIN=40; STRIDE=8; NCONF=140
KMAX=$(( KMIN + STRIDE*NCONF ))     # exclusive -> 1160 ; k=40,48,...,1152 (140 configs)
NHITS=1
DISC_TB=8

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

CONN=jj_local_ylm_conn_stoch.o
DISC=jj_local_ylm_disc_stoch.o

# ---- build once (no jobs running now; guarantees current vector+axial code) ----
echo "### compile -> $CONN, $DISC  [$(date +%F_%H:%M:%S)] ###"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_local_ylm_conn_stoch_claude.cu -o "$CONN" || { echo "### CONN BUILD FAILED ###"; exit 1; }
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_local_ylm_disc_stoch_claude.cu -o "$DISC" || { echo "### DISC BUILD FAILED ###"; exit 1; }

ENS () { echo "Nf${1}_gsq8.000000at0.200000nu01.000000nt128L1/"; }

run_conn () {   # $1=Nf  $2=GPU
  local NF=$1 GPU=$2 LOG="ylm_conn_nf${1}_claude.log"
  echo "### CONN Nf=$NF -> GPU $GPU  ($NCONF cfg, k=$KMIN..$((KMAX-STRIDE)) step $STRIDE)  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
  CUDA_VISIBLE_DEVICES=$GPU ./"$CONN" --ens-dir "$(ENS $NF)" --kmin "$KMIN" --stride "$STRIDE" --kmax "$KMAX" \
    --nhits "$NHITS" --t0 0 --spin-dilution 2>&1 | tee -a "$LOG"
  echo "### CONN Nf=$NF done (status ${PIPESTATUS[0]})  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
}

run_disc () {   # $1=Nf  $2=GPU
  local NF=$1 GPU=$2 LOG="ylm_disc_nf${1}_claude.log"
  echo "### DISC Nf=$NF -> GPU $GPU  ($NCONF cfg, tblock=$DISC_TB)  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
  CUDA_VISIBLE_DEVICES=$GPU ./"$DISC" --ens-dir "$(ENS $NF)" --kmin "$KMIN" --stride "$STRIDE" --kmax "$KMAX" \
    --nhits "$NHITS" --disc-tblock "$DISC_TB" 2>&1 | tee -a "$LOG"
  echo "### DISC Nf=$NF done (status ${PIPESTATUS[0]})  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
}

# ---- PHASE 1: CONN (the heavy tower).  Nf2 GPU0 ; Nf4, Nf6 GPU1 (<=2 procs/GPU under MPS) ----
echo "### PHASE 1: CONN  [$(date +%F_%H:%M:%S)] ###"
run_conn 2 0 &
run_conn 4 1 &
run_conn 6 1 &
wait
echo "### PHASE 1 (conn) done  [$(date +%F_%H:%M:%S)] ###"

# ---- PHASE 2: DISC (lighter).  Same distribution ----
echo "### PHASE 2: DISC  [$(date +%F_%H:%M:%S)] ###"
run_disc 2 0 &
run_disc 4 1 &
run_disc 6 1 &
wait
echo "### ALL ylm production done  [$(date +%F_%H:%M:%S)] ###"
