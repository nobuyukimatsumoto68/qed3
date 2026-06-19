#!/usr/bin/env bash
# RE-PACKED local-current Y_lm tower on the MASSIVE (Real-m / m_F) ensembles (gsq=8.0, nt128 L1, Nf=2,4,6).
# Same physics/flags as run_ylm_massive_claude.sh; ONLY the GPU packing changes.
#
# OBSERVATION (nvidia-smi): each proc uses ~380 MiB / 12 GB -> memory is NOT the limit; the runs are
# latency/compute bound.  GPU0 was at 65% with 1 proc, GPU1 at 86% with 2.  So load each GPU harder:
#
# RE-PACK: each GPU owns ONE mass at a time = that mass's 3 Nf (Nf2,Nf4,Nf6) running CONCURRENTLY on the
# same GPU (3 procs ~ 1.1 GB, co-resident under MPS), conn phase then disc phase.  Two masses run AT ONCE
# (one per GPU).  Each GPU processes TWO masses total, HEAVIER first:
#   GPU0: 0.200000 -> 0.050000      GPU1: 0.100000 -> 0.010000
# First parallel wave finishes the two heaviest (0.2, 0.1) together; then 0.05 & 0.01.
#
# Resume-safe: per-config atomic .h5 + resume-skip, SAME output dirs as run_ylm_massive_claude.sh, so this
# picks up exactly where the old (1+2 split) run left off -- kill that one, launch this, nothing is lost.
# Output: data_<ens>_vmRe<m>vmIm0.000000/{corr_ylm_conn_t00_nhits1_s1,corr_ylm_disc_tb2_nhits1}/ .
# MPS must be up: `nvidia-cuda-mps-control -d`.
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

KMIN=40
STRIDE=8
NCONF=140
KMAX=$(( KMIN + STRIDE*NCONF ))      # exclusive -> 1160 ; k=40,48,...,1152 (140 configs)
NHITS=1
DISC_TB=2                            # clean disc (interval = Nt/2 = 64 ; bias only at dt=64)

GPU0_MASSES="0.200000 0.050000"      # GPU0 owns these, heavier first
GPU1_MASSES="0.100000 0.010000"      # GPU1 owns these, heavier first

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

CONN=jj_local_ylm_conn_pack2.o       # distinct binary names (no ETXTBSY vs the running massive binaries)
DISC=jj_local_ylm_disc_pack2.o

# ---- build once (early-exit on failure) ----
echo "### compile -> $CONN, $DISC  [$(date +%F_%H:%M:%S)] ###"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_local_ylm_conn_stoch_claude.cu -o "$CONN" || { echo "### CONN BUILD FAILED ###"; exit 1; }
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_local_ylm_disc_stoch_claude.cu -o "$DISC" || { echo "### DISC BUILD FAILED ###"; exit 1; }

ENS () { echo "Nf${1}_gsq8.000000at0.200000nu01.000000mRe${2}mIm0.000000nt128L1/"; }

run_conn () {   # $1=Nf  $2=GPU  $3=mass
  local NF=$1 GPU=$2 M=$3 LOG="ylm_conn_mRe${3}_nf${1}_claude.log"
  echo "### CONN Nf=$NF mRe=$M -> GPU $GPU  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
  CUDA_VISIBLE_DEVICES=$GPU ./"$CONN" --ens-dir "$(ENS $NF $M)" --kmin "$KMIN" --stride "$STRIDE" --kmax "$KMAX" \
    --nhits "$NHITS" --t0 0 --spin-dilution --mass-re "$M" 2>&1 | tee -a "$LOG"
  echo "### CONN Nf=$NF mRe=$M done (status ${PIPESTATUS[0]})  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
}

run_disc () {   # $1=Nf  $2=GPU  $3=mass
  local NF=$1 GPU=$2 M=$3 LOG="ylm_disc_mRe${3}_nf${1}_claude.log"
  echo "### DISC Nf=$NF mRe=$M -> GPU $GPU  (tblock=$DISC_TB)  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
  CUDA_VISIBLE_DEVICES=$GPU ./"$DISC" --ens-dir "$(ENS $NF $M)" --kmin "$KMIN" --stride "$STRIDE" --kmax "$KMAX" \
    --nhits "$NHITS" --disc-tblock "$DISC_TB" --mass-re "$M" 2>&1 | tee -a "$LOG"
  echo "### DISC Nf=$NF mRe=$M done (status ${PIPESTATUS[0]})  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
}

# One GPU's whole workload: its masses in order, each mass = 3 Nf concurrent conn -> 3 Nf concurrent disc.
run_gpu () {   # $1=GPU  rest=masses (heavier first)
  local GPU=$1
  shift
  local M
  for M in "$@"; do
    echo "### [GPU $GPU] MASS BLOCK mRe=$M : conn (Nf2,4,6 concurrent)  [$(date +%F_%H:%M:%S)] ###"
    run_conn 2 "$GPU" "$M" &
    run_conn 4 "$GPU" "$M" &
    run_conn 6 "$GPU" "$M" &
    wait
    echo "### [GPU $GPU] mRe=$M conn done ; disc (Nf2,4,6 concurrent)  [$(date +%F_%H:%M:%S)] ###"
    run_disc 2 "$GPU" "$M" &
    run_disc 4 "$GPU" "$M" &
    run_disc 6 "$GPU" "$M" &
    wait
    echo "### [GPU $GPU] mRe=$M FULLY ANALYZABLE (conn + disc, all Nf)  [$(date +%F_%H:%M:%S)] ###"
  done
  echo "### [GPU $GPU] all masses done  [$(date +%F_%H:%M:%S)] ###"
}

echo "### RE-PACK: GPU0={$GPU0_MASSES} ; GPU1={$GPU1_MASSES} ; 3 Nf/mass concurrent  [$(date +%F_%H:%M:%S)] ###"
run_gpu 0 $GPU0_MASSES &
run_gpu 1 $GPU1_MASSES &
wait
echo "### ALL massive Real-m ylm runs done (re-packed)  [$(date +%F_%H:%M:%S)] ###"
