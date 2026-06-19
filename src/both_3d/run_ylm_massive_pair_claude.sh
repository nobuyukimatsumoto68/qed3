#!/usr/bin/env bash
# RE-PACKED (v2: 2 procs/GPU) local-current Y_lm tower on the MASSIVE (Real-m / m_F) ensembles
# (gsq=8.0, nt128 L1, Nf=2,4,6).  Same physics/flags as run_ylm_massive_claude.sh; ONLY packing changes.
#
# PACKING: each GPU runs TWO masses CONCURRENTLY (= 2 procs/GPU).  Each mass is its own serial stream
# that loops Nf = 2 -> 4 -> 6, and for each Nf does conn THEN disc before moving on.  So at any instant a
# GPU hosts exactly 2 procs (one per mass stream).  ~380 MiB/proc, trivial vs 12 GB; util ~86% at 2 procs.
#
# Masses paired to BALANCE the two GPUs (a GPU's wall-time = sum of its two masses; the lightest mass has
# the worst conditioning -> slowest).  Pair fastest+slowest on one GPU, the two middle on the other.  The
# two heaviest (0.2, 0.1) start immediately on separate GPUs, so they finish first (analyzable early):
#   GPU0: 0.200000 || 0.010000        GPU1: 0.100000 || 0.050000
#
# Resume-safe: per-config atomic .h5 + resume-skip, SAME output dirs as the earlier massive runs, so this
# picks up exactly where they left off.  Output:
#   data_<ens>_vmRe<m>vmIm0.000000/{corr_ylm_conn_t00_nhits1_s1,corr_ylm_disc_tb2_nhits1}/ .
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

NFS="2 4 6"                          # Nf order within each mass stream
GPU0_MASSES="0.200000 0.010000"      # run concurrently on GPU0 (2 procs)
GPU1_MASSES="0.100000 0.050000"      # run concurrently on GPU1 (2 procs)

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

CONN=jj_local_ylm_conn_pair.o        # distinct binary names (no ETXTBSY vs the running pack2 binaries)
DISC=jj_local_ylm_disc_pair.o

# ---- build once (early-exit on failure) ----
echo "### compile -> $CONN, $DISC  [$(date +%F_%H:%M:%S)] ###"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_local_ylm_conn_stoch_claude.cu -o "$CONN" || { echo "### CONN BUILD FAILED ###"; exit 1; }
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_local_ylm_disc_stoch_claude.cu -o "$DISC" || { echo "### DISC BUILD FAILED ###"; exit 1; }

ENS () { echo "Nf${1}_gsq8.000000at0.200000nu01.000000mRe${2}mIm0.000000nt128L1/"; }

run_conn () {   # $1=Nf  $2=GPU  $3=mass  (BLOCKING)
  local NF=$1 GPU=$2 M=$3 LOG="ylm_conn_mRe${3}_nf${1}_claude.log"
  echo "### CONN Nf=$NF mRe=$M -> GPU $GPU  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
  CUDA_VISIBLE_DEVICES=$GPU ./"$CONN" --ens-dir "$(ENS $NF $M)" --kmin "$KMIN" --stride "$STRIDE" --kmax "$KMAX" \
    --nhits "$NHITS" --t0 0 --spin-dilution --mass-re "$M" 2>&1 | tee -a "$LOG"
  echo "### CONN Nf=$NF mRe=$M done (status ${PIPESTATUS[0]})  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
}

run_disc () {   # $1=Nf  $2=GPU  $3=mass  (BLOCKING)
  local NF=$1 GPU=$2 M=$3 LOG="ylm_disc_mRe${3}_nf${1}_claude.log"
  echo "### DISC Nf=$NF mRe=$M -> GPU $GPU  (tblock=$DISC_TB)  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
  CUDA_VISIBLE_DEVICES=$GPU ./"$DISC" --ens-dir "$(ENS $NF $M)" --kmin "$KMIN" --stride "$STRIDE" --kmax "$KMAX" \
    --nhits "$NHITS" --disc-tblock "$DISC_TB" --mass-re "$M" 2>&1 | tee -a "$LOG"
  echo "### DISC Nf=$NF mRe=$M done (status ${PIPESTATUS[0]})  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
}

# One mass stream: serial over Nf = 2 -> 4 -> 6, conn then disc for each Nf.  Runs as one proc at a time.
run_mass_stream () {   # $1=GPU  $2=mass
  local GPU=$1 M=$2 NF
  for NF in $NFS; do
    run_conn "$NF" "$GPU" "$M"
    run_disc "$NF" "$GPU" "$M"
  done
  echo "### [GPU $GPU] mRe=$M FULLY DONE (Nf 2->4->6, conn+disc)  [$(date +%F_%H:%M:%S)] ###"
}

# One GPU: its two mass streams run CONCURRENTLY (2 procs).
run_gpu () {   # $1=GPU  $2 $3 = its two masses
  local GPU=$1
  run_mass_stream "$GPU" "$2" &
  run_mass_stream "$GPU" "$3" &
  wait
  echo "### [GPU $GPU] both masses done  [$(date +%F_%H:%M:%S)] ###"
}

echo "### PAIR PACK: GPU0={$GPU0_MASSES} ; GPU1={$GPU1_MASSES} ; 2 procs/GPU, Nf 2->4->6  [$(date +%F_%H:%M:%S)] ###"
run_gpu 0 $GPU0_MASSES &
run_gpu 1 $GPU1_MASSES &
wait
echo "### ALL massive Real-m ylm runs done (pair pack, 2/GPU)  [$(date +%F_%H:%M:%S)] ###"
