#!/bin/bash -l
# run_L4_scc_claude.sh  (_scc, 2026-07-16, NM)
# =============================================================================================
# BU SCC SGE *BATCH* script for L=4 overlap-fermion QED3 production. ONE job == ONE GPU, running
# TWO ensemble streams packed on that GPU via CUDA MPS (per NM: two streams per GPU).
#
# This script is SUBMITTED by run_wrapper_L4_scc_claude.sh (which loops over Nf/gsq/mass, builds the
# per-Nf binaries, and qsubs one instance of THIS script per GPU with a pair of ensembles via -v).
# The resource requests (-l gpus, -l gpu_c, -l h_rt, -pe omp) are passed on the qsub COMMAND LINE by
# the wrapper so a single batch file serves every GPU architecture; only account/mail/join live here.
#
# Ensemble parameters arrive as environment variables (qsub -v):
#   APPA GSQA NFA MASSA   -- slot A (required):  ./APPA GSQA NFA NU0 MASSA 0.0 MAX_SEC
#   APPB GSQB NFB MASSB   -- slot B (optional):  ./APPB GSQB NFB NU0 MASSB 0.0 MAX_SEC   (empty APPB => solo GPU)
#   NU0                   -- shared nu0 (default 1.0)
#   MAX_SEC               -- wall-time budget passed to the driver (argv[6]); the driver stops GRACEFULLY
#                            before a trajectory that would exceed it and checkpoints every trajectory
#                            (k_ckpoint=1), so the chained follow-on job resumes losslessly. 0 = unlimited.
# Nf is COMPILE-TIME in the binary (-DNF); the runtime NF arg MUST match (asserted by the driver).
# =============================================================================================

#$ -P qfe
#$ -M nmatsum@bu.edu
#$ -j y
#$ -N L4scc

set -u

echo "=========================================================="
echo "Start date : $(date)"
echo "Job name   : ${JOB_NAME:-?}   Job ID : ${JOB_ID:-?}"
echo "Host name  : ${HOSTNAME:-?}"
echo "TMPDIR     : ${TMPDIR:-?}   NSLOTS : ${NSLOTS:-?}"
echo "CUDA_VISIBLE_DEVICES (SGE-assigned) : ${CUDA_VISIBLE_DEVICES:-unset}"
echo "=========================================================="

# ---- environment (unified: cuda, gcc, repo Eigen via $QED3_INC) ---------------------------------
source /projectnb/qfe/nmatsum/qed3/env.sh

SRCDIR=/projectnb/qfe/nmatsum/qed3/src/production
cd "$SRCDIR" || { echo "ERROR: cannot cd $SRCDIR"; exit 1; }

NU0=${NU0:-1.0}
MAX_SEC=${MAX_SEC:-0.0}
: "${APPA:?APPA (slot-A binary) must be set via qsub -v}"
: "${GSQA:?GSQA must be set}"
: "${NFA:?NFA must be set}"
: "${MASSA:=0.0}"
APPB=${APPB:-}
GSQB=${GSQB:-}
NFB=${NFB:-}
MASSB=${MASSB:-0.0}
AT_TAG=${AT_TAG:-}       # optional campaign tag (e.g. at0.1) inserted into the run-log name; empty for a_t=0.2

# ---- split the assigned CPU slots across the (up to two) MPS streams -----------------------------
NPROC=1
[ -n "$APPB" ] && NPROC=2
THREADS=$(( ${NSLOTS:-2} / NPROC ))
[ "$THREADS" -lt 1 ] && THREADS=1
export OMP_NUM_THREADS=$THREADS
export OPENBLAS_NUM_THREADS=$THREADS

# ---- MPS only when packing 2 streams on one GPU. For 1 stream/GPU (NPROC=1) we run the binary DIRECTLY,
# NO MPS daemon -- the per-job MPS daemon was the source of the frequent cudaMalloc/init "CUDA error" on
# shared ece/l40s nodes (per NM 2026-07-18: switch to 1 stream/GPU). ------------------------------------
if [ "$NPROC" -ge 2 ]
then
  export CUDA_MPS_PIPE_DIRECTORY="${TMPDIR:-/tmp}/mps_pipe_${JOB_ID:-$$}"
  export CUDA_MPS_LOG_DIRECTORY="${TMPDIR:-/tmp}/mps_log_${JOB_ID:-$$}"
  mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
  echo "MPS pipe dir : $CUDA_MPS_PIPE_DIRECTORY"
  nvidia-cuda-mps-control -d
  for i in 1 2 3 4 5
  do
    pgrep -f "nvidia-cuda-mps-control" >/dev/null && break
    sleep 1
  done
  pgrep -f "nvidia-cuda-mps-control" >/dev/null \
    || { echo "ERROR: MPS daemon failed to start -- aborting"; exit 1; }
  echo "MPS daemon up; OMP_NUM_THREADS=$THREADS per stream ($NPROC streams)"
else
  echo "single stream (no MPS); OMP_NUM_THREADS=$THREADS"
fi

# ---- launch the ensemble stream(s) on the single shared GPU -------------------------------------
run_slot () {   # tag app gsq Nf mass
  local tag=$1 app=$2 gsq=$3 Nf=$4 mass=$5
  # local log="run_L4_scc_${tag}_Nf${Nf}_gsq${gsq}_m${mass}_$(date +%y%m%d%H%M).log"
  local log="run_L4_scc_${tag}_Nf${Nf}_gsq${gsq}_m${mass}${AT_TAG:+_$AT_TAG}_$(date +%y%m%d%H%M).log"
  echo "### [$tag] ./$app $gsq $Nf $NU0 $mass 0.0 $MAX_SEC   -> $log   [$(date +%F_%H:%M:%S)] ###"
  ./"$app" "$gsq" "$Nf" "$NU0" "$mass" 0.0 "$MAX_SEC" 2>&1 | tee "$log"
  echo "### [$tag] done (status ${PIPESTATUS[0]}) [$(date +%F_%H:%M:%S)] ###"
}

run_slot A "$APPA" "$GSQA" "$NFA" "$MASSA" &
PIDA=$!
if [ -n "$APPB" ]
then
  run_slot B "$APPB" "$GSQB" "$NFB" "$MASSB" &
  PIDB=$!
fi

wait "$PIDA"
[ -n "${PIDB:-}" ] && wait "$PIDB"

# ---- stop the private MPS daemon (only if we started one) ---------------------------------------
[ "$NPROC" -ge 2 ] && { echo quit | nvidia-cuda-mps-control 2>/dev/null || true; }

echo "end date : $(date)"
echo "finished"
