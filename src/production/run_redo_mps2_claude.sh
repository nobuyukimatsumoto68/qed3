#!/bin/bash
#SBATCH --account=affine.lq2_gpu
#SBATCH --qos=opp
# affine account allows only qos {opp, test}; there is NO 'normal'. The wrapper overrides --qos on the
# CLI (opp for --apply, test for --smoke); this header is the fallback for a direct sbatch.
#SBATCH --partition=lq2_gpu
#SBATCH --gpus=a100:1
#SBATCH --cpus-per-task=16
#SBATCH --output=slurm_%x_%j.out
set -u
# ---------------------------------------------------------------------------
# REDO production (corrected-gsq, bug-fixed): TWO overlap-HMC streams packed on ONE A100 under
# NVIDIA MPS. Generic client -- massless OR massive, any (L, Nf, gsq, mass) -- so one script serves
# every pair (the wrapper run_wrapper_redo_claude.sh assigns the two clients per job via --export).
# Mirrors run_massless_mps2_claude.sh: #SBATCH header + MPS daemon + backgrounded clients + wait + trap.
# Driver = the _fermilab block driver (absolute geometry paths), built per (L,Nf,KMAX) by launch_redo_claude.sh.
# Plan: /project/qed3/qed3/src/production/params_L1L2_claude.md + params_massive_claude.md.
#
# via --export (each CLIENT is BIN + its runtime args; both run in OUTDIR, each writes its OWN dir3
# encoding Nf/gsq/mRe/L -> no collision, auto-resume):
#   BIN1,GSQ1,NF1,MRE1   client 1  (massless -> MRE=0.0 ; massive -> MRE=physical mass)
#   BIN2,GSQ2,NF2,MRE2   client 2
#   OUTDIR   (= /lustre2/affine/redo)   WALL_SEC (MUST match the --time the wrapper passes; 4h -> 14400)
# - CLI: BIN <gsq> <Nf> <nu0=1.0> <mass_re> <mass_im=0.0> <MAX_SEC>   (runtime Nf MUST match the binary's -DNF)
# - Graceful wall-time stop via the binary's 6th arg MAX_SEC (stops BETWEEN trajectories, checkpoints).
# ---------------------------------------------------------------------------
BIN1=${BIN1:?set via --export}
GSQ1=${GSQ1:?set via --export}
NF1=${NF1:?set via --export}
MRE1=${MRE1:?set via --export}
BIN2=${BIN2:?set via --export}
GSQ2=${GSQ2:?set via --export}
NF2=${NF2:?set via --export}
MRE2=${MRE2:?set via --export}
OUTDIR=${OUTDIR:?set via --export}
NU0=1.0

WALL_SEC=${WALL_SEC:-28800}      # MUST match the wrapper --time (8h = 28800, opp MaxWall)
SAFETY=${SAFETY:-600}            # headroom for final checkpoint write + MPS teardown
MAX_SEC=$(( WALL_SEC - SAFETY )) # binary's 6th arg: graceful stop between trajectories
BACKSTOP=$(( WALL_SEC - 120 ))   # hard timeout backstop (hang guard)

hostname; date
source /home/nmatsum/env.sh 2>/dev/null || source /lustre2/affine/env.sh
command -v nvidia-cuda-mps-control >/dev/null || { echo "ERROR: nvidia-cuda-mps-control not found"; exit 1; }
[ -x "$BIN1" ] || { echo "ERROR: binary missing/not executable: $BIN1"; exit 1; }
[ -x "$BIN2" ] || { echo "ERROR: binary missing/not executable: $BIN2"; exit 1; }
cd "$OUTDIR" || { echo "ERROR: no OUTDIR $OUTDIR"; exit 1; }
nvidia-smi | head -15

export CUDA_MPS_LOG_DIRECTORY=$OUTDIR/mps_log_${SLURM_JOB_ID}
export CUDA_MPS_PIPE_DIRECTORY=/tmp/nvidia-mps-$USER-${SLURM_JOB_ID}
mkdir -p "$CUDA_MPS_LOG_DIRECTORY" "$CUDA_MPS_PIPE_DIRECTORY"
chmod 700 "$CUDA_MPS_PIPE_DIRECTORY"

cleanup() {
  echo quit | nvidia-cuda-mps-control 2>/dev/null
  # NO rm (global CLAUDE.md rule): the MPS pipe dir is node-local /tmp, unique per SLURM_JOB_ID,
  # harmless if left; SLURM/node epilog reclaims it.
}
trap cleanup EXIT

nvidia-cuda-mps-control -d                        # MUST launch daemon for packing

run_client() {  # $1=bin $2=gsq $3=Nf $4=mass_re
  local bin=$1 g=$2 nf=$3 mre=$4
  local jn="${SLURM_JOB_NAME:-hmc_redo}"
  local jid="${SLURM_JOB_ID:-local}"
  local tag="L$(basename "$bin" | sed -E 's/.*_L([0-9]+)_.*/\1/')_Nf${nf}_gsq${g}_mRe${mre}"
  local out="${jn}_${tag}_jid${jid}_claude.log"
  # 6th arg = MAX_SEC (graceful stop). timeout = hard hang backstop.
  timeout "${BACKSTOP}s" "$bin" "$g" "$nf" "$NU0" "$mre" 0.0 "$MAX_SEC" \
    > "$out" 2>&1
}

echo "=== REDO 2-MPS pack: [1] $(basename "$BIN1") gsq${GSQ1} Nf${NF1} mRe${MRE1} || [2] $(basename "$BIN2") gsq${GSQ2} Nf${NF2} mRe${MRE2}  MAX_SEC=${MAX_SEC}s OUTDIR=${OUTDIR} ==="
run_client "$BIN1" "$GSQ1" "$NF1" "$MRE1" &
run_client "$BIN2" "$GSQ2" "$NF2" "$MRE2" &
wait
date
echo "=== redo 2-MPS job done ==="
