#!/usr/bin/env bash
#SBATCH --account=qed3.lq2_gpu
#SBATCH --qos=normal
#SBATCH --partition=lq2_gpu
#SBATCH --gpus=a100:1
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm_%x_%j.out
# =============================================================================
# _claude -- FNAL L=1 stride-2 connected Y_lm measurement, 2 A100 GPUs (qed3).
# Chunk 2 of conn_s2_fnal_L1_impl_plan_claude.md. Runs the 4 disjoint offset
# residue classes (k = 3,5,7,9 mod 10) as 4 workers, 2/GPU under MPS, over all
# L=1 MASSLESS ensembles. Resumable (driver per-config "complete" gate). NO rm.
#
# Mirrors run_redo_mps2_claude.sh (HMC MPS pack): #SBATCH + in-job MPS daemon +
# backgrounded workers + wait + trap. Binary = the absolute-geometry _fnal copy
# so CWD can be the /lustre2 output workdir. Account/qos/time overridden by the
# wrapper's sbatch flags. GPU request (--gpus=a100:2) stays here.
#
# SMOKE hooks (Chunk 3): CONN_ENS_FILTER=<grep pat> subsets ensembles;
#   CONN_KMAX_CAP=<n> caps kmax (few configs). Submit that via qos=test.
# =============================================================================
set -u

WORKDIR=/lustre2/affine/redo/conn_s2
CFGROOT=/lustre2/affine/redo
CONN_L=${CONN_L:-1}                     # lattice refinement L (1/2/3/4); wrapper passes it via --export
BIN=$WORKDIR/jj_conn_s2_L${CONN_L}_fnal_claude.o
STRIDE=10
NHITS=1
# Offsets to run in THIS job (each = 1 worker). Default = all 4 residue classes (single 2-GPU job);
# the wrapper passes CONN_OFFSETS="2 4" / "6 8" to split into two 1-GPU jobs. Workers are mapped
# round-robin over whatever GPUs Slurm allocated (NGPU below) -> 1-GPU job packs its workers via MPS.
_raw="${CONN_OFFSETS:-2 4 6 8}"                 # wrapper passes "2_4"/"6_8" (underscore-joined); accept "2 4" too
read -ra OFFSETS <<< "${_raw//_/ }"
NW=${#OFFSETS[@]}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}

cd "$WORKDIR" || { echo "ERROR: no WORKDIR $WORKDIR"; exit 1; }
[ -x "$BIN" ] || { echo "ERROR: binary missing: $BIN (run tmp_build_conn_s2_L1_fnal_claude.sh)"; exit 1; }

# Runtime env: same as the HMC job -- `module load hdf5` puts libhdf5.so.310 on LD_LIBRARY_PATH,
# plus the gsl lib export. Without this the binary dies "libhdf5.so.310: cannot open shared object file".
source /home/nmatsum/env.sh 2>/dev/null || source /lustre2/affine/env.sh

echo "### conn_s2 L${CONN_L} job ${SLURM_JOB_ID:-nojob} on $(hostname)  $(date) ###"
nvidia-smi -L 2>/dev/null | head

# ---- L=1 MASSLESS ensembles (mRe0.000000; excludes the massive L1 condensate set) ----
mapfile -t ENS < <(cd "$CFGROOT" && ls -d Nf*mRe0.000000mIm0.000000nt128L${CONN_L}_hb*/ 2>/dev/null | sed 's#/$##')
[ -n "${CONN_ENS_FILTER:-}" ] && mapfile -t ENS < <(printf '%s\n' "${ENS[@]}" | grep -E "$CONN_ENS_FILTER")
NENS=${#ENS[@]}
echo "### ensembles: $NENS ###"; printf '  %s\n' "${ENS[@]}"
[ "$NENS" -gt 0 ] || { echo "ERROR: no L1 massless ensembles matched"; exit 1; }

# ---- MPS daemon (both GPUs served by one daemon); node-local pipe dir, unique per job ----
command -v nvidia-cuda-mps-control >/dev/null || { echo "ERROR: nvidia-cuda-mps-control not found"; exit 1; }
export CUDA_MPS_LOG_DIRECTORY=$WORKDIR/mps_log_${SLURM_JOB_ID:-manual}
export CUDA_MPS_PIPE_DIRECTORY=/tmp/nvidia-mps-$USER-${SLURM_JOB_ID:-manual}
mkdir -p "$CUDA_MPS_LOG_DIRECTORY" "$CUDA_MPS_PIPE_DIRECTORY"
chmod 700 "$CUDA_MPS_PIPE_DIRECTORY"
# NO rm (global rule): pipe dir is node-local /tmp, unique per SLURM_JOB_ID, reclaimed by the epilog.
trap 'echo quit | nvidia-cuda-mps-control 2>/dev/null' EXIT
nvidia-cuda-mps-control -d

# GPUs allocated to THIS job (1 for the split 1-GPU jobs, 2 for the legacy single job); workers
# round-robin over them. Slurm exposes the allocated GPUs as device indices 0..NGPU-1.
NGPU=$(nvidia-smi -L 2>/dev/null | wc -l); [ "${NGPU:-0}" -ge 1 ] || NGPU=1
echo "### NGPU=$NGPU  offsets=(${OFFSETS[*]})  workers=$NW ###"

# ---- one worker = one offset over all ensembles (order rotated to spread I/O) ----
run_worker() {   # $1 = worker id 0..3
  # NB: split the declarations -- under `set -u`, `local a=$1 b=$((a/2))` in ONE statement expands
  # $((a/2)) before `a` is visible -> "unbound variable". Assign, then reference on the next line.
  local wid=$1
  local off=${OFFSETS[$wid]}
  local gpu=$(( wid % NGPU ))                           # round-robin over allocated GPUs (1-GPU job -> all on 0)
  local kmin=$(( 1 + off ))
  local i j ens nf gsq d last kmax
  for (( j=0; j<NENS; j++ )); do
    i=$(( (j + wid) % NENS ))                          # rotated start per worker
    ens=${ENS[$i]}
    nf=$(printf '%s' "$ens" | grep -oE '^Nf[0-9]+' | sed 's/Nf//')
    gsq=$(printf '%s' "$ens" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//')
    d="$CFGROOT/$ens/"
    last=$(ls "$d" 2>/dev/null | grep -oE 'ckpoint_lat\.[0-9]+' | sed 's/.*\.//' | sort -n | tail -1)
    [ -n "$last" ] || { echo "  [off$off] skip $ens (no ckpoint_lat)"; continue; }
    kmax=$(( last + 1 ))                               # half-open [kmin,kmax) -> include k=last
    [ -n "${CONN_KMAX_CAP:-}" ] && [ "$kmax" -gt "$CONN_KMAX_CAP" ] && kmax=$CONN_KMAX_CAP   # smoke cap
    echo "  [off$off gpu$gpu] $ens  kmin=$kmin stride=$STRIDE kmax=$kmax"
    CUDA_VISIBLE_DEVICES=$gpu "$BIN" --gsq "$gsq" --Nf "$nf" --nu0 1.0 \
      --ens-dir "$d" --kmin "$kmin" --stride "$STRIDE" --kmax "$kmax" \
      --nhits "$NHITS" --t0 0 --spin-dilution
  done
}

for (( W=0; W<NW; W++ )); do
  run_worker "$W" > "$WORKDIR/conn_s2_L${CONN_L}_off${OFFSETS[$W]}_jid${SLURM_JOB_ID:-manual}_claude.log" 2>&1 &
done
wait
echo quit | nvidia-cuda-mps-control 2>/dev/null
echo "### conn_s2 L${CONN_L} job done  $(date) ###"
