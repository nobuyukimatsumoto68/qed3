#!/usr/bin/env bash
# run_conn_disc_ext_claude.sh -- EXTENDED-STATISTICS conn + disc jj Y_lm tower on ALL massless redo
# ensembles (incl. L4 AND the half-at at0.100000 ones), in ONE script, run SEQUENTIALLY:
#   PHASE 1 = conn (MPS ON, 2 procs/GPU, 4 workers, both GPUs)
#   PHASE 2 = disc (MPS OFF, 1 proc/GPU, 2 workers, both GPUs)  -- disc genuinely MPS-free
# STRIDE=10 (2x the old stride-20 grid), nhits 1, disc tb=2.  NO thermalization cut at measurement
# (kmin=first config; k>=20 applied in analysis).
# WHY SEQUENTIAL, not concurrent: MPS is PER-GPU, so conn(MPS)+disc(no-MPS) cannot share a GPU.
# True concurrency would force GPU0->conn / GPU1->disc, giving conn only 1 GPU (bottleneck ~220h) --
# SLOWER than these two phases on both GPUs (~198h).  So we run both phases on both GPUs in turn.
# Between phases the MPS daemon is quit (only if THIS script started it) so disc is MPS-free.
# RESUMABLE: per-config h5 "complete"-gated (old stride-20 configs are a SUBSET of the stride-10
# grid -> skipped).  No rm anywhere.
#
# Run detached:
#   nohup bash run_conn_disc_ext_claude.sh > run_conn_disc_ext_claude.log 2>&1 &
# Per-job logs: conn_*_w{0..3}_claude.log / disc_*_w{0,1}_claude.log (this directory).
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

STRIDE=10
NHITS=1
DISC_TB=2
read -ra CONN_GPUS <<< "${CONN_GPUS:-0 0 1 1}"   # conn: 2 procs/GPU (MPS)
read -ra DISC_GPUS <<< "${DISC_GPUS:-0 1}"        # disc: 1 proc/GPU (no MPS)
NWCONN=${#CONN_GPUS[@]}
NWDISC=${#DISC_GPUS[@]}

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
SRC_CONN=jj_local_ylm_scalar_conn_stoch_claude.cu
SRC_DISC=jj_local_ylm_disc_stoch_claude.cu

# ---- build per-L conn + disc binaries (skip if up-to-date; FORCE_BUILD=1 to force) ----
need_build () {
  [ -n "${FORCE_BUILD:-}" ] && return 0
  [ ! -f "$1" ] && return 0
  find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$1" -print -quit 2>/dev/null | grep -q . && return 0
  return 1
}
for L in 1 2 4
do
  for spec in "conn ${SRC_CONN} jj_local_ylm_scalar_conn_stoch_L${L}.o" "disc ${SRC_DISC} jj_local_ylm_disc_stoch_L${L}.o"
  do
    set -- $spec
    if need_build "$3"
    then
      echo "### compile $1 L${L} (-DN_REFINE_CLI=${L}) -> $3  [$(date +%F_%H:%M:%S)] ###"
      $NVCC $NVCCBASE -DN_REFINE_CLI=${L} $INCLUDES $LDFLAGS "$2" -o "$3" \
        || { echo "### $1 L${L} BUILD FAILED ###"; exit 1; }
    else
      echo "### $1 L${L} binary up-to-date, skip ###"
    fi
  done
done
echo "### build OK ###"

# ---- massless ensembles: at0.2 (L1/L2/L4) + at0.1 half-at (L1/L2); mRe0 filter ----
MASSLESS=()
for d in Nf*_*nt128L*_hb*
do
  [ -d "$d" ] || continue
  case "$d" in *mRe0.000000*) ;; *) continue;; esac
  ls "$d"/ckpoint_lat.* >/dev/null 2>&1 && MASSLESS+=("$d")
done
echo "### massless ensembles: ${#MASSLESS[@]} (incl. at0.1 half-at) ###"

get_L ()   { printf '%s' "$1" | grep -oE 'nt128L[0-9]+' | sed 's/nt128L//'; }
get_nf ()  { printf '%s' "$1" | grep -oE '^Nf[0-9]+' | sed 's/Nf//'; }
get_gsq () { printf '%s' "$1" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//'; }

substream () {   # $1=ens $2=WID $3=NW  -> sets KMINW STRIDE_EFF KMAX (return 1 if empty)
  local ens="$1" WID="$2" NW="$3" ks first last
  mapfile -t ks < <(ls "$ens"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n)
  [ "${#ks[@]}" -eq 0 ] && return 1
  first="${ks[0]}"
  last="${ks[-1]}"
  KMINW=$(( first + STRIDE * WID ))
  STRIDE_EFF=$(( STRIDE * NW ))
  KMAX=$(( last + 1 ))
  [ "$KMINW" -ge "$KMAX" ] && return 1
  return 0
}

run_conn () {   # $1=ens $2=GPU $3=WID
  local ens="$1" GPU="$2" WID="$3" L nf gsq bin LOG
  L=$(get_L "$ens"); nf=$(get_nf "$ens"); gsq=$(get_gsq "$ens")
  bin="jj_local_ylm_scalar_conn_stoch_L${L}.o"
  LOG="conn_L${L}_g${gsq}_Nf${nf}_w${WID}_claude.log"
  substream "$ens" "$WID" "$NWCONN" || { echo "### $ens w$WID conn: empty sub-stream  [$(date +%F_%H:%M:%S)] ###" >> "$LOG"; return 0; }
  {
    echo "### CONN-EXT $ens L${L} GPU${GPU} w${WID} kmin=$KMINW stride=$STRIDE_EFF kmax=$KMAX  [$(date +%F_%H:%M:%S)] ###"
    CUDA_VISIBLE_DEVICES=$GPU ./"$bin" --gsq "$gsq" --Nf "$nf" --nu0 1.0 \
      --ens-dir "$ens/" --kmin "$KMINW" --stride "$STRIDE_EFF" --kmax "$KMAX" \
      --nhits "$NHITS" --t0 0 --spin-dilution
    echo "### CONN-EXT $ens w${WID} done (status $?)  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

run_disc () {   # $1=ens $2=GPU $3=WID
  local ens="$1" GPU="$2" WID="$3" L nf gsq bin LOG
  L=$(get_L "$ens"); nf=$(get_nf "$ens"); gsq=$(get_gsq "$ens")
  bin="jj_local_ylm_disc_stoch_L${L}.o"
  LOG="disc_L${L}_g${gsq}_Nf${nf}_w${WID}_claude.log"
  substream "$ens" "$WID" "$NWDISC" || { echo "### $ens w$WID disc: empty sub-stream  [$(date +%F_%H:%M:%S)] ###" >> "$LOG"; return 0; }
  {
    echo "### DISC-EXT $ens L${L} GPU${GPU} w${WID} kmin=$KMINW stride=$STRIDE_EFF kmax=$KMAX tb=$DISC_TB  [$(date +%F_%H:%M:%S)] ###"
    CUDA_VISIBLE_DEVICES=$GPU ./"$bin" --gsq "$gsq" --Nf "$nf" --nu0 1.0 \
      --ens-dir "$ens/" --kmin "$KMINW" --stride "$STRIDE_EFF" --kmax "$KMAX" \
      --nhits "$NHITS" --disc-tblock "$DISC_TB"
    echo "### DISC-EXT $ens w${WID} done (status $?)  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

conn_worker () {   # $1=WID
  local WID="$1" GPU="${CONN_GPUS[$1]}" n=${#MASSLESS[@]} i idx
  for (( i=0; i<n; i++ ))
  do
    idx=$(( (i + WID*n/NWCONN) % n ))
    run_conn "${MASSLESS[$idx]}" "$GPU" "$WID"
  done
  echo "### conn worker $WID (GPU $GPU) done  [$(date +%F_%H:%M:%S)] ###"
}
disc_worker () {   # $1=WID
  local WID="$1" GPU="${DISC_GPUS[$1]}" n=${#MASSLESS[@]} i idx
  for (( i=0; i<n; i++ ))
  do
    idx=$(( (i + WID*n/NWDISC) % n ))
    run_disc "${MASSLESS[$idx]}" "$GPU" "$WID"
  done
  echo "### disc worker $WID (GPU $GPU) done  [$(date +%F_%H:%M:%S)] ###"
}

# ================= PHASE 1: conn (MPS ON) =================
STARTED_MPS=0
if pgrep -f nvidia-cuda-mps-control >/dev/null
then
  echo "### MPS daemon already running (conn phase) ###"
else
  echo "### starting nvidia-cuda-mps-control -d for conn phase ###"
  nvidia-cuda-mps-control -d
  for i in 1 2 3 4 5
  do
    pgrep -f nvidia-cuda-mps-control >/dev/null && { STARTED_MPS=1; break; }
    sleep 1
  done
fi
pgrep -f nvidia-cuda-mps-control >/dev/null \
  || { echo "### ERROR: MPS daemon failed to start -- aborting ###"; exit 1; }

echo "### PHASE 1 conn: ${NWCONN}w (GPUs ${CONN_GPUS[*]}), STRIDE=$STRIDE, MPS on  [$(date +%F_%H:%M:%S)] ###"
for (( W=0; W<NWCONN; W++ )); do conn_worker "$W" & done
wait
echo "### PHASE 1 conn DONE  [$(date +%F_%H:%M:%S)] ###"

# ---- quit MPS before disc IF this script started it (so disc is genuinely MPS-free) ----
if [ "$STARTED_MPS" -eq 1 ]
then
  echo "### quitting MPS daemon (this script started it) for MPS-free disc phase ###"
  echo quit | nvidia-cuda-mps-control 2>/dev/null || true
  sleep 2
else
  echo "### NOTE: MPS daemon pre-existed -- NOT quitting it. disc runs 1 proc/GPU (a single MPS ###"
  echo "###       client per GPU = no packing); stop the daemon manually if you need it fully off. ###"
fi

# ================= PHASE 2: disc (MPS OFF) =================
echo "### PHASE 2 disc: ${NWDISC}w (GPUs ${DISC_GPUS[*]}), STRIDE=$STRIDE, tb=$DISC_TB, no MPS  [$(date +%F_%H:%M:%S)] ###"
for (( W=0; W<NWDISC; W++ )); do disc_worker "$W" & done
wait
echo "### PHASE 2 disc DONE  [$(date +%F_%H:%M:%S)] ###"
echo "### ALL conn+disc-ext done  [$(date +%F_%H:%M:%S)] ###"
