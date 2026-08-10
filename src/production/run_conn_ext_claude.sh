#!/usr/bin/env bash
# run_conn_ext_claude.sh -- EXTENDED-STATISTICS connected jj+scalar Y_lm tower on ALL massless redo
# ensembles (incl. L4).  STRIDE=10 (2x the old stride-20 grid), nhits 1, --t0 0 --spin-dilution.
# NO thermalization cut at measurement (kmin = first config; k>=20 cut applied later in analysis).
# CONN ONLY (condensate already done -- run separately if needed).  Run BEFORE run_disc_ext_claude.sh.
# FOUR workers = 2 jobs per GPU under MPS (GPUs 0,1; MPS auto-started).  Each worker does 1/4 of
# EVERY ensemble (config-range split: worker W starts at first+10*W, effective stride 40; union =
# full stride-4 grid, no collision), ensemble order rotated per worker.
# RESUMABLE: driver skips per-config h5 gated on "complete" (old stride-20 configs are a SUBSET of
# the stride-10 grid -> skipped, no rework).  Re-run after each rsync top-up.  No rm anywhere.
#
# Run detached:
#   nohup bash run_conn_ext_claude.sh > run_conn_ext_claude.log 2>&1 &
# Per-job logs: conn_*_w{0,1,2,3}_claude.log (this directory) -- APPEND to the existing conn logs.
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

STRIDE=10
NHITS=1
# worker -> gpu map: 2 jobs per GPU under MPS (override with WGPU="0 1" etc.)
if [ -n "${WGPU:-}" ]
then
  read -ra GPUS <<< "$WGPU"
else
  GPUS=(0 0 1 1)
fi
NWORKERS=${#GPUS[@]}

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
SRC_CONN=jj_local_ylm_scalar_conn_stoch_claude.cu

# ---- ensure MPS daemon up (packing needs it; else workers serialize via context switch) ----
if pgrep -f nvidia-cuda-mps-control >/dev/null
then
  echo "### MPS daemon: already running ###"
else
  echo "### MPS daemon not up -- starting nvidia-cuda-mps-control -d ###"
  nvidia-cuda-mps-control -d
  for i in 1 2 3 4 5
  do
    pgrep -f nvidia-cuda-mps-control >/dev/null && break
    sleep 1
  done
fi
pgrep -f nvidia-cuda-mps-control >/dev/null \
  || { echo "### ERROR: MPS daemon failed to start -- aborting ###"; exit 1; }

# ---- build the per-L conn binaries (skip if up-to-date; FORCE_BUILD=1 to force) ----
need_build () {   # $1 = output binary
  [ -n "${FORCE_BUILD:-}" ] && return 0
  [ ! -f "$1" ] && return 0
  find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$1" -print -quit 2>/dev/null | grep -q . && return 0
  return 1
}
for L in 1 2 3 4
do
  BIN="jj_local_ylm_scalar_conn_stoch_L${L}.o"
  if need_build "$BIN"
  then
    echo "### compile conn L${L} (-DN_REFINE_CLI=${L}) -> $BIN  [$(date +%F_%H:%M:%S)] ###"
    $NVCC $NVCCBASE -DN_REFINE_CLI=${L} $INCLUDES $LDFLAGS "$SRC_CONN" -o "$BIN" \
      || { echo "### conn L${L} BUILD FAILED ###"; exit 1; }
  else
    echo "### conn L${L} binary up-to-date, skip (FORCE_BUILD=1 to rebuild) ###"
  fi
done
echo "### build OK ###"

# ---- massless ensembles from local dirs (auto-covers later arrivals, e.g. L4 Nf4/6) ----
MASSLESS=()
for d in Nf*_*nt128L*_hb*
do
  [ -d "$d" ] || continue
  case "$d" in *mRe0.000000*) ;; *) continue;; esac
  ls "$d"/ckpoint_lat.* >/dev/null 2>&1 && MASSLESS+=("$d")
done
echo "### massless: ${#MASSLESS[@]} ensembles ###"

get_L ()   { printf '%s' "$1" | grep -oE 'nt128L[0-9]+' | sed 's/nt128L//'; }
get_nf ()  { printf '%s' "$1" | grep -oE '^Nf[0-9]+' | sed 's/Nf//'; }
get_gsq () { printf '%s' "$1" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//'; }

# ---- one conn job (worker WID sub-stream of one massless ensemble) ----
run_conn () {   # $1=ens-dir  $2=GPU  $3=WID
  local ens="$1" GPU="$2" WID="$3"
  local L nf gsq bin LOG
  L=$(get_L "$ens")
  nf=$(get_nf "$ens")
  gsq=$(get_gsq "$ens")
  bin="jj_local_ylm_scalar_conn_stoch_L${L}.o"
  LOG="conn_L${L}_g${gsq}_Nf${nf}_w${WID}_claude.log"
  local ks first last kmin_w stride_eff kmax nconf
  mapfile -t ks < <(ls "$ens"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n)
  if [ "${#ks[@]}" -eq 0 ]
  then
    echo "### SKIP $ens: no ckpoint_lat  [$(date +%F_%H:%M:%S)] ###" >> "$LOG"
    return 0
  fi
  first="${ks[0]}"
  last="${ks[-1]}"
  kmin_w=$(( first + STRIDE * WID ))
  stride_eff=$(( STRIDE * NWORKERS ))
  kmax=$(( last + 1 ))
  if [ "$kmin_w" -ge "$kmax" ]
  then
    echo "### $ens w$WID: no configs in sub-stream  [$(date +%F_%H:%M:%S)] ###" >> "$LOG"
    return 0
  fi
  nconf=$(( (last - kmin_w) / stride_eff + 1 ))
  {
    echo "### CONN-EXT $ens  L${L}  GPU${GPU}  w${WID}  kmin=$kmin_w stride=$stride_eff kmax=$kmax (~$nconf cfg)  [$(date +%F_%H:%M:%S)] ###"
    CUDA_VISIBLE_DEVICES=$GPU ./"$bin" --gsq "$gsq" --Nf "$nf" --nu0 1.0 \
      --ens-dir "$ens/" --kmin "$kmin_w" --stride "$stride_eff" --kmax "$kmax" \
      --nhits "$NHITS" --t0 0 --spin-dilution
    echo "### CONN-EXT $ens w${WID} done (status $?)  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

# ---- worker: conn over all massless ensembles (order rotated so workers sit on different ens) ----
worker () {   # $1=WID
  local WID="$1"
  local GPU="${GPUS[$WID]}"
  local n=${#MASSLESS[@]}
  local i idx
  for (( i=0; i<n; i++ ))
  do
    idx=$(( (i + WID*n/NWORKERS) % n ))
    run_conn "${MASSLESS[$idx]}" "$GPU" "$WID"
  done
  echo "### worker $WID (GPU $GPU) done  [$(date +%F_%H:%M:%S)] ###"
}

echo "### START conn-ext: ${NWORKERS} workers (worker->GPU map: ${GPUS[*]}), STRIDE=$STRIDE, MPS on  [$(date +%F_%H:%M:%S)] ###"
for (( W=0; W<NWORKERS; W++ ))
do
  worker "$W" &
done
wait
echo "### ALL conn-ext done  [$(date +%F_%H:%M:%S)] ###"
