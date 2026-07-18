#!/usr/bin/env bash
# run_cond_conn_claude.sh -- GPU measurement campaign on the LOCAL redo ensembles (plan
# redo_obs_measurement_impl_plan_claude.md): PHASE 1 = scalar condensate on ALL massive ensembles,
# PHASE 2 = connected jj+scalar Y_lm tower on ALL massless ensembles.  NO thermalization cut
# (kmin = first config; cut applied later in analysis), stride 20, nhits 1.
# FOUR workers = 2 jobs per GPU under MPS (GPUs 0,1; MPS auto-started, same as the gsq8 conn
# campaign -- NM 2026-07-17).  Each worker does 1/4 of EVERY ensemble (config-range split:
# worker W starts at kmin + 20*W with effective stride 80; union = full stride-20 grid, no
# collision), ensemble order rotated per worker so workers sit on different ensembles.
# Valence mass (condensate) = sea physical mRe (R=1, same at every L -- same convention as
# run_condensate_eo_heavy_claude.sh).  Massless conn: --t0 0 --spin-dilution (gsq8-campaign
# settings).
# RESUMABLE: both drivers skip per-config h5 gated on "complete" -> re-run this script after each
# rsync top-up to extend to the new configs.
#
# Run detached:
#   nohup bash run_cond_conn_claude.sh > run_cond_conn_claude.log 2>&1 &
# Per-job logs: cond_*_w{0,1}_claude.log / conn_*_w{0,1}_claude.log (this directory).
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

STRIDE=20
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
SRC_COND=jj_scalar_condensate_eo_stoch_claude.cu
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

# ---- build the per-L binaries for both drivers (skip if up-to-date; FORCE_BUILD=1 to force) ----
need_build () {   # $1 = output binary
  [ -n "${FORCE_BUILD:-}" ] && return 0
  [ ! -f "$1" ] && return 0
  find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$1" -print -quit 2>/dev/null | grep -q . && return 0
  return 1
}
for L in 1 2 4
do
  for pre in cond conn
  do
    if [ "$pre" = cond ]
    then
      SRC="$SRC_COND"
      BIN="jj_scalar_condensate_eo_L${L}.o"
    else
      SRC="$SRC_CONN"
      BIN="jj_local_ylm_scalar_conn_stoch_L${L}.o"
    fi
    if need_build "$BIN"
    then
      echo "### compile ${pre} L${L} (-DN_REFINE_CLI=${L}) -> $BIN  [$(date +%F_%H:%M:%S)] ###"
      $NVCC $NVCCBASE -DN_REFINE_CLI=${L} $INCLUDES $LDFLAGS "$SRC" -o "$BIN" \
        || { echo "### ${pre} L${L} BUILD FAILED ###"; exit 1; }
    else
      echo "### ${pre} L${L} binary up-to-date, skip (FORCE_BUILD=1 to rebuild) ###"
    fi
  done
done
echo "### build OK ###"

# ---- ensemble discovery from the local dirs (auto-covers later arrivals, e.g. L4 Nf4/6) ----
MASSIVE=()
MASSLESS=()
for d in Nf*_*nt128L*_hb*
do
  [ -d "$d" ] || continue
  ls "$d"/ckpoint_lat.* >/dev/null 2>&1 || continue
  if [[ "$d" == *mRe0.000000* ]]
  then
    MASSLESS+=("$d")
  else
    MASSIVE+=("$d")
  fi
done
echo "### massive: ${#MASSIVE[@]} ensembles ; massless: ${#MASSLESS[@]} ensembles ###"

# ---- helpers to parse the dir name ----
get_L ()   { printf '%s' "$1" | grep -oE 'nt128L[0-9]+' | sed 's/nt128L//'; }
get_mre () { printf '%s' "$1" | grep -oE 'mRe[0-9.]+mIm' | sed 's/mRe//;s/mIm//'; }
get_nf ()  { printf '%s' "$1" | grep -oE '^Nf[0-9]+' | sed 's/Nf//'; }
get_gsq () { printf '%s' "$1" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//'; }

# ---- one condensate job (worker WID sub-stream of one massive ensemble) ----
run_cond () {   # $1=ens-dir  $2=GPU  $3=WID
  local ens="$1" GPU="$2" WID="$3"
  local L mre nf gsq bin LOG
  L=$(get_L "$ens")
  mre=$(get_mre "$ens")
  nf=$(get_nf "$ens")
  gsq=$(get_gsq "$ens")
  bin="jj_scalar_condensate_eo_L${L}.o"
  LOG="cond_L${L}_g${gsq}_Nf${nf}_mRe${mre}_w${WID}_claude.log"
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
    echo "### COND $ens  L${L}  mass=$mre  GPU${GPU}  w${WID}  kmin=$kmin_w stride=$stride_eff kmax=$kmax (~$nconf cfg)  [$(date +%F_%H:%M:%S)] ###"
    CUDA_VISIBLE_DEVICES=$GPU ./"$bin" --ens-dir "$ens/" --kmin "$kmin_w" --stride "$stride_eff" --kmax "$kmax" \
      --nhits "$NHITS" --mass-re "$mre"
    echo "### COND $ens w${WID} done (status $?)  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

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
    echo "### CONN $ens  L${L}  GPU${GPU}  w${WID}  kmin=$kmin_w stride=$stride_eff kmax=$kmax (~$nconf cfg)  [$(date +%F_%H:%M:%S)] ###"
    CUDA_VISIBLE_DEVICES=$GPU ./"$bin" --gsq "$gsq" --Nf "$nf" --nu0 1.0 \
      --ens-dir "$ens/" --kmin "$kmin_w" --stride "$stride_eff" --kmax "$kmax" \
      --nhits "$NHITS" --t0 0 --spin-dilution
    echo "### CONN $ens w${WID} done (status $?)  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

# ---- worker: PHASE 1 condensate (all massive), then PHASE 2 conn (all massless); order rotated ----
worker () {   # $1=WID
  local WID="$1"
  local GPU="${GPUS[$WID]}"
  local n i idx
  echo "### worker $WID on GPU $GPU: phase 1 condensate (${#MASSIVE[@]})  [$(date +%F_%H:%M:%S)] ###"
  n=${#MASSIVE[@]}
  for (( i=0; i<n; i++ ))
  do
    idx=$(( (i + WID*n/NWORKERS) % n ))
    run_cond "${MASSIVE[$idx]}" "$GPU" "$WID"
  done
  echo "### worker $WID on GPU $GPU: phase 2 conn (${#MASSLESS[@]})  [$(date +%F_%H:%M:%S)] ###"
  n=${#MASSLESS[@]}
  for (( i=0; i<n; i++ ))
  do
    idx=$(( (i + WID*n/NWORKERS) % n ))
    run_conn "${MASSLESS[$idx]}" "$GPU" "$WID"
  done
  echo "### worker $WID done  [$(date +%F_%H:%M:%S)] ###"
}

echo "### START: ${NWORKERS} workers (worker->GPU map: ${GPUS[*]}), STRIDE=$STRIDE, MPS on  [$(date +%F_%H:%M:%S)] ###"
for (( W=0; W<NWORKERS; W++ ))
do
  worker "$W" &
done
wait
echo "### ALL done  [$(date +%F_%H:%M:%S)] ###"
