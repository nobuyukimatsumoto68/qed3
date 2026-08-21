#!/usr/bin/env bash
# run_disc_stride2_L2_claude.sh -- STRIDE=2 disconnected jj Y_lm tower, nhits=2, L=2 ONLY.
# Boosts the L2 disc statistics from the stride-10 grid (~398 matched cfg, the disc bottleneck) to the
# full stride-2 grid (~5x), so the L2 vector conn-disc stops being disc-limited.
#
# PACKING: 4 workers = 2 procs/GPU on 2 GPUs under MPS (L2 util ~74% -> packing pays ~1.35x, NM).
# Each worker does a config sub-stream (worker W starts at first+STRIDE*W, effective stride STRIDE*NW);
# the union of the 4 workers is the full stride-2 grid.
#
# RESUMABLE / SAFE: each (config,hit) is its own file corr.<k>.h<h>.h5 gated on "complete", RNG seed is
# (esnid,k,h) independent of stride/nhits, so re-running skips finished work and the stride-10 h0/h1
# files already present stay bit-identical.  NO rm and NO kill anywhere (per project policy).  The
# stride-2 run just FILLS the gaps (configs 22,24,26,... between the existing 20,30,40,...).
#
# NOTE (analysis): while in flight the dir is inhomogeneous (some cfg 2 hits, some 1); hit-0 readers
# ('corr.*.h0.h5') are unaffected; a hit-averaging reader must require equal hit count per config.
#
# The disc binary REBUILDS automatically if the .cu/.h are newer than the .o.
#
# RUN DETACHED (this is a long GPU job -- you launch it, not Claude):
#   nohup bash run_disc_stride2_L2_claude.sh > run_disc_stride2_L2_claude.log 2>&1 &
# Per-job logs: disc_s2_L2_*_w{0..3}_claude.log (this directory).
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

STRIDE=2
NHITS=2
DISC_TB=2
TARGET_L=2
read -ra GPUS <<< "${GPUS:-0 0 1 1}"     # 4 workers, 2 per GPU (MPS packed)

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
SRC=jj_local_ylm_disc_stoch_claude.cu

# ---- L2 binary (rebuild when .cu/.h newer; FORCE_BUILD=1 to force) ----
need_build () {
  [ -n "${FORCE_BUILD:-}" ] && return 0
  [ ! -f "$1" ] && return 0
  find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$1" -print -quit 2>/dev/null | grep -q . && return 0
  return 1
}
BIN="jj_local_ylm_disc_stoch_L${TARGET_L}.o"
if need_build "$BIN"
then
  echo "### compile L${TARGET_L} (-DN_REFINE_CLI=${TARGET_L}) -> $BIN  [$(date +%F_%H:%M:%S)] ###"
  $NVCC $NVCCBASE -DN_REFINE_CLI=${TARGET_L} $INCLUDES $LDFLAGS "$SRC" -o "$BIN" \
    || { echo "### DISC L${TARGET_L} BUILD FAILED ###"; exit 1; }
else
  echo "### L${TARGET_L} binary up-to-date, skip ###"
fi
echo "### build OK ###"

get_L ()   { printf '%s' "$1" | grep -oE 'nt128L[0-9]+' | sed 's/nt128L//'; }
get_nf ()  { printf '%s' "$1" | grep -oE '^Nf[0-9]+' | sed 's/Nf//'; }
get_gsq () { printf '%s' "$1" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//'; }

# ---- collect L2 massless ensembles ----
ENS=()
for d in Nf*_*nt128L*_hb*
do
  [ -d "$d" ] || continue
  case "$d" in *mRe0.000000*) ;; *) continue;; esac
  ls "$d"/ckpoint_lat.* >/dev/null 2>&1 || continue
  [ "$(get_L "$d")" = "$TARGET_L" ] || continue
  ENS+=("$d")
done
echo "### STRIDE=$STRIDE  NHITS=$NHITS  L=$TARGET_L  tb=$DISC_TB ###"
echo "### ${#ENS[@]} L2 ensembles, ${#GPUS[@]} workers (GPUs ${GPUS[*]}), MPS packed ###"
for e in "${ENS[@]}"; do echo "###   $e ###"; done

run_disc () {   # $1=ens-dir  $2=GPU  $3=WID  $4=NWORKERS
  local ens="$1" GPU="$2" WID="$3" NW="$4"
  local L nf gsq bin LOG ks first last kmin_w stride_eff kmax nconf
  L=$(get_L "$ens")
  nf=$(get_nf "$ens")
  gsq=$(get_gsq "$ens")
  bin="jj_local_ylm_disc_stoch_L${L}.o"
  LOG="disc_s2_L${L}_g${gsq}_Nf${nf}_w${WID}_claude.log"
  mapfile -t ks < <(ls "$ens"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n)
  if [ "${#ks[@]}" -eq 0 ]
  then
    echo "### SKIP $ens: no ckpoint_lat  [$(date +%F_%H:%M:%S)] ###" >> "$LOG"
    return 0
  fi
  first="${ks[0]}"
  last="${ks[-1]}"
  kmin_w=$(( first + STRIDE * WID ))
  stride_eff=$(( STRIDE * NW ))
  kmax=$(( last + 1 ))
  if [ "$kmin_w" -ge "$kmax" ]
  then
    echo "### $ens w$WID: no configs in sub-stream  [$(date +%F_%H:%M:%S)] ###" >> "$LOG"
    return 0
  fi
  nconf=$(( (last - kmin_w) / stride_eff + 1 ))
  {
    echo "### DISC-s2 $ens  L${L}  GPU${GPU}  w${WID}  kmin=$kmin_w stride=$stride_eff kmax=$kmax tb=$DISC_TB nhits=$NHITS (~$nconf cfg)  [$(date +%F_%H:%M:%S)] ###"
    CUDA_VISIBLE_DEVICES=$GPU ./"$bin" --gsq "$gsq" --Nf "$nf" --nu0 1.0 \
      --ens-dir "$ens/" --kmin "$kmin_w" --stride "$stride_eff" --kmax "$kmax" \
      --nhits "$NHITS" --disc-tblock "$DISC_TB"
    echo "### DISC-s2 $ens w${WID} done (status $?)  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

worker () {   # $1=WID
  local WID="$1"
  local GPU="${GPUS[$WID]}"
  local n=${#ENS[@]} NW=${#GPUS[@]}
  local i idx
  for (( i=0; i<n; i++ ))
  do
    idx=$(( (i + WID*n/NW) % n ))
    run_disc "${ENS[$idx]}" "$GPU" "$WID" "$NW"
  done
  echo "### worker $WID (GPU $GPU) done  [$(date +%F_%H:%M:%S)] ###"
}

# ---- MPS up (packed phase) ----
if pgrep -f nvidia-cuda-mps-control >/dev/null
then
  echo "### MPS daemon: already running ###"
else
  echo "### starting nvidia-cuda-mps-control -d ###"
  nvidia-cuda-mps-control -d
  for i in 1 2 3 4 5
  do
    pgrep -f nvidia-cuda-mps-control >/dev/null && break
    sleep 1
  done
fi
pgrep -f nvidia-cuda-mps-control >/dev/null \
  || { echo "### ERROR: MPS daemon failed to start -- aborting ###"; exit 1; }

if [ "${#ENS[@]}" -eq 0 ]
then
  echo "### no L2 ensembles found -- nothing to do ###"
  exit 0
fi
echo "### START stride-2 L2 disc: ${#GPUS[@]} workers (GPUs ${GPUS[*]}), MPS  [$(date +%F_%H:%M:%S)] ###"
for (( W=0; W<${#GPUS[@]}; W++ ))
do
  worker "$W" &
done
wait
echo "### ALL stride-2 L2 disc done  [$(date +%F_%H:%M:%S)] ###"
