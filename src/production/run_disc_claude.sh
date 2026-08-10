#!/usr/bin/env bash
# run_disc_claude.sh -- DISCONNECTED jj Y_lm tower (jj_local_ylm_disc_stoch_claude.cu) on ALL
# massless redo ensembles.  stride 20, nhits 1, disc-tblock=2 (tb2, self-contraction-bias fixed).
# NO MPS: ONE process per GPU (2 workers, GPU 0 and GPU 1).  Each worker does half of EVERY
# ensemble (config-range split: worker W starts at first+20*W, effective stride 40; union = full
# stride-20 grid).  Output -> data_<ens>_vmRe0.../corr_ylm_disc_tb2/corr.<k>.h<h>.h5.
# RESUMABLE: per-config h5 "complete"-gated -> re-run after each rsync top-up.  No rm anywhere.
#
# Run detached:
#   nohup bash run_disc_claude.sh > run_disc_claude.log 2>&1 &
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

STRIDE=20
NHITS=1
DISC_TB=2
GPUS=(0 1)              # one process per GPU, NO MPS
NWORKERS=${#GPUS[@]}

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
SRC=jj_local_ylm_disc_stoch_claude.cu

# ---- per-L binaries (skip if up-to-date; FORCE_BUILD=1 to force) ----
need_build () {
  [ -n "${FORCE_BUILD:-}" ] && return 0
  [ ! -f "$1" ] && return 0
  find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$1" -print -quit 2>/dev/null | grep -q . && return 0
  return 1
}
for L in 1 2 4
do
  BIN="jj_local_ylm_disc_stoch_L${L}.o"
  if need_build "$BIN"
  then
    echo "### compile L${L} (-DN_REFINE_CLI=${L}) -> $BIN  [$(date +%F_%H:%M:%S)] ###"
    $NVCC $NVCCBASE -DN_REFINE_CLI=${L} $INCLUDES $LDFLAGS "$SRC" -o "$BIN" \
      || { echo "### DISC L${L} BUILD FAILED ###"; exit 1; }
  else
    echo "### L${L} binary up-to-date, skip ###"
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

run_disc () {   # $1=ens-dir  $2=GPU  $3=WID
  local ens="$1" GPU="$2" WID="$3"
  local L nf gsq bin LOG ks first last kmin_w stride_eff kmax nconf
  L=$(get_L "$ens")
  nf=$(get_nf "$ens")
  gsq=$(get_gsq "$ens")
  bin="jj_local_ylm_disc_stoch_L${L}.o"
  LOG="disc_L${L}_g${gsq}_Nf${nf}_w${WID}_claude.log"
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
    echo "### DISC $ens  L${L}  GPU${GPU}  w${WID}  kmin=$kmin_w stride=$stride_eff kmax=$kmax tb=$DISC_TB (~$nconf cfg)  [$(date +%F_%H:%M:%S)] ###"
    CUDA_VISIBLE_DEVICES=$GPU ./"$bin" --gsq "$gsq" --Nf "$nf" --nu0 1.0 \
      --ens-dir "$ens/" --kmin "$kmin_w" --stride "$stride_eff" --kmax "$kmax" \
      --nhits "$NHITS" --disc-tblock "$DISC_TB"
    echo "### DISC $ens w${WID} done (status $?)  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

# ---- worker: 1 process per GPU; each does 1/NWORKERS of every ensemble (order rotated) ----
worker () {   # $1=WID
  local WID="$1"
  local GPU="${GPUS[$WID]}"
  local n=${#MASSLESS[@]}
  local i idx
  for (( i=0; i<n; i++ ))
  do
    idx=$(( (i + WID*n/NWORKERS) % n ))
    run_disc "${MASSLESS[$idx]}" "$GPU" "$WID"
  done
  echo "### worker $WID (GPU $GPU) done  [$(date +%F_%H:%M:%S)] ###"
}

echo "### START disc: ${NWORKERS} workers (GPUs ${GPUS[*]}), STRIDE=$STRIDE, tb=$DISC_TB, NO MPS  [$(date +%F_%H:%M:%S)] ###"
for (( W=0; W<NWORKERS; W++ ))
do
  worker "$W" &
done
wait
echo "### ALL disc done  [$(date +%F_%H:%M:%S)] ###"
