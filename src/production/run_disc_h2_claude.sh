#!/usr/bin/env bash
# run_disc_h2_claude.sh -- SECOND stochastic hit (nhits=2) for the disconnected jj Y_lm tower on ALL
# massless redo ensembles (42: at0.2 L1/L2/L3/L4 + at0.1 L1/L2).  STRIDE=10, disc-tblock=2 (tb2).
#
# WHY THIS IS CHEAP: the output dir no longer encodes nhits (renamed 2026-08-10 to
# data_<ens>/corr_ylm_disc_tb2/), each hit is its own file corr.<k>.h<h>.h5 gated on "complete", and
# the RNG seed is (esnid, k, h) -- INDEPENDENT of nhits.  So --nhits 2 SKIPS every existing h0 and
# computes ONLY the new h1: no duplicated work, and h0 stays bit-identical to what is already there.
#
# NOTE (analysis): the dir is no longer homogeneous by construction.  While this run is in flight some
# configs have 2 hits and some have 1.  The current readers glob 'corr.*.h0.h5' (hit 0 ONLY), so they
# are unaffected -- but a hit-averaging reader MUST require an equal hit count per config and drop the
# stragglers, else configs get unevenly weighted.
#
# PACKING, split by L (NM 2026-08-10).  The earlier "1 disc proc already saturates the GPU" number
# (80-92% util) was measured on L3 and L4 ONLY -- it does NOT transfer to L1, whose 12-site lattice is
# latency-bound like conn.  So:
#   PHASE A: L in $PACK_L (default "1")  -> 2 procs/GPU under MPS, 4 workers  (packed)
#   PHASE B: all remaining L             -> 1 proc/GPU,            2 workers  (unpacked)
# Add L2 to the pack list with  PACK_L="1 2"  if the util check below says it is worth it.
# Each worker does 1/NWORKERS of EVERY ensemble in its phase (worker W starts at first+STRIDE*W,
# effective stride STRIDE*NWORKERS; the union is the full stride-10 grid in BOTH phases).
# MEASURED util: L3/L4 disc 80-92% (saturated -> left solo); L1 disc ~75% (NM) -> packed here.
# 75% implies ~1.33x from time-filling alone, but nvidia-smi sm% only says A kernel was resident, NOT
# that the SMs were full -- on a 12-site lattice occupancy is far lower than util, which is how conn
# got ~1.8x while showing 50-60%.  So expect 1.3x or a little better for L1.
# L2 util is still UNMEASURED and L2 is >half the total work (139 of 265 GPU-h): if it also sits ~75%,
# `PACK_L="1 2"` would take the run from ~133 h to ~100 h.  Check with `nvidia-smi dmon -c 5` once
# phase B reaches an L2 ensemble.
# RESUMABLE: re-run anytime; finished (config,hit) pairs are skipped.  No rm anywhere.
#
# The disc binaries REBUILD automatically on launch: jj_local_ylm_disc_stoch_claude.cu is newer than
# the .o files (the dir_out edit), so need_build() triggers.
#
# Run detached:
#   nohup bash run_disc_h2_claude.sh > run_disc_h2_claude.log 2>&1 &
# Per-job logs: disc_h2_*_w{0,1}_claude.log (this directory).
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

STRIDE=10
NHITS=2                 # h0 already exists and is skipped -> this computes h1 only
DISC_TB=2
PACK_L="${PACK_L:-1}"                             # which L run PACKED (2 procs/GPU under MPS)
read -ra PACK_GPUS <<< "${PACK_GPUS:-0 0 1 1}"    # phase A: 2 per GPU
read -ra SOLO_GPUS <<< "${SOLO_GPUS:-0 1}"        # phase B: 1 per GPU

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
SRC=jj_local_ylm_disc_stoch_claude.cu

# ---- per-L binaries (rebuild when the .cu/.h are newer; FORCE_BUILD=1 to force) ----
need_build () {
  [ -n "${FORCE_BUILD:-}" ] && return 0
  [ ! -f "$1" ] && return 0
  find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$1" -print -quit 2>/dev/null | grep -q . && return 0
  return 1
}
for L in 1 2 3 4
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

get_L ()   { printf '%s' "$1" | grep -oE 'nt128L[0-9]+' | sed 's/nt128L//'; }
get_nf ()  { printf '%s' "$1" | grep -oE '^Nf[0-9]+' | sed 's/Nf//'; }
get_gsq () { printf '%s' "$1" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//'; }

# ---- massless ensembles, split into the PACKED and SOLO sets by L ----
PACKED=()
SOLO=()
for d in Nf*_*nt128L*_hb*
do
  [ -d "$d" ] || continue
  case "$d" in *mRe0.000000*) ;; *) continue;; esac
  ls "$d"/ckpoint_lat.* >/dev/null 2>&1 || continue
  L=$(get_L "$d")
  is_packed=0
  for pl in $PACK_L
  do
    [ "$L" = "$pl" ] && is_packed=1
  done
  if [ "$is_packed" = "1" ]
  then
    PACKED+=("$d")
  else
    SOLO+=("$d")
  fi
done
echo "### NHITS=$NHITS (h0 skipped, h1 computed);  PACK_L='$PACK_L' ###"
echo "### phase A PACKED: ${#PACKED[@]} ensembles, ${#PACK_GPUS[@]} workers (GPUs ${PACK_GPUS[*]}), MPS ###"
echo "### phase B SOLO  : ${#SOLO[@]} ensembles, ${#SOLO_GPUS[@]} workers (GPUs ${SOLO_GPUS[*]}) ###"

run_disc () {   # $1=ens-dir  $2=GPU  $3=WID  $4=NWORKERS (for the config-range split)
  local ens="$1" GPU="$2" WID="$3" NW="$4"
  local L nf gsq bin LOG ks first last kmin_w stride_eff kmax nconf
  L=$(get_L "$ens")
  nf=$(get_nf "$ens")
  gsq=$(get_gsq "$ens")
  bin="jj_local_ylm_disc_stoch_L${L}.o"
  LOG="disc_h2_L${L}_g${gsq}_Nf${nf}_w${WID}_claude.log"
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
    echo "### DISC-h2 $ens  L${L}  GPU${GPU}  w${WID}  kmin=$kmin_w stride=$stride_eff kmax=$kmax tb=$DISC_TB nhits=$NHITS (~$nconf cfg)  [$(date +%F_%H:%M:%S)] ###"
    CUDA_VISIBLE_DEVICES=$GPU ./"$bin" --gsq "$gsq" --Nf "$nf" --nu0 1.0 \
      --ens-dir "$ens/" --kmin "$kmin_w" --stride "$stride_eff" --kmax "$kmax" \
      --nhits "$NHITS" --disc-tblock "$DISC_TB"
    echo "### DISC-h2 $ens w${WID} done (status $?)  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

packed_worker () {   # $1=WID
  local WID="$1"
  local GPU="${PACK_GPUS[$WID]}"
  local n=${#PACKED[@]} NW=${#PACK_GPUS[@]}
  local i idx
  for (( i=0; i<n; i++ ))
  do
    idx=$(( (i + WID*n/NW) % n ))
    run_disc "${PACKED[$idx]}" "$GPU" "$WID" "$NW"
  done
  echo "### packed worker $WID (GPU $GPU) done  [$(date +%F_%H:%M:%S)] ###"
}

solo_worker () {   # $1=WID
  local WID="$1"
  local GPU="${SOLO_GPUS[$WID]}"
  local n=${#SOLO[@]} NW=${#SOLO_GPUS[@]}
  local i idx
  for (( i=0; i<n; i++ ))
  do
    idx=$(( (i + WID*n/NW) % n ))
    run_disc "${SOLO[$idx]}" "$GPU" "$WID" "$NW"
  done
  echo "### solo worker $WID (GPU $GPU) done  [$(date +%F_%H:%M:%S)] ###"
}

# ---- MPS up for the packed phase (kept UP afterwards, per NM: the solo phase is 1 client/GPU
#      = no packing anyway, so leaving the daemon running changes nothing) ----
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

# ================= PHASE A: packed (default L1) =================
if [ "${#PACKED[@]}" -gt 0 ]
then
  echo "### PHASE A disc-h2 PACKED: ${#PACK_GPUS[@]} workers (GPUs ${PACK_GPUS[*]}), MPS on  [$(date +%F_%H:%M:%S)] ###"
  for (( W=0; W<${#PACK_GPUS[@]}; W++ ))
  do
    packed_worker "$W" &
  done
  wait
  echo "### PHASE A done  [$(date +%F_%H:%M:%S)] ###"
fi

# ================= PHASE B: solo (L2/L3/L4) =================
if [ "${#SOLO[@]}" -gt 0 ]
then
  echo "### PHASE B disc-h2 SOLO: ${#SOLO_GPUS[@]} workers (GPUs ${SOLO_GPUS[*]})  [$(date +%F_%H:%M:%S)] ###"
  for (( W=0; W<${#SOLO_GPUS[@]}; W++ ))
  do
    solo_worker "$W" &
  done
  wait
  echo "### PHASE B done  [$(date +%F_%H:%M:%S)] ###"
fi

echo "### ALL disc-h2 done  [$(date +%F_%H:%M:%S)] ###"
