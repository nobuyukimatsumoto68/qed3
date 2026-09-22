#!/usr/bin/env bash
# run_disc_remainder_claude.sh -- disconnected jj Y_lm tower (nhits=2) for the OUTSTANDING disc remainders:
#   (a) L1 at=0.1 -- disc entirely UN-RUN on all 9 ensembles (both hits h0+h1 fresh).
#   (b) L4 Nf4 gsq2.0 at=0.2 -- 2nd-hit (h1) topup (was h0=49 / h1=30) + any off-by-1.
# Same driver/grid as run_disc_ext_claude.sh (jj_local_ylm_disc_stoch_claude.cu, STRIDE=10, tb2) but
# NHITS=2 and SCOPED to the two remainder sets above.  a_t auto-derived from the ens-dir (at bug fix),
# so L1 at=0.1 uses a_t=0.1 correctly.  Per-config h5 "complete"-gated (skips done h0/h1 -> for L4 Nf4 g2.0
# only the missing h1 is computed; for L1 at=0.1 both hits are computed fresh).  NO rm anywhere.
#
# GPU0 ONLY, MPS 2-pack by default (GPU1 is running the distillation sweep).  L1 disc ~75% util -> packing
# pays ~1.33x (NM).  Override WGPU="0 1" to spread across both GPUs if GPU1 is free.
#
# LAUNCH AFTER the conn stride-10 run finishes:
#   nohup bash run_disc_remainder_claude.sh > run_disc_remainder_claude.log 2>&1 &
# Per-job logs: disc_rem_*_w{0,1}_claude.log (this directory).
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

STRIDE=10
NHITS=2
DISC_TB=2
read -ra GPUS <<< "${WGPU:-0 0}"   # 2 workers on GPU0 (MPS 2-pack); WGPU="0 1" to spread if GPU1 free
NWORKERS=${#GPUS[@]}

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
SRC=jj_local_ylm_disc_stoch_claude.cu

# ---- per-L binaries (only L1 and L4 needed) ----
need_build () {
  [ -n "${FORCE_BUILD:-}" ] && return 0
  [ ! -f "$1" ] && return 0
  find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$1" -print -quit 2>/dev/null | grep -q . && return 0
  return 1
}
for L in 1 4
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

# ---- REMAINDER ensembles: (a) L1 at=0.1 all 9 ; (b) L4 Nf4 gsq2.0 at=0.2 ----
MASSLESS=()
for d in Nf*_*nt128L*_hb*
do
  [ -d "$d" ] || continue
  case "$d" in *_vmRe*) continue;; esac
  case "$d" in *mRe0.000000*) ;; *) continue;; esac
  ls "$d"/ckpoint_lat.* >/dev/null 2>&1 || continue
  keep=0
  case "$d" in
    Nf*at0.100000*nt128L1_hb*)             keep=1;;   # (a) L1 at=0.1
    Nf4_gsq2.000000at0.200000*nt128L4_hb*) keep=1;;   # (b) L4 Nf4 g2.0 h1 topup
  esac
  [ "$keep" = 1 ] && MASSLESS+=("$d")
done
echo "### remainder ensembles: ${#MASSLESS[@]} ###"
printf '###   %s\n' "${MASSLESS[@]}"

# ---- MPS daemon (2-pack on GPU0 needs it) ----
if pgrep -f nvidia-cuda-mps-control >/dev/null
then
  echo "### MPS daemon: already running ###"
else
  echo "### starting MPS daemon ###"
  nvidia-cuda-mps-control -d
  for i in 1 2 3 4 5
  do
    pgrep -f nvidia-cuda-mps-control >/dev/null && break
    sleep 1
  done
fi
pgrep -f nvidia-cuda-mps-control >/dev/null || { echo "### ERROR: MPS daemon failed to start -- aborting ###"; exit 1; }

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
  LOG="disc_rem_L${L}_g${gsq}_Nf${nf}_w${WID}_claude.log"
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
    echo "### DISC-REM $ens  L${L}  GPU${GPU}  w${WID}  kmin=$kmin_w stride=$stride_eff kmax=$kmax tb=$DISC_TB nhits=$NHITS (~$nconf cfg)  [$(date +%F_%H:%M:%S)] ###"
    CUDA_VISIBLE_DEVICES=$GPU ./"$bin" --gsq "$gsq" --Nf "$nf" --nu0 1.0 \
      --ens-dir "$ens/" --kmin "$kmin_w" --stride "$stride_eff" --kmax "$kmax" \
      --nhits "$NHITS" --disc-tblock "$DISC_TB"
    echo "### DISC-REM $ens w${WID} done (status $?)  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

# ---- worker: each does 1/NWORKERS of every remainder ensemble (order rotated) ----
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

echo "### START disc-remainder: ${NWORKERS} workers (GPUs ${GPUS[*]}), STRIDE=$STRIDE, tb=$DISC_TB, nhits=$NHITS, MPS  [$(date +%F_%H:%M:%S)] ###"
for (( W=0; W<NWORKERS; W++ ))
do
  worker "$W" &
done
wait
echo "### ALL disc-remainder done  [$(date +%F_%H:%M:%S)] ###"
