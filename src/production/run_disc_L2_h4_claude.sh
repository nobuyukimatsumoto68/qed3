#!/usr/bin/env bash
# run_disc_L2_h4_claude.sh -- L2 at0.2 disconnected jj Y_lm tower: FOURTH HIT h3 on every disc config (NM 2026-10-23).
# Variant copy of run_disc_L2_h3_claude.sh (kept unchanged; it did the extra-grid fill + h2, finished 2026-10-20).
#   Same config set: GRID A k = 1 mod 5 (800) U GRID B k = r mod 8 (r per ens: Nf2 g3.0 r=3, Nf4 g2.0 r=5, Nf6 g2.0 r=7,
#   else 1) -> 1200 cfg/ens, all with h0,h1,h2 complete.  NHITS=4: the driver skips every complete (k,h) and seeds hit h
#   from (esnid,k,h), so only h3 is computed; h0-h2 stay bit-identical.
#   After the run: 12 disc cfg x 4 hits per M=40 bin on all 9 ensembles (analysis: hit-averaged loop).
# Scope: L=2, a_t=0.2 only (the at0.1 L2 dirs are excluded).
# PACKING: 4 workers = 2 procs/GPU on 2 GPUs under MPS, as in the stride-5 run.  GPUS="0 0" to leave GPU1 free.
# TOOLCHAIN (2026-09-30): CUDA 12.6 is gone -> /usr/local/cuda-12.9/bin/nvcc -ccbin g++-12; cuda-12.9 lib64
#   first on LD_LIBRARY_PATH.  Rebuild ONLY if the disc .cu or ./includes/*.h is newer than the binary
#   (the old any-*.cu check would rebuild for unrelated files); FORCE_BUILD=1 forces it.
# Cost estimate: 10800 (config,hit) solves at ~105/155/230 s (g1/2/3) per proc (measured in the h3 run), 4 procs -> ~5.1 days.
# RESUMABLE: re-launching skips finished work.  NO rm and NO kill anywhere (per project policy).
# Per-job logs: disc_h4_L2_*_{A,B}_w{0..3}_claude.log (this directory).
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
export LD_LIBRARY_PATH=/usr/local/cuda-12.9/lib64:${LD_LIBRARY_PATH:-}

NHITS=4
DISC_TB=2
TARGET_L=2
read -ra GPUS <<< "${GPUS:-0 0 1 1}"     # 4 workers, 2 per GPU (MPS packed)

# NVCC=nvcc
NVCC=/usr/local/cuda-12.9/bin/nvcc
NVCCBASE="-ccbin g++-12 -arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
SRC=jj_local_ylm_disc_stoch_claude.cu

# ---- L2 binary (rebuild when .cu/.h newer; FORCE_BUILD=1 to force) ----
need_build () {
  [ -n "${FORCE_BUILD:-}" ] && return 0
  [ ! -f "$1" ] && return 0
  # find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$1" -print -quit 2>/dev/null | grep -q . && return 0
  [ "$SRC" -nt "$1" ] && return 0
  find ./includes -maxdepth 1 -name '*.h' -newer "$1" -print -quit 2>/dev/null | grep -q . && return 0
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
  case "$d" in *at0.200000*) ;; *) continue;; esac
  ENS+=("$d")
done
echo "### GRID A: k=1 mod 5 ; GRID B: k=r mod 8 (r per ens) ; NHITS=$NHITS  L=$TARGET_L at0.2  tb=$DISC_TB ###"
echo "### ${#ENS[@]} L2 ensembles, ${#GPUS[@]} workers (GPUs ${GPUS[*]}), MPS packed ###"
for e in "${ENS[@]}"; do echo "###   $e ###"; done

grid_b_offset () {   # $1=ens-dir -> residue r of grid B (k = r mod 8), matching the existing extra configs
  case "$1" in
    Nf2_gsq3.000000at0.200000*) echo 3 ;;
    Nf4_gsq2.000000at0.200000*) echo 5 ;;
    Nf6_gsq2.000000at0.200000*) echo 7 ;;
    *) echo 1 ;;
  esac
}

run_disc () {   # $1=ens-dir  $2=GPU  $3=WID  $4=NWORKERS  $5=grid name  $6=grid stride  $7=first k of the grid
  local ens="$1" GPU="$2" WID="$3" NW="$4" GRID="$5" STRIDE="$6" K0="$7"
  local L nf gsq bin LOG last kmin_w stride_eff kmax nconf
  L=$(get_L "$ens")
  nf=$(get_nf "$ens")
  gsq=$(get_gsq "$ens")
  bin="jj_local_ylm_disc_stoch_L${L}.o"
  LOG="disc_h4_L${L}_g${gsq}_Nf${nf}_${GRID}_w${WID}_claude.log"
  last=$(find "$ens" -maxdepth 1 -name 'ckpoint_lat.*' | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n | tail -1)
  if [ -z "$last" ]
  then
    echo "### SKIP $ens: no ckpoint_lat  [$(date +%F_%H:%M:%S)] ###" >> "$LOG"
    return 0
  fi
  kmin_w=$(( K0 + STRIDE * WID ))
  stride_eff=$(( STRIDE * NW ))
  kmax=$(( last + 1 ))
  if [ "$kmin_w" -ge "$kmax" ]
  then
    echo "### $ens $GRID w$WID: no configs in sub-stream  [$(date +%F_%H:%M:%S)] ###" >> "$LOG"
    return 0
  fi
  nconf=$(( (last - kmin_w) / stride_eff + 1 ))
  {
    echo "### DISC-h4 $GRID $ens  L${L}  GPU${GPU}  w${WID}  kmin=$kmin_w stride=$stride_eff kmax=$kmax tb=$DISC_TB nhits=$NHITS (~$nconf cfg)  [$(date +%F_%H:%M:%S)] ###"
    CUDA_VISIBLE_DEVICES=$GPU ./"$bin" --gsq "$gsq" --Nf "$nf" --nu0 1.0 \
      --ens-dir "$ens/" --kmin "$kmin_w" --stride "$stride_eff" --kmax "$kmax" \
      --nhits "$NHITS" --disc-tblock "$DISC_TB"
    echo "### DISC-h4 $GRID $ens w${WID} done (status $?)  [$(date +%F_%H:%M:%S)] ###"
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
    run_disc "${ENS[$idx]}" "$GPU" "$WID" "$NW" B 8 "$(grid_b_offset "${ENS[$idx]}")"
    run_disc "${ENS[$idx]}" "$GPU" "$WID" "$NW" A 5 1
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
echo "### START L2 disc h4 (grid B, then grid A): ${#GPUS[@]} workers (GPUs ${GPUS[*]}), MPS  [$(date +%F_%H:%M:%S)] ###"
for (( W=0; W<${#GPUS[@]}; W++ ))
do
  worker "$W" &
done
wait
echo "### ALL L2 disc h4 done  [$(date +%F_%H:%M:%S)] ###"
