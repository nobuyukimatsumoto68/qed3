#!/usr/bin/env bash
# run_glue_therm_flow_claude.sh -- PRODUCTION flow-trajectory sweep (glue_therm_flow_claude.cu)
# over ALL local redo ensembles (massless AND massive), EVERY config (stride 1, kmin 1).
# Protocol v3 (2026-07-18): 400 steps to tmax (eps=tmax/400), save every 2 (201 points);
# tmax PER ENSEMBLE = 20*gsq/L (scaling curve reaches r/t=0.05) -> h5 per ensemble
# data_<basename>/therm_flow_claude.h5 (tlist/E[t,cfg]/klist; plan Sec 4d).
# NOTE: grid change vs the 2026-07-17 files -- the driver ABORTS on old-grid h5; remove the
# old therm_flow_claude.h5 files YOURSELF before the first run of the new protocol.
# 12 single-threaded workers, ONE ensemble per worker at a time (NM 2026-07-17), largest
# ensemble first for load balance. CPU-only.
# RESUMABLE: the driver reads the h5 and flows only missing configs -> re-run after each
# rsync top-up. No rm anywhere in this script.
#
# Run detached:
#   nohup bash run_glue_therm_flow_claude.sh > run_glue_therm_flow_claude.log 2>&1 &
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS=1
NWORK="${NWORK:-12}"

KMAX=100000
KMIN=1
STRIDE=1

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lm'
SRC=glue_therm_flow_claude.cu

LOG=run_glue_therm_flow_claude.log
echo "================ therm-flow sweep START $(date) ================" | tee -a "$LOG"

# ---- build the three per-L binaries (skip if up-to-date; FORCE_BUILD=1 to force) ----
need_build () {   # $1 = output binary
  [ -n "${FORCE_BUILD:-}" ] && return 0
  [ ! -f "$1" ] && return 0
  find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$1" -print -quit 2>/dev/null | grep -q . && return 0
  return 1
}
for L in 1 2 3 4
do
  BIN=glue_therm_flow_L${L}_claude.o
  if need_build "$BIN"
  then
    echo "### compile L${L} (-DN_REFINE_CLI=${L}) -> $BIN  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
    $NVCC $NVCCBASE $INCLUDES $LDFLAGS -DN_REFINE_CLI=$L -o "$BIN" "$SRC" 2>&1 | tee -a "$LOG"
    rc=${PIPESTATUS[0]}
    if [ "$rc" -ne 0 ]
    then
      echo "### ERROR: build failed for L=$L -- aborting ###" | tee -a "$LOG"
      exit "$rc"
    fi
  else
    echo "### L${L} binary up-to-date, skip (FORCE_BUILD=1 to rebuild) ###" | tee -a "$LOG"
  fi
done

# ---- ALL local ensembles (massless + massive), LARGEST ncfg FIRST; auto-covers new arrivals ----
ENS_SORTED=$(
  for d in Nf*_*nt128L*_hb*
  do
    [ -d "$d" ] || continue
    n=$(ls "$d"/ckpoint_lat.* 2>/dev/null | wc -l)
    [ "$n" -gt 0 ] || continue
    echo "$n $d"
  done | sort -rn | awk '{print $2}'
)
NENS=$(printf '%s\n' "$ENS_SORTED" | grep -c . || true)
echo "---- ensembles: $NENS (massless + massive) ----" | tee -a "$LOG"

# ---- one worker job: flow one ensemble ----
run_one () {
  local ens="$1"
  local L nf gsq mre tmax tag elog
  L=$(printf '%s' "$ens" | grep -oE 'nt128L[0-9]+' | sed 's/nt128L//')
  nf=$(printf '%s' "$ens" | grep -oE '^Nf[0-9]+' | sed 's/Nf//')
  gsq=$(printf '%s' "$ens" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//')
  mre=$(printf '%s' "$ens" | grep -oE 'mRe[0-9.]+mIm' | sed 's/mRe//;s/mIm//')
  tmax=$(awk -v g="$gsq" -v l="$L" 'BEGIN{printf "%.6f", 20.0*g/l}')
  tag="Nf${nf}_gsq${gsq}_mRe${mre}_L${L}"
  elog="therm_flow_${tag}_claude.log"
  echo "[start $(date '+%F %T')] $tag  tmax=$tmax" >> "$LOG"
  ./glue_therm_flow_L${L}_claude.o "$gsq" "$nf" 1.0 "$KMAX" "$KMIN" "$STRIDE" "$ens"/ "$tmax" >> "$elog" 2>&1
  echo "[done  $(date '+%F %T')] $tag (status $?, log: $elog)" >> "$LOG"
}
export -f run_one
export KMAX KMIN STRIDE LOG

echo "---- launch $NWORK-worker pool $(date) ----" | tee -a "$LOG"
printf '%s\n' "$ENS_SORTED" | xargs -P "$NWORK" -I{} bash -c 'run_one "$1"' _ {}

echo "================ therm-flow sweep DONE $(date) ================" | tee -a "$LOG"
