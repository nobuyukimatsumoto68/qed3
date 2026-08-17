#!/usr/bin/env bash
# run_glue_f2_v2_sweep_claude.sh -- production sweep of the POWER-EXTENDED glue basis (p = 2 and 4,
# i.e. F^2 AND F^4, on the 7 production shapes) over ALL massless redo ensembles.
# Driver: glue_f2_v2_shapes_claude.cu.  Plan: glue_f2_v2_shapes_impl_plan_claude.md (verdict in 6b).
#
#   output : data_<ens>/glue_f2_v2_shapes.<k>.h5   (n_shapes = 14, nops = 56)
#            shapes 0..6 = p=2 (bit-identical to the production glue_f2_shapes block),
#            shapes 7..13 = p=4.
#   grid   : EVERY config (kmin 1, stride 1) -- same as the production glue sweep; the k<20
#            thermalization cut is applied later in analysis.
#   scope  : 42 massless ensembles, 93707 configs, ~9 GB of h5.
#
# RESUMABLE: the driver skips per-config h5 gated on "complete", so the already-finished
# L3 Nf2 gsq1.5 ensemble (799 configs) is skipped in seconds, and the sweep can be re-run after any
# rsync top-up.  NO rm anywhere in this script.
#
# CPU-only (gradient flow + host Wilson-loop holonomies), single-threaded workers, one ensemble per
# worker.  This does NOT touch the GPUs, so it can run alongside the fermionic jj jobs.
# ORDERING: by the cost proxy ncfg * L^2 (per-config cost grows with the lattice, config counts
# shrink with it), largest first, so the pool drains evenly.
#
# ETA IS UNCERTAIN.  The only measured rate is L3 = 4.5 s/config/worker (799 configs, 7.5 min on 8
# workers), and L3+L4 are just 10% of the configs.  Scaling that by lattice size suggests ~45
# worker-hours total (~6 h on 8 workers), but the previous production glue sweep (F + F^2, Jul 25 to
# Aug 7) implies a rate ~10x slower, so treat 6 h as a floor and watch the first per-ensemble
# timings in the log.  L1 and L2 carry 90% of the configs and neither has been timed with this
# driver.
#
# Run detached:
#   nohup bash run_glue_f2_v2_sweep_claude.sh > run_glue_f2_v2_sweep_claude.log 2>&1 &
# Knobs: NWORK (default 8), FORCE_BUILD=1
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS=1
NWORK="${NWORK:-8}"

KMAX=100000
KMIN=1
STRIDE=1

SRC=glue_f2_v2_shapes_claude.cu
NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O2 -lcusolver -std=c++17 -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"

LOG=run_glue_f2_v2_sweep_claude.log
echo "================ glue F^2-v2 SWEEP START $(date) ================" | tee -a "$LOG"
echo "NWORK=$NWORK  grid: kmin=$KMIN stride=$STRIDE (every config)" | tee -a "$LOG"

# ---- build one binary per refinement level ----
echo "---- build $(date) ----" | tee -a "$LOG"
for L in 1 2 3 4
do
  BIN=glue_f2_v2_shapes_L${L}_claude.o
  DO=1
  if [ -z "${FORCE_BUILD:-}" ] && [ -f "$BIN" ]
  then
    NEWER=$(find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$BIN" -print -quit 2>/dev/null)
    if [ -z "$NEWER" ]
    then
      DO=0
      echo "L$L binary up-to-date, skip" | tee -a "$LOG"
    fi
  fi
  if [ "$DO" = "1" ]
  then
    if [ "$L" = "1" ]
    then
      DEF=""
    else
      DEF="-DN_REFINE_CLI=$L"
    fi
    echo "compile L$L -> $BIN" | tee -a "$LOG"
    "$NVCC" "$SRC" $DEF $NVCCFLAGS $INCLUDES $LDFLAGS -o "$BIN" 2>&1 | tee -a "$LOG"
    if [ "${PIPESTATUS[0]}" -ne 0 ]
    then
      echo "L$L BUILD FAILED -- stopping" | tee -a "$LOG"
      exit 1
    fi
  fi
done
echo "build OK $(date)" | tee -a "$LOG"

# ---- enumerate massless ensembles, most EXPENSIVE first (cost proxy = ncfg * L^2) ----
ENS_SORTED=$(
  for d in Nf*_*nt128L*_hb*
  do
    [ -d "$d" ] || continue
    case "$d" in *mRe0.000000*) ;; *) continue;; esac
    n=$(ls "$d"/ckpoint_lat.* 2>/dev/null | wc -l)
    [ "$n" -gt 0 ] || continue
    L=$(printf '%s' "$d" | grep -oE 'nt128L[0-9]+' | sed 's/nt128L//')
    cost=$(( n * L * L ))
    echo "$cost $d"
  done | sort -rn | awk '{print $2}'
)
NENS=$(printf '%s\n' "$ENS_SORTED" | grep -c . || true)
NCFG=$(
  for d in $ENS_SORTED
  do
    ls "$d"/ckpoint_lat.* 2>/dev/null | wc -l
  done | awk '{s+=$1} END{print s}'
)
NDONE=$(find . -maxdepth 2 -name "glue_f2_v2_shapes.*.h5" 2>/dev/null | wc -l)
echo "---- $NENS massless ensembles, $NCFG configs total, $NDONE already done ----" | tee -a "$LOG"

# ---- one worker job = one ensemble ----
run_one () {
  local ens="$1"
  local L nf gsq tag elog t0 t1
  L=$(printf '%s' "$ens" | grep -oE 'nt128L[0-9]+' | sed 's/nt128L//')
  nf=$(printf '%s' "$ens" | grep -oE '^Nf[0-9]+' | sed 's/Nf//')
  gsq=$(printf '%s' "$ens" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//')
  tag="Nf${nf}_gsq${gsq}_L${L}"
  elog="glue_f2_v2_${tag}_claude.log"
  t0=$(date +%s)
  echo "[start $(date '+%F %T')] $tag" >> "$LOG"
  {
    echo "==== $tag  F^2-v2 (p=2,4)  $(date) ===="
    ./glue_f2_v2_shapes_L${L}_claude.o "$gsq" "$nf" 1.0 "$KMAX" "$KMIN" "$STRIDE" "$ens"/
    echo "==== $tag done (status $?)  $(date) ===="
  } >> "$elog" 2>&1
  t1=$(date +%s)
  echo "[done  $(date '+%F %T')] $tag  ($(( t1 - t0 ))s)  log: $elog" >> "$LOG"
}
export -f run_one
export KMAX KMIN STRIDE LOG

echo "---- launch $NWORK-worker pool $(date) ----" | tee -a "$LOG"
printf '%s\n' "$ENS_SORTED" | xargs -P "$NWORK" -I{} bash -c 'run_one "$1"' _ {}

NV2=$(find . -maxdepth 2 -name "glue_f2_v2_shapes.*.h5" 2>/dev/null | wc -l)
echo "---- summary $(date) ----" | tee -a "$LOG"
echo "glue_f2_v2_shapes h5 written : $NV2  (configs available: $NCFG)" | tee -a "$LOG"
if [ "$NV2" -lt "$NCFG" ]
then
  echo "NOTE: $(( NCFG - NV2 )) still missing -- re-run to top up (resumable)." | tee -a "$LOG"
fi
echo "================ glue F^2-v2 SWEEP DONE $(date) ================" | tee -a "$LOG"
