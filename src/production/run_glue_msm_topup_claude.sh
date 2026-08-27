#!/usr/bin/env bash
# run_glue_msm_topup_claude.sh -- TOPUP of the LINEAR F_12 shape basis (glue_msm_shapes) over all
# massless redo ensembles.  Sibling of run_glue_f2_v2_sweep_claude.sh, which does the SQUARED/QUARTIC
# basis (glue_f2_v2_shapes).  Together they cover every glue observable used by the figures:
#   glue_msm_shapes    -> F l=1 and F l=2 GEVP masses      (this script)
#   glue_f2_v2_shapes  -> F^2 0++ and the F^4 operators    (run_glue_f2_v2_sweep_claude.sh)
#
# Driver: glue2_msm_shapes_claude.cu, REVERTED 2026-08-17 to the production configuration --
# N_FLOW = 1 at FLOW_T = 2.0 and H5PREFIX = "glue_msm_shapes".  It had been left at N_FLOW = 4 with
# the "_mf" test prefix from the NULL multi-flow experiment; running it in that state would have cost
# 4x the shape evaluation AND written a separate test set instead of topping up production.
# VERIFY before launching:  grep -n 'N_FLOW = \|H5PREFIX = ' glue2_msm_shapes_claude.cu
#
#   output : data_<ens>/glue_msm_shapes.<k>.h5
#   grid   : EVERY config (kmin 1, stride 1); the k<20 thermalization cut is applied in analysis.
#
# RESUMABLE / CHEAP TOPUP: the driver's "complete" gate is checked BEFORE U.read and before the flow
# integration, so already-measured configs cost only an h5 open (~ms).  Only the deficit is computed.
# NO rm anywhere in this script.
#
# CPU-only (gradient flow + host Wilson-loop holonomies), single-threaded workers, one ensemble per
# worker.  Does NOT touch the GPUs -- safe to run alongside the fermionic jj jobs.
# ORDERING: cost proxy ncfg * L^2, largest first, so the pool drains evenly.
#
# Run detached:
#   NWORK=2 nohup bash run_glue_msm_topup_claude.sh > run_glue_msm_topup_claude.log 2>&1 &
# Knobs: NWORK (default 2), FORCE_BUILD=1
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS=1
NWORK="${NWORK:-2}"

KMAX=100000
KMIN=1
STRIDE=1

SRC=glue2_msm_shapes_claude.cu
NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O2 -lcusolver -std=c++17 -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"

LOG=run_glue_msm_topup_claude.log
echo "================ glue linear-F (msm) TOPUP START $(date) ================" | tee -a "$LOG"
echo "NWORK=$NWORK  grid: kmin=$KMIN stride=$STRIDE (every config)" | tee -a "$LOG"

# ---- guard: refuse to run if the driver is still in the multi-flow test state ----
NF_LINE=$(grep -E "^\s*constexpr int N_FLOW = " "$SRC" | head -1)
PF_LINE=$(grep -E "^\s*const std::string H5PREFIX = " "$SRC" | tail -1)
echo "driver state: $NF_LINE | $PF_LINE" | tee -a "$LOG"
if ! printf '%s' "$NF_LINE" | grep -q "N_FLOW = 1"
then
  echo "ABORT: $SRC is not at N_FLOW = 1 -- it would write the multi-flow test set." | tee -a "$LOG"
  exit 1
fi
if ! printf '%s' "$PF_LINE" | grep -q '"glue_msm_shapes"'
then
  echo "ABORT: $SRC H5PREFIX is not \"glue_msm_shapes\" -- it would not top up production." | tee -a "$LOG"
  exit 1
fi

# ---- build one binary per refinement level ----
echo "---- build $(date) ----" | tee -a "$LOG"
for L in 1 2 3 4
do
  BIN=glue2_msm_shapes_L${L}_claude.o
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
NDONE=$(find . -maxdepth 2 -name "glue_msm_shapes.*.h5" 2>/dev/null | wc -l)
echo "---- $NENS massless ensembles, $NCFG configs total, $NDONE done, $(( NCFG - NDONE )) to compute ----" | tee -a "$LOG"

# ---- one worker job = one ensemble ----
run_one () {
  local ens="$1"
  local L nf gsq tag elog t0 t1
  L=$(printf '%s' "$ens" | grep -oE 'nt128L[0-9]+' | sed 's/nt128L//')
  nf=$(printf '%s' "$ens" | grep -oE '^Nf[0-9]+' | sed 's/Nf//')
  gsq=$(printf '%s' "$ens" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//')
  tag="Nf${nf}_gsq${gsq}_L${L}"
  elog="glue_msm_topup_${tag}_claude.log"
  t0=$(date +%s)
  echo "[start $(date '+%F %T')] $tag" >> "$LOG"
  {
    echo "==== $tag  linear F (msm) topup  $(date) ===="
    ./glue2_msm_shapes_L${L}_claude.o "$gsq" "$nf" 1.0 "$KMAX" "$KMIN" "$STRIDE" "$ens"/
    echo "==== $tag done (status $?)  $(date) ===="
  } >> "$elog" 2>&1
  t1=$(date +%s)
  echo "[done  $(date '+%F %T')] $tag  ($(( t1 - t0 ))s)  log: $elog" >> "$LOG"
}
export -f run_one
export KMAX KMIN STRIDE LOG

echo "---- launch $NWORK-worker pool $(date) ----" | tee -a "$LOG"
printf '%s\n' "$ENS_SORTED" | xargs -P "$NWORK" -I{} bash -c 'run_one "$1"' _ {}

NMSM=$(find . -maxdepth 2 -name "glue_msm_shapes.*.h5" 2>/dev/null | wc -l)
echo "---- summary $(date) ----" | tee -a "$LOG"
echo "glue_msm_shapes h5 : $NMSM  (configs available: $NCFG)" | tee -a "$LOG"
if [ "$NMSM" -lt "$NCFG" ]
then
  echo "NOTE: $(( NCFG - NMSM )) still missing -- re-run to top up (resumable)." | tee -a "$LOG"
fi
echo "================ glue linear-F (msm) TOPUP DONE $(date) ================" | tee -a "$LOG"
