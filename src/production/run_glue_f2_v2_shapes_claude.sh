#!/usr/bin/env bash
# run_glue_f2_v2_shapes_claude.sh -- Chunk 3 of glue_f2_v2_shapes_impl_plan_claude.md.
#
# Measures the POWER-EXTENDED shape basis (p = 2 AND p = 4 on the 7 production shapes) on ONE test
# ensemble, so the F^4 operators can be added to the 0++ GEVP.  Driver: glue_f2_v2_shapes_claude.cu.
#
#   output   : data_<ens>/glue_f2_v2_shapes.<k>.h5     (n_shapes = 14, nops = 56)
#   shapes   : 0..6  = p=2 -> the PRODUCTION F^2 basis, must reproduce glue_f2_shapes (Chunk 4 gate)
#              7..13 = p=4 -> F^4
#   default  : L3 Nf2 gsq1.5, ALL configs (stride 1), which is the same ensemble+grid as the
#              existing glue_f2_shapes production data, so the two are directly comparable.
#
# The prefix is DISTINCT from glue_f2_shapes / glue_f2_shapes_s9, so this neither clobbers nor is
# skipped by the "complete" gate of any existing data.  RESUMABLE: finished configs are skipped.
# CPU-only (gradient flow + host Wilson-loop holonomies), single-threaded workers.
# NO rm anywhere in this script.
#
# Parallelism: NWORK processes over ONE ensemble, split by disjoint (kmin, stride) -- worker w takes
# k = 1+w, 1+w+NWORK, ...  This is the driver's documented collision-free packing (each k is its own
# file).  Default NWORK=8 matches the previous glue sweep (run_glue_shapes_redo_claude.sh); override
# with NWORK=4 etc. if the machine is busy.
#
# Run detached:
#   nohup bash run_glue_f2_v2_shapes_claude.sh > run_glue_f2_v2_shapes_claude.log 2>&1 &
# Knobs: NWORK, LEVEL, ENS, KMAX, FORCE_BUILD=1
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS=1

NWORK="${NWORK:-8}"
LEVEL="${LEVEL:-3}"
ENS="${ENS:-Nf2_gsq1.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L3_hb0.400000-1.000000}"
GSQ="${GSQ:-1.500000}"
NF="${NF:-2}"
NU0="${NU0:-1.0}"
KMAX="${KMAX:-100000}"

SRC=glue_f2_v2_shapes_claude.cu
BIN=glue_f2_v2_shapes_L${LEVEL}_claude.o
LOG=run_glue_f2_v2_shapes_claude.log

NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O2 -lcusolver -std=c++17 -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"

echo "================ glue F^2-v2 (p=2,4) START $(date) ================" | tee -a "$LOG"
echo "L=$LEVEL  Nf=$NF  gsq=$GSQ  NWORK=$NWORK" | tee -a "$LOG"
echo "ens = $ENS" | tee -a "$LOG"

# ---- preflight ----
echo "---- preflight $(date) ----" | tee -a "$LOG"
if [ ! -d "$ENS" ]
then
  echo "ERROR: ensemble dir not found: $ENS" | tee -a "$LOG"
  exit 1
fi
NLAT=$(ls "$ENS"/ckpoint_lat.* 2>/dev/null | wc -l)
if [ "$NLAT" -eq 0 ]
then
  echo "ERROR: no ckpoint_lat.* in $ENS" | tee -a "$LOG"
  exit 1
fi
NPROD=$(ls "data_${ENS}"/glue_f2_shapes.*.h5 2>/dev/null | wc -l)
echo "configs available      : $NLAT" | tee -a "$LOG"
echo "production glue_f2_shapes h5 (Chunk 4 comparison set): $NPROD" | tee -a "$LOG"

# ---- phase 1: build ----
echo "---- phase 1 build $(date) ----" | tee -a "$LOG"
DO_BUILD=1
if [ -z "${FORCE_BUILD:-}" ] && [ -f "$BIN" ]
then
  NEWER=$(find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$BIN" -print -quit 2>/dev/null)
  if [ -z "$NEWER" ]
  then
    DO_BUILD=0
    echo "binary $BIN up-to-date, skip build (FORCE_BUILD=1 to force)" | tee -a "$LOG"
  fi
fi
if [ "$DO_BUILD" = "1" ]
then
  if [ "$LEVEL" = "1" ]
  then
    DEF=""
  else
    DEF="-DN_REFINE_CLI=$LEVEL"
  fi
  echo "$NVCC $SRC $DEF ... -o $BIN" | tee -a "$LOG"
  "$NVCC" "$SRC" $DEF $NVCCFLAGS $INCLUDES $LDFLAGS -o "$BIN" 2>&1 | tee -a "$LOG"
  if [ "${PIPESTATUS[0]}" -ne 0 ]
  then
    echo "BUILD FAILED -- stopping" | tee -a "$LOG"
    exit 1
  fi
fi
if [ ! -x "$BIN" ]
then
  echo "ERROR: $BIN missing after build" | tee -a "$LOG"
  exit 1
fi
echo "build OK $(date)" | tee -a "$LOG"

# ---- phase 2: measure ----
# worker w: kmin = 1+w, stride = NWORK  -> the union over w is every config, no two workers share a k.
run_worker () {
  local w="$1"
  local kmin=$(( 1 + w ))
  local wlog="glue_f2_v2_${ENS_TAG}_w${w}_claude.log"
  {
    echo "==== v2 worker $w  kmin=$kmin stride=$NWORK kmax=$KMAX  $(date) ===="
    ./"$BIN" "$GSQ" "$NF" "$NU0" "$KMAX" "$kmin" "$NWORK" "$ENS"/
    echo "==== v2 worker $w done (status $?)  $(date) ===="
  } >> "$wlog" 2>&1
  echo "[worker $w done $(date '+%F %T')] log: $wlog" >> "$LOG"
}
export -f run_worker
ENS_TAG="Nf${NF}_gsq${GSQ}_L${LEVEL}"
export ENS_TAG BIN GSQ NF NU0 KMAX NWORK ENS LOG

echo "---- phase 2 measure: $NWORK workers $(date) ----" | tee -a "$LOG"
seq 0 $(( NWORK - 1 )) | xargs -P "$NWORK" -I{} bash -c 'run_worker "$1"' _ {}

# ---- phase 3: summary ----
echo "---- phase 3 summary $(date) ----" | tee -a "$LOG"
NV2=$(ls "data_${ENS}"/glue_f2_v2_shapes.*.h5 2>/dev/null | wc -l)
echo "glue_f2_v2_shapes h5 written : $NV2" | tee -a "$LOG"
echo "glue_f2_shapes    h5 (production, for the Chunk 4 gate) : $NPROD" | tee -a "$LOG"
echo "configs available            : $NLAT" | tee -a "$LOG"
if [ "$NV2" -lt "$NLAT" ]
then
  echo "NOTE: $(( NLAT - NV2 )) configs still missing -- re-run this script to top up (resumable)." | tee -a "$LOG"
fi
echo "================ glue F^2-v2 DONE $(date) ================" | tee -a "$LOG"
