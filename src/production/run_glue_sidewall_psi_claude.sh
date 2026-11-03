#!/usr/bin/env bash
# run_glue_sidewall_psi_claude.sh -- chunk 3 of final/analysis_gluonic/jtop_psi_impl_plan_claude.md:
# MEASUREMENT of the j_top Psi component (side-wall flux, ops A/B x 5 smearing levels, l=1,2) on every config.
# Driver: glue_sidewall_psi_claude.cu (test passed 2026-10-12: Gauss law 5e-15, gauge invariance 1e-14, flow stable).
#   output : data_<ens>/glue_sidewall_psi.<k>.h5   (same layout as glue_msm_shapes)
#   grid   : EVERY config (kmin 1, stride 1); the k<20 thermalization cut and M=40 binning are applied in analysis.
# RESUMABLE: the driver skips configs whose h5 already has the "complete" flag.  NO rm, no kill.
# CPU-only (gradient flows + temporal plaquettes), single-threaded workers, one ensemble per worker; no GPU.
# Default scope = NM order: L1 at0.2 (9 ensembles) first.  Then L2:  L_ONLY=2 bash run_glue_sidewall_psi_claude.sh
# Run detached:
#   NWORK=3 nohup bash run_glue_sidewall_psi_claude.sh > run_glue_sidewall_psi_claude.out 2>&1 &
# Knobs: NWORK (default 3), L_ONLY (default 1), AT_ONLY (default 0.200000), FORCE_BUILD=1
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=/usr/local/cuda-12.9/lib64:/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib:${LD_LIBRARY_PATH:-}
NWORK="${NWORK:-3}"
L_ONLY="${L_ONLY:-1}"
AT_ONLY="${AT_ONLY:-0.200000}"

SRC=glue_sidewall_psi_claude.cu
NVCC=/usr/local/cuda-12.9/bin/nvcc
NVCCFLAGS="-ccbin g++-12 -w -arch=sm_70 -O2 -std=c++17 -lcusolver -Xcompiler -fopenmp"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"

LOG=run_glue_sidewall_psi_claude.log
echo "================ glue side-wall Psi MEASUREMENT START $(date) ================" | tee -a "$LOG"
echo "NWORK=$NWORK L_ONLY=$L_ONLY AT_ONLY=$AT_ONLY grid: every config" | tee -a "$LOG"

# ---- build one production binary per requested L (no -DSIDEWALL_TEST) ----
for L in $L_ONLY
do
  BIN=glue_sidewall_psi_L${L}_claude.o
  DO=1
  if [ -z "${FORCE_BUILD:-}" ] && [ -f "$BIN" ]
  then
    NEWER=$(find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$BIN" -print -quit 2>/dev/null)
    if [ -z "$NEWER" ]
    then
      DO=0
      echo "L$L binary up-to-date, skip build" | tee -a "$LOG"
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
    if [ "${PIPESTATUS[0]}" -ne 0 ] || [ ! -x "$BIN" ]
    then
      echo "L$L BUILD FAILED -- stopping" | tee -a "$LOG"
      exit 1
    fi
  fi
done

# ---- enumerate massless ensembles for the requested L and a_t ----
ENS_LIST=$(
  for d in Nf*_*at${AT_ONLY}*nt128L*_hb*
  do
    [ -d "$d" ] || continue
    case "$d" in *mRe0.000000*) ;; *) continue;; esac
    L=$(printf '%s' "$d" | grep -oE 'nt128L[0-9]+' | sed 's/nt128L//')
    case " $L_ONLY " in
      *" $L "*) ;;
      *) continue;;
    esac
    [ -e "$d/ckpoint_lat.20" ] || continue
    echo "$d"
  done
)
NENS=$(printf '%s\n' "$ENS_LIST" | grep -c . || true)
echo "---- $NENS ensembles ----" | tee -a "$LOG"
printf '%s\n' "$ENS_LIST" | tee -a "$LOG"

run_one () {
  local ens="$1"
  local L nf gsq tag elog t0 t1
  L=$(printf '%s' "$ens" | grep -oE 'nt128L[0-9]+' | sed 's/nt128L//')
  nf=$(printf '%s' "$ens" | sed -E 's/^Nf([0-9]+)_.*/\1/')
  gsq=$(printf '%s' "$ens" | sed -E 's/.*_gsq([0-9.]+)at.*/\1/')
  tag="Nf${nf}_g${gsq}_L${L}"
  elog="sidewall_psi_${tag}_at${AT_ONLY}_claude.log"
  t0=$(date +%s)
  # args: gsq Nf nu0 kmax_run kmin stride ens_dir
  ./glue_sidewall_psi_L${L}_claude.o "$gsq" "$nf" 1.0 100000 1 1 "$ens" > "$elog" 2>&1
  t1=$(date +%s)
  echo "done $tag  $((t1 - t0)) s  ($(grep -c '# read from k' "$elog") configs measured this run)"
}
export -f run_one
export AT_ONLY

printf '%s\n' "$ENS_LIST" | xargs -P "$NWORK" -I{} bash -c 'run_one "$1"' _ {} 2>&1 | tee -a "$LOG"
NDONE=$(find . -maxdepth 2 -name "glue_sidewall_psi.*.h5" 2>/dev/null | wc -l)
echo "================ DONE $(date) : $NDONE glue_sidewall_psi h5 files in total ================" | tee -a "$LOG"
