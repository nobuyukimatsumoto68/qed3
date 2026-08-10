#!/usr/bin/env bash
# run_glue_shapes_redo_claude.sh -- glueball correlator sweep (F linear + F^2 squared, 7-shape
# consolidated basis) on ALL local MASSLESS redo ensembles, EVERY config (stride 1, kmin 1; no
# thermalization cut -- applied later in analysis).  CPU-only (gradient flow + host Wilson-loop
# holonomies); NWORK single-threaded workers (default 8, override NWORK=4 etc.), one ensemble per
# worker at a time (F then F^2), largest ensemble first for load balance.
# Drivers = the production copies with the ens_dir arg7 patch -> output h5 lands in
# data_<ensemble-basename>/glue_{msm,f2}_shapes.<k>.h5 (shared with the fermionic obs).
# RESUMABLE: drivers skip per-config h5 gated on "complete" -> re-run after each rsync top-up.
# No rm anywhere in this script.
#
# Run detached:
#   nohup bash run_glue_shapes_redo_claude.sh > run_glue_shapes_redo_claude.log 2>&1 &
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS=1
NWORK="${NWORK:-8}"

KMAX=100000
KMIN=1
STRIDE=1

NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O2 -lcusolver -std=c++17 -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"

LOG=run_glue_shapes_redo_claude.log
echo "================ REDO glue shapes sweep START $(date) ================" | tee -a "$LOG"

# ---- build six binaries -----------------------------------------------------
echo "---- build $(date) ----" | tee -a "$LOG"
for L in 1 2 3 4
do
  if [ "$L" = "1" ]; then DEF=""; else DEF="-DN_REFINE_CLI=$L"; fi
  "$NVCC" glue2_msm_shapes_claude.cu $DEF $NVCCFLAGS $INCLUDES $LDFLAGS -o glue2_msm_shapes_L${L}_claude.o 2>&1 | tee -a "$LOG"
  [ "${PIPESTATUS[0]}" -ne 0 ] && { echo "F  L$L BUILD FAILED" | tee -a "$LOG"; exit 1; }
  "$NVCC" glue_f2_shapes_claude.cu   $DEF $NVCCFLAGS $INCLUDES $LDFLAGS -o glue_f2_shapes_L${L}_claude.o   2>&1 | tee -a "$LOG"
  [ "${PIPESTATUS[0]}" -ne 0 ] && { echo "F^2 L$L BUILD FAILED" | tee -a "$LOG"; exit 1; }
done

# ---- massless ensembles from the local dirs, LARGEST ncfg FIRST (load balance);
#      auto-covers later arrivals (e.g. L4 Nf4/6) on a re-run ----
ENS_SORTED=$(
  for d in Nf*_*nt128L*_hb*
  do
    [ -d "$d" ] || continue
    case "$d" in *mRe0.000000*) ;; *) continue;; esac
    n=$(ls "$d"/ckpoint_lat.* 2>/dev/null | wc -l)
    [ "$n" -gt 0 ] || continue
    echo "$n $d"
  done | sort -rn | awk '{print $2}'
)
NENS=$(printf '%s\n' "$ENS_SORTED" | grep -c . || true)
echo "---- massless ensembles: $NENS ----" | tee -a "$LOG"

# ---- one worker job: F then F^2 for one ensemble ----------------------------
run_one () {
  local ens="$1"
  local L nf gsq tag elog
  L=$(printf '%s' "$ens" | grep -oE 'nt128L[0-9]+' | sed 's/nt128L//')
  nf=$(printf '%s' "$ens" | grep -oE '^Nf[0-9]+' | sed 's/Nf//')
  gsq=$(printf '%s' "$ens" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//')
  tag="Nf${nf}_gsq${gsq}_L${L}"
  elog="glue_shapes_${tag}_claude.log"
  echo "[start $(date '+%F %T')] $tag" >> "$LOG"
  {
    echo "==== $tag  F (linear)  $(date) ===="
    ./glue2_msm_shapes_L${L}_claude.o "$gsq" "$nf" 1.0 "$KMAX" "$KMIN" "$STRIDE" "$ens"/
    echo "==== $tag  F^2 (squared)  $(date) ===="
    ./glue_f2_shapes_L${L}_claude.o   "$gsq" "$nf" 1.0 "$KMAX" "$KMIN" "$STRIDE" "$ens"/
  } >> "$elog" 2>&1
  echo "[done  $(date '+%F %T')] $tag  (log: $elog)" >> "$LOG"
}
export -f run_one
export KMAX KMIN STRIDE LOG

echo "---- launch $NWORK-worker pool $(date) ----" | tee -a "$LOG"
printf '%s\n' "$ENS_SORTED" | xargs -P "$NWORK" -I{} bash -c 'run_one "$1"' _ {}

echo "================ REDO glue shapes sweep DONE $(date) ================" | tee -a "$LOG"
