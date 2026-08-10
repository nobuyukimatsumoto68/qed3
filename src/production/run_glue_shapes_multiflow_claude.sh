#!/usr/bin/env bash
# run_glue_shapes_multiflow_claude.sh -- MULTI-FLOW glueball correlator sweep (F linear + F^2 squared,
# 7-shape basis measured at N_FLOW=4 cumulative Wilson-flow times t = 0.5, 1.0, 2.0, 3.0 at the
# UNCHANGED dt = 0.01) on ALL massless redo ensembles (42: at0.2 L1/L2/L3/L4 + at0.1 L1/L2).
# Basis: 7 -> 28 operators per (l,m) block, the flow level folded into the shape axis
# (ishape = iflow*7 + is), so the GEVP analysis binary needs NO change -- and flow subsets can be
# selected at ANALYSIS time via its orbits_arg operator subset.
# Method: multi-smearing variational basis, C. Morningstar & M. Peardon, arXiv:hep-lat/9901004
# (Phys. Rev. D 60, 034509 (1999)); Wilson flow as the smearing, M. Luescher, arXiv:1006.4518.
#
# Plan: glue_multiflow_impl_plan_claude.md
# Output -> data_<ens>/glue_{msm,f2}_shapes_mf.<k>.h5   (DISTINCT "_mf" prefix: the existing
# single-flow glue_{msm,f2}_shapes.*.h5 production data is neither read, clobbered, nor skipped).
# CPU-only (gradient flow + host Wilson-loop holonomies) -> safe to run CONCURRENTLY with the GPU
# conn/disc job.  NWORK single-threaded workers (default 16), one ensemble per worker at a time
# (F then F^2), largest ensemble first for load balance.
# RESUMABLE: drivers skip per-config h5 gated on "complete" -> re-run after each rsync top-up.
# Binaries are the *_mf_L<L>_claude.o set, SEPARATE from the single-flow production binaries.
# No rm anywhere in this script.
#
# Run detached:
#   NWORK=16 nohup bash run_glue_shapes_multiflow_claude.sh > run_glue_shapes_multiflow_claude.log 2>&1 &
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS=1
NWORK="${NWORK:-16}"

KMAX=100000
KMIN=1
STRIDE=1

NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O2 -lcusolver -std=c++17 -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"

LOG=run_glue_shapes_multiflow_claude.log
echo "================ MULTI-FLOW glue shapes sweep START $(date) ================" | tee -a "$LOG"

# ---- build the eight _mf binaries (L1..L4 x {F, F^2}) ----
echo "---- build $(date) ----" | tee -a "$LOG"
for L in 1 2 3 4
do
  if [ "$L" = "1" ]
  then
    DEF=""
  else
    DEF="-DN_REFINE_CLI=$L"
  fi
  "$NVCC" glue2_msm_shapes_claude.cu $DEF $NVCCFLAGS $INCLUDES $LDFLAGS -o glue2_msm_shapes_mf_L${L}_claude.o 2>&1 | tee -a "$LOG"
  [ "${PIPESTATUS[0]}" -ne 0 ] && { echo "F  L$L BUILD FAILED" | tee -a "$LOG"; exit 1; }
  "$NVCC" glue_f2_shapes_claude.cu   $DEF $NVCCFLAGS $INCLUDES $LDFLAGS -o glue_f2_shapes_mf_L${L}_claude.o   2>&1 | tee -a "$LOG"
  [ "${PIPESTATUS[0]}" -ne 0 ] && { echo "F^2 L$L BUILD FAILED" | tee -a "$LOG"; exit 1; }
done
echo "---- build OK $(date) ----" | tee -a "$LOG"

# ---- massless ensembles from the local dirs, LARGEST ncfg FIRST (load balance);
#      auto-covers later arrivals (L3 still filling, L4 still growing) on a re-run ----
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

# ---- one worker job: F then F^2 for one ensemble ----
run_one () {
  local ens="$1"
  local L nf gsq tag elog
  L=$(printf '%s' "$ens" | grep -oE 'nt128L[0-9]+' | sed 's/nt128L//')
  nf=$(printf '%s' "$ens" | grep -oE '^Nf[0-9]+' | sed 's/Nf//')
  gsq=$(printf '%s' "$ens" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//')
  tag="Nf${nf}_gsq${gsq}_L${L}"
  elog="glue_mf_shapes_${tag}_claude.log"
  echo "[start $(date '+%F %T')] $tag" >> "$LOG"
  {
    echo "==== $tag  F (linear) MULTI-FLOW  $(date) ===="
    ./glue2_msm_shapes_mf_L${L}_claude.o "$gsq" "$nf" 1.0 "$KMAX" "$KMIN" "$STRIDE" "$ens"/
    echo "==== $tag  F^2 (squared) MULTI-FLOW  $(date) ===="
    ./glue_f2_shapes_mf_L${L}_claude.o   "$gsq" "$nf" 1.0 "$KMAX" "$KMIN" "$STRIDE" "$ens"/
  } >> "$elog" 2>&1
  echo "[done  $(date '+%F %T')] $tag  (log: $elog)" >> "$LOG"
}
export -f run_one
export KMAX KMIN STRIDE LOG

echo "---- launch $NWORK-worker pool $(date) ----" | tee -a "$LOG"
printf '%s\n' "$ENS_SORTED" | xargs -P "$NWORK" -I{} bash -c 'run_one "$1"' _ {}

echo "================ MULTI-FLOW glue shapes sweep DONE $(date) ================" | tee -a "$LOG"
