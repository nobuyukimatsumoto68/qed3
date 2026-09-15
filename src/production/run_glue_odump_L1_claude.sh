#!/usr/bin/env bash
# run_glue_odump_L1_claude.sh -- ADD the per-timeslice operator series O_i(t) ("O" dataset) to the L1
# glue h5, for the F^2 - sigma\sigma (0++) mixing cross-correlator <O_F(t) sigma\sigma(0)>_c.
#
# Driver glue_f2_v2_shapes_claude.cu now writes "O" (nops x Nt) and, with GLUE_REQUIRE_O=1, RE-MEASURES a
# config unless it already has "O" (the Truncate write recomputes F_corr_blk BIT-IDENTICALLY via the
# deterministic Wilson flow, so no data loss -- it only ADDS O to existing complete files).
#
# Matched to the scalar-loop run (run_sigma_L1_claude.sh): L=1, at0.2 massless (Nf2/4/6 x gsq0.5/1/1.5),
# KMIN=1 STRIDE=10 -> exactly the loop configs k=1,11,21,...  CPU-only (coexists with the GPU loop job).
# NO rm.  NO kill.  Launch DETACHED:
#   setsid nohup bash run_glue_odump_L1_claude.sh > glue_odump_L1_claude.log 2>&1 < /dev/null &
set -u

# ---- knobs ----
NWORK="${NWORK:-4}"          # CPU workers (bump with NWORK=8 if the box is idle)
KMIN=1
STRIDE=10                    # match the loop-run stride (configs k=1,11,21,...)
KMAX="${KMAX:-100000000}"
L=1
AT_TOKEN="${AT_TOKEN:-at0.200000}"
ENS_GLOB="${ENS_GLOB:-Nf*_gsq*${AT_TOKEN}*mRe0.000000mIm0.000000nt128L${L}_hb*}"

export GLUE_REQUIRE_O=1       # <-- forces re-measure of complete-but-no-O configs (adds O)

SRC=glue_f2_v2_shapes_claude.cu
BIN=glue_f2_v2_shapes_L${L}_claude.o
NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O2 -lcusolver -std=c++17 -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"
LOG=run_glue_odump_L1_claude.log

echo "================ glue O-dump L1 START $(date) ================" | tee -a "$LOG"
echo "NWORK=$NWORK  KMIN=$KMIN STRIDE=$STRIDE  GLUE_REQUIRE_O=$GLUE_REQUIRE_O  AT=$AT_TOKEN" | tee -a "$LOG"

# ---- build the L1 binary (FORCE, since the .cu was edited to add O) ----
echo "---- compile L${L} -> $BIN  $(date) ----" | tee -a "$LOG"
"$NVCC" "$SRC" $NVCCFLAGS $INCLUDES $LDFLAGS -o "$BIN" 2>&1 | tee -a "$LOG"
if [ "${PIPESTATUS[0]}" -ne 0 ]
then
  echo "L${L} BUILD FAILED -- stopping" | tee -a "$LOG"
  exit 1
fi
echo "build OK $(date)" | tee -a "$LOG"

# ---- enumerate matched L1 ensembles ----
ENS_SORTED=$(
  for d in $ENS_GLOB
  do
    [ -d "$d" ] || continue
    ls "$d"/ckpoint_lat.* >/dev/null 2>&1 || continue
    echo "$d"
  done
)
NENS=$(printf '%s\n' "$ENS_SORTED" | grep -c . || true)
echo "---- $NENS L${L} ${AT_TOKEN} massless ensembles ----" | tee -a "$LOG"
printf '  %s\n' $ENS_SORTED | tee -a "$LOG"

# ---- one worker job = one ensemble ----
run_one () {
  local ens="$1"
  local nf gsq elog
  nf=$(printf '%s' "$ens" | grep -oE '^Nf[0-9]+' | sed 's/Nf//')
  gsq=$(printf '%s' "$ens" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//')
  elog="glue_odump_Nf${nf}_gsq${gsq}_L1_claude.log"
  {
    echo "==== Nf${nf} g${gsq} L1  O-dump (kmin=$KMIN stride=$STRIDE)  $(date) ===="
    ./glue_f2_v2_shapes_L1_claude.o "$gsq" "$nf" 1.0 "$KMAX" "$KMIN" "$STRIDE" "$ens"/
    echo "==== Nf${nf} g${gsq} L1 done (status $?)  $(date) ===="
  } >> "$elog" 2>&1
}
export -f run_one
export KMIN STRIDE KMAX

echo "---- launch $NWORK-worker pool $(date) ----" | tee -a "$LOG"
printf '%s\n' "$ENS_SORTED" | xargs -P "$NWORK" -I{} bash -c 'run_one "$1"' _ {}
echo "================ glue O-dump L1 DONE $(date) ================" | tee -a "$LOG"
