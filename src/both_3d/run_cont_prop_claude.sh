#!/usr/bin/env bash
# run_cont_prop_claude.sh
#
# Build and run the continuum free fermion propagator (Eq. C.28 of qed3_v2-6.pdf) and write the
# all-to-all matrix in the other project's Dinv.h5 schema (jj_propagator_exact_claude.cu), so the
# SAME contraction / notebooks can read it.  Pure CPU + GSL + HighFive; single-threaded.
#
# Usage:
#   bash run_cont_prop_claude.sh                # build, run --test, then L=1 and L=2
#   bash run_cont_prop_claude.sh 1 2 4          # build, then those L (4 is ~27 GB / heavier)
#   bash run_cont_prop_claude.sh --no-build 1   # skip the build step
#
# Output: cont_prop_L<L>/Dinv.0.h5 .  Log: cont_prop_run_claude.log .
set -u

SRC=cont_prop_eigbasis_claude.cpp
BIN=cont_prop_eigbasis_claude
LOG=cont_prop_run_claude.log
GEO=../../geometry/data
NT=128; AT=0.2; NMAX=40

# match the hdf5/highfive paths used by run_reweighting_R_claude.sh; boost for the stable Jacobi xi
INC='-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/ -I../../qfe_mod/include/boost_1_86_0'
LIB='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5 -lgsl -lgslcblas -lm'

DO_BUILD=1
LS=()
for a in "$@"; do
  case "$a" in
    --no-build) DO_BUILD=0 ;;
    *) LS+=("$a") ;;
  esac
done
[ ${#LS[@]} -eq 0 ] && LS=(1 2)

export OMP_NUM_THREADS=1
: > "$LOG"

if [ "$DO_BUILD" -eq 1 ]; then
  echo "==== BUILD ====" | tee -a "$LOG"
  g++ -O2 -std=c++17 $INC "$SRC" -o "$BIN" $LIB 2>&1 | tee -a "$LOG"
  rc=${PIPESTATUS[0]}
  if [ "$rc" -ne 0 ]; then echo "BUILD FAILED (rc=$rc)" | tee -a "$LOG"; exit 1; fi
fi
[ -x "./$BIN" ] || { echo "no binary ./$BIN (run without --no-build)" | tee -a "$LOG"; exit 1; }

echo "==== VALIDATION SUITE (Chunks 1-3) ====" | tee -a "$LOG"
"./$BIN" --test 2>&1 | tee -a "$LOG"

for L in "${LS[@]}"; do
  echo "==== PRODUCTION L=$L ====" | tee -a "$LOG"
  "./$BIN" --L "$L" --nt "$NT" --at "$AT" --nmax "$NMAX" --geo "$GEO" --out "cont_prop_L$L" \
      2>&1 | tee -a "$LOG"
  rc=${PIPESTATUS[0]}
  if [ "$rc" -ne 0 ]; then echo "L=$L FAILED (rc=$rc)" | tee -a "$LOG"; exit 1; fi
done
echo "==== DONE ====" | tee -a "$LOG"
