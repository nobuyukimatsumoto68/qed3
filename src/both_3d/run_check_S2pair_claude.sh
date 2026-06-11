#!/usr/bin/env bash
# run_check_S2pair_claude.sh
#
# Generic off-pole S^2 two-point check (off-diagonal sign machinery) vs the CFT
# gamma-structure.  PURE HOST MATH (GSL) -- no CUDA, no GPU.  The source is named
# .cu only by project convention; we compile it as plain C++ with g++ -x c++.
# Single-thread.
#
# Usage:  bash run_check_S2pair_claude.sh
# Output tee'd to check_S2pair_claude.log .
set -u

SRC=check_S2_generic_pair_claude.cu
BIN=check_S2_generic_pair_claude
LOG=check_S2pair_claude.log
BOOST=-I../../qfe_mod/include/boost_1_86_0

export OMP_NUM_THREADS=1
: > "$LOG"

echo "==== BUILD (g++, no CUDA) ====" | tee -a "$LOG"
g++ -O2 -std=c++17 -x c++ $BOOST "$SRC" -o "$BIN" -lgsl -lgslcblas -lm 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]; then echo "BUILD FAILED (rc=$rc)" | tee -a "$LOG"; exit 1; fi

echo "==== RUN (timed) ====" | tee -a "$LOG"
t0=$(date +%s.%N)
"./$BIN" 2>&1 | tee -a "$LOG"
t1=$(date +%s.%N)
printf "==== DONE -- wall time: %.2f s ====\n" "$(echo "$t1 - $t0" | bc -l)" | tee -a "$LOG"
