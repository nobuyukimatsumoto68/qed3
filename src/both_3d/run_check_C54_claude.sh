#!/usr/bin/env bash
# run_check_C54_claude.sh
#
# Build and run the Appendix C.3 analytic check: the pure-S^2 ("N_t = 1") free fermion
# propagator, derived eigenfunction sum (C.55) vs the Ising-CFT closed form (C.54).
# Pure host math (GSL + Boost Jacobi for the stable xi) -- NO CUDA; g++ -x c++.  Single-thread.
#
# Usage:  bash run_check_C54_claude.sh
# Output to stdout and tee'd to check_C54_claude.log .
set -u

SRC=check_C54_s2_prop_claude.cu
BIN=check_C54_s2_prop_claude
LOG=check_C54_claude.log
BOOST=-I../../qfe_mod/include/boost_1_86_0

export OMP_NUM_THREADS=1
: > "$LOG"

echo "==== BUILD (g++, no CUDA) ====" | tee -a "$LOG"
g++ -O2 -std=c++17 -x c++ $BOOST "$SRC" -o "$BIN" -lgsl -lgslcblas -lm 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]; then echo "BUILD FAILED (rc=$rc)" | tee -a "$LOG"; exit 1; fi

echo "==== RUN ====" | tee -a "$LOG"
"./$BIN" 2>&1 | tee -a "$LOG"
echo "==== DONE ====" | tee -a "$LOG"
