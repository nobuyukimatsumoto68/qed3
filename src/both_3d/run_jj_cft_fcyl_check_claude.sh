#!/usr/bin/env bash
# run_jj_cft_fcyl_check_claude.sh
#
# Build and run the Eq. (4.26) check: general n1 != n2 cylinder conserved-current correlator,
# from contracting the continuum eigenbasis propagator into the vector current.  Frame-free
# checks: f^{ab} has three equal singular values = |C_j| e^{-2t} D^{-2}, det f < 0, and the
# temporal component G_t = f^33 matches its closed form.
# Pure host math (GSL + Boost Jacobi) -- NO CUDA; g++ -x c++.  Single-thread.
#
# Usage:  bash run_jj_cft_fcyl_check_claude.sh [nmax]    (default nmax=40, overflow-safe / IR)
# Output tee'd to jj_cft_fcyl_check_claude.log .
set -u

SRC=jj_cft_fcyl_check_claude.cc
BIN=jj_cft_fcyl_check_claude
LOG=jj_cft_fcyl_check_claude.log
BOOST=-I../../qfe_mod/include/boost_1_86_0
NMAX=${1:-40}

export OMP_NUM_THREADS=1
: > "$LOG"

echo "==== BUILD (g++, no CUDA) ====" | tee -a "$LOG"
g++ -O2 -std=c++17 -x c++ $BOOST "$SRC" -o "$BIN" -lgsl -lgslcblas -lm 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]; then echo "BUILD FAILED (rc=$rc)" | tee -a "$LOG"; exit 1; fi

echo "==== RUN (nmax=$NMAX) ====" | tee -a "$LOG"
"./$BIN" "$NMAX" 2>&1 | tee -a "$LOG"
echo "==== DONE ====" | tee -a "$LOG"
