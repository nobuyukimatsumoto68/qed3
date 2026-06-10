#!/usr/bin/env bash
# run_jj_cft_sptp_check_claude.sh
#
# sp/tp projections of f^{ab}(t;n1,n2) for n1 != n2 (random IR t) vs the canonical-frame (4.26):
#   tp: G_t = f^33  vs  nhat1.T.nhat2 ;   sp: G_s = f^11+f^22  vs  theta1.T.theta2 + phi1.T.phi2.
# Plus an isotropy test (matched n12, different absolute orientation).
# Pure host math (GSL + Boost Jacobi) -- NO CUDA; g++ -x c++.  Single-thread.
#
# Usage:  bash run_jj_cft_sptp_check_claude.sh [nmax]    (default nmax=40, IR / overflow-safe)
# Output tee'd to jj_cft_sptp_check_claude.log .
set -u

SRC=jj_cft_sptp_check_claude.cc
BIN=jj_cft_sptp_check_claude
LOG=jj_cft_sptp_check_claude.log
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
