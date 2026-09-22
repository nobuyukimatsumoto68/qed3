#!/bin/bash
# Build + run the primal-link-table dump (config-independent geometry) for the interacting T_00 vertex.
# Self-contained; tees all output to t00_dump_links_claude.log.  CPU-only (no GPU kernels).
set -u

LOG=t00_dump_links_claude.log
: > "$LOG"

echo "=== [1/2] module load + compile t00_dump_links_claude.cu ===" | tee -a "$LOG"
module load cuda/12.8 2>/dev/null
module load gcc/13.2.0 2>/dev/null

# no CUDA kernels -> compile the .cu as C++ with g++ (Eigen headers only)
g++ -x c++ -std=c++17 -O2 \
    -I../../qfe_mod/include -Iincludes \
    t00_dump_links_claude.cu -o t00_dump_links_claude.o 2>&1 | tee -a "$LOG"

if [ ! -x ./t00_dump_links_claude.o ]; then
  echo "!! BUILD FAILED (no executable produced) -- see $LOG" | tee -a "$LOG"
  exit 1
fi

echo "" | tee -a "$LOG"
echo "=== [2/2] run: dump primal_links_n1_claude.dat ===" | tee -a "$LOG"
./t00_dump_links_claude.o 2>&1 | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "=== head of primal_links_n1_claude.dat ===" | tee -a "$LOG"
head -8 primal_links_n1_claude.dat 2>&1 | tee -a "$LOG"
echo "=== done ===" | tee -a "$LOG"
