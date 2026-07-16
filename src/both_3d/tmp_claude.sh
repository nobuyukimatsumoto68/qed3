#!/usr/bin/env bash
# BUILD + RUN the block-solve microbenchmark (BlockedMat mrhs vs serial op_Dmsq) for L1/L2/L4 on a real
# massless config.  Single process on GPU 0 (no MPS -- measuring one-process block throughput).
# Reads back: block_solve_bench_claude.log .
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

GPU=0
NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
SRC=test_block_solve_bench_claude.cu

LOG=block_solve_bench_claude.log
# (log appends; timestamps + '### BUILD' headers delimit each run.  If you want a clean log, rename the old one.)

# ens-dir | config k  (real massless configs; base=1 contiguous so these exist)
ENS_L1="Nf2_gsq8.000000at0.200000nu01.000000nt128L1"; K_L1=2000
ENS_L2="Nf2_gsq8.000000at0.200000nu01.000000nt128L2"; K_L2=1000
ENS_L4="Nf2_gsq8.000000at0.200000nu01.000000nt128L4"; K_L4=100

run_L () {   # $1=L  $2=ens  $3=k
  local L="$1" ens="$2" k="$3"
  local BIN="test_block_solve_bench_L${L}.o"
  echo "### BUILD L${L} (-DN_REFINE_CLI=${L}) -> $BIN  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
  $NVCC $NVCCBASE -DN_REFINE_CLI=${L} $INCLUDES $LDFLAGS "$SRC" -o "$BIN" 2>>"$LOG" \
    || { echo "### L${L} BUILD FAILED (see $LOG) ###" | tee -a "$LOG"; return 1; }
  echo "### RUN L${L} on $ens k=$k  GPU$GPU  (NB<=8)  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
  CUDA_VISIBLE_DEVICES=$GPU ./"$BIN" --ens-dir "${ens}/" --k "$k" --maxnb 8 2>&1 | tee -a "$LOG"
  echo "### L${L} done  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
}

# L1 settled (~2.1x @ NB=16, tiny lattice).  Re-check L2 + L4 with the added NB=2,4 to fill the gap
# (L2 best was NB=8=1.66x; L4 was NB=1=1.36x then NB>=8 hurts).  L2 is cheap; for L4 the results print
# per-NB, so you can stop (Ctrl-C) after NB=8 -- NB=16/32/64 are already known worse for L4.
# run_L 1 "$ENS_L1" "$K_L1"
run_L 2 "$ENS_L2" "$K_L2"
run_L 4 "$ENS_L4" "$K_L4"
echo "### L2+L4 block-solve re-check done  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
