#!/usr/bin/env bash
# Chunk 4 -- MPS worth-it benchmark for L in {1,4} x Nf in {2,6}, MASSLESS.
# For each config: measure 1 solo process (T1) vs 2 concurrent MPS-packed processes (T2 each), on ONE
# GPU. MPS throughput gain = 2*T1 / avg(T2) : >1 => packing 2 clients beats a single process (MPS worth
# it); ~1 or <1 => MPS does not help (the process already fills the GPU, or contention dominates).
#   Nf=2 = production SERIAL path (-DSERIAL_REF, NSTACK=1); Nf=6 = block ACTION+FORCE (-DBLOCK_FORCE).
# 1 trajectory per measurement (KMAX=2). gsq=8, nu0=1, mass=0. GPU set by GPU= below (default 1).
#
# ASSUMES the target GPU has NO other COMPUTE jobs (graphics/`G` jobs are fine). Starts the MPS daemon
# on that GPU and LEAVES IT UP (does not tear it down -- that is global). Fresh scratch dirs one level
# under src/ (so ../../geometry resolves; production dirs untouched). SINGLE-USE: delete the mps_*_claude
# dirs before a re-run. No deletes here.  Run after: module load cuda/12.8 ; module load gcc/13.2.0.
# ~2 h wall (the L4 pairs dominate). Output tee'd to bench_mps_Nfblock_claude.log.

set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
VROOT=/mnt/barracuda22/qed3/qed3/src
GPU=1

cd "$SRCDIR"
LOG="$SRCDIR/bench_mps_Nfblock_claude.log"
: > "$LOG"

NVCC=nvcc
NVCCFLAGS="-t 4 -arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
SRC=hmc_Nfblocked_claude.cu

build () {
  local out="$1"; shift
  echo "==================== BUILD $out ====================" | tee -a "$LOG"
  if $NVCC $NVCCFLAGS $INCLUDES "$@" "$SRC" -o "$out" 2>&1 | tee -a "$LOG" ; then
    test -f "$out" || { echo "FAIL build $out" | tee -a "$LOG"; exit 1; }
    echo "OK build $out" | tee -a "$LOG"
  else
    echo "FAIL build $out (nvcc)" | tee -a "$LOG"; exit 1
  fi
}

build bench_L1_Nf2.out -DSERIAL_REF  -DLREFINE=1 -DNFPF=2 -DDIR_NO_MASS -DKMAX=2
build bench_L1_Nf6.out -DBLOCK_FORCE -DLREFINE=1 -DNFPF=6 -DDIR_NO_MASS -DKMAX=2
build bench_L4_Nf2.out -DSERIAL_REF  -DLREFINE=4 -DNFPF=2 -DDIR_NO_MASS -DKMAX=2
build bench_L4_Nf6.out -DBLOCK_FORCE -DLREFINE=4 -DNFPF=6 -DDIR_NO_MASS -DKMAX=2

# --- use the AMBIENT MPS (do NOT start a second daemon / custom pipe dir -- that conflicts with an
# existing default-pipe daemon and crashes CUDA init). All GPU-$GPU processes route through whatever
# MPS server is already up (the user's run scripts keep one running). Just pin the GPU. ---
export CUDA_VISIBLE_DEVICES=$GPU
if pgrep -f nvidia-cuda-mps-server > /dev/null 2>&1 ; then
  echo "# MPS server detected -> packed runs will share the GPU via MPS" | tee -a "$LOG"
else
  echo "# WARNING: no MPS server running -> packed runs will TIME-SLICE (MPS gain will read ~1x or <1x)" | tee -a "$LOG"
fi

first_hmc () { grep '# HMC :' "$1" 2>/dev/null | awk 'NR==1{print $4}'; }

# $1=tag $2=binary $3=gsq $4=Nf
bench_one () {
  local tag="$1"; local bin="$2"; local gsq="$3"; local nf="$4"
  local sd="$VROOT/mps_${tag}_solo_claude"
  local p0="$VROOT/mps_${tag}_pk0_claude"
  local p1="$VROOT/mps_${tag}_pk1_claude"
  mkdir -p "$sd" "$p0" "$p1"

  echo "-------------------- $tag : SOLO (1 client) --------------------" | tee -a "$LOG"
  ( cd "$sd" && "$SRCDIR/$bin" "$gsq" "$nf" 1.0 > run.out 2>&1 )
  local t1=$(first_hmc "$sd/run.out")

  echo "-------------------- $tag : PACKED (2 MPS clients) --------------------" | tee -a "$LOG"
  ( cd "$p0" && "$SRCDIR/$bin" "$gsq" "$nf" 1.0 > run.out 2>&1 ) &
  local j0=$!
  ( cd "$p1" && "$SRCDIR/$bin" "$gsq" "$nf" 1.0 > run.out 2>&1 ) &
  local j1=$!
  wait "$j0"
  wait "$j1"
  local t2a=$(first_hmc "$p0/run.out")
  local t2b=$(first_hmc "$p1/run.out")

  echo "RESULT $tag : solo=$t1  packed=($t2a,$t2b)" | tee -a "$LOG"
  awk -v tag="$tag" -v t1="$t1" -v a="$t2a" -v b="$t2b" 'BEGIN{
    if(t1+0<=0 || a+0<=0 || b+0<=0){
      printf "  %-8s : incomplete (a run produced no # HMC line -- see run.out)\n", tag;
    } else {
      t2=(a+b)/2.0; g=2.0*t1/t2;
      verdict = (g>1.05 ? "MPS WORTH IT" : (g<0.95 ? "MPS HURTS" : "~neutral"));
      printf "  %-8s : solo %8.1f s | packed(2) avg %8.1f s | throughput gain %.2fx | %s\n", tag, t1, t2, g, verdict;
    }
  }' | tee -a "$LOG"
}

echo "" | tee -a "$LOG"
echo "==================== MPS BENCH (gsq=8, massless, GPU $GPU) ====================" | tee -a "$LOG"
bench_one L1_Nf2 bench_L1_Nf2.out 8 2
bench_one L1_Nf6 bench_L1_Nf6.out 8 6
bench_one L4_Nf2 bench_L4_Nf2.out 8 2
bench_one L4_Nf6 bench_L4_Nf6.out 8 6

echo "" | tee -a "$LOG"
echo "==================== MPS BENCH DONE ====================" | tee -a "$LOG"
echo "# gain = 2*T_solo/T_packed ; >1.05 => MPS packing (2 clients) raises throughput." | tee -a "$LOG"
