#!/usr/bin/env bash
# Chunk 2 (force) -- validate the Nf-block FORCE against the action-only block (serial force).
# Both builds use the SAME (already-validated) block ACTION solve; they differ ONLY in the force path
# (block BlockedForce vs serial per-flavor grad_l4). From an identical cold start + seeded RNG the
# momentum, phi and eta are identical, so the FIRST-trajectory dH must agree to ~1e-6 -> that isolates
# and validates the block force. massless, GPU1, local geometry.
#
# L1/Nf6: a few trajectories (fast). L4/Nf6: 1 trajectory (KMAX=2) -- L4 is slow (per NM).
# Runs in FRESH sibling scratch dirs one level under src/ (so ../../geometry/ resolves; production dirs
# in src/both_3d are untouched). SINGLE-USE: to re-run, delete the scratch dirs first. No deletes here.
#
# Run after: module load cuda/12.8 ; module load gcc/13.2.0

set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
REF=/mnt/barracuda22/qed3/qed3/src/valNfblk_fref_claude    # action-only block (serial force) = reference
BLK=/mnt/barracuda22/qed3/qed3/src/valNfblk_fblk_claude    # block force (-DBLOCK_FORCE)
mkdir -p "$REF"
mkdir -p "$BLK"

LOG="$SRCDIR/run_Nfblock_forceval_claude.log"
: > "$LOG"

NVCC=nvcc
NVCCFLAGS="-t 4 -arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
SRC=hmc_Nfblocked_claude.cu

cd "$SRCDIR"

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

# reference = action-only block (serial force) ; test = block force. KMAX=2 => 1 trajectory (the
# first-traj dH is the isolating correctness check; L1 timing is NOT the metric -- see L4).
build fref_L1.out              -DLREFINE=1 -DNFPF=6 -DDIR_NO_MASS -DKMAX=2
build fblk_L1.out -DBLOCK_FORCE -DLREFINE=1 -DNFPF=6 -DDIR_NO_MASS -DKMAX=2
build fref_L4.out              -DLREFINE=4 -DNFPF=6 -DDIR_NO_MASS -DKMAX=2
build fblk_L4.out -DBLOCK_FORCE -DLREFINE=4 -DNFPF=6 -DDIR_NO_MASS -DKMAX=2

run_one () {
  local pdir="$1"; local bin="$2"; local ofile="$3"; local tmo="$4"
  echo "-------------------- RUN $bin (GPU1) --------------------" | tee -a "$LOG"
  ( cd "$pdir" && CUDA_VISIBLE_DEVICES=1 timeout "$tmo" "$SRCDIR/$bin" 8 6 1.0 > "$ofile" 2>&1 )
  echo "exit=$? -> $pdir/$ofile" | tee -a "$LOG"
}

run_one "$REF" fref_L1.out ref_L1.out 1800
run_one "$BLK" fblk_L1.out blk_L1.out 1800
run_one "$REF" fref_L4.out ref_L4.out 14400
run_one "$BLK" fblk_L4.out blk_L4.out 14400

cmp_pair () {
  local tag="$1"; local serf="$2"; local blkf="$3"
  echo "" | tee -a "$LOG"
  echo "==================== COMPARE $tag  (action-only serial-force | block-force) ====================" | tee -a "$LOG"
  echo "-- dH per trajectory --" | tee -a "$LOG"
  paste <(grep '# dH :' "$serf" | awk '{print $4}') <(grep '# dH :' "$blkf" | awk '{print $4}') | head -8 | tee -a "$LOG"
  local s1=$(grep '# dH :' "$serf" | awk 'NR==1{print $4}')
  local b1=$(grep '# dH :' "$blkf" | awk 'NR==1{print $4}')
  awk -v s="$s1" -v b="$b1" 'BEGIN{ d=s-b; if(d<0)d=-d; printf "first-traj |dH_serialforce - dH_blockforce| = %.3e   (ref=%s blk=%s)\n", d, s, b }' | tee -a "$LOG"
  echo "-- per-traj sec (force-block speedup) --" | tee -a "$LOG"
  local sa=$(grep '# HMC :' "$serf" | awk '{s+=$4;n++} END{if(n>0)printf "%.4f",s/n; else print "0"}')
  local ba=$(grep '# HMC :' "$blkf" | awk '{s+=$4;n++} END{if(n>0)printf "%.4f",s/n; else print "0"}')
  awk -v s="$sa" -v b="$ba" 'BEGIN{ if(b>0) printf "serial-force %.4f s | block-force %.4f s | speedup %.2fx\n", s, b, s/b; else print "timing NA" }' | tee -a "$LOG"
}

cmp_pair L1_Nf6 "$REF/ref_L1.out" "$BLK/blk_L1.out"
cmp_pair L4_Nf6 "$REF/ref_L4.out" "$BLK/blk_L4.out"

echo "" | tee -a "$LOG"
echo "==================== FORCE VALIDATION DONE ====================" | tee -a "$LOG"
