#!/usr/bin/env bash
# Chunk 2 -- validate the Nf-block action-inversion path against the SERIAL reference.
# Same source (hmc_Nfblocked_claude.cu): block build vs -DSERIAL_REF build, IDENTICAL params, so a
# cold start + seeded RNG gives the SAME momentum + same phi -> the FIRST-trajectory dH must agree to
# ~solver tolerance (block-CG shared-stop only over-converges); later trajectories diverge chaotically
# but track statistically. Also reports serial-vs-block per-trajectory wall time (the speedup).
#
# Real target: L4/Nf6. A fast L1/Nf6 sanity runs first. Massless (mass=0), GPU1, local geometry.
# Runs in SIBLING scratch dirs one level under src/ (so the drivers' relative ../../geometry/ resolves)
# and NEVER writes into src/both_3d (protects production output dirs). Compile+run; USER runs this.
# SINGLE-USE: to re-run, delete the two scratch dirs first (checkpoints would otherwise auto-resume).
# No deletes anywhere in this script.
#
# Run after: module load cuda/12.8 ; module load gcc/13.2.0

set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
SER=/mnt/barracuda22/qed3/qed3/src/valNfblk_serial_claude
BLK=/mnt/barracuda22/qed3/qed3/src/valNfblk_block_claude
mkdir -p "$SER"
mkdir -p "$BLK"

LOG="$SRCDIR/run_Nfblock_validate_claude.log"
: > "$LOG"

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
SRC=hmc_Nfblocked_claude.cu

cd "$SRCDIR"

build () {
  # $1 = out binary ; $2... = -D macros
  local out="$1"; shift
  echo "==================== BUILD $out ====================" | tee -a "$LOG"
  if $NVCC $NVCCFLAGS $INCLUDES "$@" "$SRC" -o "$out" 2>&1 | tee -a "$LOG" ; then
    test -f "$out" || { echo "FAIL build $out" | tee -a "$LOG"; exit 1; }
    echo "OK build $out" | tee -a "$LOG"
  else
    echo "FAIL build $out (nvcc)" | tee -a "$LOG"; exit 1
  fi
}

# --- builds (block vs serial-ref), Nf6, massless dir convention ---
build valblk_L1.out              -DLREFINE=1 -DNFPF=6 -DDIR_NO_MASS -DKMAX=9
build valser_L1.out  -DSERIAL_REF -DLREFINE=1 -DNFPF=6 -DDIR_NO_MASS -DKMAX=9
build valblk_L4.out              -DLREFINE=4 -DNFPF=6 -DDIR_NO_MASS -DKMAX=6
build valser_L4.out  -DSERIAL_REF -DLREFINE=4 -DNFPF=6 -DDIR_NO_MASS -DKMAX=6

# --- runs : each in its own scratch parent so ../../geometry resolves and dir3 is isolated ---
run_one () {
  # $1 = parent dir ; $2 = binary ; $3 = stdout file ; $4 = timeout sec
  local pdir="$1"; local bin="$2"; local ofile="$3"; local tmo="$4"
  echo "-------------------- RUN $bin (GPU1) --------------------" | tee -a "$LOG"
  ( cd "$pdir" && CUDA_VISIBLE_DEVICES=1 timeout "$tmo" "$SRCDIR/$bin" 8 6 1.0 > "$ofile" 2>&1 )
  echo "exit=$? -> $pdir/$ofile" | tee -a "$LOG"
}

run_one "$SER" valser_L1.out ser_L1.out 900
run_one "$BLK" valblk_L1.out blk_L1.out 900
run_one "$SER" valser_L4.out ser_L4.out 10800
run_one "$BLK" valblk_L4.out blk_L4.out 10800

# --- compare ---
cmp_pair () {
  local tag="$1"; local serf="$2"; local blkf="$3"
  echo "" | tee -a "$LOG"
  echo "==================== COMPARE $tag ====================" | tee -a "$LOG"
  echo "-- dH per trajectory :  serial   block --" | tee -a "$LOG"
  paste <(grep '# dH :' "$serf" | awk '{print $4}') <(grep '# dH :' "$blkf" | awk '{print $4}') | head -12 | tee -a "$LOG"
  local s1=$(grep '# dH :' "$serf" | awk 'NR==1{print $4}')
  local b1=$(grep '# dH :' "$blkf" | awk 'NR==1{print $4}')
  awk -v s="$s1" -v b="$b1" 'BEGIN{ d=s-b; if(d<0)d=-d; printf "first-traj |dH_ser - dH_blk| = %.3e   (ser=%s blk=%s)\n", d, s, b }' | tee -a "$LOG"
  local savg=$(grep '# HMC :' "$serf" | awk '{s+=$4; n++} END{ if(n>0) printf "%.4f", s/n; else print "0" }')
  local bavg=$(grep '# HMC :' "$blkf" | awk '{s+=$4; n++} END{ if(n>0) printf "%.4f", s/n; else print "0" }')
  awk -v s="$savg" -v b="$bavg" 'BEGIN{ if(b>0) printf "per-traj: serial %.4f s | block %.4f s | speedup %.2fx\n", s, b, s/b; else print "timing NA" }' | tee -a "$LOG"
}

cmp_pair L1_Nf6 "$SER/ser_L1.out" "$BLK/blk_L1.out"
cmp_pair L4_Nf6 "$SER/ser_L4.out" "$BLK/blk_L4.out"

echo "" | tee -a "$LOG"
echo "==================== VALIDATION DONE ====================" | tee -a "$LOG"
