#!/usr/bin/env bash
# Chunk 1 compile-verify for the unified Nf-block HMC driver (hmc_Nfblocked_claude.cu).
# Builds three macro combos to exercise both LREFINE branches, DIR_NO_MASS, and NSTACK=2/3.
# Run from src/both_3d after: module load cuda/12.8 ; module load gcc/13.2.0
# Compile-only (no GPU needed to build). Output tee'd to hmc_Nfblocked_build_claude.log. No deletes.

set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=hmc_Nfblocked_build_claude.log
: > "$LOG"

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
SRC=hmc_Nfblocked_claude.cu

build () {
  # $1 = tag (-> output binary + log header), $2... = extra -D macros
  local tag="$1"; shift
  echo "==================== BUILD $tag ====================" | tee -a "$LOG"
  echo "$NVCC $NVCCFLAGS $INCLUDES $* $SRC -o hmc_block_${tag}.out" | tee -a "$LOG"
  if $NVCC $NVCCFLAGS $INCLUDES "$@" "$SRC" -o "hmc_block_${tag}.out" 2>&1 | tee -a "$LOG" ; then
    if [ -f "hmc_block_${tag}.out" ]; then
      echo "OK   $tag" | tee -a "$LOG"
    else
      echo "FAIL $tag (no binary produced)" | tee -a "$LOG"; exit 1
    fi
  else
    echo "FAIL $tag (nvcc error)" | tee -a "$LOG"; exit 1
  fi
  echo | tee -a "$LOG"
}

# L2, Nf6 (fermilab convention: dir includes mRe/mIm) -- the validation target
build L2_Nf6 -DLREFINE=2 -DNFPF=6

# L1, Nf6, local massless convention (DIR_NO_MASS + local geometry path default) -- exercises LREFINE==1
build L1_Nf6_local -DLREFINE=1 -DNFPF=6 -DDIR_NO_MASS

# L4, Nf4 -- exercises NSTACK=2
build L4_Nf4 -DLREFINE=4 -DNFPF=4

# Phase 2: FORCE also blocked (-DBLOCK_FORCE) -- exercises BlockedForce + MinimumNorm2BlockF
build L2_Nf6_bf -DLREFINE=2 -DNFPF=6 -DBLOCK_FORCE
build L4_Nf4_bf -DLREFINE=4 -DNFPF=4 -DBLOCK_FORCE

echo "==================== ALL BUILDS OK ====================" | tee -a "$LOG"
