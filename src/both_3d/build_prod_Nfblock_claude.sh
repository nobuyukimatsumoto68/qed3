#!/usr/bin/env bash
# Chunk 4 -- LOCAL production Nf-block binaries (action + FORCE blocked, -DBLOCK_FORCE).
# L1 x {Nf4, Nf6}, massless, local geometry, DIR_NO_MASS (matches the existing hmc_L1 output-dir
# convention -> AUTO-RESUMES the existing L1 Nf4/6 ensembles). Nf2 stays on the serial hmc_L1 + MPS.
# Production caps mirror hmc_L1_claude.cu: KMAX=320, KRNG=10.
# L2/L4 production runs on FNAL (remote builds with the fermilab geometry paths -- see the blackboard
# memo); a LOCAL L2/L4 build for testing is the two commented lines at the bottom.
#
# Run from src/both_3d after: module load cuda/12.8 ; module load gcc/13.2.0
# Compile-only (no GPU). Output tee'd to build_prod_Nfblock_claude.log. No deletes.

set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=build_prod_Nfblock_claude.log
: > "$LOG"

NVCC=nvcc
NVCCFLAGS="-t 4 -arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
SRC=hmc_Nfblocked_claude.cu

build () {
  local out="$1"; shift
  echo "==================== BUILD $out ====================" | tee -a "$LOG"
  echo "flags: $*" | tee -a "$LOG"
  if $NVCC $NVCCFLAGS $INCLUDES "$@" "$SRC" -o "$out" 2>&1 | tee -a "$LOG" ; then
    test -f "$out" || { echo "FAIL build $out" | tee -a "$LOG"; exit 1; }
    echo "OK build $out" | tee -a "$LOG"
  else
    echo "FAIL build $out (nvcc)" | tee -a "$LOG"; exit 1
  fi
}

# LOCAL production: L1, block force, massless local convention, production caps.
build hmc_Nfblk_L1_Nf4.out -DBLOCK_FORCE -DLREFINE=1 -DNFPF=4 -DDIR_NO_MASS -DKMAX=320 -DKRNG=10
build hmc_Nfblk_L1_Nf6.out -DBLOCK_FORCE -DLREFINE=1 -DNFPF=6 -DDIR_NO_MASS -DKMAX=320 -DKRNG=10

echo "==================== PROD BUILDS OK ====================" | tee -a "$LOG"
echo "Run e.g.: CUDA_VISIBLE_DEVICES=<gpu> ./hmc_Nfblk_L1_Nf6.out <gsq> 6 1.0   (resumes Nf6_gsq*...L1/)" | tee -a "$LOG"

# --- LOCAL L2/L4 test builds (production L2/L4 runs on FNAL with fermilab -DGEODESIC_H/-DGEOM_DATA):
# build hmc_Nfblk_L2_Nf6.out -DBLOCK_FORCE -DLREFINE=2 -DNFPF=6 -DKMAX=620 -DKRNG=5
# build hmc_Nfblk_L4_Nf6.out -DBLOCK_FORCE -DLREFINE=4 -DNFPF=6 -DKMAX=320 -DKRNG=5
