#!/usr/bin/env bash
# run_exactjj_free_a100_claude.sh -- STOCHASTIC exact conserved-current jj correlators in the FREE case
# (U=1), on an A100-80. Tests whether the EXACT sp sign flips with L (the cutoff-effect / sp-sign question).
# Matrix-free (multishift K on O(N) vectors) BUT the mrhs/BlockedMat scratch (~3*N*Nt*npole*16 B per pool)
# OOMs a 12 GB GPU at L=4 -> needs the A100-80.
#   L=1 = jj_corr_block_t_claude.cu (orig); L=2/4 = jj_corr_block_t_L{2,4}_claude.cu (only N_REFINE differs,
#   + output dir gets an _L<N_REFINE> suffix so the L's don't collide).
# Free mode = OMIT --ens-dir (U=1). Output: data_free_vmRe0.000000vmIm0.000000/corr_nt0<n>_nhits<H>_L<L>/
# Usage:  bash run_exactjj_free_a100_claude.sh "4" 8 2>&1 | tee exactjj_free_a100.log   (LS, NHITS; def "2 4", 8)
set -u
cd "$(dirname "$0")" || exit 1
module load cuda/12.8 2>/dev/null; module load gcc/13.2.0 2>/dev/null
export CUDA_VISIBLE_DEVICES=0                       # <-- the A100
LS="${1:-2 4}"; NHITS="${2:-8}"; NT0=2

NVCC=nvcc
# -arch=sm_80 for A100 (repo default is sm_70 for Titan V).
NVCCFLAGS="-arch=sm_80 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
# ADJUST INCLUDES/LDFLAGS to this server's Eigen / HighFive / HDF5 locations.

src_for_L(){ [ "$1" = 1 ] && echo jj_corr_block_t_claude.cu || echo jj_corr_block_t_L$1_claude.cu; }

for L in $LS; do
  SRC=$(src_for_L "$L"); BIN=jj_corr_block_t_L${L}.o
  echo "===== [L=$L] compile $SRC (sm_80) ====="
  $NVCC $NVCCFLAGS $INCLUDES $LDFLAGS "$SRC" -o "$BIN" || { echo "BUILD FAILED L=$L"; exit 1; }
  echo "===== [L=$L] run FREE stochastic jj (nhits=$NHITS, n_t0=$NT0) [$(date +%H:%M:%S)] ====="
  ./"$BIN" --nhits "$NHITS" --n-t0 "$NT0" || { echo "RUN FAILED L=$L"; exit 1; }
  echo "===== [L=$L] DONE [$(date +%H:%M:%S)] ====="
done
echo "ALL DONE -> data_free_vmRe0.000000vmIm0.000000/corr_nt0${NT0}_nhits${NHITS}_L<L>/"
