#!/usr/bin/env bash
# run_exactjj_free_claude.sh -- STOCHASTIC exact conserved-current jj correlators in the FREE case (U=1),
# at L=1,2,4. Tests whether the EXACT sp sign flips with L (the cutoff-effect question). Matrix-free
# (multishift K on O(N) vectors) -> fits the 12 GB Titan V even at L=4 (NO dense matrices, no A100).
#   L=1 = jj_corr_block_t_claude.cu (orig, N_REFINE=1); L=2/4 = jj_corr_block_t_L{2,4}_claude.cu copies.
# Free mode = OMIT --ens-dir (U=1). Output: data_free_vmRe0.000000vmIm0.000000/corr_nt0<n>_nhits<H>/corr.<k>.h5
# Usage:  bash run_exactjj_free_claude.sh "1 2 4" 8 2>&1 | tee exactjj_free_claude.log
#         (args: LS="1 2 4"  NHITS=8 ; defaults "1 2" and 8)
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
module load cuda/12.8 2>/dev/null; module load gcc/13.2.0 2>/dev/null
LS="${1:-1 2}"; NHITS="${2:-8}"; GPU="${3:-0}"; NT0=2
export CUDA_VISIBLE_DEVICES=$GPU

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

src_for_L(){ [ "$1" = 1 ] && echo jj_corr_block_t_claude.cu || echo jj_corr_block_t_L$1_claude.cu; }

for L in $LS; do
  SRC=$(src_for_L "$L"); BIN=jj_corr_block_t_L${L}.o
  echo "===== [L=$L] compile $SRC ====="
  $NVCC $NVCCFLAGS $INCLUDES $LDFLAGS "$SRC" -o "$BIN" || { echo "BUILD FAILED L=$L"; exit 1; }
  echo "===== [L=$L] run FREE stochastic jj (nhits=$NHITS, n_t0=$NT0) [$(date +%H:%M:%S)] ====="
  ./"$BIN" --nhits "$NHITS" --n-t0 "$NT0" || { echo "RUN FAILED L=$L"; exit 1; }
  echo "===== [L=$L] DONE [$(date +%H:%M:%S)] ====="
done
echo "ALL DONE -> data_free_vmRe0.000000vmIm0.000000/corr_nt0${NT0}_nhits${NHITS}/"
