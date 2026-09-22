#!/usr/bin/env bash
# tmp_peram_mrhs_build_claude.sh -- BUILD ONLY (timed) of the mrhs + multi-source distillation perambulator
# variant distill_peram_mrhs_claude.cu.  Does NOT run the binary and does NOT touch the currently-running
# distill_peram_L1_claude.o / distill_peram_claude.cu.  NO rm, NO kill.
#
# Flags/includes match run_distill_L1_claude.sh (the proven distill build); the only new dependency,
# blocked_mat_claude.h, already lives under includes/ so -Iincludes covers it.
#
# Run yourself (foreground is fine; build is quick):
#   bash tmp_peram_mrhs_build_claude.sh
# then hand me back tmp_peram_mrhs_build_claude.log.
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1

L=1
SRC=distill_peram_mrhs_claude.cu
BIN=distill_peram_mrhs_L${L}_claude.o
NVCC=/usr/local/cuda-12.6/bin/nvcc
# link set matched to the proven distill/jj_sigma build (overlap/MatPoly + BlockedMat need cublas/cusolver/
# cusparse + GSL).
NVCCFLAGS="-w -arch=sm_70 -O2 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
LOG=tmp_peram_mrhs_build_claude.log

echo "================ peram-mrhs BUILD START $(date) ================" | tee "$LOG"
echo "---- source: $SRC -> $BIN ----" | tee -a "$LOG"
echo "---- NVCCFLAGS: $NVCCFLAGS ----" | tee -a "$LOG"

BUILD_T0=$(date +%s)
"$NVCC" "$SRC" $NVCCFLAGS $INCLUDES $LDFLAGS -o "$BIN" 2>&1 | tee -a "$LOG"
RC="${PIPESTATUS[0]}"
BUILD_T1=$(date +%s)

if [ "$RC" -ne 0 ]
then
  echo "BUILD FAILED (rc=$RC) after $((BUILD_T1-BUILD_T0)) s -- stopping (binary NOT produced)" | tee -a "$LOG"
  exit 1
fi

echo "build OK in $((BUILD_T1-BUILD_T0)) s -> $BIN  $(date)" | tee -a "$LOG"
ls -l "$BIN" 2>&1 | tee -a "$LOG"
echo "NOTE: run separately with e.g.  ./$BIN --Nf 2 --gsq 0.5 --nv 24 --nsrc 2   (or --tsrc-list 0,64)" | tee -a "$LOG"
echo "================ peram-mrhs BUILD DONE $(date) ================" | tee -a "$LOG"
