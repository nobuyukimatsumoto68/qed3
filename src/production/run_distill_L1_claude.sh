#!/usr/bin/env bash
# run_distill_L1_claude.sh -- CHUNK 1+2 build + validation of exact distillation: BASIS V(t) = low modes of
# the timeslice Wilson normal operator D_{W,2}^dag D_{W,2}, AND the PERAMBULATOR tau = V^dag D_ov^{-1} V over
# the source window [tsrc0, tsrc0+twin).  Plan: distillation_impl_plan_claude.md.
# Driver: distill_peram_claude.cu.  L1 (N_SITES=12 -> 2 N_s = 24 = full/exact basis).  Overlap solves on GPU1.
#
# Two phases:
#   (1) FREE FIELD (U=1): nv=24, twin=32 -> T1a-T1d (basis) + T2a/T2b/T2c (peram: exactness, IV.17, tau').
#   (2) ONE real L1 config (first discovered at0.2 massless ensemble): same, on a gauge background.
#
# Basis math is Eigen/CPU; the perambulator overlap solves run on GPU1 (CUDA_VISIBLE_DEVICES=1).  NO rm.  NO
# kill.  Launch detached:
#   setsid nohup bash run_distill_L1_claude.sh > run_distill_L1_claude.log 2>&1 < /dev/null &
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export CUDA_VISIBLE_DEVICES=1        # use GPU1 (chunk-2 overlap solves); device 0 in-code = physical GPU1

L=1
NV=24                          # full rank at L1 (= 2 * N_SITES = exact local)
AT_TOKEN="${AT_TOKEN:-at0.200000}"
ENS_GLOB="${ENS_GLOB:-Nf*_gsq*${AT_TOKEN}*mRe0.000000mIm0.000000nt128L${L}_hb*}"
KMIN="${KMIN:-1}"
STRIDE="${STRIDE:-10}"
NCFG_TEST="${NCFG_TEST:-3}"    # number of real configs to validate in phase 2

SRC=distill_peram_claude.cu
BIN=distill_peram_L${L}_claude.o
NVCC=/usr/local/cuda-12.6/bin/nvcc
# link set matched to the proven jj_sigma build (overlap/MatPoly need cublas/cusolver/cusparse + GSL).
NVCCFLAGS="-w -arch=sm_70 -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
LOG=run_distill_L1_claude.log

echo "================ distill chunk1 START $(date) ================" | tee -a "$LOG"

# ---- build ----
echo "---- compile L${L} -> $BIN  $(date) ----" | tee -a "$LOG"
"$NVCC" "$SRC" $NVCCFLAGS $INCLUDES $LDFLAGS -o "$BIN" 2>&1 | tee -a "$LOG"
if [ "${PIPESTATUS[0]}" -ne 0 ]
then
  echo "L${L} BUILD FAILED -- stopping" | tee -a "$LOG"
  exit 1
fi
echo "build OK $(date)" | tee -a "$LOG"

# ---- phase 1: free field (U=1), full validation ----
echo "---- phase 1: FREE FIELD nv=$NV + T1a/T1b/T1c/T1d  $(date) ----" | tee -a "$LOG"
./"$BIN" --Nf 2 --gsq 0.5 --nu0 1.0 --nv "$NV" --gauge-check 2>&1 | tee -a "$LOG"

# ---- phase 2: first real L1 ensemble, a few configs ----
ENS=""
for d in $ENS_GLOB
do
  [ -d "$d" ] || continue
  ls "$d"/ckpoint_lat.* >/dev/null 2>&1 || continue
  ENS="$d"
  break
done
if [ -z "$ENS" ]
then
  echo "---- phase 2 SKIPPED: no L${L} ${AT_TOKEN} massless ensemble with configs found ----" | tee -a "$LOG"
else
  nf=$(printf '%s' "$ENS" | grep -oE '^Nf[0-9]+' | sed 's/Nf//')
  gsq=$(printf '%s' "$ENS" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//')
  KMAX=$(( KMIN + STRIDE*NCFG_TEST ))
  echo "---- phase 2: real config basis  ens=$ENS  Nf=$nf gsq=$gsq  k in [$KMIN,$KMAX) stride $STRIDE  $(date) ----" | tee -a "$LOG"
  ./"$BIN" --Nf "$nf" --gsq "$gsq" --nu0 1.0 --nv "$NV" --ens-dir "$ENS"/ \
           --kmin "$KMIN" --kmax "$KMAX" --stride "$STRIDE" --gauge-check 2>&1 | tee -a "$LOG"
fi

echo "================ distill chunk1 DONE $(date) ================" | tee -a "$LOG"
