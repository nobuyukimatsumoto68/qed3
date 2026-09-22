#!/usr/bin/env bash
# run_distill_sym_L1L2_Nf2_claude.sh -- regenerate perambulators with the SYMMETRIZED distillation basis
#   (-DBASIS_SYM=1: basis = low modes of D_W^H D_W + D_W D_W^H, which respects the theory symmetry so a
#   TRUNCATED basis no longer leaks the single meson m_PS into sigma^2).  Nf2, at0.2, smallest-gsq each L:
#     L1 gsq=0.5  (N_REFINE=1, Nv=24 = 2*12 = complete)
#     L2 gsq=1.0  (N_REFINE=2, Nv=24 = truncated 24-of-84 -- the case the fix targets)
#   Output -> data_<ENS>/distill_Nv24_sym/  (--out-suffix _sym; NEVER touches the existing _v2 / production dirs).
#   STRIDE=4 from k=first.  nsrc=2 fused, twin=32, tsrc0=0.  RESUMABLE, skip-protected.  NO rm.  NO kill.
#
# MPS-packed, 2 workers/GPU, STAGGERED BY STRIDE START so the workers interleave the stride-4 grid with NO
#   overwrite: worker0 = k=first, first+8, ... ; worker1 = k=first+4, first+12, ... (each stride=2*STRIDE).
# L1 runs FIRST (both workers), then L2 (both workers).
# Launch ONCE.  DETACHED:
#   GPU=1 setsid nohup bash run_distill_sym_L1L2_Nf2_claude.sh > run_distill_sym_L1L2_Nf2_claude.log 2>&1 < /dev/null &
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1

ENS_L1="Nf2_gsq0.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000"
ENS_L2="Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000"

GPU="${GPU:-1}"
NPACK=2
STRIDE="${STRIDE:-4}"
NV=24
TSRC0=0
TWIN=32
NSRC=2
OUT_SUFFIX=_sym

NVCC=/usr/local/cuda-12.6/bin/nvcc
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
export LD_LIBRARY_PATH=/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib:${LD_LIBRARY_PATH:-}
SRC=distill_peram_mrhs_claude.cu

build_bin () {   # $1 = L (N_REFINE) -> echoes binary name; builds with BASIS_SYM=1
  local L="$1"
  local BIN="distill_peram_mrhs_L${L}_Nv${NV}_sym_claude.o"
  local FLAGS="-w -arch=sm_70 -O3 -std=c++20 -DN_REFINE_CLI=${L} -DNSTACK_CLI=${NV} -DBASIS_SYM=1 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
  if [ ! -f "$BIN" ] || find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$BIN" -print -quit 2>/dev/null | grep -q .
  then
    echo "### compile SYM L${L} (N_REFINE=${L}, NSTACK=${NV}, BASIS_SYM=1) -> $BIN  [$(date +%F_%H:%M:%S)] ###" >&2
    "$NVCC" "$SRC" $FLAGS $INCLUDES $LDFLAGS -o "$BIN" >&2 || { echo "### BUILD FAILED L${L} ###" >&2; return 1; }
  else
    echo "### SYM L${L} Nv${NV} binary up-to-date ###" >&2
  fi
  echo "$BIN"
}

start_mps () {
  if [ "$NPACK" -gt 1 ] && ! pgrep -f nvidia-cuda-mps-control >/dev/null; then
    echo "### starting MPS daemon ###"
    nvidia-cuda-mps-control -d
    for i in 1 2 3 4 5; do pgrep -f nvidia-cuda-mps-control >/dev/null && break; sleep 1; done
  fi
}

run_stagger () {   # $1 = ENS, $2 = L, $3 = BIN, $4 = WID (0 or 1)
  local ENS="$1" L="$2" BIN="$3" WID="$4"
  if [ ! -d "$ENS" ] || ! ls "$ENS"/ckpoint_lat.* >/dev/null 2>&1; then echo "### SKIP (no ckpoint): $ENS ###"; return; fi
  local nf gsq at first last kmax wstride wstart
  nf=$(printf '%s' "$ENS" | grep -oE 'Nf[0-9]+' | head -1 | sed 's/Nf//')
  gsq=$(printf '%s' "$ENS" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//')
  at=$(printf '%s' "$ENS" | grep -oE 'at[0-9.]+nu0' | sed 's/at//;s/nu0//')
  first=$(ls "$ENS"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n | head -1)
  last=$( ls "$ENS"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n | tail -1)
  kmax=$(( last + 1 ))
  wstride=$(( STRIDE * NPACK ))                 # each worker advances by NPACK*STRIDE
  wstart=$(( first + WID * STRIDE ))            # staggered start: w0=first, w1=first+STRIDE, ...
  echo "### SYM ENS $ENS  L${L} w${WID} kmin=${wstart} stride=${wstride} kmax=${kmax} (overall stride ${STRIDE}) nsrc=${NSRC} Nv=${NV} FUSED out${OUT_SUFFIX}  [$(date +%F_%H:%M:%S)] ###"
  CUDA_VISIBLE_DEVICES=$GPU ./"$BIN" --gsq "$gsq" --Nf "$nf" --nu0 1.0 --at "$at" \
    --nv "$NV" --tsrc0 "$TSRC0" --twin "$TWIN" --nsrc "$NSRC" --fused-sources --out-suffix "$OUT_SUFFIX" \
    --ens-dir "$ENS/" --kmin "$wstart" --stride "$wstride" --kmax "$kmax"
  echo "### SYM ENS $ENS L${L} w${WID} DONE (status $?)  [$(date +%F_%H:%M:%S)] ###"
}

run_L () {   # $1 = ENS, $2 = L : build binary, launch NPACK staggered workers, wait
  local ENS="$1" L="$2" BIN w
  BIN=$(build_bin "$L") || { echo "### build failed L${L}, skipping ###"; return 1; }
  for (( w=0; w<NPACK; w++ )); do
    ( run_stagger "$ENS" "$L" "$BIN" "$w" ) >> "distill_sym_L${L}_w${w}_claude.log" 2>&1 &
  done
  wait
  echo "### SYM L${L} ALL WORKERS DONE  [$(date +%F_%H:%M:%S)] ###"
}

echo "### START SYM regen: GPU${GPU} x ${NPACK} MPS, STRIDE=${STRIDE}, BASIS_SYM=1, out${OUT_SUFFIX}  [$(date +%F_%H:%M:%S)] ###"
start_mps
echo "### === L1 (Nf2 gsq0.5, N_REFINE=1, Nv24 complete) FIRST ==="
run_L "$ENS_L1" 1
echo "### === L2 (Nf2 gsq1.0, N_REFINE=2, Nv24 truncated) NEXT ==="
run_L "$ENS_L2" 2
echo "### SYM regen ALL DONE  [$(date +%F_%H:%M:%S)] ###"
