#!/usr/bin/env bash
# run_distill_sym_L3_Nf2_claude.sh -- SYMMETRIZED-basis (-DBASIS_SYM=1) perambulators for L3 Nf2 gsq1.5 (at0.2, massless).
#   Variant of run_distill_sym_L1L2_Nf2_claude.sh (that file is untouched): ONE ensemble, N_REFINE=3, and MULTI-GPU:
#   NPACK MPS-packed workers per GPU on every GPU in GPU_LIST, all STAGGERED on the stride grid so no two workers
#   ever touch the same k (worker WID: kmin = first + WID*STRIDE, stride = STRIDE*NWORKERS).
#   Output -> data_<ENS>/distill_Nv24_sym/  (--out-suffix _sym).  RESUMABLE: binary skips k whose peram exists.
#   NO rm.  NO kill.  Never touches _v2 / production dirs.
# Launch ONCE, DETACHED (NM runs this):
#   setsid nohup bash run_distill_sym_L3_Nf2_claude.sh > run_distill_sym_L3_Nf2_claude.log 2>&1 < /dev/null &
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1

ENS="Nf2_gsq1.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L3_hb0.400000-1.000000"
L=3
GPU_LIST="${GPU_LIST:-1}"          # NM 2026-09-22: GPU1 is the free one; GPU_LIST="0 1" to use both
NPACK="${NPACK:-2}"
STRIDE="${STRIDE:-2}"          # L3 has 999 configs -> stride 2 = ~500 perams (L1/L2 used stride 4 on 4000)
NV="${NV:-24}"
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
BIN="distill_peram_mrhs_L${L}_Nv${NV}_sym_claude.o"

build_bin () {
  local FLAGS="-w -arch=sm_70 -O3 -std=c++20 -DN_REFINE_CLI=${L} -DNSTACK_CLI=${NV} -DBASIS_SYM=1 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
  if [ ! -f "$BIN" ] || find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$BIN" -print -quit 2>/dev/null | grep -q .
  then
    echo "### compile SYM L${L} (N_REFINE=${L}, NSTACK=${NV}, BASIS_SYM=1) -> $BIN  [$(date +%F_%H:%M:%S)] ###"
    "$NVCC" "$SRC" $FLAGS $INCLUDES $LDFLAGS -o "$BIN" || { echo "### BUILD FAILED L${L} ###"; return 1; }
  else
    echo "### SYM L${L} Nv${NV} binary up-to-date ###"
  fi
}

start_mps () {
  if [ "$NPACK" -gt 1 ] && ! pgrep -f nvidia-cuda-mps-control >/dev/null; then
    echo "### starting MPS daemon ###"
    nvidia-cuda-mps-control -d
    for i in 1 2 3 4 5; do pgrep -f nvidia-cuda-mps-control >/dev/null && break; sleep 1; done
  fi
}

run_worker () {   # $1 = WID (global), $2 = GPU, $3 = NWORKERS
  local WID="$1" GPU="$2" NW="$3"
  local nf gsq at first last kmax wstride wstart
  nf=$(printf '%s' "$ENS" | grep -oE 'Nf[0-9]+' | head -1 | sed 's/Nf//')
  gsq=$(printf '%s' "$ENS" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//')
  at=$(printf '%s' "$ENS" | grep -oE 'at[0-9.]+nu0' | sed 's/at//;s/nu0//')
  first=$(ls "$ENS"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n | head -1)
  last=$( ls "$ENS"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n | tail -1)
  kmax=$(( last + 1 ))
  wstride=$(( STRIDE * NW ))
  wstart=$(( first + WID * STRIDE ))
  echo "### SYM ENS $ENS L${L} w${WID} GPU${GPU} kmin=${wstart} stride=${wstride} kmax=${kmax} (overall stride ${STRIDE}) nsrc=${NSRC} Nv=${NV} FUSED out${OUT_SUFFIX}  [$(date +%F_%H:%M:%S)] ###"
  CUDA_VISIBLE_DEVICES=$GPU ./"$BIN" --gsq "$gsq" --Nf "$nf" --nu0 1.0 --at "$at" \
    --nv "$NV" --tsrc0 "$TSRC0" --twin "$TWIN" --nsrc "$NSRC" --fused-sources --out-suffix "$OUT_SUFFIX" \
    --ens-dir "$ENS/" --kmin "$wstart" --stride "$wstride" --kmax "$kmax"
  echo "### SYM ENS $ENS L${L} w${WID} DONE (status $?)  [$(date +%F_%H:%M:%S)] ###"
}

if [ ! -d "$ENS" ] || ! ls "$ENS"/ckpoint_lat.* >/dev/null 2>&1; then echo "### ABORT: no ckpoint in $ENS ###"; exit 1; fi
echo "### START SYM L${L} regen: GPUs [${GPU_LIST}] x ${NPACK} MPS each, STRIDE=${STRIDE}, Nv=${NV}, BASIS_SYM=1, out${OUT_SUFFIX}  [$(date +%F_%H:%M:%S)] ###"
build_bin || exit 1
start_mps
NGPU=$(echo $GPU_LIST | wc -w)
NW=$(( NGPU * NPACK ))
WID=0
for GPU in $GPU_LIST; do
  for (( p=0; p<NPACK; p++ )); do
    ( run_worker "$WID" "$GPU" "$NW" ) >> "distill_sym_L${L}_w${WID}_claude.log" 2>&1 &
    WID=$(( WID + 1 ))
  done
done
wait
echo "### SYM L${L} ALL ${NW} WORKERS DONE  [$(date +%F_%H:%M:%S)] ###"
