#!/usr/bin/env bash
# run_distill_L2_Nf2_finestride_claude.sh -- MORE statistics on the Nf2 L2 ensembles (instead of Nf4).
# Narrows the stride 10 -> STRIDE (default 5): doubles each Nf2 ensemble to ~800 configs by ADDING the
# fill-in points (k=6,16,26,...); the existing stride-10 perams (k=1,11,...) are skip-protected (kept).
# Nv=24 truncated + NSTACK=24 + nsrc=2 fused (same as the grid).  Nf4 is DROPPED for now (its partial
# perams stay on disk, resumable later).  RESUMABLE, skip-protected.  NO rm.  NO kill.
#
# BALANCED 2-worker split: both MPS workers process EVERY ensemble together, each taking one HALF of the
# stride-5 config range (worker0 = lower half [first,mid), worker1 = upper half [mid,last]).  Because the
# already-done configs are spread across both halves, each worker computes ~half the NEW configs -> no idle
# worker (the offset-0/offset-5 tiling would have made worker0 skip everything for g1/g2).
# Launch ONCE.  DETACHED:
#   GPU=1 setsid nohup bash run_distill_L2_Nf2_finestride_claude.sh > run_distill_L2_Nf2_finestride_claude.log 2>&1 < /dev/null &
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1

ENS_LIST=(
  "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000"
  "Nf2_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000"
  "Nf2_gsq3.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000"
)

GPU="${GPU:-1}"
NPACK=2
STRIDE="${STRIDE:-5}"           # narrowed from 10 -> ~800 cfg/ensemble (2x); set STRIDE=4 etc. for finer
NV=24
TSRC0=0
TWIN=32
NSRC=2
L=2

NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O3 -std=c++20 -DN_REFINE_CLI=2 -DNSTACK_CLI=${NV} -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
export LD_LIBRARY_PATH=/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib:${LD_LIBRARY_PATH:-}
SRC=distill_peram_mrhs_claude.cu
BIN="distill_peram_mrhs_L${L}_Nv${NV}_claude.o"

if [ ! -f "$BIN" ] || find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$BIN" -print -quit 2>/dev/null | grep -q .
then
  echo "### compile mrhs L${L} (N_REFINE=2, NSTACK=${NV}) -> $BIN  [$(date +%F_%H:%M:%S)] ###"
  "$NVCC" "$SRC" $NVCCFLAGS $INCLUDES $LDFLAGS -o "$BIN" || { echo "### BUILD FAILED ###"; exit 1; }
else
  echo "### mrhs L${L} Nv${NV} binary up-to-date ###"
fi

if [ "$NPACK" -gt 1 ] && ! pgrep -f nvidia-cuda-mps-control >/dev/null
then
  echo "### starting MPS daemon ###"
  nvidia-cuda-mps-control -d
  for i in 1 2 3 4 5; do pgrep -f nvidia-cuda-mps-control >/dev/null && break; sleep 1; done
fi

run_half () {   # $1 = ENS, $2 = WID (0 lower half, 1 upper half)
  local ENS="$1" WID="$2"
  if [ ! -d "$ENS" ] || ! ls "$ENS"/ckpoint_lat.* >/dev/null 2>&1; then echo "### SKIP (no ckpoint): $ENS ###"; return; fi
  local nf gsq at first last kmax nstr mid
  nf=$(printf '%s' "$ENS" | grep -oE 'Nf[0-9]+' | head -1 | sed 's/Nf//')
  gsq=$(printf '%s' "$ENS" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//')
  at=$(printf '%s' "$ENS" | grep -oE 'at[0-9.]+nu0' | sed 's/at//;s/nu0//')
  first=$(ls "$ENS"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n | head -1)
  last=$( ls "$ENS"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n | tail -1)
  kmax=$(( last + 1 ))
  nstr=$(( (last - first) / STRIDE + 1 ))            # total stride-STRIDE configs
  mid=$(( first + STRIDE * (nstr / 2) ))             # split index (aligned to the stride grid)
  local wk_kmin wk_kmax
  if [ "$WID" -eq 0 ]; then wk_kmin=$first; wk_kmax=$mid; else wk_kmin=$mid; wk_kmax=$kmax; fi
  echo "### ENS $ENS  w${WID} half [${wk_kmin},${wk_kmax}) stride=${STRIDE} nsrc=${NSRC} Nv=${NV} FUSED  [$(date +%F_%H:%M:%S)] ###"
  CUDA_VISIBLE_DEVICES=$GPU ./"$BIN" --gsq "$gsq" --Nf "$nf" --nu0 1.0 --at "$at" \
    --nv "$NV" --tsrc0 "$TSRC0" --twin "$TWIN" --nsrc "$NSRC" --fused-sources \
    --ens-dir "$ENS/" --kmin "$wk_kmin" --stride "$STRIDE" --kmax "$wk_kmax"
  echo "### ENS $ENS w${WID} DONE (status $?)  [$(date +%F_%H:%M:%S)] ###"
}

run_worker () {   # $1 = WID ; do every ensemble (its half), in order
  local WID="$1"
  local LOG="distill_L2_Nf2_finestride_w${WID}_claude.log"
  {
    echo "### START L2 Nf2 finestride worker ${WID}/${NPACK}  GPU${GPU}  STRIDE=${STRIDE}  [$(date +%F_%H:%M:%S)] ###"
    local e
    for e in "${ENS_LIST[@]}"; do run_half "$e" "$WID"; done
    echo "### worker ${WID} ALL DONE  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

echo "### START L2 Nf2 finestride: GPU${GPU} x ${NPACK} MPS, STRIDE=${STRIDE}, 3 Nf2 ensembles, halved per worker  [$(date +%F_%H:%M:%S)] ###"
for (( w=0; w<NPACK; w++ )); do run_worker "$w" & done
wait
echo "### L2 Nf2 finestride ALL DONE  [$(date +%F_%H:%M:%S)] ###"
