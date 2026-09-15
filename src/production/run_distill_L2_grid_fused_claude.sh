#!/usr/bin/env bash
# run_distill_L2_grid_fused_claude.sh -- L2 perambulators over the Nf{2,4} x gsq{1,2,3} at0.2 grid.
# Nv=24 (TRUNCATED L2 basis -> distillation smearing, NOT exact; matches the L1 basis SIZE), NSTACK=24 for the
# real solve speedup, --nsrc 2 (windows {0,64}) + --fused-sources (nsrc=2 stats at nsrc=1 solve cost).
# STRIDE=10 -> ~400 configs/ensemble (k=1,11,...,3991), like the L1 400-config sets.  Writes each to
# data_<ens>/distill_Nv24_v2/peram.<k>.h5.  RESUMABLE, skip-protected.  NO rm.  NO kill.
#
# SCHEDULE: ONE worker owns a WHOLE ensemble (no config-tiling); the 2 MPS workers run TWO ensembles
# concurrently on GPU1.  Ascending order 1->6 via static interleave: worker0 = {1,3,5}, worker1 = {2,4,6},
# so the concurrent pairs are (1,2) -> (3,4) -> (5,6).  Each worker does its ensembles sequentially.
# Launch ONCE.  DETACHED:
#   GPU=1 setsid nohup bash run_distill_L2_grid_fused_claude.sh > run_distill_L2_grid_fused_claude.log 2>&1 < /dev/null &
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1

# ---- ensembles, ASCENDING order 1..6 ----
ENS_LIST=(
  "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000"   # 1
  "Nf2_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000"   # 2
  "Nf2_gsq3.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000"   # 3
  "Nf4_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000"   # 4
  "Nf4_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000"   # 5
  "Nf4_gsq3.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000"   # 6
)

# ---- knobs ----
GPU="${GPU:-1}"
NPACK=2                         # two MPS workers = two ensembles concurrent
STRIDE="${STRIDE:-10}"          # ~400 configs/ensemble
NV=24                           # truncated basis (smeared L2); NSTACK also 24 (build flag) for the solve speedup
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

run_one_ensemble () {   # $1 = ens config-dir name
  local ENS="$1"
  if [ ! -d "$ENS" ] || ! ls "$ENS"/ckpoint_lat.* >/dev/null 2>&1; then
    echo "### SKIP (no ckpoint_lat): $ENS ###"
    return
  fi
  local nf gsq at first last kmax
  nf=$(printf '%s' "$ENS" | grep -oE 'Nf[0-9]+' | head -1 | sed 's/Nf//')
  gsq=$(printf '%s' "$ENS" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//')
  at=$(printf '%s' "$ENS" | grep -oE 'at[0-9.]+nu0' | sed 's/at//;s/nu0//')
  first=$(ls "$ENS"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n | head -1)
  last=$( ls "$ENS"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n | tail -1)
  kmax=$(( last + 1 ))
  echo "### ENS $ENS  Nf${nf} g${gsq} at${at}  stride=${STRIDE} nsrc=${NSRC} Nv=${NV} FUSED  ~$(( (last-first)/STRIDE + 1 )) cfg  [$(date +%F_%H:%M:%S)] ###"
  CUDA_VISIBLE_DEVICES=$GPU ./"$BIN" --gsq "$gsq" --Nf "$nf" --nu0 1.0 --at "$at" \
    --nv "$NV" --tsrc0 "$TSRC0" --twin "$TWIN" --nsrc "$NSRC" --fused-sources \
    --ens-dir "$ENS/" --kmin "$first" --stride "$STRIDE" --kmax "$kmax"
  echo "### ENS $ENS DONE (status $?)  [$(date +%F_%H:%M:%S)] ###"
}

run_worker () {   # $1 = WID ; owns ensembles i with i%NPACK==WID (0-based), in ascending order
  local WID="$1"
  local LOG="distill_L${L}_grid_w${WID}_claude.log"
  {
    echo "### START L${L} GRID worker ${WID}/${NPACK}  GPU${GPU}  Nv=${NV}  [$(date +%F_%H:%M:%S)] ###"
    local i
    for (( i=WID; i<${#ENS_LIST[@]}; i+=NPACK )); do
      run_one_ensemble "${ENS_LIST[$i]}"
    done
    echo "### worker ${WID} ALL DONE  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

echo "### START L${L} GRID: GPU${GPU} x ${NPACK} MPS  worker0={1,3,5} worker1={2,4,6}  stride=${STRIDE} Nv=${NV} --nsrc ${NSRC} --fused-sources  [$(date +%F_%H:%M:%S)] ###"
for (( w=0; w<NPACK; w++ )); do run_worker "$w" & done
wait
echo "### L${L} GRID ALL DONE  [$(date +%F_%H:%M:%S)] ###"
