#!/usr/bin/env bash
# run_distill_L1_sweep_claude.sh -- PRODUCTION perambulator sweep for exact distillation of the sigma\sigma
# four-point / coupled GEVP (Chester-Pufu 0++ mixing).  Driver: distill_peram_claude.cu.  Writes, per config,
# data_<ens>/distill_Nv<NV>/peram.<k>.h5 = /V, /evals, /peram/{tau,tau_gw} (windowed W x W x Nv x Nv).
#
# Basis = timeslice Wilson D_{W,2}^dag D_{W,2} low modes (Nv=2N_s=24 at L1 = EXACT local).  Forward peram tau
# + furnished tau'.  Source window [tsrc0, tsrc0+twin), twin=32 brackets the L1 scalar plateau dt[8,24].
# Plans: distillation_impl_plan_claude.md + distill_fourpoint_diagrams_claude.md.
#
# GPU1 only, MPS NPACK=2 -> 2 workers.  Sub-stream tiling: worker w does kmin=first+STRIDE*w,
# stride_eff=STRIDE*NPACK, so the 2 workers together cover the full stride-10 grid.  RESUMABLE: the driver
# skips per-config h5 that already have /peram/tau.  NO rm.  NO kill.
#
# Launch DETACHED:  setsid nohup bash run_distill_L1_sweep_claude.sh > /dev/null 2>&1 < /dev/null &
# Per-worker output -> distill_L1_sweep_w<WID>_claude.log .  (~285 s/config; ~360 configs / 2 packs.)
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1

# ---- knobs ----
GPU=1
NPACK=2
STRIDE=10
NV=24                  # full 2*N_SITES at L1 (exact local)
TSRC0=0
TWIN=32
AT_TOKEN="${AT_TOKEN:-at0.200000}"   # override AT_TOKEN=at0.100000 for the half-a_t set
L=1
ENS_GLOB="${ENS_GLOB:-Nf*_gsq*${AT_TOKEN}*mRe0.000000mIm0.000000nt128L${L}_hb*}"

NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=distill_peram_claude.cu
BIN="distill_peram_L${L}_claude.o"

# ---- build (skip if up-to-date; FORCE_BUILD=1 to force) ----
need_build () {
  [ -n "${FORCE_BUILD:-}" ] && return 0
  [ ! -f "$BIN" ] && return 0
  find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$BIN" -print -quit 2>/dev/null | grep -q . && return 0
  return 1
}
if need_build
then
  echo "### compile distill L${L} -> $BIN  [$(date +%F_%H:%M:%S)] ###"
  "$NVCC" "$SRC" $NVCCFLAGS $INCLUDES $LDFLAGS -o "$BIN" \
    || { echo "### distill L${L} BUILD FAILED ###"; exit 1; }
else
  echo "### distill L${L} binary up-to-date, skip (FORCE_BUILD=1 to rebuild) ###"
fi

# ---- ensure MPS daemon up (2-pack needs it) ----
if pgrep -f nvidia-cuda-mps-control >/dev/null
then
  echo "### MPS daemon: already running ###"
else
  echo "### MPS daemon not up -- starting nvidia-cuda-mps-control -d ###"
  nvidia-cuda-mps-control -d
  for i in 1 2 3 4 5
  do
    pgrep -f nvidia-cuda-mps-control >/dev/null && break
    sleep 1
  done
fi
pgrep -f nvidia-cuda-mps-control >/dev/null \
  || { echo "### ERROR: MPS daemon failed to start -- aborting ###"; exit 1; }

# ---- collect ensembles ----
ENS=()
for d in $ENS_GLOB
do
  [ -d "$d" ] || continue
  ls "$d"/ckpoint_lat.* >/dev/null 2>&1 || continue
  ENS+=( "$d" )
done
echo "### ${#ENS[@]} L${L} massless ensembles matched (${AT_TOKEN}) ###"
printf '  %s\n' "${ENS[@]}"

# ---- one worker: sub-stream WID of NPACK over ALL matched ensembles ----
run_worker () {   # $1=WID
  local WID="$1"
  local stride_eff=$(( STRIDE * NPACK ))
  local LOG="distill_L${L}_sweep_w${WID}_claude.log"
  local ens nf gsq at first last kmin_w kmax
  {
    echo "### START distill sweep L${L} worker ${WID}/${NPACK}  GPU${GPU}  stride_eff=${stride_eff}  Nv=${NV} twin=${TWIN}  [$(date +%F_%H:%M:%S)] ###"
    for ens in "${ENS[@]}"
    do
      nf=$(printf '%s' "$ens" | grep -oE 'Nf[0-9]+' | head -1 | sed 's/Nf//')
      gsq=$(printf '%s' "$ens" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//')
      at=$(printf '%s' "$ens" | grep -oE 'at[0-9.]+nu0' | sed 's/at//;s/nu0//')
      first=$(ls "$ens"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n | head -1)
      last=$( ls "$ens"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n | tail -1)
      [ -z "$first" ] && { echo "### SKIP $ens (no ckpoint_lat) ###"; continue; }
      kmin_w=$(( first + STRIDE * WID ))
      kmax=$(( last + 1 ))
      [ "$kmin_w" -ge "$kmax" ] && { echo "### SKIP $ens w${WID} (kmin>=kmax) ###"; continue; }
      echo "### $ens  Nf${nf} g${gsq} at${at}  kmin=${kmin_w} stride=${stride_eff} kmax=${kmax}  [$(date +%F_%H:%M:%S)] ###"
      CUDA_VISIBLE_DEVICES=$GPU ./"$BIN" --gsq "$gsq" --Nf "$nf" --nu0 1.0 --at "$at" \
        --nv "$NV" --tsrc0 "$TSRC0" --twin "$TWIN" \
        --ens-dir "$ens/" --kmin "$kmin_w" --stride "$stride_eff" --kmax "$kmax"
      echo "### $ens w${WID} done (status $?)  [$(date +%F_%H:%M:%S)] ###"
    done
    echo "### distill sweep L${L} worker ${WID} ALL DONE  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

# ---- launch NPACK workers on GPU1 (MPS-packed), wait for all ----
echo "### START distill sweep L${L}: GPU${GPU} x MPS ${NPACK}-pack, STRIDE=${STRIDE} Nv=${NV} twin=${TWIN} tsrc0=${TSRC0}  [$(date +%F_%H:%M:%S)] ###"
for (( w=0; w<NPACK; w++ ))
do
  run_worker "$w" &
done
wait
echo "### ALL WORKERS DONE  [$(date +%F_%H:%M:%S)] ###"
