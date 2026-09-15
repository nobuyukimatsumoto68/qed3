#!/usr/bin/env bash
# run_sigma_L1_claude.sh -- build + launch the MASSLESS scalar-loop driver (D_S + D'_S) on the L=1
# massless ensembles, for the sigma^2 - F^2 (0++) mixing check (arXiv:1603.05582).  Driver:
# jj_sigma_loops_stoch_claude.cu.  Measures h0/sigma_loops/{DS,DS_1mD,Dp,Dp_1mD} per (config,hit).
#
# GPU1 only, MPS NPACK=2 -> 2 workers.  Sub-stream tiling: worker w does configs kmin=first+STRIDE*w,
# stride_eff=STRIDE*NPACK, so the 2 workers together cover the full stride-10 grid.
# RESUMABLE: the driver skips per-(config,hit) files that already have the DS group.  NO rm.  NO kill.
#
# Launch DETACHED, e.g.:   nohup bash run_sigma_L1_claude.sh > /dev/null 2>&1 &
# All output is tee'd into sigma_L1_*_claude.log (one per worker).
set -u

# ---- knobs ----
GPU=1
NPACK=2
STRIDE=10
NHITS=2
DISC_TBLOCK=2          # interval = Nt/DISC_TBLOCK = Nt/2 time-dilution classes
SEED_TAG=sigmaloops    # INDEPENDENT of disc-loop / FNAL-conn / D'_S streams
AT_TOKEN="${AT_TOKEN:-at0.200000}"   # override AT_TOKEN=at0.100000 for the half-a_t set
L=1

# massless L1 ensembles (mRe0 mIm0); override ENS_GLOB to narrow (e.g. only Nf2)
ENS_GLOB="${ENS_GLOB:-Nf*_gsq*${AT_TOKEN}*mRe0.000000mIm0.000000nt128L${L}_hb*}"

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
SRC=jj_sigma_loops_stoch_claude.cu
BIN="jj_sigma_loops_L${L}.o"

# ---- build (skip if up-to-date; FORCE_BUILD=1 to force) ----
need_build () {
  [ -n "${FORCE_BUILD:-}" ] && return 0
  [ ! -f "$BIN" ] && return 0
  find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$BIN" -print -quit 2>/dev/null | grep -q . && return 0
  return 1
}
if need_build
then
  echo "### compile sigma-loops L${L} (-DN_REFINE_CLI=${L}) -> $BIN  [$(date +%F_%H:%M:%S)] ###"
  $NVCC $NVCCBASE -DN_REFINE_CLI=${L} $INCLUDES $LDFLAGS "$SRC" -o "$BIN" \
    || { echo "### sigma-loops L${L} BUILD FAILED ###"; exit 1; }
else
  echo "### sigma-loops L${L} binary up-to-date, skip (FORCE_BUILD=1 to rebuild) ###"
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
  local LOG="sigma_L${L}_w${WID}_claude.log"
  local ens nf gsq at first last kmin_w kmax
  {
    echo "### START sigma-loops L${L} worker ${WID}/${NPACK}  GPU${GPU}  stride_eff=${stride_eff}  [$(date +%F_%H:%M:%S)] ###"
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
        --ens-dir "$ens/" --kmin "$kmin_w" --stride "$stride_eff" --kmax "$kmax" \
        --nhits "$NHITS" --disc-tblock "$DISC_TBLOCK" --seed-tag "$SEED_TAG"
      echo "### $ens w${WID} done (status $?)  [$(date +%F_%H:%M:%S)] ###"
    done
    echo "### sigma-loops L${L} worker ${WID} ALL DONE  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

# ---- launch NPACK workers on GPU1 (MPS-packed), wait for all ----
echo "### START sigma-loops L${L}: GPU${GPU} x MPS ${NPACK}-pack, STRIDE=${STRIDE} nhits=${NHITS} tblock=${DISC_TBLOCK}  [$(date +%F_%H:%M:%S)] ###"
for (( w=0; w<NPACK; w++ ))
do
  run_worker "$w" &
done
wait
echo "### ALL WORKERS DONE  [$(date +%F_%H:%M:%S)] ###"
