#!/usr/bin/env bash
# run_distill_L1_nsrc2_claude.sh -- HIGH-STATISTICS L1 perambulators for the sigma^2 mixing / two-meson study.
# MRHS driver, STRIDE=4 (~1000 configs from 3999 ckpoint) x --nsrc 2 (two source windows {0,64}) = ~5x the
# current 400 stride-10 nsrc-1 set.  Cheap at L1 (nv=24).  Writes data_<ens>/distill_Nv24_v2/peram.<k>.h5 with
# /peram/{tau,tau_gw} shape (nsrc=2, twin, twin, Nv, Nv) + meta/tsrc_list (SEPARATE _v2 dir -- does NOT touch the
# stride-10 distill_Nv24/).  Analysis: distill_contract_claude.load_peram_windows -> 2 windows as independent
# source sets.  Goal: enough statistics to lift C_FF's light-0++ tail above noise -> resolve the {F^2,sigma^2}
# mixing GEVP, and sharpen m1 - 2 m_PS (the sub-threshold two-meson).
#
# Ensemble = Nf2 gsq1.0 at0.2 L1 (the primary study ensemble).  GPU1, MPS NPACK=2, sub-stream tiled.  RESUMABLE,
# skip-protected.  NO rm.  NO kill.  Launch ONCE (each setsid duplicates the worker!).
# Launch DETACHED: NPACK=2 setsid nohup bash run_distill_L1_nsrc2_claude.sh > run_distill_L1_nsrc2_claude.log 2>&1 < /dev/null &
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1

# ---- knobs ----
GPU="${GPU:-1}"
NPACK="${NPACK:-2}"
STRIDE="${STRIDE:-4}"
NV=24
TSRC0=0
TWIN=32
NSRC=2
L=1
ENS="${ENS:-Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000}"

NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
export LD_LIBRARY_PATH=/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib:${LD_LIBRARY_PATH:-}
SRC=distill_peram_mrhs_claude.cu
BIN="distill_peram_mrhs_L${L}_claude.o"

if [ ! -f "$BIN" ] || find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$BIN" -print -quit 2>/dev/null | grep -q .
then
  echo "### compile mrhs L${L} -> $BIN  [$(date +%F_%H:%M:%S)] ###"
  "$NVCC" "$SRC" $NVCCFLAGS $INCLUDES $LDFLAGS -o "$BIN" || { echo "### BUILD FAILED ###"; exit 1; }
else
  echo "### mrhs L${L} binary up-to-date ###"
fi

[ -d "$ENS" ] && ls "$ENS"/ckpoint_lat.* >/dev/null 2>&1 || { echo "### ERROR: no ckpoint_lat in $ENS ###"; exit 1; }
nf=$(printf '%s' "$ENS" | grep -oE 'Nf[0-9]+' | head -1 | sed 's/Nf//')
gsq=$(printf '%s' "$ENS" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//')
at=$(printf '%s' "$ENS" | grep -oE 'at[0-9.]+nu0' | sed 's/at//;s/nu0//')
first=$(ls "$ENS"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n | head -1)
last=$( ls "$ENS"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n | tail -1)
kmax=$(( last + 1 ))
echo "### ens=$ENS  Nf${nf} g${gsq} at${at}  stride=${STRIDE} nsrc=${NSRC}  ~$(( (last-first)/STRIDE + 1 )) configs  [$(date +%F_%H:%M:%S)] ###"

if [ "$NPACK" -gt 1 ] && ! pgrep -f nvidia-cuda-mps-control >/dev/null
then
  echo "### starting MPS daemon (NPACK=$NPACK) ###"
  nvidia-cuda-mps-control -d
  for i in 1 2 3 4 5; do pgrep -f nvidia-cuda-mps-control >/dev/null && break; sleep 1; done
fi

run_worker () {   # $1=WID ; sub-stream: kmin=first+STRIDE*WID, stride_eff=STRIDE*NPACK
  local WID="$1"
  local stride_eff=$(( STRIDE * NPACK ))
  local kmin_w=$(( first + STRIDE * WID ))
  local LOG="distill_L${L}_nsrc2_w${WID}_claude.log"
  {
    echo "### START L${L} nsrc2 worker ${WID}/${NPACK}  GPU${GPU}  Nv=${NV} nsrc=${NSRC} kmin=${kmin_w} stride=${stride_eff} kmax=${kmax}  [$(date +%F_%H:%M:%S)] ###"
    CUDA_VISIBLE_DEVICES=$GPU ./"$BIN" --gsq "$gsq" --Nf "$nf" --nu0 1.0 --at "$at" \
      --nv "$NV" --tsrc0 "$TSRC0" --twin "$TWIN" --nsrc "$NSRC" \
      --ens-dir "$ENS/" --kmin "$kmin_w" --stride "$stride_eff" --kmax "$kmax"
    echo "### worker ${WID} done (status $?)  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

echo "### START L${L} nsrc2: GPU${GPU} x ${NPACK} MPS, stride=${STRIDE} Nv=${NV} --nsrc ${NSRC} -> distill_Nv24_v2/  [$(date +%F_%H:%M:%S)] ###"
for (( w=0; w<NPACK; w++ )); do run_worker "$w" & done
wait
echo "### L${L} nsrc2 DONE  [$(date +%F_%H:%M:%S)] ###"
