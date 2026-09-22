#!/usr/bin/env bash
# run_distill_L2_mrhs_claude.sh -- L2 perambulator generation (NEW) for the finite-volume test of the
# sub-threshold two-meson (bound state vs Luscher shift): the L2 partner of the L1 Nf2 gsq1.0 bound-state study.
# MRHS driver, compiled -DN_REFINE_CLI=2 (L2: N_SITES=42, NSTACK=Nx=84 = full exact basis).  --nsrc 2 (TWO
# source windows {0, Nt/2=64}) -> ~2x independent statistics per config.  Writes data_<ens>/distill_Nv84_v2/
# peram.<k>.h5 with /peram/{tau,tau_gw} shape (nsrc=2, twin, twin, Nv, Nv) + meta/tsrc_list (SEPARATE _v2 dir,
# never touches production).  DOWNSTREAM: distill_contract_claude.py must read the leading nsrc axis (chunk 3).
#
# Ensemble = Nf2 gsq1.0 at0.2 L2 (3999 configs).  GPU1, single worker by default (nv=84 mrhs block is heavy;
# set NPACK=2 to try MPS 2-pack if memory allows).  RESUMABLE, skip-protected.  NO rm.  NO kill.
# Launch DETACHED:  setsid nohup bash run_distill_L2_mrhs_claude.sh > /dev/null 2>&1 < /dev/null &
# Log -> distill_L2_mrhs_w<WID>_claude.log .  RUN after the _v1 sweep is stopped (frees GPU1).
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1

# ---- knobs ----
GPU="${GPU:-1}"
NPACK="${NPACK:-1}"          # nv=84 block is memory-heavy; start single-worker, try NPACK=2 if it fits
STRIDE="${STRIDE:-10}"
NV=84                        # full 2*N_SITES at L2 (exact)
TSRC0=0
TWIN=32
NSRC=2                       # two source windows {0, 64}
L=2
ENS="${ENS:-Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000}"

NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O3 -std=c++20 -DN_REFINE_CLI=2 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
export LD_LIBRARY_PATH=/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib:${LD_LIBRARY_PATH:-}
SRC=distill_peram_mrhs_claude.cu
BIN="distill_peram_mrhs_L${L}_claude.o"

# ---- build (L2 = separate binary; ALWAYS rebuild if source newer) ----
if [ ! -f "$BIN" ] || find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$BIN" -print -quit 2>/dev/null | grep -q .
then
  echo "### compile mrhs L${L} (-DN_REFINE_CLI=2, NSTACK=84) -> $BIN  [$(date +%F_%H:%M:%S)] ###"
  "$NVCC" "$SRC" $NVCCFLAGS $INCLUDES $LDFLAGS -o "$BIN" \
    || { echo "### mrhs L${L} BUILD FAILED ###"; exit 1; }
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

# ---- MPS daemon only if 2-packing ----
if [ "$NPACK" -gt 1 ] && ! pgrep -f nvidia-cuda-mps-control >/dev/null
then
  echo "### starting MPS daemon (NPACK=$NPACK) ###"
  nvidia-cuda-mps-control -d
  for i in 1 2 3 4 5; do pgrep -f nvidia-cuda-mps-control >/dev/null && break; sleep 1; done
fi

run_worker () {   # $1=WID
  local WID="$1"
  local stride_eff=$(( STRIDE * NPACK ))
  local kmin_w=$(( first + STRIDE * WID ))
  local LOG="distill_L${L}_mrhs_w${WID}_claude.log"
  {
    echo "### START mrhs L${L} worker ${WID}/${NPACK}  GPU${GPU}  Nv=${NV} nsrc=${NSRC}  ens=$ENS  kmin=${kmin_w} stride=${stride_eff} kmax=${kmax}  [$(date +%F_%H:%M:%S)] ###"
    CUDA_VISIBLE_DEVICES=$GPU ./"$BIN" --gsq "$gsq" --Nf "$nf" --nu0 1.0 --at "$at" \
      --nv "$NV" --tsrc0 "$TSRC0" --twin "$TWIN" --nsrc "$NSRC" \
      --ens-dir "$ENS/" --kmin "$kmin_w" --stride "$stride_eff" --kmax "$kmax"
    echo "### mrhs L${L} worker ${WID} done (status $?)  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

echo "### START mrhs L${L}: GPU${GPU} x ${NPACK}, Nv=${NV} --nsrc ${NSRC} -> distill_Nv84_v2/  ens=$ENS  [$(date +%F_%H:%M:%S)] ###"
for (( w=0; w<NPACK; w++ )); do run_worker "$w" & done
wait
echo "### L${L} DONE  [$(date +%F_%H:%M:%S)] ###"
