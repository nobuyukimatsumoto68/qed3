#!/usr/bin/env bash
# run_distill_L2_fused_claude.sh -- L2 perambulators for the sigma^2 / two-meson bounding-energy study.
# MRHS driver at L2 (N_REFINE=2, N_SITES=42, Nv=84 = complete local basis), --nsrc 2 (windows {0,64}) with
# --fused-sources: ONE combined-source block solve per offset feeds BOTH windows, so the solve count is twin
# (not nsrc*twin) -- nsrc=2 statistics at nsrc=1 solve cost (the L2 wall-time win; largest here since L2 solves
# dominate).  Writes data_<ens>/distill_Nv24_v2/peram.<k>.h5, shape (nsrc=2,twin,twin,24,24)+meta/tsrc_list
# (Nv=24 = TRUNCATED L2 basis -> distillation smearing, NOT exact; NSTACK=24 for the real solve speedup).
# Fused correctness+contamination validated at L1 (analysis-range dev ~2.5e-5 << gauge noise); the binary's
# per-config [fused-check] prints the L2 contamination on the first config -- WATCH it (larger mass gap or
# smaller L2 gap changes the tail).  SEPARATE _v2 dir; RESUMABLE, skip-protected.  NO rm.  NO kill.
#
# Ensemble = Nf2 gsq1.0 at0.2 L2 (larger-volume partner of the L1 primary).  GPU1, MPS NPACK=2, sub-stream tiled.
# Launch ONCE (each setsid duplicates the worker!).
# Launch DETACHED:
#   NPACK=2 setsid nohup bash run_distill_L2_fused_claude.sh > run_distill_L2_fused_claude.log 2>&1 < /dev/null &
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1

# ---- knobs ----
GPU="${GPU:-1}"
NPACK="${NPACK:-2}"
STRIDE="${STRIDE:-20}"          # ~200 configs from 3999 ckpoint (L2 is costly); lower for more, env-overridable
NV=24                           # TRUNCATED to 24 modes to match the L1 basis SIZE (NOT the complete L2 basis of
                                # 84) -> distillation SMEARING at L2, no longer exact.  Block width NSTACK also 24
                                # (-DNSTACK_CLI below) so the block matvec is 24-wide -> the real ~3.5x solve win.
TSRC0=0
TWIN=32
NSRC=2
L=2
ENS="${ENS:-Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000}"

NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O3 -std=c++20 -DN_REFINE_CLI=2 -DNSTACK_CLI=${NV} -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
export LD_LIBRARY_PATH=/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib:${LD_LIBRARY_PATH:-}
SRC=distill_peram_mrhs_claude.cu
BIN="distill_peram_mrhs_L${L}_Nv${NV}_claude.o"  # L2 Nv24 binary (N_REFINE=2, NSTACK=24); distinct name so it
                                                 # never reuses a full-NSTACK L2 build (the build guard is mtime-only)

if [ ! -f "$BIN" ] || find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$BIN" -print -quit 2>/dev/null | grep -q .
then
  echo "### compile mrhs L${L} (N_REFINE=2) -> $BIN  [$(date +%F_%H:%M:%S)] ###"
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
echo "### ens=$ENS  Nf${nf} g${gsq} at${at}  stride=${STRIDE} nsrc=${NSRC} FUSED  ~$(( (last-first)/STRIDE + 1 )) configs  [$(date +%F_%H:%M:%S)] ###"

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
  local LOG="distill_L${L}_fused_w${WID}_claude.log"
  {
    echo "### START L${L} FUSED worker ${WID}/${NPACK}  GPU${GPU}  Nv=${NV} nsrc=${NSRC} kmin=${kmin_w} stride=${stride_eff} kmax=${kmax}  [$(date +%F_%H:%M:%S)] ###"
    CUDA_VISIBLE_DEVICES=$GPU ./"$BIN" --gsq "$gsq" --Nf "$nf" --nu0 1.0 --at "$at" \
      --nv "$NV" --tsrc0 "$TSRC0" --twin "$TWIN" --nsrc "$NSRC" --fused-sources \
      --ens-dir "$ENS/" --kmin "$kmin_w" --stride "$stride_eff" --kmax "$kmax"
    echo "### worker ${WID} done (status $?)  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

echo "### START L${L} FUSED: GPU${GPU} x ${NPACK} MPS, stride=${STRIDE} Nv=${NV} NSTACK=${NV} --nsrc ${NSRC} --fused-sources -> distill_Nv${NV}_v2/  [$(date +%F_%H:%M:%S)] ###"
for (( w=0; w<NPACK; w++ )); do run_worker "$w" & done
wait
echo "### L${L} FUSED DONE  [$(date +%F_%H:%M:%S)] ###"
