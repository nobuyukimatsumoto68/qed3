#!/bin/bash
# SUPERSEDED 2026-10-02 (NM): superjackknife handles mismatched conn/disc config sets, so the odd-k
# restriction is not needed. Kept for reference only; not launched.
# run_disc_scalar_oddk_claude.sh -- scalar FS disconnected loop, ODD-k config sets (all paired with conn).
# COPY of run_disc_scalar_L1L2_claude.sh (kept unchanged as the reference).  Change vs the original:
#   configs = two stride-10 passes per ensemble, k = K0 + 10 j with K0 = 1 and K0 = 5, i.e.
#   k = 1, 5, 11, 15, 21, 25, ...  -> 800 configs per ensemble at BOTH L1 and L2, all at ODD k, so every
#   disc config has a connected config (conn exists at odd k only).  The original used L2 stride 5
#   (k = 1, 6, 11, 16, ...): only its odd half (400) pairs with conn.
#   Nf loop: NFS (default "2 4") runs Nf2 first, then Nf4.
#   Already-finished hits are skipped (per-hit "complete" gating), so for Nf2 the K0=1 pass is a fast
#   skip and only the K0=5 pass computes.
# Plan: final/analysis_scalar/scalar_fs_disc_impl_plan_claude.md.
#
# Scope (NM 2026-10-02): massless, at0.2, L1 + L2, nhits=2 (one mrhs block, NB_CLI=2), Nf2 then Nf4.
# PACKING: 4 workers = 2 per GPU under MPS; each worker takes an interleaved sub-stream of the current
# (ensemble, K0) pass (worker W starts at K0 + 10 W, effective stride 40).
# RESUME-SAFE; no rm, no kill.
#
# RUN DETACHED (long GPU job -- you launch it, not Claude):
#     nohup bash run_disc_scalar_oddk_claude.sh > /dev/null 2>&1 &
#     NFS="6" nohup bash run_disc_scalar_oddk_claude.sh > /dev/null 2>&1 &          # later
# Master log: run_disc_scalar_oddk_claude.log ; per worker: disc_scalar_L<L>_g<gsq>_Nf<nf>_k0<K0>_w<W>_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

NFS="${NFS:-2 4}"
LS="${LS:-1 2}"
read -ra GPUS <<< "${GPUS:-0 0 1 1}"     # 4 workers, 2 per GPU (MPS packed)
NHITS=2
DISC_TB=2
NB=2                                     # hits per block solve (must equal NHITS here)
# OLD (L1 stride 10, L2 stride 5 from k=1):
# stride_of_L () { case "$1" in 1) echo 10 ;; 2) echo 5 ;; *) echo 10 ;; esac }
STRIDE=10
K0S="${K0S:-1 5}"                        # two stride-10 passes: k = 1,11,21,... and 5,15,25,...

MASTER="run_disc_scalar_oddk_claude.log""
exec > >(tee -a "$MASTER") 2>&1

# ---- toolchain (CUDA 12.9: last release with sm_70; host compiler g++-12, system gcc 15 is rejected) ----
NVCC="${NVCC:-/usr/local/cuda-12.9/bin/nvcc}"
export LD_LIBRARY_PATH="/usr/local/cuda-12.9/lib64:${LD_LIBRARY_PATH:-}"
CCBIN="${CCBIN:-g++-12}"
NVCCBASE="-ccbin $CCBIN -arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
SRC=jj_local_ylm_scalar_disc_stoch_claude.cu

echo "############################################################"
echo "### scalar FS disc ODD-k  NFS=$NFS  K0S=$K0S  LS=$LS  GPUS=${GPUS[*]}  nhits=$NHITS tb=$DISC_TB NB=$NB  [$(date +%F_%H:%M:%S)] ###"
echo "############################################################"

# ---- binaries (rebuild when .cu/.h newer; FORCE_BUILD=1 to force) ----
need_build () {
  [ -n "${FORCE_BUILD:-}" ] && return 0
  [ ! -f "$1" ] && return 0
  find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$1" -print -quit 2>/dev/null | grep -q . && return 0
  return 1
}
for L in $LS
do
  BIN="jj_local_ylm_scalar_disc_stoch_L${L}_NB${NB}.o"
  if need_build "$BIN"
  then
    echo "### compile L${L} NB${NB} -> $BIN  [$(date +%F_%H:%M:%S)] ###"
    $NVCC $NVCCBASE -DN_REFINE_CLI=${L} -DNB_CLI=${NB} $INCLUDES $LDFLAGS "$SRC" -o "$BIN"
    if [ $? -ne 0 ]
    then
      echo "### BUILD FAILED: L${L} NB${NB} -- aborting ###"
      exit 1
    fi
  else
    echo "### L${L} NB${NB} binary up-to-date, skip ###"
  fi
done
echo "### build OK ###"

get_L ()   { printf '%s' "$1" | grep -oE 'nt128L[0-9]+' | sed 's/nt128L//'; }
get_nf ()  { printf '%s' "$1" | grep -oE '^Nf[0-9]+' | sed 's/Nf//'; }
get_gsq () { printf '%s' "$1" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//'; }

# ---- collect (ensemble, K0) passes: Nf in $NFS order, then L ascending, gsq ascending, then K0 ----
ENS=()
PASSK0=()
for NF in $NFS
do
  for L in $LS
  do
    for d in $(ls -d Nf${NF}_gsq*at0.200000*mRe0.000000mIm0.000000nt128L${L}_hb* 2>/dev/null | sort -t q -k2,2)
    do
      [ -d "$d" ] || continue
      ls "$d"/ckpoint_lat.* >/dev/null 2>&1 || continue
      for K0 in $K0S
      do
        ENS+=("$d")
        PASSK0+=("$K0")
      done
    done
  done
done
if [ "${#ENS[@]}" -eq 0 ]
then
  echo "### no ensembles found for NFS=$NFS LS=$LS -- nothing to do ###"
  exit 0
fi
echo "### ${#ENS[@]} (ensemble, K0) passes, ${#GPUS[@]} workers (GPUs ${GPUS[*]}), MPS packed ###"
for (( i=0; i<${#ENS[@]}; i++ ))
do
  echo "###   ${ENS[$i]}  K0=${PASSK0[$i]} stride $STRIDE"
done

run_disc () {   # $1=ens-dir  $2=GPU  $3=WID  $4=NWORKERS  $5=K0
  local ens="$1"
  local GPU="$2"
  local WID="$3"
  local NW="$4"
  local K0="$5"
  local L nf gsq bin LOG ks first last kmin_w stride stride_eff kmax nconf
  L=$(get_L "$ens")
  nf=$(get_nf "$ens")
  gsq=$(get_gsq "$ens")
  stride=$STRIDE
  bin="jj_local_ylm_scalar_disc_stoch_L${L}_NB${NB}.o"
  LOG="disc_scalar_L${L}_g${gsq}_Nf${nf}_k0${K0}_w${WID}_claude.log"
  mapfile -t ks < <(ls "$ens"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n)
  if [ "${#ks[@]}" -eq 0 ]
  then
    echo "### SKIP $ens: no ckpoint_lat  [$(date +%F_%H:%M:%S)] ###" >> "$LOG"
    return 0
  fi
  first="$K0"
  last="${ks[-1]}"
  # OLD: kmin_w=$(( ${ks[0]} + stride * WID ))
  kmin_w=$(( first + stride * WID ))
  stride_eff=$(( stride * NW ))
  kmax=$(( last + 1 ))
  if [ "$kmin_w" -ge "$kmax" ]
  then
    echo "### $ens w$WID: no configs in sub-stream  [$(date +%F_%H:%M:%S)] ###" >> "$LOG"
    return 0
  fi
  nconf=$(( (last - kmin_w) / stride_eff + 1 ))
  {
    echo "### DISC-scalar $ens  L${L}  K0=${K0}  GPU${GPU}  w${WID}  kmin=$kmin_w stride=$stride_eff kmax=$kmax tb=$DISC_TB nhits=$NHITS NB=$NB (~$nconf cfg)  [$(date +%F_%H:%M:%S)] ###"
    CUDA_VISIBLE_DEVICES=$GPU ./"$bin" --gsq "$gsq" --Nf "$nf" --nu0 1.0 \
      --ens-dir "$ens/" --kmin "$kmin_w" --stride "$stride_eff" --kmax "$kmax" \
      --nhits "$NHITS" --disc-tblock "$DISC_TB"
    echo "### DISC-scalar $ens w${WID} done (status $?)  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

worker () {   # $1=WID
  local WID="$1"
  local GPU="${GPUS[$WID]}"
  local NW=${#GPUS[@]}
  local i
  for (( i=0; i<${#ENS[@]}; i++ ))
  do
    echo "### worker $WID (GPU $GPU) -> ${ENS[$i]} K0=${PASSK0[$i]}  [$(date +%F_%H:%M:%S)] ###"
    run_disc "${ENS[$i]}" "$GPU" "$WID" "$NW" "${PASSK0[$i]}"
  done
  echo "### worker $WID (GPU $GPU) done  [$(date +%F_%H:%M:%S)] ###"
}

# ---- MPS up (packed phase) ----
if pgrep -f nvidia-cuda-mps-control >/dev/null
then
  echo "### MPS daemon: already running ###"
else
  echo "### starting nvidia-cuda-mps-control -d ###"
  nvidia-cuda-mps-control -d
  for i in 1 2 3 4 5
  do
    pgrep -f nvidia-cuda-mps-control >/dev/null && break
    sleep 1
  done
fi
pgrep -f nvidia-cuda-mps-control >/dev/null \
  || { echo "### ERROR: MPS daemon failed to start -- aborting ###"; exit 1; }

echo "### START scalar FS disc ODD-k NFS=$NFS: ${#GPUS[@]} workers (GPUs ${GPUS[*]}), MPS  [$(date +%F_%H:%M:%S)] ###"
for (( W=0; W<${#GPUS[@]}; W++ ))
do
  worker "$W" &
done
wait
echo "### ALL WORKERS DONE ODD-k NFS=$NFS  [$(date +%F_%H:%M:%S)] ###"
