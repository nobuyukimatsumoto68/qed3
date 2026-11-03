#!/bin/bash
# run_disc_scalar_L1L2_claude.sh -- PRODUCTION run of the scalar FS disconnected loop
# (jj_local_ylm_scalar_disc_stoch_claude.cu: forward-only J_1 = tr[W_00 D_ov^{-1}], l=0, massless).
# Plan: final/analysis_scalar/scalar_fs_disc_impl_plan_claude.md.  Smoke: scalar_disc_smoke_claude.log.
#
# Scope (NM 2026-09-29): massless, at0.2 only, L1 (400 cfg, stride 10) + L2 (800 cfg, stride 5),
# nhits=2 solved as ONE mrhs block (NB_CLI=2 binaries), Nf2 first then Nf4, Nf6.
# Config sets = the same k's as the existing vector disc (k = 1, 11, 21, ... at L1; 1, 6, 11, ... at L2).
#
# PACKING: 4 workers = 2 processes per GPU on 2 GPUs under MPS.  Each worker takes an interleaved
# config sub-stream of the CURRENT ensemble (worker W starts at first+STRIDE*W, effective stride
# STRIDE*NW), so all workers finish an ensemble together and the load stays balanced.  Measured in the
# smoke run (L2): blocked + 2-packed = 54 s per hit per GPU (unblocked solo 122 s).
# RESUME-SAFE: per-hit "complete" gating + config-level early skip in the driver; re-running the same
# command continues where it stopped.  Output: data_<ens>_vm.../corr_ylm_disc_scalar_tb2/corr.<k>.h<h>.h5
# (the vector disc dir corr_ylm_disc_tb2/ is never opened).  No rm, no kill anywhere in this script.
#
# RUN DETACHED (long GPU job -- you launch it, not Claude):
#     nohup bash run_disc_scalar_L1L2_claude.sh > /dev/null 2>&1 &
#     NF=4 nohup bash run_disc_scalar_L1L2_claude.sh > /dev/null 2>&1 &      # later: Nf4, then NF=6
#     LS="2" GPUS="0 0" NF=2 bash run_disc_scalar_L1L2_claude.sh              # L2 only, 2 workers on GPU0
# Master log: run_disc_scalar_Nf<NF>_claude.log ; per-worker logs: disc_scalar_L<L>_g<gsq>_Nf<nf>_w<W>_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

NF="${NF:-2}"
LS="${LS:-1 2}"
read -ra GPUS <<< "${GPUS:-0 0 1 1}"     # 4 workers, 2 per GPU (MPS packed)
NHITS=2
DISC_TB=2
NB=2                                     # hits per block solve (must equal NHITS here)
stride_of_L () {                         # L1 -> 400 cfg of 4000 traj; L2 -> 800 cfg
  case "$1" in
    1) echo 10 ;;
    2) echo 5 ;;
    *) echo 10 ;;
  esac
}

MASTER="run_disc_scalar_Nf${NF}_claude.log"
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
echo "### scalar FS disc production  NF=$NF  LS=$LS  GPUS=${GPUS[*]}  nhits=$NHITS tb=$DISC_TB NB=$NB  [$(date +%F_%H:%M:%S)] ###"
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

# ---- collect ensembles: massless, at0.2, Nf=$NF, L in $LS; order = L ascending, gsq ascending ----
ENS=()
for L in $LS
do
  for d in $(ls -d Nf${NF}_gsq*at0.200000*mRe0.000000mIm0.000000nt128L${L}_hb* 2>/dev/null | sort -t q -k2,2)
  do
    [ -d "$d" ] || continue
    ls "$d"/ckpoint_lat.* >/dev/null 2>&1 || continue
    ENS+=("$d")
  done
done
if [ "${#ENS[@]}" -eq 0 ]
then
  echo "### no ensembles found for NF=$NF LS=$LS -- nothing to do ###"
  exit 0
fi
echo "### ${#ENS[@]} ensembles, ${#GPUS[@]} workers (GPUs ${GPUS[*]}), MPS packed ###"
for e in "${ENS[@]}"
do
  echo "###   $e  (stride $(stride_of_L "$(get_L "$e")"))"
done

run_disc () {   # $1=ens-dir  $2=GPU  $3=WID  $4=NWORKERS
  local ens="$1"
  local GPU="$2"
  local WID="$3"
  local NW="$4"
  local L nf gsq bin LOG ks first last kmin_w stride stride_eff kmax nconf
  L=$(get_L "$ens")
  nf=$(get_nf "$ens")
  gsq=$(get_gsq "$ens")
  stride=$(stride_of_L "$L")
  bin="jj_local_ylm_scalar_disc_stoch_L${L}_NB${NB}.o"
  LOG="disc_scalar_L${L}_g${gsq}_Nf${nf}_w${WID}_claude.log"
  mapfile -t ks < <(ls "$ens"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n)
  if [ "${#ks[@]}" -eq 0 ]
  then
    echo "### SKIP $ens: no ckpoint_lat  [$(date +%F_%H:%M:%S)] ###" >> "$LOG"
    return 0
  fi
  first="${ks[0]}"
  last="${ks[-1]}"
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
    echo "### DISC-scalar $ens  L${L}  GPU${GPU}  w${WID}  kmin=$kmin_w stride=$stride_eff kmax=$kmax tb=$DISC_TB nhits=$NHITS NB=$NB (~$nconf cfg)  [$(date +%F_%H:%M:%S)] ###"
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
  local e
  for e in "${ENS[@]}"
  do
    echo "### worker $WID (GPU $GPU) -> $e  [$(date +%F_%H:%M:%S)] ###"
    run_disc "$e" "$GPU" "$WID" "$NW"
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

echo "### START scalar FS disc NF=$NF: ${#GPUS[@]} workers (GPUs ${GPUS[*]}), MPS  [$(date +%F_%H:%M:%S)] ###"
for (( W=0; W<${#GPUS[@]}; W++ ))
do
  worker "$W" &
done
wait
echo "### ALL WORKERS DONE NF=$NF  [$(date +%F_%H:%M:%S)] ###"
