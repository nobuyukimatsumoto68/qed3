#!/usr/bin/env bash
# run_conn_s2_local_claude.sh -- LOCAL (barracuda22) stride-2 connected Y_lm current-correlator tower.
# Plan: conn_s2_local_impl_plan_claude.md.  Blackboard: redo_ensembles_claude.txt UPDATE 2026-08-28.
#
# Fills the odd-k stride-2 grid (k = 1,3,5,7,9 mod 10) using the LOCAL relative-geometry driver
# jj_local_ylm_scalar_conn_stoch_claude.cu (NOT the _fnal copy, which hardcodes /project/... geometry).
# The existing conn tower already has k=1 mod10 (from run_conn_ext_claude.sh); this job adds the four
# new residue classes k=3,5,7,9 mod10 (driver --kmin 3|5|7|9 --stride 10), whose union with k=1 is the
# stride-2 grid.  Residue-class scheme mirrors run_conn_s2_fnal_claude.sh.  jj block-avg: arXiv:0804.1501.
#
# GPU0 only, MPS NPACK=2 -> 2 workers (override GPU_LIST="0 1" for both TITAN V = 4 workers).  PHASES:
#   PHASE A : L3 at=0.2 -- offsets 1 3 5 7 9  (INCLUDES the k=1 mod10 topup: L3 grew 799->999).
#   PHASE B : L1 at=0.1 -- offsets 1 3 5 7 9  (FNAL is OUT of conn -> local owns the FULL odd-k grid,
#             incl k=1 which was entirely absent locally: k=1 count was 0 on ALL 9 at0.1 ensembles).
# RESUMABLE: the driver skips per-config h5 gated on "complete" -> re-run is safe, no rework.  NO rm.
#
# Run detached:
#   nohup bash run_conn_s2_local_claude.sh > run_conn_s2_local_claude.log 2>&1 &
# Per-unit logs: conn_s2L<L>_g<gsq>_Nf<nf>_off<off>_w<wid>_claude.log (this dir).
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

# GPUs to use and MPS packing per GPU.  GPU0 only by default; override GPU_LIST="0 1" for both TITAN V.
read -ra GPUS <<< "${GPU_LIST:-0}"
NGPU=${#GPUS[@]}
NPACK=2
NWORKERS=$(( NGPU * NPACK ))
STRIDE=10
NHITS=1

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
SRC_CONN=jj_local_ylm_scalar_conn_stoch_claude.cu

# ---- build the per-L conn binary (skip if up-to-date; FORCE_BUILD=1 to force) ----
need_build () {
  [ -n "${FORCE_BUILD:-}" ] && return 0
  [ ! -f "$1" ] && return 0
  find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$1" -print -quit 2>/dev/null | grep -q . && return 0
  return 1
}
build_L () {
  local L="$1"
  local BIN="jj_local_ylm_scalar_conn_stoch_L${L}.o"
  if need_build "$BIN"
  then
    echo "### compile conn L${L} (-DN_REFINE_CLI=${L}) -> $BIN  [$(date +%F_%H:%M:%S)] ###"
    $NVCC $NVCCBASE -DN_REFINE_CLI=${L} $INCLUDES $LDFLAGS "$SRC_CONN" -o "$BIN" \
      || { echo "### conn L${L} BUILD FAILED ###"; exit 1; }
  else
    echo "### conn L${L} binary up-to-date, skip (FORCE_BUILD=1 to rebuild) ###"
  fi
}

# ---- ensure MPS daemon up (2-pack needs it; else the two workers serialize via context switch) ----
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

# ---- one driver call for one (ensemble, offset) work unit ----
run_unit () {   # $1=ens-dir  $2=offset(kmin)  $3=at  $4=L  $5=WID  $6=gpu
  local ens="$1" off="$2" at="$3" L="$4" WID="$5" gpu="$6"
  local nf gsq bin last kmax LOG
  nf=$(printf '%s' "$ens" | grep -oE '^Nf[0-9]+' | sed 's/Nf//')
  gsq=$(printf '%s' "$ens" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//')
  bin="jj_local_ylm_scalar_conn_stoch_L${L}.o"
  last=$(ls "$ens"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n | tail -1)
  if [ -z "$last" ]
  then
    echo "### SKIP $ens off$off (no ckpoint_lat) ###"
    return 0
  fi
  kmax=$(( last + 1 ))
  LOG="conn_s2L${L}_g${gsq}_Nf${nf}_off${off}_w${WID}_claude.log"
  {
    echo "### CONN-S2 $ens L${L} off$off  kmin=$off stride=$STRIDE kmax=$kmax at=$at gpu=$gpu  [$(date +%F_%H:%M:%S)] ###"
    CUDA_VISIBLE_DEVICES=$gpu ./"$bin" --gsq "$gsq" --Nf "$nf" --nu0 1.0 --at "$at" \
      --ens-dir "$ens/" --kmin "$off" --stride "$STRIDE" --kmax "$kmax" \
      --nhits "$NHITS" --t0 0 --spin-dilution
    echo "### CONN-S2 $ens off$off done (status $?)  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

# ---- one phase: (ensembles x offsets) work units, round-robin over NWORKERS (NGPU x NPACK, MPS) ----
run_phase () {   # $1=L  $2=at  $3=glob  $4..=offsets
  local L="$1" at="$2" glob="$3"
  shift 3
  local OFFS=( "$@" )
  local ENS=()
  local d
  for d in $glob
  do
    [ -d "$d" ] || continue
    case "$d" in *_vmRe*) continue;; esac
    case "$d" in *mRe0.000000*) ;; *) continue;; esac
    ls "$d"/ckpoint_lat.* >/dev/null 2>&1 && ENS+=("$d")
  done
  echo "### PHASE L${L} at=$at : ${#ENS[@]} ensembles x offsets(${OFFS[*]})  [$(date +%F_%H:%M:%S)] ###"
  local UNITS=()
  local e o
  for e in "${ENS[@]}"
  do
    for o in "${OFFS[@]}"
    do
      UNITS+=("$e|$o")
    done
  done
  local w
  for (( w=0; w<NWORKERS; w++ ))
  do
    (
      gpu=${GPUS[ w % NGPU ]}
      idx=$w
      while [ "$idx" -lt "${#UNITS[@]}" ]
      do
        u="${UNITS[$idx]}"
        ens="${u%|*}"
        off="${u#*|}"
        run_unit "$ens" "$off" "$at" "$L" "$w" "$gpu"
        idx=$(( idx + NWORKERS ))
      done
    ) &
  done
  wait
  echo "### PHASE L${L} at=$at done  [$(date +%F_%H:%M:%S)] ###"
}

build_L 3
build_L 1
echo "### build OK ###"

echo "### START conn-s2 local: GPUs=(${GPUS[*]}) x MPS ${NPACK}-pack = ${NWORKERS} workers, STRIDE=$STRIDE  [$(date +%F_%H:%M:%S)] ###"

# PHASE A : L3 at=0.2 -- k=1 topup + 4 new classes = full odd-k stride-2
run_phase 3 0.200000 'Nf*at0.200000*nt128L3_hb*' 1 3 5 7 9

# PHASE B : L1 at=0.1 -- FULL odd-k (FNAL is out of conn -> local owns k=1 too ; k=1 was 0 on all 9)
run_phase 1 0.100000 'Nf*at0.100000*nt128L1_hb*' 1 3 5 7 9

echo "### ALL conn-s2 local done  [$(date +%F_%H:%M:%S)] ###"
# ---- glueball REMAINDER (CPU, runs concurrently -- launch separately in its OWN terminal): ----
#   nohup bash run_glue_f2_v2_sweep_claude.sh > run_glue_f2_v2_sweep_claude.log 2>&1 &
# complete-gated; picks up L1 at=0.1 (Nf4/Nf6 g0.5/g1.5, un-run) + L4 partial. CPU-only -> no GPU contention.
