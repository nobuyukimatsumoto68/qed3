#!/usr/bin/env bash
# CONNECTED local-current + scalar Y_lm tower (FULL: vector s1/s2/s3 + axial + sigma_PS/sigma_FS, l=0..3) for the
# 9 MASSLESS gsq=8 ensembles (Nf 2/4/6 x L1/L2/L4), at STRIDE=2, nhit=1, t0=0, spin-dilution.  FRESH mode
# (NOT --IsScalarOnly) -> writes h0/ylm/s{a}, h0/ylm_axial/s{a}, h0/scalar, h0/scalar_fs + "complete".
# Resumable: the driver skips a config whose file is "complete", so re-running continues.  Plan:
# conn_gsq8_massless_impl_plan_claude.md.
#
# Run detached:
#   nohup bash run_conn_gsq8_massless_claude.sh > run_conn_gsq8_massless_claude.log 2>&1 &
# Env: GPUS="0 1" (space-sep list), JOBS_PER_GPU=2, NF=all|2|4|6, FILTER=all|L1|L2|L4, STRIDE=2, KMIN_TRAJ=1.
# Scheduling: each worker does 1/NWORKERS of EVERY selected ensemble (config-range split), workers start on
# rotated ensembles so all L run in parallel + load balanced.  Per-worker logs conn_massless_L{L}_Nf{nf}_w{W}.
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

GPUS="${GPUS:-0 1}"
JOBS_PER_GPU="${JOBS_PER_GPU:-2}"
FILTER="${FILTER:-all}"
STRIDE="${STRIDE:-2}"
KMIN_TRAJ="${KMIN_TRAJ:-1}"     # skip trajectories numbered below this (thermalization; massless already long)
NHITS=1

# worker -> gpu map (each GPU repeated JOBS_PER_GPU times), or explicit WGPU="0 0 1 1"
if [ -n "${WGPU:-}" ]; then
  read -ra WG <<< "$WGPU"
else
  WG=()
  for g in $GPUS; do for ((j=0; j<JOBS_PER_GPU; j++)); do WG+=("$g"); done; done
fi
NWORKERS=${#WG[@]}
WGPU_STR="${WG[*]}"

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
SRC=jj_local_ylm_scalar_conn_stoch_claude.cu

# ---- ensure MPS daemon up (packing needs it; else workers serialize via context switch) ----
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

# ---- build the three per-L binaries (skip if up-to-date; FORCE_BUILD=1 to force) ----
need_build () {   # $1 = output binary
  [ -n "${FORCE_BUILD:-}" ] && return 0
  [ ! -f "$1" ] && return 0
  find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$1" -print -quit 2>/dev/null | grep -q . && return 0
  return 1
}
declare -A BIN
for L in 1 2 4
do
  B="jj_local_ylm_scalar_conn_stoch_L${L}.o"
  if need_build "$B"
  then
    echo "### compile L${L} (-DN_REFINE_CLI=${L}) -> $B  [$(date +%F_%H:%M:%S)] ###"
    $NVCC $NVCCBASE -DN_REFINE_CLI=${L} $INCLUDES $LDFLAGS "$SRC" -o "$B" \
      || { echo "### CONN L${L} BUILD FAILED ###"; exit 1; }
  else
    echo "### L${L} binary up-to-date, skip (FORCE_BUILD=1 to rebuild) ###"
  fi
  BIN[$L]="$B"
done
echo "### build OK.  $NWORKERS workers, map=[$WGPU_STR], FILTER=$FILTER STRIDE=$STRIDE ###"

# ---- 9 massless ensembles as "ens|L" (valence mass = 0) ----
ENS_ALL=(
  "Nf2_gsq8.000000at0.200000nu01.000000nt128L1|1"
  "Nf4_gsq8.000000at0.200000nu01.000000nt128L1|1"
  "Nf6_gsq8.000000at0.200000nu01.000000nt128L1|1"
  "Nf2_gsq8.000000at0.200000nu01.000000nt128L2|2"
  "Nf4_gsq8.000000at0.200000nu01.000000nt128L2|2"
  "Nf6_gsq8.000000at0.200000nu01.000000nt128L2|2"
  "Nf2_gsq8.000000at0.200000nu01.000000nt128L4|4"
  "Nf4_gsq8.000000at0.200000nu01.000000nt128L4|4"
  "Nf6_gsq8.000000at0.200000nu01.000000nt128L4|4"
)
# apply FILTER (by L) and NF (by flavor; NF=all|2|4|6)
NF="${NF:-all}"
ENS_LIST=()
for spec in "${ENS_ALL[@]}"
do
  L="${spec##*|}"
  nf="${spec#Nf}"
  nf="${nf%%_*}"
  case "$FILTER" in
    all) ;;
    L1)  [ "$L" = 1 ] || continue;;
    L2)  [ "$L" = 2 ] || continue;;
    L4)  [ "$L" = 4 ] || continue;;
  esac
  if [ "$NF" != "all" ] && [ "$nf" != "$NF" ]; then continue; fi
  ENS_LIST+=("$spec")
done
NENS=${#ENS_LIST[@]}
echo "### ensembles (NF=$NF FILTER=$FILTER): $NENS -> ${ENS_LIST[*]} ###"

# Each worker does 1/NWORKERS of EVERY ensemble's stride-STRIDE grid (config-range split): worker WID takes the
# sub-stream starting at kmin + STRIDE*WID with effective stride STRIDE*NWORKERS.  Union over workers = full
# stride-STRIDE grid, no collision (base spacing is 1 here, so kmin+STRIDE*WID is always a real config).  The
# ensemble ORDER is ROTATED by WID so at t=0 the workers sit on DIFFERENT ensembles -> L1/L2/L4 all progress in
# PARALLEL from the start (and load stays balanced: each worker = 1/N of every ensemble).
run_one () {   # $1="ens|L"  $2=GPU  $3=WID
  local spec="$1" GPU="$2" WID="$3" ens L Nf
  ens="${spec%|*}"
  L="${spec##*|}"
  Nf=${ens#Nf}
  Nf=${Nf%%_*}
  local LOG="conn_massless_L${L}_Nf${Nf}_w${WID}_claude.log"
  local ks first last kmin kmin_w stride_eff kmax nconf k
  mapfile -t ks < <(ls "$ens"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n)
  if [ "${#ks[@]}" -eq 0 ]; then echo "### SKIP $ens: no ckpoint_lat  [$(date +%F_%H:%M:%S)] ###" >> "$LOG"; return 0; fi
  first="${ks[0]}"
  last="${ks[-1]}"
  kmin="$first"
  for k in "${ks[@]}"; do if [ "$k" -ge "$KMIN_TRAJ" ]; then kmin="$k"; break; fi; done
  kmin_w=$(( kmin + STRIDE * WID ))
  stride_eff=$(( STRIDE * NWORKERS ))
  kmax=$(( last + 1 ))
  if [ "$kmin_w" -ge "$kmax" ]; then echo "### $ens w$WID: no configs in sub-stream  [$(date +%F_%H:%M:%S)] ###" >> "$LOG"; return 0; fi
  nconf=$(( (last - kmin_w) / stride_eff + 1 ))
  {
    echo "### CONN $ens  L${L}  GPU${GPU}  w${WID}  kmin=$kmin_w stride=$stride_eff kmax=$kmax (~$nconf cfg)  [$(date +%F_%H:%M:%S)] ###"
    CUDA_VISIBLE_DEVICES=$GPU ./"${BIN[$L]}" --ens-dir "$ens/" --kmin "$kmin_w" --stride "$stride_eff" --kmax "$kmax" \
      --nhits "$NHITS" --t0 0 --spin-dilution
    echo "### CONN $ens w${WID} done (status $?)  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

run_worker () {   # $1=WID  $2=GPU  -- rotated ensemble order, config-split sub-stream WID
  local WID="$1" GPU="$2" j idx spec
  for ((j=0; j<NENS; j++))
  do
    idx=$(( (WID + j) % NENS ))
    spec="${ENS_LIST[$idx]}"
    run_one "$spec" "$GPU" "$WID"
  done
}

echo "### START conn-massless ($NWORKERS workers, config-split + rotated order)  [$(date +%F_%H:%M:%S)] ###"
for ((w=0; w<NWORKERS; w++))
do
  echo "### launch worker $w on GPU ${WG[$w]} ###"
  run_worker "$w" "${WG[$w]}" &
done
wait
echo "### ALL conn-massless done (NF=$NF FILTER=$FILTER STRIDE=$STRIDE)  [$(date +%F_%H:%M:%S)] ###"
