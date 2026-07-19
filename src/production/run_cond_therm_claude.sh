#!/usr/bin/env bash
# run_cond_therm_claude.sh -- DENSE condensate over the first 100 configs (stride 1) of every
# massive ensemble, to set the per-ensemble thermalization cut from the sigma_PS/sigma_FS MC
# series. Reuses the existing per-L binaries + the stride-20 h5 (complete-gated -> only new k's
# flow). Output -> data_<ens>_vmRe<m>.../corr_condensate_eo_nhits1/corr.<k>.h0.h5.
# 4 workers, 2/GPU MPS (same as run_cond_conn). Run detached:
#   nohup bash run_cond_therm_claude.sh > run_cond_therm_claude.log 2>&1 &
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
KMAX=100
NHITS=1
if [ -n "${WGPU:-}" ]; then read -ra GPUS <<< "$WGPU"; else GPUS=(0 0 1 1); fi
NW=${#GPUS[@]}

pgrep -f nvidia-cuda-mps-control >/dev/null || { nvidia-cuda-mps-control -d; sleep 2; }
pgrep -f nvidia-cuda-mps-control >/dev/null || { echo "### MPS failed ###"; exit 1; }

get () { printf '%s' "$1" | grep -oE "$2" | sed "$3"; }

run_one () {   # $1=ens  $2=GPU  $3=WID
  local ens="$1" GPU="$2" WID="$3" L mre gsq bin LOG kmin_w seff
  L=$(get "$ens" 'nt128L[0-9]+' 's/nt128L//')
  mre=$(get "$ens" 'mRe[0-9.]+mIm' 's/mRe//;s/mIm//')
  gsq=$(get "$ens" 'gsq[0-9.]+at' 's/gsq//;s/at//')
  bin="jj_scalar_condensate_eo_L${L}.o"
  LOG="condtherm_L${L}_g${gsq}_mRe${mre}_w${WID}_claude.log"
  kmin_w=$(( 1 + WID ))
  seff=$NW
  echo "### $ens GPU$GPU w$WID kmin=$kmin_w stride=$seff kmax=$KMAX [$(date +%T)] ###" >> "$LOG"
  CUDA_VISIBLE_DEVICES=$GPU ./"$bin" --ens-dir "$ens/" --kmin "$kmin_w" --stride "$seff" \
    --kmax "$KMAX" --nhits "$NHITS" --mass-re "$mre" >> "$LOG" 2>&1
  echo "### done (status $?) [$(date +%T)] ###" >> "$LOG"
}

MASSIVE=()
for d in Nf2_*_hb*; do
  [ -d "$d" ] || continue
  case "$d" in *mRe0.000000*) continue;; esac
  ls "$d"/ckpoint_lat.* >/dev/null 2>&1 && MASSIVE+=("$d")
done
echo "### ${#MASSIVE[@]} massive ensembles, first $KMAX cfgs, $NW workers [$(date +%T)] ###"

worker () {
  local WID="$1" GPU="${GPUS[$WID]}" n=${#MASSIVE[@]} i idx
  for (( i=0; i<n; i++ )); do
    idx=$(( (i + WID*n/NW) % n ))
    run_one "${MASSIVE[$idx]}" "$GPU" "$WID"
  done
}
for (( W=0; W<NW; W++ )); do worker "$W" & done
wait
echo "### ALL done [$(date +%T)] ###"
