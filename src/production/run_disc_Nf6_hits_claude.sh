#!/usr/bin/env bash
# run_disc_Nf6_hits_claude.sh -- extra stochastic hits for the vector disc on the Nf6 strong-coupling ensembles
# (NM 2026-10-23): L1 Nf6 g1.0, g1.5 -> 4 hits ; L2 Nf6 g2.0, g3.0 -> 6 hits.  a_t = 0.2, massless.
# Same driver and output dir as before (jj_local_ylm_disc_stoch_claude.cu -> data_<ens>/corr_ylm_disc_tb2/), on each
# ensemble's EXISTING config set, so every config ends with the same number of hits:
#   L1 Nf6 g1.0 : k = 1 mod 10                       (400 cfg, now 2 hits) -> NHITS=4
#   L1 Nf6 g1.5 : k = 1 mod 5                        (800 cfg, now 2 hits) -> NHITS=4
#   L2 Nf6 g2.0 : k = 1 mod 5  U  k = 7 mod 8        (1200 cfg, now 3 hits) -> NHITS=6
#   L2 Nf6 g3.0 : k = 1 mod 5  U  k = 1 mod 8        (1200 cfg, now 3 hits) -> NHITS=6
# The driver skips every (k,h) already "complete" and seeds hit h from (esnid,k,h): only the new hits are computed,
# existing hit files stay bit-identical.  A stopped run resumes cleanly on relaunch.
# Cost estimate: L1 g1.0 800 + L1 g1.5 1600 solves at ~80-100 s; L2 g2.0 3600 at ~155 s, L2 g3.0 3600 at ~230 s
#   (per proc, measured) -> ~1.6M proc-s / 4 procs ~ 4.7 days.
# PACKING: 4 workers = 2 procs/GPU on 2 GPUs under MPS (as the earlier disc runs).  GPUS="0 0" leaves GPU1 free.
# TOOLCHAIN: /usr/local/cuda-12.9/bin/nvcc -ccbin g++-12; cuda-12.9 lib64 first.  Rebuild only if the disc .cu or
#   ./includes/*.h is newer than the binary (FORCE_BUILD=1 forces).
# NO rm and NO kill anywhere (project policy).  Per-job logs: disc_nf6h_L*_g*_<grid>_w<W>_claude.log (this directory).
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
export LD_LIBRARY_PATH=/usr/local/cuda-12.9/lib64:${LD_LIBRARY_PATH:-}

DISC_TB=2
read -ra GPUS <<< "${GPUS:-0 0 1 1}"     # 4 workers, 2 per GPU (MPS packed)

NVCC=/usr/local/cuda-12.9/bin/nvcc
NVCCBASE="-ccbin g++-12 -arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
SRC=jj_local_ylm_disc_stoch_claude.cu

need_build () {
  [ -n "${FORCE_BUILD:-}" ] && return 0
  [ ! -f "$1" ] && return 0
  [ "$SRC" -nt "$1" ] && return 0
  find ./includes -maxdepth 1 -name '*.h' -newer "$1" -print -quit 2>/dev/null | grep -q . && return 0
  return 1
}
for TL in 1 2
do
  BIN="jj_local_ylm_disc_stoch_L${TL}.o"
  if need_build "$BIN"
  then
    echo "### compile L${TL} (-DN_REFINE_CLI=${TL}) -> $BIN  [$(date +%F_%H:%M:%S)] ###"
    $NVCC $NVCCBASE -DN_REFINE_CLI=${TL} $INCLUDES $LDFLAGS "$SRC" -o "$BIN" \
      || { echo "### DISC L${TL} BUILD FAILED ###"; exit 1; }
  else
    echo "### L${TL} binary up-to-date, skip ###"
  fi
done
echo "### build OK ###"

# ---- job list: ens-dir | L | Nf | gsq | grid name | first k | stride | NHITS ----
JOBS=(
  "Nf6_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000|1|6|1.000000|A|1|10|4"
  "Nf6_gsq1.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000|1|6|1.500000|A|1|5|4"
  "Nf6_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000|2|6|2.000000|A|1|5|6"
  "Nf6_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000|2|6|2.000000|B|7|8|6"
  "Nf6_gsq3.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000|2|6|3.000000|A|1|5|6"
  "Nf6_gsq3.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000|2|6|3.000000|B|1|8|6"
)
for j in "${JOBS[@]}"
do
  ens="${j%%|*}"
  [ -d "$ens" ] || { echo "### ERROR: missing ensemble dir $ens ###"; exit 1; }
  echo "###   $j ###"
done

run_job () {   # $1=job string  $2=GPU  $3=WID  $4=NWORKERS
  local job="$1" GPU="$2" WID="$3" NW="$4"
  local ens L nf gsq GRID K0 STRIDE NH last kmin_w stride_eff kmax nconf LOG
  IFS='|' read -r ens L nf gsq GRID K0 STRIDE NH <<< "$job"
  LOG="disc_nf6h_L${L}_g${gsq}_Nf${nf}_${GRID}_w${WID}_claude.log"
  last=$(find "$ens" -maxdepth 1 -name 'ckpoint_lat.*' | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n | tail -1)
  if [ -z "$last" ]
  then
    echo "### SKIP $ens: no ckpoint_lat  [$(date +%F_%H:%M:%S)] ###" >> "$LOG"
    return 0
  fi
  kmin_w=$(( K0 + STRIDE * WID ))
  stride_eff=$(( STRIDE * NW ))
  kmax=$(( last + 1 ))
  if [ "$kmin_w" -ge "$kmax" ]
  then
    echo "### $ens $GRID w$WID: no configs in sub-stream  [$(date +%F_%H:%M:%S)] ###" >> "$LOG"
    return 0
  fi
  nconf=$(( (last - kmin_w) / stride_eff + 1 ))
  {
    echo "### DISC-nf6h $GRID $ens  L${L}  GPU${GPU}  w${WID}  kmin=$kmin_w stride=$stride_eff kmax=$kmax tb=$DISC_TB nhits=$NH (~$nconf cfg)  [$(date +%F_%H:%M:%S)] ###"
    CUDA_VISIBLE_DEVICES=$GPU ./jj_local_ylm_disc_stoch_L${L}.o --gsq "$gsq" --Nf "$nf" --nu0 1.0 \
      --ens-dir "$ens/" --kmin "$kmin_w" --stride "$stride_eff" --kmax "$kmax" \
      --nhits "$NH" --disc-tblock "$DISC_TB"
    echo "### DISC-nf6h $GRID $ens w${WID} done (status $?)  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

worker () {   # $1=WID
  local WID="$1"
  local GPU="${GPUS[$WID]}"
  local n=${#JOBS[@]} NW=${#GPUS[@]}
  local i idx
  for (( i=0; i<n; i++ ))
  do
    idx=$(( (i + WID*n/NW) % n ))
    run_job "${JOBS[$idx]}" "$GPU" "$WID" "$NW"
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

echo "### START Nf6 extra hits: ${#GPUS[@]} workers (GPUs ${GPUS[*]}), MPS  [$(date +%F_%H:%M:%S)] ###"
for (( W=0; W<${#GPUS[@]}; W++ ))
do
  worker "$W" &
done
wait
echo "### ALL Nf6 extra hits done  [$(date +%F_%H:%M:%S)] ###"
