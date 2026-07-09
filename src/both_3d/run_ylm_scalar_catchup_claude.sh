#!/usr/bin/env bash
# GPU CATCH-UP (MPS-packed): append the scalar (sigma_PS/sigma_FS) Y_lm tower onto every ensemble that ALREADY
# has the jj ylm tower (conn + disc tb2).  Runs the FURNISHED scalar drivers in --IsScalarOnly mode, which
# re-solves (needs the gauge fields) and APPENDS h0/scalar + h0/scalar_fs (conn) and h0/disc/ylm/s0 + s0_1mD
# (disc) into the existing per-config .h5, leaving the vector/axial keys untouched.
#
# Covers all 30 interacting ensembles auto-discovered from data_*/corr_ylm_conn_t00_nhits1_s1/ :
#   L1 (NREF=1): Nf{2,4,6} x {massless, mRe 0.01,0.05,0.10,0.20}                    (140 cfg each)
#   L2 (NREF=2): Nf{2,4,6} x {massless, mRe 0.010572,0.052862,0.105725,0.211450}    (62-140 cfg)
# Per ensemble the config range (kmin/stride/kmax) is auto-derived from the existing conn files (L1 step 8,
# L2 step 4), so we hit EXACTLY the jj-processed configs.  Resume-safe: the driver skips a config whose scalar
# keys already exist, so re-running continues where it left off.
#
# PACKING: this script now AUTO-STARTS CUDA MPS (and aborts if it can't) -- REQUIRED for packing, else two
# procs on one GPU serialize via context-switch (the ~8x slowdown we hit was a DEAD MPS daemon, not a packing
# pathology; with MPS up, 2/GPU gives ~2x aggregate).  Default = 2 jobs on GPU 0 (WGPU="0 0").  Workers split
# the ensemble list by cost-weighted LPT.  Each worker logs to ylm_scalar_catchup_w<w>_claude.log ;
# build+orchestration -> ylm_scalar_catchup_claude.log .
#
# Env: GPUS="0" (space-sep list), JOBS_PER_GPU=2, FILTER=all|L1|L2, PHASE=both|conn|disc.
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

GPUS="${GPUS:-0}"
JOBS_PER_GPU="${JOBS_PER_GPU:-2}"
FILTER="${FILTER:-all}"
PHASE="${PHASE:-both}"
# WGPU = explicit GPU id per worker (space-sep), e.g. "0 0 1" = 2 workers on GPU0 + 1 on GPU1.
# If unset, expand GPUS x JOBS_PER_GPU (each GPU repeated JOBS_PER_GPU times).
if [ -n "${WGPU:-}" ]; then
  read -ra WG <<< "$WGPU"
else
  WG=()
  for g in $GPUS; do for ((j=0; j<JOBS_PER_GPU; j++)); do WG+=("$g"); done; done
fi
NWORKERS=${#WG[@]}
WGPU_STR="${WG[*]}"

LOG=ylm_scalar_catchup_claude.log
: > "$LOG"
exec > >(tee -a "$LOG") 2>&1

# ---- ensure CUDA MPS is up BEFORE launching workers (else 2-packing serializes via context-switch = the
#      ~8x slowdown we hit; a dead MPS daemon, NOT a packing pathology).  Mirrors the HMC run scripts. ----
if pgrep -f nvidia-cuda-mps-control >/dev/null; then
  echo "### MPS daemon: already running ###"
else
  echo "### MPS daemon not up -- starting nvidia-cuda-mps-control -d ###"
  nvidia-cuda-mps-control -d
  for i in 1 2 3 4 5; do pgrep -f nvidia-cuda-mps-control >/dev/null && break; sleep 1; done
fi
pgrep -f nvidia-cuda-mps-control >/dev/null \
  || { echo "### ERROR: MPS daemon failed to start -- aborting (would run un-packed/serialized) ###"; exit 1; }

NVCC=nvcc
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
flags () { echo "-arch=sm_70 -g -O3 -std=c++20 -DN_REFINE_CLI=$1 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"; }

# Build only when needed: skip if the binary exists and NO source (.cu/.h under ./ or includes/) is newer.
# Set FORCE_BUILD=1 to always rebuild.  (On a restart with unchanged source this avoids the redundant nvcc.)
need_build () {   # $1 = output binary
  [ -n "${FORCE_BUILD:-}" ] && return 0
  [ ! -f "$1" ] && return 0
  find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$1" -print -quit 2>/dev/null | grep -q . && return 0
  return 1
}
build_one () {   # $1=src  $2=out  $3=nref  $4=label
  if need_build "$2"; then
    echo "### [build] $4 -> $2 ###"
    $NVCC $(flags $3) $INCLUDES $LDFLAGS "$1" -o "$2" || { echo "### $4 BUILD FAILED ###"; exit 1; }
  else
    echo "### [build] $2 up-to-date, skip (FORCE_BUILD=1 to rebuild) ###"
  fi
}
echo "### [build] check 4 binaries (conn/disc x NREF 1/2)  [$(date +%F_%H:%M:%S)] ###"
declare -A CONN DISC
for NREF in 1 2; do
  CB="jj_local_ylm_scalar_conn_stoch_L${NREF}.o"
  DB="jj_local_ylm_scalar_disc_stoch_L${NREF}.o"
  build_one jj_local_ylm_scalar_conn_stoch_claude.cu "$CB" "$NREF" "CONN L$NREF"
  build_one jj_local_ylm_scalar_disc_stoch_claude.cu "$DB" "$NREF" "DISC L$NREF"
  CONN[$NREF]="$CB"
  DISC[$NREF]="$DB"
done
echo "### build OK.  $NWORKERS workers, worker->gpu map = [$WGPU_STR], FILTER=$FILTER PHASE=$PHASE ###"

# emit this worker's share of ensembles (co-residency-weighted LPT): ENSDIR|MASS|NREF|KMIN|STRIDE|KMAX|NCFG|UNIFORM|SEA
# Balances by COMPLETION TIME: worker w minimizes (load_w + cost)/speed_w, where speed_w reflects MPS sharing
# (a GPU running k packed workers => each ~ solo[k]).  So an un-contended GPU (e.g. GPU1 solo) gets more work.
emit () {   # $1=WGPU string (gpu id per worker)  $2=WID
  python3 - "$FILTER" "$1" "$2" <<'PY'
import glob, os, re, sys
from collections import Counter
FILTER = sys.argv[1]
wgpu = sys.argv[2].split()
WID = int(sys.argv[3])
NW = len(wgpu)
cnt = Counter(wgpu)
solo = {1: 1.0, 2: 0.85, 3: 0.72}                      # fraction of solo speed per packed worker
speed = [solo.get(cnt[g], max(0.5, 2.2/cnt[g])) for g in wgpu]
ens = []
for cd in sorted(glob.glob('data_*/corr_ylm_conn_t00_nhits1_s1')):
    esnid = cd.split('/')[0][len('data_'):]
    if esnid.startswith('free'):
        continue
    nref = 2 if 'nt128L2' in esnid else 1
    if FILTER == 'L1' and nref != 1:
        continue
    if FILTER == 'L2' and nref != 2:
        continue
    ks = sorted(int(f.rsplit('/', 1)[-1].split('.')[1]) for f in glob.glob(cd + '/corr.*.h0.h5'))
    if not ks:
        continue
    ens_base = re.sub(r'_vmRe[0-9.]+vmIm[0-9.]+$', '', esnid)
    mass = re.search(r'_vmRe([0-9.]+)vmIm', esnid).group(1)
    diffs = sorted(set(ks[i+1]-ks[i] for i in range(len(ks)-1))) if len(ks) > 1 else [1]
    stride = diffs[0]
    uni = 'UNIFORM' if len(diffs) == 1 else 'NON-UNIF' + str(diffs)
    sea = os.path.isdir(ens_base + '/') and bool(glob.glob(ens_base + '/ckpoint_lat.*'))
    cost = (4 if nref == 2 else 1) * len(ks)           # L2 ~4x per-config (N=42 vs 12)
    line = f"{ens_base}|{mass}|{nref}|{ks[0]}|{stride}|{ks[-1]+1}|{len(ks)}|{uni}|{sea}"
    ens.append((cost, line))
ens.sort(key=lambda e: -e[0])                          # heaviest first (LPT); deterministic -> all workers agree
loads = [0.0]*NW
assign = [[] for _ in range(NW)]
for cost, line in ens:
    w = min(range(NW), key=lambda k: (loads[k]+cost)/speed[k])
    assign[w].append(line)
    loads[w] += cost
for line in assign[WID]:
    print(line)
PY
}

run_worker () {   # $1=WID  $2=GPU
  local WID=$1 GPU=$2
  local wlog="ylm_scalar_catchup_w${WID}_claude.log"
  local P
  P=$(emit "$WGPU_STR" "$WID")
  {
    local n
    n=$(printf '%s\n' "$P" | grep -c .)
    echo "### worker $WID -> GPU $GPU : $n ensembles  [$(date +%F_%H:%M:%S)] ###"
    local i=0
    printf '%s\n' "$P" | while IFS='|' read -r EB MASS NREF KMIN STRIDE KMAX NC UNI SEA; do
      i=$((i+1))
      if [ "$SEA" != "True" ]; then
        echo "### w$WID [$i/$n] SKIP $EB -- sea config dir missing ###"
        continue
      fi
      echo "### w$WID [$i/$n] $EB  mass=$MASS L$NREF  k=$KMIN..$((KMAX-1)) step $STRIDE ($NC cfg)  [$(date +%F_%H:%M:%S)] ###"
      COMMON="--ens-dir ${EB}/ --mass-re $MASS --kmin $KMIN --stride $STRIDE --kmax $KMAX --nhits 1 --IsScalarOnly"
      if [ "$PHASE" = "both" ] || [ "$PHASE" = "conn" ]; then
        echo "### w$WID   conn append ###"
        CUDA_VISIBLE_DEVICES=$GPU ./"${CONN[$NREF]}" $COMMON --t0 0 --spin-dilution || echo "### w$WID CONN FAILED $EB (status $?) -- continuing ###"
      fi
      if [ "$PHASE" = "both" ] || [ "$PHASE" = "disc" ]; then
        echo "### w$WID   disc append (tb2) ###"
        CUDA_VISIBLE_DEVICES=$GPU ./"${DISC[$NREF]}" $COMMON --disc-tblock 2 || echo "### w$WID DISC FAILED $EB (status $?) -- continuing ###"
      fi
    done
    echo "### worker $WID done  [$(date +%F_%H:%M:%S)] ###"
  } > "$wlog" 2>&1
}

# launch NWORKERS per the worker->gpu map WG (MPS-packed where a GPU repeats)
for ((w=0; w<NWORKERS; w++)); do
  echo "### launch worker $w on GPU ${WG[$w]} -> ylm_scalar_catchup_w${w}_claude.log ###"
  run_worker "$w" "${WG[$w]}" &
done
wait
echo "### ALL scalar catch-up done ($NWORKERS workers, FILTER=$FILTER PHASE=$PHASE)  [$(date +%F_%H:%M:%S)] ###"
