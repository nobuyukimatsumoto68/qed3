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
# PACKING: CUDA MPS must be up (`nvidia-cuda-mps-control -d`).  Default = 2 jobs on GPU 0.  Workers split the
# ensemble list round-robin (worker w takes ensembles with index % NWORKERS == w).  Each worker logs to
# ylm_scalar_catchup_w<w>_claude.log ; build+orchestration -> ylm_scalar_catchup_claude.log .
#
# Env: GPUS="0" (space-sep list), JOBS_PER_GPU=2, FILTER=all|L1|L2, PHASE=both|conn|disc.
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

GPUS="${GPUS:-0}"
JOBS_PER_GPU="${JOBS_PER_GPU:-2}"
FILTER="${FILTER:-all}"
PHASE="${PHASE:-both}"
NGPU=$(echo $GPUS | wc -w)
NWORKERS=$(( NGPU * JOBS_PER_GPU ))

LOG=ylm_scalar_catchup_claude.log
: > "$LOG"
exec > >(tee -a "$LOG") 2>&1

NVCC=nvcc
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
flags () { echo "-arch=sm_70 -g -O3 -std=c++20 -DN_REFINE_CLI=$1 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"; }

echo "### [build] 4 binaries (conn/disc x NREF 1/2)  [$(date +%F_%H:%M:%S)] ###"
declare -A CONN DISC
for NREF in 1 2; do
  CB="jj_local_ylm_scalar_conn_stoch_L${NREF}.o"
  DB="jj_local_ylm_scalar_disc_stoch_L${NREF}.o"
  $NVCC $(flags $NREF) $INCLUDES $LDFLAGS jj_local_ylm_scalar_conn_stoch_claude.cu -o "$CB" || { echo "### CONN L$NREF BUILD FAILED ###"; exit 1; }
  $NVCC $(flags $NREF) $INCLUDES $LDFLAGS jj_local_ylm_scalar_disc_stoch_claude.cu -o "$DB" || { echo "### DISC L$NREF BUILD FAILED ###"; exit 1; }
  CONN[$NREF]="$CB"
  DISC[$NREF]="$DB"
done
echo "### build OK.  $NWORKERS workers (GPUS='$GPUS' x $JOBS_PER_GPU/gpu), FILTER=$FILTER PHASE=$PHASE ###"

# emit this worker's share of ensembles (LPT cost-balanced): ENSDIR|MASS|NREF|KMIN|STRIDE|KMAX|NCFG|UNIFORM|SEA
emit () {   # $1=NWORKERS  $2=WID
  python3 - "$FILTER" "$1" "$2" <<'PY'
import glob, os, re, sys
FILTER, NW, WID = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
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
    cost = (4 if nref == 2 else 1) * len(ks)   # L2 ~4x per-config (N=42 vs 12); rough balance weight
    line = f"{ens_base}|{mass}|{nref}|{ks[0]}|{stride}|{ks[-1]+1}|{len(ks)}|{uni}|{sea}"
    ens.append((cost, line))
# LPT: heaviest first, each to the currently-lightest worker (deterministic -> all workers agree)
ens.sort(key=lambda e: -e[0])
loads = [0]*NW
assign = [[] for _ in range(NW)]
for cost, line in ens:
    w = min(range(NW), key=lambda k: loads[k])
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
  P=$(emit "$NWORKERS" "$WID")
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

# launch NWORKERS: JOBS_PER_GPU on each GPU in GPUS (MPS-packed)
wid=0
for GPU in $GPUS; do
  for ((j=0; j<JOBS_PER_GPU; j++)); do
    echo "### launch worker $wid on GPU $GPU -> ylm_scalar_catchup_w${wid}_claude.log ###"
    run_worker "$wid" "$GPU" &
    wid=$((wid+1))
  done
done
wait
echo "### ALL scalar catch-up done ($NWORKERS workers, FILTER=$FILTER PHASE=$PHASE)  [$(date +%F_%H:%M:%S)] ###"
