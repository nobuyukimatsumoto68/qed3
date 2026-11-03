#!/bin/bash
# tmp_scalar_disc_smoke_claude.sh -- BUILD + SMOKE for the scalar FS disconnected loop driver
# jj_local_ylm_scalar_disc_stoch_claude.cu  (plan: final/analysis_scalar/scalar_fs_disc_impl_plan_claude.md).
#
# You run this (GPU job), Claude reads the log back:
#     bash tmp_scalar_disc_smoke_claude.sh
#     GPU=1 bash tmp_scalar_disc_smoke_claude.sh          # use GPU 1 instead of 0
#     TAG=b bash tmp_scalar_disc_smoke_claude.sh          # RE-RUN: finished hits are skipped, so a repeat
#                                                         # needs a new TAG (new scratch output dirs)
# All output -> scalar_disc_smoke_claude.log (appended).  Per-case logs: scalar_disc_smoke_<case>_claude.log.
#
# Phases:
#   1  build 4 binaries: L1, L2  x  NB_CLI=1 (unblocked reference), NB_CLI=2 (2 hits in one block)
#   2  L1 Nf2 g1.0, config k=1001, 2 hits: unblocked and blocked, SOLO
#   3  L2 Nf2 g2.0, config k=1001, 2 hits: unblocked and blocked, SOLO
#   4  L2 Nf2 g2.0, 2 processes on one GPU under MPS (k=1011 and k=2001): blocked, then unblocked
#   5  compare J_1 blocked vs unblocked; Re J_1 vs (1/2) tr W; size of Im J_1; timing table
# Output goes to SCRATCH dirs  data_<ens>/corr_ylm_disc_scalar_tb2_smoke<TAG>NB<n>[p<A|B>]/ ,
# never to the production dir corr_ylm_disc_scalar_tb2/ and never to the vector disc dir.
# Expected wall time: about 40 minutes.
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
LOG=scalar_disc_smoke_claude.log
exec > >(tee -a "$LOG") 2>&1

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
GPU="${GPU:-0}"
TAG="${TAG:-a}"
DISC_TB=2
NHITS=2

# OLD (nvcc from PATH):
# NVCC=nvcc
# CUDA 12.9 = last release that supports sm_70 (TITAN V); CUDA 13 dropped it.  Override with NVCC=<path>.
NVCC="${NVCC:-/usr/local/cuda-12.9/bin/nvcc}"
export LD_LIBRARY_PATH="/usr/local/cuda-12.9/lib64:${LD_LIBRARY_PATH:-}"
# Host compiler: the system default is gcc 15.2, which nvcc 12.9 rejects (supports gcc <= 14).
# g++-12 is installed (g++-13 is not), so pass it explicitly.  Override with CCBIN=<compiler>.
CCBIN="${CCBIN:-g++-12}"
# OLD (system default host compiler):
# NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
NVCCBASE="-ccbin $CCBIN -arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
SRC=jj_local_ylm_scalar_disc_stoch_claude.cu

ENS_L1="Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000"
ENS_L2="Nf2_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000"

echo "############################################################"
echo "### scalar disc smoke  GPU=$GPU  TAG=$TAG  [$(date +%F_%H:%M:%S)] ###"
echo "############################################################"

# ---------------- PHASE 1: build ----------------
echo "### PHASE 1: build ###"
for L in 1 2
do
  for NB in 1 2
  do
    BIN="jj_local_ylm_scalar_disc_stoch_L${L}_NB${NB}.o"
    echo "### compile L${L} NB${NB} -> $BIN  [$(date +%F_%H:%M:%S)] ###"
    $NVCC $NVCCBASE -DN_REFINE_CLI=${L} -DNB_CLI=${NB} $INCLUDES $LDFLAGS "$SRC" -o "$BIN"
    if [ $? -ne 0 ]
    then
      echo "### BUILD FAILED: L${L} NB${NB} ###"
      exit 1
    fi
  done
done
echo "### PHASE 1 build OK ###"

# run one case.  $1=L  $2=NB  $3=ens  $4=gsq  $5=k  $6=case-name
run_case () {
  local L="$1"
  local NB="$2"
  local ens="$3"
  local gsq="$4"
  local k="$5"
  local name="$6"
  local bin="jj_local_ylm_scalar_disc_stoch_L${L}_NB${NB}.o"
  local clog="scalar_disc_smoke_${name}_${TAG}_claude.log"
  local t0
  local t1
  t0=$(date +%s)
  echo "### CASE $name: L$L NB$NB k=$k  [$(date +%F_%H:%M:%S)] ###" > "$clog"
  CUDA_VISIBLE_DEVICES=$GPU ./"$bin" --gsq "$gsq" --Nf 2 --nu0 1.0 \
    --ens-dir "$ens/" --kmin "$k" --kmax "$(( k + 1 ))" --stride 1 \
    --nhits "$NHITS" --disc-tblock "$DISC_TB" --out-suffix "_smoke${TAG}${name}" >> "$clog" 2>&1
  local st=$?
  t1=$(date +%s)
  echo "### CASE $name done: status $st, wall $(( t1 - t0 )) s ###" >> "$clog"
}

# ---------------- PHASE 2: L1 solo ----------------
echo "### PHASE 2: L1 solo (unblocked, then blocked) ###"
run_case 1 1 "$ENS_L1" 1.000000 1001 L1NB1
cat "scalar_disc_smoke_L1NB1_${TAG}_claude.log"
run_case 1 2 "$ENS_L1" 1.000000 1001 L1NB2
cat "scalar_disc_smoke_L1NB2_${TAG}_claude.log"

# ---------------- PHASE 3: L2 solo ----------------
echo "### PHASE 3: L2 solo (unblocked, then blocked) ###"
run_case 2 1 "$ENS_L2" 2.000000 1001 L2NB1
cat "scalar_disc_smoke_L2NB1_${TAG}_claude.log"
run_case 2 2 "$ENS_L2" 2.000000 1001 L2NB2
cat "scalar_disc_smoke_L2NB2_${TAG}_claude.log"

# ---------------- PHASE 4: L2, 2 processes packed under MPS ----------------
echo "### PHASE 4: L2 2-packed under MPS ###"
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
if ! pgrep -f nvidia-cuda-mps-control >/dev/null
then
  echo "### ERROR: MPS daemon failed to start -- skipping PHASE 4 ###"
else
  echo "### packed, blocked (NB2): 2 processes ###"
  run_case 2 2 "$ENS_L2" 2.000000 1011 L2NB2pA &
  run_case 2 2 "$ENS_L2" 2.000000 2001 L2NB2pB &
  wait
  cat "scalar_disc_smoke_L2NB2pA_${TAG}_claude.log"
  cat "scalar_disc_smoke_L2NB2pB_${TAG}_claude.log"
  echo "### packed, unblocked (NB1): 2 processes ###"
  run_case 2 1 "$ENS_L2" 2.000000 1011 L2NB1pA &
  run_case 2 1 "$ENS_L2" 2.000000 2001 L2NB1pB &
  wait
  cat "scalar_disc_smoke_L2NB1pA_${TAG}_claude.log"
  cat "scalar_disc_smoke_L2NB1pB_${TAG}_claude.log"
fi

# ---------------- PHASE 5: checks + timing table ----------------
echo "### PHASE 5: checks ###"
TAG="$TAG" ENS_L1="$ENS_L1" ENS_L2="$ENS_L2" python3 - <<'PY'
import os
import numpy as np
import h5py

tag = os.environ["TAG"]
ens = {1: os.environ["ENS_L1"], 2: os.environ["ENS_L2"]}
vm = "_vmRe0.000000vmIm0.000000"


def load(L, name, k, h):
    d = "data_%s%s/corr_ylm_disc_scalar_tb2_smoke%s%s/" % (ens[L], vm, tag, name)
    fn = d + "corr.%d.h%d.h5" % (k, h)
    if not os.path.exists(fn):
        return None
    with h5py.File(fn, "r") as f:
        re = f["h0/disc/ylm/s0/l0/m0/J/real"][()]
        im = f["h0/disc/ylm/s0/l0/m0/J/imag"][()]
        half = float(f["h0/disc/ylm/s0/half_trW"][()][0])
    return re + 1j * im, half


for L in (1, 2):
    for h in (0, 1):
        a = load(L, "L%dNB1" % L, 1001, h)
        b = load(L, "L%dNB2" % L, 1001, h)
        if a is None or b is None:
            print("L%d hit %d: MISSING output (NB1 %s, NB2 %s)" % (L, h, a is not None, b is not None))
            continue
        ja = a[0]
        jb = b[0]
        half = a[1]
        dmax = np.max(np.abs(ja - jb))
        scale = np.max(np.abs(ja))
        print("L%d hit %d: max|J_NB2 - J_NB1| = %.3e   (max|J| = %.3e, relative %.3e)"
              % (L, h, dmax, scale, dmax / scale))
        dre = ja.real - half
        print("L%d hit %d: half_trW = %.6e ; Re J_1 - half_trW: mean over t = %.3e, rms over t = %.3e"
              % (L, h, half, dre.mean(), np.sqrt(np.mean(dre ** 2))))
        print("L%d hit %d: Im J_1: mean over t = %.3e, rms over t = %.3e"
              % (L, h, ja.imag.mean(), np.sqrt(np.mean(ja.imag ** 2))))
PY

echo "### timing table (seconds) ###"
for name in L1NB1 L1NB2 L2NB1 L2NB2 L2NB2pA L2NB2pB L2NB1pA L2NB1pB
do
  clog="scalar_disc_smoke_${name}_${TAG}_claude.log"
  [ -f "$clog" ] || continue
  echo "--- $name"
  grep -E 's per hit|wall|status' "$clog"
done
echo "### scalar disc smoke DONE  [$(date +%F_%H:%M:%S)] ###"
