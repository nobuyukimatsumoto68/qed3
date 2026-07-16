#!/bin/bash
# _claude: HMC PARAM TUNING (massless) for the CORRECTED gsq range. Sandbox = tuning/ (own includes/
# copy + driver copy; production/ symlinks both_3d and is untouched). Tune steps/tau/ladder in
# tuning/includes/hasenbusch_ladder_claude.h, rebuild, re-run, read dH/acceptance from the logs.
#
#   L1: gsq = 0.5, 1.0, 1.5      L2: gsq = 1.0, 2.0, 3.0      L4: gsq = 2.0, 4.0, 6.0   (Nf2, massless)
#
# MPS packing (per NM): ALL L packed 2/GPU across both GPUs (L4 MPS OK).
#   gpu0 = {L1 gsq0.5,1.0,1.5 slotA}  ||  {L4 gsq2.0,4.0 slotB}
#   gpu1 = {L2 gsq1.0,2.0,3.0 slotA}  ||  {L4 gsq6.0     slotB}
# -> 2 procs/GPU (4 total); weights ~balanced (gpu0: 3xL1+2xL4, gpu1: 3xL2+1xL4).
#
# KMAX bounds trajectory count for fast iteration (tuning copy of the driver, -DKMAX). Each (L,gsq)
# writes its own auto-resuming dir; per-traj "# dH : .. is_accept : .. rate : .." lines drive the tune.
#
# Run detached:  nohup bash run_tune_massless_claude.sh > run_tune_massless_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/tuning
cd "$SRCDIR" || exit 1
source ../../env.sh

NF=2
NU0=1.0
KMAX=25           # ~20 configs to see thermalization (per NM); bump if needed
KRNG=1            # FULL checkpointing (keep every ckpoint_rng; -DKRNG=1)

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=hmc_hasenbusch_claude.cu

# ---- 1. ensure the ambient MPS daemon is up (default pipe; do NOT start a custom pipe dir) --------
if pgrep -f nvidia-cuda-mps-control >/dev/null
then
  echo "MPS daemon: already running"
else
  echo "MPS daemon not up -- starting nvidia-cuda-mps-control -d"
  nvidia-cuda-mps-control -d
  for i in 1 2 3 4 5
  do
    pgrep -f nvidia-cuda-mps-control >/dev/null && break
    sleep 1
  done
fi
pgrep -f nvidia-cuda-mps-control >/dev/null \
  || { echo "ERROR: MPS daemon failed to start -- aborting"; exit 1; }

# ---- 2. build one binary per L (gsq is a runtime arg) --------------------------------------------
build () {
  local L=$1
  local out=hmc_tune_L${L}.out
  echo "===== build $out (LREF=$L KMAX=$KMAX)  [$(date +%F_%H:%M:%S)] ====="
  $NVCC $NVCCFLAGS $INCLUDES -DLREF=$L -DKMAX=$KMAX -DKRNG=$KRNG "$SRC" $LDFLAGS -o "$out" 2>&1 | tee build_tune_L${L}_claude.log
  test -f "$out" || { echo "BUILD FAILED: $out"; exit 1; }
}
build 1
build 2
build 4

# ---- 3. launch ----------------------------------------------------------------------------------
run_one () {   # gpu L gsq : ONE run (massless), tee to its own log
  local gpu=$1 L=$2 gsq=$3
  echo "### [gpu$gpu L$L gsq$gsq] start [$(date +%F_%H:%M:%S)] ###"
  CUDA_VISIBLE_DEVICES=$gpu ./hmc_tune_L${L}.out $gsq $NF $NU0 0.0 0.0 \
    2>&1 | tee run_tune_L${L}_gsq${gsq}_claude.log
  echo "### [gpu$gpu L$L gsq$gsq] done (status ${PIPESTATUS[0]}) [$(date +%F_%H:%M:%S)] ###"
}
run_seq () {   # gpu L gsq... : runs sequentially on one GPU (a packing slot / solo queue)
  local gpu=$1 L=$2; shift 2
  for g in "$@"; do run_one $gpu $L "$g"; done
}

echo "===== launch (all MPS, 2/GPU): gpu0={L1 || L4 2.0,4.0}, gpu1={L2 || L4 6.0}  [$(date +%F_%H:%M:%S)] ====="

run_seq 0 1 0.5 1.0 1.5 &   # gpu0 slot A: L1 all gsq
run_seq 0 4 2.0 4.0     &   # gpu0 slot B: L4 gsq2, gsq4
run_seq 1 2 1.0 2.0 3.0 &   # gpu1 slot A: L2 all gsq
run_seq 1 4 6.0         &   # gpu1 slot B: L4 gsq6

wait
echo "===== tuning runs finished  [$(date +%F_%H:%M:%S)] ====="
