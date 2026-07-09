#!/bin/bash
# L=1 massless gsq=16.0, Nf=2/4/6, kmax=320 -- QUEUED behind gsq=4 on GPU1.
# Waits for the gsq=4 streams to finish (GPU1 free), then launches with the same 2-slot
# grouping as gsq=4: slot1 = {Nf2 -> Nf4 sequential}, slot2 = {Nf6} (2 procs/GPU under MPS).
# All fresh (k=0). Uses current hmc_L1_claude.o (npole=21/window0.001, nsteps=10 all Nf).
# gsq=16 is stronger coupling than 12 -> slower + more freeze-prone; watch via /check-configs.
#
# Run detached NOW (it self-triggers when GPU1 frees):
#   nohup bash run_L1_gsq16_claude.sh > run_L1_gsq16_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
cd "$SRCDIR" || exit 1
source ../../env.sh

GSQ=16.0
NU0=1.0
GPU=1

# ---- 1. queue: wait for the gsq=4 streams on GPU1 to finish -----------------
echo "===== gsq=16 queued -- waiting for gsq=4 (GPU1) to finish  [$(date +%F_%H:%M:%S)] ====="
while pgrep -u "$USER" -f "hmc_L1_claude.o 4.0 " >/dev/null
do
  sleep 300
done
echo "===== GPU1 free -- launching gsq=16  [$(date +%F_%H:%M:%S)] ====="

# ---- 2. ensure the CUDA MPS daemon is up (SM-share, not time-slice) --------
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
  || { echo "ERROR: MPS daemon failed to start -- aborting (would run un-packed)"; exit 1; }

# ---- 3. build (up to date if unchanged) -----------------------------------
echo "===== build hmc_L1_claude.o  [$(date +%F_%H:%M:%S)] ====="
make hmc_L1_claude.o 2>&1 | tee build_hmc_L1_claude.log
test "${PIPESTATUS[0]}" -eq 0 || { echo "L1 BUILD FAILED"; exit 1; }

# ---- 4. launch on GPU1: slot1={Nf2->Nf4 seq}, slot2={Nf6} (2 procs, MPS) ---
echo "===== launch gsq=16 on GPU${GPU}: slot1={Nf2->Nf4}, slot2={Nf6}  [$(date +%F_%H:%M:%S)] ====="

(
  for nf in 2 4
  do
    echo "### [slot1] gsq16 Nf${nf} start [$(date +%F_%H:%M:%S)] ###"
    CUDA_VISIBLE_DEVICES=${GPU} ./hmc_L1_claude.o ${GSQ} ${nf} ${NU0} \
      2>&1 | tee run_L1_gsq16_Nf${nf}_claude.log
    echo "### [slot1] gsq16 Nf${nf} done (status ${PIPESTATUS[0]}) [$(date +%F_%H:%M:%S)] ###"
  done
) &
SLOT1=$!

CUDA_VISIBLE_DEVICES=${GPU} ./hmc_L1_claude.o ${GSQ} 6 ${NU0} \
  2>&1 | tee run_L1_gsq16_Nf6_claude.log &
SLOT2=$!

wait $SLOT1
wait $SLOT2
echo "===== gsq=16 (Nf2/4/6) finished  [$(date +%F_%H:%M:%S)] ====="
