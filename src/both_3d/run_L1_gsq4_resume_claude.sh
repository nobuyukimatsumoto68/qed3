#!/bin/bash
# Resume the paused L=1 gsq=4.0 scan on GPU1, 2 jobs/GPU under MPS, in two concurrent slots:
#   slot 1 (group {2,4}): Nf=2 to kmax, THEN Nf=4   (sequential within the group)
#   slot 2 (group {6}) : Nf=6
# => starts running Nf2 + Nf6; when Nf2 finishes (~20h) Nf4 takes its slot. Always 2 procs.
# Uses current hmc_L1_claude.o (npole=21/window0.001, nsteps=10 all Nf, kmax=320).
# Resume points: Nf2 k=252, Nf4 k=135, Nf6 k=93.
#
# Run detached:
#   nohup bash run_L1_gsq4_resume_claude.sh > run_L1_gsq4_resume_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
cd "$SRCDIR" || exit 1
source ../../env.sh

GSQ=4.0
NU0=1.0
GPU=1

# ---- ensure the CUDA MPS daemon is up BEFORE launching (so the co-resident procs
# ---- actually SM-share instead of time-slicing). Clients must start AFTER the daemon.
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
echo "MPS server pipe: $(ls -d /tmp/nvidia-mps 2>/dev/null || echo default-not-yet-created)"

echo "===== build hmc_L1_claude.o  [$(date +%F_%H:%M:%S)] ====="
make hmc_L1_claude.o 2>&1 | tee build_hmc_L1_claude.log
test "${PIPESTATUS[0]}" -eq 0 || { echo "L1 BUILD FAILED"; exit 1; }

echo "===== resume gsq=4.0 on GPU${GPU}: slot1={Nf2->Nf4}, slot2={Nf6}  [$(date +%F_%H:%M:%S)] ====="

# slot 1: group {2,4} sequential
(
  for nf in 2 4
  do
    echo "### [slot1] gsq4 Nf${nf} start [$(date +%F_%H:%M:%S)] ###"
    CUDA_VISIBLE_DEVICES=${GPU} ./hmc_L1_claude.o ${GSQ} ${nf} ${NU0} \
      2>&1 | tee run_L1_gsq4_Nf${nf}_claude.log
    echo "### [slot1] gsq4 Nf${nf} done (status ${PIPESTATUS[0]}) [$(date +%F_%H:%M:%S)] ###"
  done
) &
SLOT1=$!

# slot 2: group {6}
CUDA_VISIBLE_DEVICES=${GPU} ./hmc_L1_claude.o ${GSQ} 6 ${NU0} \
  2>&1 | tee run_L1_gsq4_Nf6_claude.log &
SLOT2=$!

wait $SLOT1
wait $SLOT2
echo "===== gsq=4.0 resume (all of Nf2/4/6) finished  [$(date +%F_%H:%M:%S)] ====="
