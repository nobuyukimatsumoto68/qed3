#!/bin/bash
# Rebuild hmc_L4_claude.o (now nsteps=8, matching hmc_fermilab_wmass_L4_claude.cu)
# and restart BOTH L=4 streams (Nf=4 and Nf=6) on GPU1 under MPS.
# The L=2 streams on GPU0 are left running untouched (their nsteps=5 already
# matches the fermilab L2 code).
#
# Both L4 streams resume from their existing checkpoints (Nf4 ~k=83, Nf6 ~k=0).
#
# Run detached, e.g.:
#   nohup bash restart_L4_nsteps8_claude.sh > restart_L4_nsteps8_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
cd "$SRCDIR" || exit 1
source ../../env.sh

GSQ=8.0
NU0=1.0
PAT="hmc_L4_claude.o"

# ---- 1. stop BOTH L=4 streams (Nf=4 and Nf=6); leave L=2 alone ------------
echo "===== stopping L=4 streams  [$(date +%F_%H:%M:%S)] ====="
ps -eo pid,etime,args | grep -E "$PAT" | grep -v grep
pkill -u "$USER" -f "$PAT" 2>/dev/null
for i in 1 2 3 4 5
do
  pgrep -u "$USER" -f "$PAT" >/dev/null || break
  sleep 2
done
if pgrep -u "$USER" -f "$PAT" >/dev/null
then
  echo "WARNING: L4 still alive after SIGTERM, sending SIGKILL"
  pkill -9 -u "$USER" -f "$PAT" 2>/dev/null
  sleep 2
fi
echo "L4 remaining: $(pgrep -u "$USER" -f "$PAT" | tr '\n' ' ')(none if blank)"
echo "L2 still running (should be 2): $(pgrep -u "$USER" -f "hmc_L2_claude.o" | tr '\n' ' ')"

# ---- 2. MPS sanity --------------------------------------------------------
# ---- ensure the CUDA MPS daemon is up BEFORE launching (co-resident procs SM-share,
# ---- not time-slice). Clients must start AFTER the daemon.
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

# ---- 3. rebuild (source changed -> will recompile) ------------------------
echo "===== build hmc_L4_claude.o  [$(date +%F_%H:%M:%S)] ====="
make hmc_L4_claude.o 2>&1 | tee build_hmc_L4_claude.log
test "${PIPESTATUS[0]}" -eq 0 || { echo "L4 BUILD FAILED"; exit 1; }

# ---- 4. relaunch L4 Nf=4 and Nf=6 on GPU1 (co-resident under MPS) ---------
echo "===== launch L4 Nf=4 and Nf=6 on GPU1 (MPS)  [$(date +%F_%H:%M:%S)] ====="
for nf in 4 6
do
  CUDA_VISIBLE_DEVICES=1 ./hmc_L4_claude.o ${GSQ} ${nf} ${NU0} \
    2>&1 | tee run_L4_Nf${nf}_gsq${GSQ}_claude.log &
done
wait
echo "===== L4 streams finished  [$(date +%F_%H:%M:%S)] ====="
