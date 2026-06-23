#!/bin/bash
# Kill the L=8 job (GPU0) and resume the two L=2 streams on GPU0 under MPS.
# The L=4 streams on GPU1 are left running untouched.
#
# L2 resumes from its existing checkpoints (Nf2 ~k=1311, Nf4 ~k=658) with the
# already-built hmc_L2_claude.o (tmax=1.0, nsteps=5, k_ckpoint_rng=1000).
#
# Run detached, e.g.:
#   nohup bash kill_L8_restart_L2_claude.sh > kill_L8_restart_L2_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
cd "$SRCDIR" || exit 1
source ../../env.sh

GSQ=8.0
NU0=1.0

# ---- 1. kill the L=8 job (leave L=4 alone) --------------------------------
echo "===== killing L=8  [$(date +%F_%H:%M:%S)] ====="
ps -eo pid,etime,args | grep -E "hmc_L8_claude.o" | grep -v grep
pkill -u "$USER" -f "hmc_L8_claude.o" 2>/dev/null
for i in 1 2 3 4 5
do
  pgrep -u "$USER" -f "hmc_L8_claude.o" >/dev/null || break
  sleep 2
done
if pgrep -u "$USER" -f "hmc_L8_claude.o" >/dev/null
then
  echo "WARNING: L8 still alive after SIGTERM, sending SIGKILL"
  pkill -9 -u "$USER" -f "hmc_L8_claude.o" 2>/dev/null
  sleep 2
fi
echo "L8 remaining: $(pgrep -u "$USER" -f "hmc_L8_claude.o" | tr '\n' ' ')(none if blank)"
echo "L4 still running (should be 2): $(pgrep -u "$USER" -f "hmc_L4_claude.o" | tr '\n' ' ')"

# ---- 2. MPS sanity (two L2 procs co-resident on GPU0) ---------------------
pgrep -f nvidia-cuda-mps-control >/dev/null && echo "MPS daemon: running" \
  || echo "WARNING: MPS daemon NOT running -- start it:  nvidia-cuda-mps-control -d"

# ---- 3. build (up to date if unchanged) + relaunch L2 Nf2,Nf4 on GPU0 -----
echo "===== build hmc_L2_claude.o  [$(date +%F_%H:%M:%S)] ====="
make hmc_L2_claude.o 2>&1 | tee build_hmc_L2_claude.log
test "${PIPESTATUS[0]}" -eq 0 || { echo "L2 BUILD FAILED"; exit 1; }

echo "===== launch L2 Nf=2 and Nf=4 on GPU0 (MPS)  [$(date +%F_%H:%M:%S)] ====="
for nf in 2 4
do
  CUDA_VISIBLE_DEVICES=0 ./hmc_L2_claude.o ${GSQ} ${nf} ${NU0} \
    2>&1 | tee run_L2_Nf${nf}_gsq${GSQ}_claude.log &
done
wait
echo "===== L2 streams finished  [$(date +%F_%H:%M:%S)] ====="
