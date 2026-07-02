#!/bin/bash
# Stop both L=4 streams, rebuild with the updated Zolotarev setup, and relaunch.
# Matches hmc_fermilab_wmass_L4_claude.cu (2026-06-26):
#   overlap  : Fermion D(DW, mass, 21, 0.001)  -- npole 11->21 + FIXED window k=0.001
#              (10x wider; adaptive re-fit removed in includes/overlap_wmass_claude.h)
#   nsteps   : Nf6=10, else 8   (tmax=1.0)
# L4 resumes from its latest checkpoints (not frozen) -- the more-accurate sign fn
# samples the same overlap distribution, so no re-thermalization needed.
#
# Run detached:
#   nohup bash restart_L4_zolo_claude.sh > restart_L4_zolo_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
cd "$SRCDIR" || exit 1
source ../../env.sh

GSQ=8.0
NU0=1.0
PAT="hmc_L4_claude.o"

# ---- 1. stop both L=4 streams ---------------------------------------------
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

# ---- 2. MPS + rebuild -----------------------------------------------------
pgrep -f nvidia-cuda-mps-control >/dev/null && echo "MPS daemon: running" \
  || echo "WARNING: MPS daemon NOT running -- start it:  nvidia-cuda-mps-control -d"

echo "===== rebuild hmc_L4_claude.o  [$(date +%F_%H:%M:%S)] ====="
make hmc_L4_claude.o 2>&1 | tee build_hmc_L4_claude.log
test "${PIPESTATUS[0]}" -eq 0 || { echo "L4 BUILD FAILED"; exit 1; }

# ---- 3. relaunch L4 Nf4 + Nf6 on GPU1 (MPS) -------------------------------
echo "===== launch L4 Nf=4 and Nf=6 on GPU1 (MPS)  [$(date +%F_%H:%M:%S)] ====="
for nf in 4 6
do
  CUDA_VISIBLE_DEVICES=1 ./hmc_L4_claude.o ${GSQ} ${nf} ${NU0} \
    2>&1 | tee run_L4_Nf${nf}_gsq${GSQ}_claude.log &
done
wait
echo "===== L4 streams finished  [$(date +%F_%H:%M:%S)] ====="
