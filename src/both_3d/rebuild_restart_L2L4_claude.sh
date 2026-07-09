#!/bin/bash
# Rebuild the L=2,4 massless HMC binaries (after the tmax / nsteps / k_ckpoint_rng
# edits) and relaunch the four streams under MPS.
#
# Safe to run even if the old processes are already stopped (the stop step is a
# no-op then).  Does NOT touch any checkpoint files.
#
# Run detached, e.g.:
#   nohup bash rebuild_restart_L2L4_claude.sh > rebuild_restart_L2L4_claude.log 2>&1 &
# or inside tmux/screen.
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
cd "$SRCDIR" || exit 1

PAT='hmc_L[24]_claude.o'

# ---- 1. stop any lingering L2/L4 HMC processes (no-op if already stopped) ----
echo "===== stopping any running ${PAT} processes  [$(date +%F_%H:%M:%S)] ====="
ps -eo pid,etime,args | grep -E "$PAT" | grep -v grep
pkill -u "$USER" -f "$PAT" 2>/dev/null

# wait up to ~10 s for a clean exit (SIGTERM); checkpoints are written per-traj
# so a graceful stop avoids truncating the latest checkpoint
for i in 1 2 3 4 5
do
  pgrep -u "$USER" -f "$PAT" >/dev/null || break
  sleep 2
done
if pgrep -u "$USER" -f "$PAT" >/dev/null
then
  echo "WARNING: still alive after SIGTERM, sending SIGKILL"
  pkill -9 -u "$USER" -f "$PAT" 2>/dev/null
  sleep 2
fi
echo "remaining matching pids: $(pgrep -u "$USER" -f "$PAT" | tr '\n' ' ')(none if blank)"

# ---- 2. ensure the CUDA MPS daemon is up (co-resident procs SM-share, not time-slice).
# ---- Clients must start AFTER the daemon; abort if it can't come up.
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

# ---- 3. rebuild (make detects the edited sources) + relaunch via wrapper -----
# run_nf_L2L4_massless_claude.sh rebuilds hmc_L{2,4}_claude.o and launches the
# four streams (GPU0=L2{Nf2,Nf4}, GPU1=L4{Nf2,Nf4}) under MPS, resume-safe.
echo "===== rebuild + relaunch via run_nf_L2L4_massless_claude.sh  [$(date +%F_%H:%M:%S)] ====="
exec bash run_nf_L2L4_massless_claude.sh
