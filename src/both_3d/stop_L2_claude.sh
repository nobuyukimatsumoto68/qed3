#!/bin/bash
# Stop BOTH L=2 streams (Nf4 and Nf6) entirely -- no relaunch. Frees GPU0.
# The L=4 streams on GPU1 are left running untouched.
#
# Intended to be run AFTER restart_L2Nf6_nsteps12_claude.sh has recovered L2/Nf6
# (so its dir ends at the clean rolled-back k=289 state, not frozen).
# Checkpoints/logs are left in place for analysis.
#
# Run:  bash stop_L2_claude.sh
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
cd "$SRCDIR" || exit 1

PAT="hmc_L2_claude.o"

echo "===== stopping ALL L=2 streams  [$(date +%F_%H:%M:%S)] ====="
ps -eo pid,etime,args | grep -E "$PAT" | grep -v grep
pkill -u "$USER" -f "$PAT" 2>/dev/null
for i in 1 2 3 4 5
do
  pgrep -u "$USER" -f "$PAT" >/dev/null || break
  sleep 2
done
if pgrep -u "$USER" -f "$PAT" >/dev/null
then
  echo "WARNING: L2 still alive after SIGTERM, sending SIGKILL"
  pkill -9 -u "$USER" -f "$PAT" 2>/dev/null
  sleep 2
fi
echo "L2 remaining: $(pgrep -u "$USER" -f "$PAT" | tr '\n' ' ')(none if blank)"
echo "L4 still running (should be 2): $(pgrep -u "$USER" -f "hmc_L4_claude.o" | tr '\n' ' ')"
echo "===== L2 stopped; GPU0 now free  [$(date +%F_%H:%M:%S)] ====="
