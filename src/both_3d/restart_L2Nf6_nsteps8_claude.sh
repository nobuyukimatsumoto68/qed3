#!/bin/bash
# Recover the frozen L2/Nf6 (trap at k=155): roll back to k=154 and restart with
# the finer step (nsteps=8, step 0.125; L2 nsteps is now Nf-dependent, Nf4 stays 6).
# Only L2/Nf6 is touched -- L2/Nf4 (GPU0) and both L4 streams (GPU1) keep running.
#
# TWO-STEP (rollback uses rm, forbidden inside scripts):
#   1) run this once -> kills L2/Nf6, then ABORTS at the guard (rng>154 still present).
#   2) roll back L2/Nf6 yourself:
#        d=Nf6_gsq8.000000at0.200000nu01.000000nt128L2
#        for k in $(seq 155 300); do rm -f "$d/ckpoint_lat.$k" "$d/ckpoint_rng.$k"; done
#      (preview: bash prune_ckpoints_claude.sh "$d" 155 300)
#   3) re-run this -> guard passes, rebuilds hmc_L2_claude.o, relaunches L2/Nf6 on GPU0.
#
# Run detached:
#   nohup bash restart_L2Nf6_nsteps8_claude.sh > restart_L2Nf6_nsteps8_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
cd "$SRCDIR" || exit 1
source ../../env.sh

GSQ=8.0
NU0=1.0
PAT="hmc_L2_claude.o ${GSQ} 6 ${NU0}"
L2NF6=Nf6_gsq8.000000at0.200000nu01.000000nt128L2

# ---- 1. stop ONLY L2/Nf6 --------------------------------------------------
echo "===== stopping L2/Nf6  [$(date +%F_%H:%M:%S)] ====="
ps -eo pid,etime,args | grep -E "$PAT" | grep -v grep
pkill -u "$USER" -f "$PAT" 2>/dev/null
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
echo "L2/Nf6 remaining: $(pgrep -u "$USER" -f "$PAT" | tr '\n' ' ')(none if blank)"
echo "still running (L2/Nf4 + both L4 = 3): $(pgrep -u "$USER" -f "hmc_L[24]_claude.o ${GSQ} 4 ${NU0}|hmc_L4_claude.o ${GSQ} 6 ${NU0}" | tr '\n' ' ')"

# ---- 2. rollback guard ----------------------------------------------------
hi_rng=$(ls "$L2NF6"/ckpoint_rng.* 2>/dev/null | sed 's/.*\.//' | sort -n | tail -1)
if [ -n "$hi_rng" ] && [ "$hi_rng" -gt 154 ]
then
  echo "ABORT: L2/Nf6 not rolled back yet (highest rng=$hi_rng > 154)."
  echo "  d=$L2NF6"
  echo "  for k in \$(seq 155 300); do rm -f \"\$d/ckpoint_lat.\$k\" \"\$d/ckpoint_rng.\$k\"; done"
  exit 1
fi
echo "L2/Nf6 rollback OK (highest rng=${hi_rng:-none} <= 154)"

# ---- 3. rebuild hmc_L2_claude.o (Nf-dependent nsteps now) ------------------
echo "===== rebuild hmc_L2_claude.o  [$(date +%F_%H:%M:%S)] ====="
make hmc_L2_claude.o 2>&1 | tee build_hmc_L2_claude.log
test "${PIPESTATUS[0]}" -eq 0 || { echo "L2 BUILD FAILED"; exit 1; }

# ---- 4. relaunch L2/Nf6 on GPU0 (co-resident with L2/Nf4 under MPS) --------
echo "===== relaunch L2/Nf6 on GPU0  [$(date +%F_%H:%M:%S)] ====="
CUDA_VISIBLE_DEVICES=0 ./hmc_L2_claude.o ${GSQ} 6 ${NU0} \
  2>&1 | tee run_L2_Nf6_gsq${GSQ}_claude.log &
wait $!
echo "===== L2/Nf6 finished  [$(date +%F_%H:%M:%S)] ====="
