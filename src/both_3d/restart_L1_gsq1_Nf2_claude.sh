#!/bin/bash
# Recover the frozen L1 gsq=1.0 Nf2 (trap at k=113): roll back to k=109 and rerun
# to kmax=320 at nsteps=10 (step 0.19; was 8/step 0.2375). npole=21 unchanged.
# Nf4/Nf6 are already complete+clean -- only Nf2 reruns.
#
# TWO-STEP (rollback uses rm, forbidden inside scripts):
#   1) run this once -> kills any L1/Nf2 proc, then ABORTS at the guard (rng>109).
#   2) roll back Nf2 yourself:
#        d=Nf2_gsq1.000000at0.200000nu01.000000nt128L1
#        for k in $(seq 110 320); do rm -f "$d/ckpoint_lat.$k" "$d/ckpoint_rng.$k"; done
#      (preview: bash prune_ckpoints_claude.sh "$d" 110 320)
#   3) re-run this -> guard passes, rebuilds hmc_L1_claude.o, reruns Nf2 on GPU1.
#
# Run detached:
#   nohup bash restart_L1_gsq1_Nf2_claude.sh > restart_L1_gsq1_Nf2_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
cd "$SRCDIR" || exit 1
source ../../env.sh

GSQ=1.0
NU0=1.0
PAT="hmc_L1_claude.o ${GSQ} 2 ${NU0}"
NF2DIR=Nf2_gsq1.000000at0.200000nu01.000000nt128L1

# ---- 1. stop any L1/Nf2 proc (likely already exited) ----------------------
echo "===== stopping L1/Nf2  [$(date +%F_%H:%M:%S)] ====="
pkill -u "$USER" -f "$PAT" 2>/dev/null
for i in 1 2 3 4 5
do
  pgrep -u "$USER" -f "$PAT" >/dev/null || break
  sleep 2
done
pgrep -u "$USER" -f "$PAT" >/dev/null && { echo "SIGKILL"; pkill -9 -u "$USER" -f "$PAT" 2>/dev/null; sleep 2; }
echo "L1/Nf2 remaining: $(pgrep -u "$USER" -f "$PAT" | tr '\n' ' ')(none if blank)"

# ---- 2. rollback guard (k=109) --------------------------------------------
hi_rng=$(ls "$NF2DIR"/ckpoint_rng.* 2>/dev/null | sed 's/.*\.//' | sort -n | tail -1)
if [ -n "$hi_rng" ] && [ "$hi_rng" -gt 109 ]
then
  echo "ABORT: Nf2 not rolled back yet (highest rng=$hi_rng > 109)."
  echo "  d=$NF2DIR"
  echo "  for k in \$(seq 110 320); do rm -f \"\$d/ckpoint_lat.\$k\" \"\$d/ckpoint_rng.\$k\"; done"
  exit 1
fi
echo "Nf2 rollback OK (highest rng=${hi_rng:-none} <= 109)"

# ---- 3. rebuild (Nf2 now nsteps=10) ---------------------------------------
echo "===== rebuild hmc_L1_claude.o  [$(date +%F_%H:%M:%S)] ====="
make hmc_L1_claude.o 2>&1 | tee build_hmc_L1_claude.log
test "${PIPESTATUS[0]}" -eq 0 || { echo "L1 BUILD FAILED"; exit 1; }

# ---- 4. rerun Nf2 at gsq=1.0 on GPU1 --------------------------------------
echo "===== rerun L1 gsq=1.0 Nf2 on GPU1  [$(date +%F_%H:%M:%S)] ====="
CUDA_VISIBLE_DEVICES=1 ./hmc_L1_claude.o ${GSQ} 2 ${NU0} \
  2>&1 | tee run_L1_gsq1_Nf2_claude.log &
wait $!
echo "===== L1 gsq=1.0 Nf2 finished  [$(date +%F_%H:%M:%S)] ====="
