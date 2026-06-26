#!/bin/bash
# Stop ALL L2/L4 streams and restart everything with the settings now in the
# sources: rng interval L2=5 / L4=1 (fermilab); nsteps L2=6 / L4=8, tmax=1.0
# (L4=8 is fermilab; L2=6 is one finer than fermilab's 5 to curb the massless
# exceptional-config freeze that hit L2/Nf4 at nsteps=5).
# Packing unchanged: GPU0 = L2{Nf4, Nf6}, GPU1 = L4{Nf4, Nf6} (MPS).
#
# Streams resume from their latest checkpoints EXCEPT the frozen L2/Nf4, which must
# be rolled back to k=149 first (the trap is at k=799..847; see guard below).
#
# TWO-STEP because rollback uses rm (forbidden inside scripts):
#   1) run this script once -> it kills all streams, then ABORTS at the guard
#      because L2/Nf4 still has its frozen checkpoints (rng>149).
#   2) roll back L2/Nf4 yourself:
#        d=Nf4_gsq8.000000at0.200000nu01.000000nt128L2
#        for k in $(seq 150 1000); do rm -f "$d/ckpoint_lat.$k" "$d/ckpoint_rng.$k"; done
#      (preview first: bash prune_ckpoints_claude.sh "$d" 150 1000)
#   3) re-run this script -> guard passes, it rebuilds + relaunches all 4.
#
# Run detached, e.g.:
#   nohup bash restart_all_fermrng_L2L4_claude.sh > restart_all_fermrng_L2L4_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
cd "$SRCDIR" || exit 1
source ../../env.sh

GSQ=8.0
NU0=1.0
L2NF4=Nf4_gsq8.000000at0.200000nu01.000000nt128L2

# ---- 1. stop ALL L2/L4 streams --------------------------------------------
echo "===== stopping ALL L2/L4 streams  [$(date +%F_%H:%M:%S)] ====="
ps -eo pid,etime,args | grep -E "hmc_L[24]_claude.o" | grep -v grep
pkill -u "$USER" -f "hmc_L[24]_claude.o" 2>/dev/null
for i in 1 2 3 4 5
do
  pgrep -u "$USER" -f "hmc_L[24]_claude.o" >/dev/null || break
  sleep 2
done
if pgrep -u "$USER" -f "hmc_L[24]_claude.o" >/dev/null
then
  echo "WARNING: still alive after SIGTERM, sending SIGKILL"
  pkill -9 -u "$USER" -f "hmc_L[24]_claude.o" 2>/dev/null
  sleep 2
fi
echo "remaining: $(pgrep -u "$USER" -f "hmc_L[24]_claude.o" | tr '\n' ' ')(none if blank)"

# ---- 2. rollback guard for the frozen L2/Nf4 ------------------------------
hi_rng=$(ls "$L2NF4"/ckpoint_rng.* 2>/dev/null | sed 's/.*\.//' | sort -n | tail -1)
if [ -n "$hi_rng" ] && [ "$hi_rng" -gt 149 ]
then
  echo "ABORT: L2/Nf4 not rolled back yet (highest rng=$hi_rng > 149)."
  echo "Roll it back, then re-run this script:"
  echo "  d=$L2NF4"
  echo "  for k in \$(seq 150 1000); do rm -f \"\$d/ckpoint_lat.\$k\" \"\$d/ckpoint_rng.\$k\"; done"
  exit 1
fi
echo "L2/Nf4 rollback OK (highest rng=${hi_rng:-none} <= 149)"

# ---- 3. MPS + rebuild both binaries ---------------------------------------
pgrep -f nvidia-cuda-mps-control >/dev/null && echo "MPS daemon: running" \
  || echo "WARNING: MPS daemon NOT running -- start it:  nvidia-cuda-mps-control -d"

echo "===== rebuild hmc_L2/L4  [$(date +%F_%H:%M:%S)] ====="
make hmc_L2_claude.o 2>&1 | tee build_hmc_L2_claude.log
test "${PIPESTATUS[0]}" -eq 0 || { echo "L2 BUILD FAILED"; exit 1; }
make hmc_L4_claude.o 2>&1 | tee build_hmc_L4_claude.log
test "${PIPESTATUS[0]}" -eq 0 || { echo "L4 BUILD FAILED"; exit 1; }

# ---- 4. relaunch all 4 : GPU0 = L2{Nf4,Nf6}, GPU1 = L4{Nf4,Nf6} (MPS) -----
echo "===== relaunch all 4  [$(date +%F_%H:%M:%S)] ====="
for nf in 4 6
do
  CUDA_VISIBLE_DEVICES=0 ./hmc_L2_claude.o ${GSQ} ${nf} ${NU0} \
    2>&1 | tee run_L2_Nf${nf}_gsq${GSQ}_claude.log &
  CUDA_VISIBLE_DEVICES=1 ./hmc_L4_claude.o ${GSQ} ${nf} ${NU0} \
    2>&1 | tee run_L4_Nf${nf}_gsq${GSQ}_claude.log &
done
wait
echo "===== all L2/L4 streams finished  [$(date +%F_%H:%M:%S)] ====="
