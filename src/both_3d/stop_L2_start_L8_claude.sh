#!/bin/bash
# Stop the two L=2 streams (GPU0), build the L=8 binary, and launch L=8 Nf=2 on
# GPU0.  The L=4 streams on GPU1 are left running untouched.
#
# L8 massless run, numbers mimic hmc_fermilab_wmass_L8_claude.cu:
#   N_REFINE=8, Nt=128, at=0.1, M5=-1.0, npole=13, tmax=1.0, nsteps=16,
#   kmax=80, k_ckpoint=1, k_ckpoint_rng=2.   Only Nf=2.
#
# Run detached, e.g.:
#   nohup bash stop_L2_start_L8_claude.sh > stop_L2_start_L8_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
cd "$SRCDIR" || exit 1
source ../../env.sh

GSQ=8.0
NU0=1.0

# ---- 1. stop ONLY the L=2 streams (leave L=4 alone) -----------------------
echo "===== stopping L=2 streams  [$(date +%F_%H:%M:%S)] ====="
ps -eo pid,etime,args | grep -E "hmc_L2_claude.o" | grep -v grep
pkill -u "$USER" -f "hmc_L2_claude.o" 2>/dev/null
for i in 1 2 3 4 5
do
  pgrep -u "$USER" -f "hmc_L2_claude.o" >/dev/null || break
  sleep 2
done
if pgrep -u "$USER" -f "hmc_L2_claude.o" >/dev/null
then
  echo "WARNING: L2 still alive after SIGTERM, sending SIGKILL"
  pkill -9 -u "$USER" -f "hmc_L2_claude.o" 2>/dev/null
  sleep 2
fi
echo "L2 remaining: $(pgrep -u "$USER" -f "hmc_L2_claude.o" | tr '\n' ' ')(none if blank)"
echo "L4 still running (should be 2): $(pgrep -u "$USER" -f "hmc_L4_claude.o" | tr '\n' ' ')"

# ---- 2. build the L=8 binary ----------------------------------------------
echo "===== build hmc_L8_claude.o  [$(date +%F_%H:%M:%S)] ====="
make hmc_L8_claude.o 2>&1 | tee build_hmc_L8_claude.log
test "${PIPESTATUS[0]}" -eq 0 || { echo "L8 BUILD FAILED"; exit 1; }

# ---- 3. launch L=8 Nf=2 on GPU0 (single process; no MPS packing needed) ---
echo "===== launch L8 Nf=2 on GPU0  [$(date +%F_%H:%M:%S)] ====="
CUDA_VISIBLE_DEVICES=0 ./hmc_L8_claude.o ${GSQ} 2 ${NU0} \
  2>&1 | tee run_L8_Nf2_gsq${GSQ}_claude.log &
wait $!
echo "===== L8 Nf=2 finished  [$(date +%F_%H:%M:%S)] ====="
