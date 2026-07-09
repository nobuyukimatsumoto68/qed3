#!/bin/bash
# Repackage: stop the Nf=2 streams (L=2 and L=4) and start Nf=6 (L=2 and L=4).
# The Nf=4 streams keep running untouched.  New packing under MPS:
#   GPU 0 = L=2 : Nf=4 (running) + Nf=6 (new)
#   GPU 1 = L=4 : Nf=4 (running) + Nf=6 (new)
#
# Nf=6 starts fresh (no prior Nf6 L2/L4 checkpoints -> thermalize from k=0) and
# inherits the per-L nsteps already in the sources (L2=5, L4=6; same as Nf=4).
#
# Run detached, e.g.:
#   nohup bash repackage_Nf6_L2L4_claude.sh > repackage_Nf6_L2L4_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
cd "$SRCDIR" || exit 1
source ../../env.sh

GSQ=8.0
NU0=1.0
NF2_PAT="hmc_L[24]_claude.o ${GSQ} 2 ${NU0}"

# ---- 1. stop ONLY the Nf=2 streams (match the exact arg string) ------------
echo "===== stopping Nf=2 (L2 + L4)  [$(date +%F_%H:%M:%S)] ====="
ps -eo pid,args | grep -E "$NF2_PAT" | grep -v grep
pkill -u "$USER" -f "$NF2_PAT" 2>/dev/null
for i in 1 2 3 4 5
do
  pgrep -u "$USER" -f "$NF2_PAT" >/dev/null || break
  sleep 2
done
if pgrep -u "$USER" -f "$NF2_PAT" >/dev/null
then
  echo "WARNING: Nf=2 still alive after SIGTERM, sending SIGKILL"
  pkill -9 -u "$USER" -f "$NF2_PAT" 2>/dev/null
  sleep 2
fi
echo "Nf=2 remaining: $(pgrep -u "$USER" -f "$NF2_PAT" | tr '\n' ' ')(none if blank)"
echo "Nf=4 still running (should be 2): $(pgrep -u "$USER" -f "hmc_L[24]_claude.o ${GSQ} 4 ${NU0}" | tr '\n' ' ')"

# ---- 2. MPS sanity (two procs co-resident per GPU) ------------------------
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

# ---- 3. build (up to date if unchanged) -----------------------------------
echo "===== build hmc_L2/L4  [$(date +%F_%H:%M:%S)] ====="
make hmc_L2_claude.o 2>&1 | tee build_hmc_L2_claude.log
test "${PIPESTATUS[0]}" -eq 0 || { echo "L2 BUILD FAILED"; exit 1; }
make hmc_L4_claude.o 2>&1 | tee build_hmc_L4_claude.log
test "${PIPESTATUS[0]}" -eq 0 || { echo "L4 BUILD FAILED"; exit 1; }

# ---- 4. launch Nf=6 : L2 on GPU0, L4 on GPU1 (co-resident with Nf=4) -------
echo "===== launch Nf=6 : L2 GPU0, L4 GPU1  [$(date +%F_%H:%M:%S)] ====="
CUDA_VISIBLE_DEVICES=0 ./hmc_L2_claude.o ${GSQ} 6 ${NU0} \
  2>&1 | tee run_L2_Nf6_gsq${GSQ}_claude.log &
CUDA_VISIBLE_DEVICES=1 ./hmc_L4_claude.o ${GSQ} 6 ${NU0} \
  2>&1 | tee run_L4_Nf6_gsq${GSQ}_claude.log &
wait
echo "===== Nf=6 streams finished  [$(date +%F_%H:%M:%S)] ====="
