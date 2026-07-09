#!/bin/bash
# L=1 massless overlap HMC, gsq=1.0, Nf=2/4/6 -- three jobs packed on GPU1 under MPS.
# kmax=320 cap; Zolotarev npole=21 + FIXED window k=0.001 (matches the L4 fix);
# other settings = L1 production (tmax=1.9, nsteps Nf2=8/Nf4=10/Nf6=10, at=0.2, m=0).
# Fresh dirs Nf{2,4,6}_gsq1.000000at0.200000nu01.000000nt128L1/.
#
# GPU1 must be free first -- kill L4 yourself:  pkill -u "$USER" -f "hmc_L4_claude.o"
# (GPU0 is busy with the ylm L2 sweep.)
#
# Run detached:
#   nohup bash run_L1_gsq1_claude.sh > run_L1_gsq1_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
cd "$SRCDIR" || exit 1
source ../../env.sh

GSQ=1.0
NU0=1.0

# ---- guard: GPU1 should be free (L4 stopped) ------------------------------
if pgrep -u "$USER" -f "hmc_L4_claude.o" >/dev/null
then
  echo "WARNING: hmc_L4_claude.o still running on GPU1 -- kill it first:"
  echo "  pkill -u \"$USER\" -f \"hmc_L4_claude.o\""
  echo "(continuing anyway would oversubscribe GPU1)"
  exit 1
fi
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

# ---- build ----------------------------------------------------------------
echo "===== build hmc_L1_claude.o  [$(date +%F_%H:%M:%S)] ====="
make hmc_L1_claude.o 2>&1 | tee build_hmc_L1_claude.log
test "${PIPESTATUS[0]}" -eq 0 || { echo "L1 BUILD FAILED"; exit 1; }

# ---- launch Nf=2,4,6 at gsq=1.0 on GPU1 (3 procs co-resident, MPS) --------
echo "===== launch L1 gsq=1.0 Nf=2,4,6 on GPU1 (MPS)  [$(date +%F_%H:%M:%S)] ====="
for nf in 2 4 6
do
  CUDA_VISIBLE_DEVICES=1 ./hmc_L1_claude.o ${GSQ} ${nf} ${NU0} \
    2>&1 | tee run_L1_gsq1_Nf${nf}_claude.log &
done
wait
echo "===== L1 gsq=1.0 runs finished  [$(date +%F_%H:%M:%S)] ====="
