#!/bin/bash
# L=1 massless overlap HMC, gsq=12.0, Nf=2/4/6 -- three jobs packed on GPU1 under MPS.
# Same recipe as run_L1_gsq{1,4}_claude.sh: kmax=320, Zolotarev npole=21 + window k=0.001,
# nsteps=10 all Nf (step 0.19), tmax=1.9, at=0.2, m=0.
# Fresh dirs Nf{2,4,6}_gsq12.000000at0.200000nu01.000000nt128L1/ (the old data_-prefixed
# gsq12 dir is a different name -- no collision).
#
# GPU1 is free (gsq=1.0 Nf2 recovery finished). GPU0 is running the gsq=4.0 3-pack.
#
# Run detached:
#   nohup bash run_L1_gsq12_claude.sh > run_L1_gsq12_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
cd "$SRCDIR" || exit 1
source ../../env.sh

GSQ=12.0
NU0=1.0

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

# ---- build (up to date if unchanged) --------------------------------------
echo "===== build hmc_L1_claude.o  [$(date +%F_%H:%M:%S)] ====="
make hmc_L1_claude.o 2>&1 | tee build_hmc_L1_claude.log
test "${PIPESTATUS[0]}" -eq 0 || { echo "L1 BUILD FAILED"; exit 1; }

# ---- launch Nf=2,4,6 at gsq=12.0 on GPU1 (3 procs co-resident, MPS) -------
echo "===== launch L1 gsq=12.0 Nf=2,4,6 on GPU1 (MPS)  [$(date +%F_%H:%M:%S)] ====="
for nf in 2 4 6
do
  CUDA_VISIBLE_DEVICES=1 ./hmc_L1_claude.o ${GSQ} ${nf} ${NU0} \
    2>&1 | tee run_L1_gsq12_Nf${nf}_claude.log &
done
wait
echo "===== L1 gsq=12.0 runs finished  [$(date +%F_%H:%M:%S)] ====="
