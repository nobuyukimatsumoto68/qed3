#!/bin/bash
# L=1 massless overlap HMC, gsq=4.0, Nf=2/4/6 -- three jobs packed on GPU0 under MPS.
# Same recipe as run_L1_gsq1_claude.sh: kmax=320, Zolotarev npole=21 + window k=0.001,
# nsteps=10 all Nf (step 0.19), tmax=1.9, at=0.2, m=0.
# Fresh dirs Nf{2,4,6}_gsq4.000000at0.200000nu01.000000nt128L1/.
#
# GPU0 should be free (the ylm L2 sweep is done). GPU1 is running the gsq=1.0 Nf2 recovery.
#
# Run detached:
#   nohup bash run_L1_gsq4_claude.sh > run_L1_gsq4_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
cd "$SRCDIR" || exit 1
source ../../env.sh

GSQ=4.0
NU0=1.0

pgrep -f nvidia-cuda-mps-control >/dev/null && echo "MPS daemon: running" \
  || echo "WARNING: MPS daemon NOT running -- start it:  nvidia-cuda-mps-control -d"

# ---- build (up to date if unchanged) --------------------------------------
echo "===== build hmc_L1_claude.o  [$(date +%F_%H:%M:%S)] ====="
make hmc_L1_claude.o 2>&1 | tee build_hmc_L1_claude.log
test "${PIPESTATUS[0]}" -eq 0 || { echo "L1 BUILD FAILED"; exit 1; }

# ---- launch Nf=2,4,6 at gsq=4.0 on GPU0 (3 procs co-resident, MPS) --------
echo "===== launch L1 gsq=4.0 Nf=2,4,6 on GPU0 (MPS)  [$(date +%F_%H:%M:%S)] ====="
for nf in 2 4 6
do
  CUDA_VISIBLE_DEVICES=0 ./hmc_L1_claude.o ${GSQ} ${nf} ${NU0} \
    2>&1 | tee run_L1_gsq4_Nf${nf}_claude.log &
done
wait
echo "===== L1 gsq=4.0 runs finished  [$(date +%F_%H:%M:%S)] ====="
