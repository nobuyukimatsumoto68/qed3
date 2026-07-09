#!/bin/bash
# Resume the paused L=1 gsq=4.0 scan: Nf=2 and Nf=4 packed on GPU1 under MPS (2 jobs/GPU).
# Nf=6 is the deferred group (run separately later -- see run_L1_gsq4_resume_Nf6_claude.sh).
# Uses the current hmc_L1_claude.o (npole=21/window0.001, nsteps=10 all Nf, kmax=320).
# Each resumes from its last checkpoint: Nf2 k=252, Nf4 k=135.
#
# Run detached:
#   nohup bash run_L1_gsq4_resume_Nf24_claude.sh > run_L1_gsq4_resume_Nf24_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
cd "$SRCDIR" || exit 1
source ../../env.sh

GSQ=4.0
NU0=1.0

pgrep -f nvidia-cuda-mps-control >/dev/null && echo "MPS daemon: running" \
  || echo "WARNING: MPS daemon NOT running -- start it:  nvidia-cuda-mps-control -d"

echo "===== build hmc_L1_claude.o  [$(date +%F_%H:%M:%S)] ====="
make hmc_L1_claude.o 2>&1 | tee build_hmc_L1_claude.log
test "${PIPESTATUS[0]}" -eq 0 || { echo "L1 BUILD FAILED"; exit 1; }

echo "===== resume gsq=4.0 Nf=2,4 on GPU1 (MPS, 2 jobs)  [$(date +%F_%H:%M:%S)] ====="
for nf in 2 4
do
  CUDA_VISIBLE_DEVICES=1 ./hmc_L1_claude.o ${GSQ} ${nf} ${NU0} \
    2>&1 | tee run_L1_gsq4_Nf${nf}_claude.log &
done
wait
echo "===== gsq=4.0 Nf=2,4 finished  [$(date +%F_%H:%M:%S)] ====="
