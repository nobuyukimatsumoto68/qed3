#!/bin/bash
# Resume the deferred L=1 gsq=4.0 Nf=6 (single, long-pole ~183h to 320). Resumes from k=93.
# GPU defaults to 1; pass a different index as arg1 (e.g. 0 once GPU0 frees).
#   bash run_L1_gsq4_resume_Nf6_claude.sh        # -> GPU1
#   bash run_L1_gsq4_resume_Nf6_claude.sh 0      # -> GPU0
# Uses current hmc_L1_claude.o (npole=21/window0.001, nsteps=10, kmax=320).
#
# Run detached:
#   nohup bash run_L1_gsq4_resume_Nf6_claude.sh > run_L1_gsq4_resume_Nf6_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
cd "$SRCDIR" || exit 1
source ../../env.sh

GSQ=4.0
NU0=1.0
GPU=${1:-1}

pgrep -f nvidia-cuda-mps-control >/dev/null && echo "MPS daemon: running" \
  || echo "WARNING: MPS daemon NOT running -- start it:  nvidia-cuda-mps-control -d"

echo "===== build hmc_L1_claude.o  [$(date +%F_%H:%M:%S)] ====="
make hmc_L1_claude.o 2>&1 | tee build_hmc_L1_claude.log
test "${PIPESTATUS[0]}" -eq 0 || { echo "L1 BUILD FAILED"; exit 1; }

echo "===== resume gsq=4.0 Nf=6 on GPU${GPU}  [$(date +%F_%H:%M:%S)] ====="
CUDA_VISIBLE_DEVICES=${GPU} ./hmc_L1_claude.o ${GSQ} 6 ${NU0} \
  2>&1 | tee run_L1_gsq4_Nf6_claude.log &
wait $!
echo "===== gsq=4.0 Nf=6 finished  [$(date +%F_%H:%M:%S)] ====="
