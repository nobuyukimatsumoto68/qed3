#!/bin/bash
# L=2,4 massless overlap HMC SEA runs on the local 2x TITAN V box.
#
# Packing under CUDA MPS (two procs co-resident per GPU, SM-shared ~1.6-1.9x;
# each binary runs kmax=4000 trajectories then exits, resume-safe via
# ckpoint_lat/rng):
#   GPU 0 : hmc_L2_claude.o  ->  Nf=2  AND  Nf=4   (concurrent under MPS)
#   GPU 1 : hmc_L4_claude.o  ->  Nf=2  AND  Nf=4   (concurrent under MPS)
# Nf=6 is deferred (run later by adding it to the NF_LIST below).
#
# MPS must be up first:  nvidia-cuda-mps-control -d
# (12 GB TITAN V easily holds two of these small L2/L4 HMC procs.)
#
# Physics params match the L=1 production set: gsq=8.0, Nt=128, at=0.2, nu0=1.0,
# pole count n=11 (all baked into the sources; only N_REFINE differs).
#
# Run from this directory (the binary uses relative paths ../../geometry/data/).
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
cd "$SRCDIR" || exit 1
source ../../env.sh

GSQ=8.0
NU0=1.0
NF_LIST="2 4"

# ---- build both binaries --------------------------------------------------
echo "===== build hmc_L2_claude.o ====="
date
make hmc_L2_claude.o 2>&1 | tee build_hmc_L2_claude.log
test "${PIPESTATUS[0]}" -eq 0 || { echo "L2 BUILD FAILED"; exit 1; }

echo "===== build hmc_L4_claude.o ====="
date
make hmc_L4_claude.o 2>&1 | tee build_hmc_L4_claude.log
test "${PIPESTATUS[0]}" -eq 0 || { echo "L4 BUILD FAILED"; exit 1; }

# ---- launch one HMC per (GPU,Nf); two co-resident per GPU under MPS --------
run_one() {
  local gpu=$1
  local app=$2
  local tag=$3
  local nf=$4
  echo "### [GPU ${gpu}] ${tag} Nf=${nf} start  [$(date +%F_%H:%M:%S)] ###"
  CUDA_VISIBLE_DEVICES=${gpu} ./${app} ${GSQ} ${nf} ${NU0} \
    2>&1 | tee run_${tag}_Nf${nf}_gsq${GSQ}_claude.log
  echo "### [GPU ${gpu}] ${tag} Nf=${nf} done (status ${PIPESTATUS[0]})  [$(date +%F_%H:%M:%S)] ###"
}

echo "===== launching L2 (GPU0) + L4 (GPU1), two procs/GPU under MPS ====="
date
for nf in $NF_LIST
do
  run_one 0 hmc_L2_claude.o L2 ${nf} &
  run_one 1 hmc_L4_claude.o L4 ${nf} &
done
wait
echo "===== all L2/L4 massless runs finished ====="
date
