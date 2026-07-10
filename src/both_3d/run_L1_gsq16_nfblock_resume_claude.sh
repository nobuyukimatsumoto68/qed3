#!/bin/bash
# Resume the L=1 gsq=16.0 production with the FASTER Nf-block (block-force) HMC code
# (hmc_Nfblocked_claude.cu). GPU1, 2 procs/GPU under MPS, grouping {Nf2->Nf4}+{Nf6}.
#   Nf2 -> serial hmc_L1_claude.o          (NSTACK=1, block is a no-op for Nf=2)
#   Nf4 -> hmc_block_L1_Nf4.out  (-DNFPF=4 -DBLOCK_FORCE)
#   Nf6 -> hmc_block_L1_Nf6.out  (-DNFPF=6 -DBLOCK_FORCE)
# Block code is dH-validated (serial==block to ~1e-9) -> valid to resume serial-generated dirs.
# NOTE: at L1 the block gain is only ~1.1x (force is a minority of the L1 trajectory); the big
# win (2.5x) is at L4/Nf6. Physics identical: Zolotarev n=21/k=0.001, L1 tmax=1.9/nsteps=10, kmax=320.
# Resume points: Nf2 k=170, Nf6 k=69, Nf4 fresh (k=0). Dirs Nf{2,4,6}_gsq16...L1 (-DDIR_NO_MASS).
#
# Run detached:
#   nohup bash run_L1_gsq16_nfblock_resume_claude.sh > run_L1_gsq16_nfblock_resume_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
cd "$SRCDIR" || exit 1
source ../../env.sh

GSQ=16.0
NU0=1.0
GPU=1

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
SRC=hmc_Nfblocked_claude.cu
BLKDEF="-DLREFINE=1 -DBLOCK_FORCE -DDIR_NO_MASS -DKMAX=320 -DKRNG=10 -DNPGAUGE=1 -DNPSORT=1"

# ---- 1. ensure MPS daemon is up (ambient default pipe -- do NOT start a custom pipe dir) ----
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

# ---- 2. build: serial Nf2 binary + the two block-force binaries -------------
echo "===== build hmc_L1_claude.o (Nf2 serial)  [$(date +%F_%H:%M:%S)] ====="
make hmc_L1_claude.o 2>&1 | tee build_hmc_L1_claude.log
test "${PIPESTATUS[0]}" -eq 0 || { echo "L1 serial BUILD FAILED"; exit 1; }

echo "===== build hmc_block_L1_Nf4.out  [$(date +%F_%H:%M:%S)] ====="
$NVCC $NVCCFLAGS $INCLUDES $BLKDEF -DNFPF=4 "$SRC" -o hmc_block_L1_Nf4.out 2>&1 | tee build_block_L1_Nf4_claude.log
test -f hmc_block_L1_Nf4.out || { echo "L1/Nf4 block BUILD FAILED"; exit 1; }

echo "===== build hmc_block_L1_Nf6.out  [$(date +%F_%H:%M:%S)] ====="
$NVCC $NVCCFLAGS $INCLUDES $BLKDEF -DNFPF=6 "$SRC" -o hmc_block_L1_Nf6.out 2>&1 | tee build_block_L1_Nf6_claude.log
test -f hmc_block_L1_Nf6.out || { echo "L1/Nf6 block BUILD FAILED"; exit 1; }

# ---- 3. launch on GPU1: slot1={Nf2 serial -> Nf4 block}, slot2={Nf6 block} --
echo "===== launch gsq=16 (Nf-block) on GPU${GPU}  [$(date +%F_%H:%M:%S)] ====="

(
  echo "### [slot1] gsq16 Nf2 (serial) start [$(date +%F_%H:%M:%S)] ###"
  CUDA_VISIBLE_DEVICES=${GPU} ./hmc_L1_claude.o ${GSQ} 2 ${NU0} \
    2>&1 | tee run_L1_gsq16_Nf2_claude.log
  echo "### [slot1] gsq16 Nf2 done -> starting Nf4 (block) [$(date +%F_%H:%M:%S)] ###"
  CUDA_VISIBLE_DEVICES=${GPU} ./hmc_block_L1_Nf4.out ${GSQ} 4 ${NU0} \
    2>&1 | tee run_L1_gsq16_Nf4_claude.log
  echo "### [slot1] gsq16 Nf4 done [$(date +%F_%H:%M:%S)] ###"
) &
SLOT1=$!

CUDA_VISIBLE_DEVICES=${GPU} ./hmc_block_L1_Nf6.out ${GSQ} 6 ${NU0} \
  2>&1 | tee run_L1_gsq16_Nf6_claude.log &
SLOT2=$!

wait $SLOT1
wait $SLOT2
echo "===== gsq=16 (Nf-block) finished  [$(date +%F_%H:%M:%S)] ====="
