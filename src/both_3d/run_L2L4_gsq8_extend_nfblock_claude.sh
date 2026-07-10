#!/bin/bash
# Extend the LOCAL massless gsq8 SEA ensembles with the faster Nf-block HMC code.
#   L2 -> kmax 2000  (Nf2 k1601, Nf4 k518, Nf6 k289)
#   L4 -> kmax 600   (Nf2 k168,  Nf4 k190, Nf6 k68)
# BOTH GPUs, 2 procs/GPU under MPS:
#   GPU0 = L4 : slot {Nf6}  ||  slot {Nf4 -> Nf2}
#   GPU1 = L2 : slot {Nf6}  ||  slot {Nf4 -> Nf2}
# (Nf6 heaviest gets its own slot; Nf4->Nf2 share the other -> the two slots roughly balance.
#  L4/Nf6 is the ~2.5x block-force win. 4 procs total, 2/GPU.)
#
# Binaries: hmc_Nfblocked_claude.cu, per (L,Nf) via -D. Nf4/Nf6 = BLOCK_FORCE (the win);
# Nf2 = NSTACK=1 (no block, effectively serial) but built from the same source so -DKMAX applies.
# Physics = current standard (Zolotarev n=21/k=0.001, L2/L4 tmax=1.0/nsteps=10 [-DNSTEPS=10;
# bumped from 8 after the L2/Nf4 k=541 near-zero-mode freeze], at=0.2, massless);
# dir naming byte-matches (-DDIR_NO_MASS) so it AUTO-RESUMES the existing Nf{2,4,6}_gsq8...L{2,4} dirs.
# NPGAUGE=NPSORT=1 (local: 1 host thread/proc, no CPU oversubscription when packing).
#
# Run detached:
#   nohup bash run_L2L4_gsq8_extend_nfblock_claude.sh > run_L2L4_gsq8_extend_nfblock_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/both_3d
cd "$SRCDIR" || exit 1
source ../../env.sh

GSQ=8.0
NU0=1.0
KRNG=5

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
SRC=hmc_Nfblocked_claude.cu

# ---- 1. ensure MPS daemon up (ambient default pipe; do NOT start a custom pipe dir) --------
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
  || { echo "ERROR: MPS daemon failed to start -- aborting"; exit 1; }

# ---- 2. build the six binaries (Nf4/Nf6 block-force; Nf2 NSTACK=1 no-block) ----------------
build () {
  local L=$1 nf=$2 kmax=$3 bf=$4
  local out=hmc_block_L${L}_Nf${nf}.out
  echo "===== build $out (LREFINE=$L NFPF=$nf KMAX=$kmax $bf)  [$(date +%F_%H:%M:%S)] ====="
  $NVCC $NVCCFLAGS $INCLUDES \
    -DLREFINE=$L -DNFPF=$nf -DDIR_NO_MASS -DKMAX=$kmax -DKRNG=$KRNG -DNSTEPS=10 -DNPGAUGE=1 -DNPSORT=1 $bf \
    "$SRC" -o "$out" 2>&1 | tee build_${out%.out}_claude.log
  test -f "$out" || { echo "BUILD FAILED: $out"; exit 1; }
}

build 2 6 2000 "-DBLOCK_FORCE"
build 2 4 2000 "-DBLOCK_FORCE"
build 2 2 2000 ""
build 4 6 600  "-DBLOCK_FORCE"
build 4 4 600  "-DBLOCK_FORCE"
build 4 2 600  ""

# ---- 3. launch: GPU0 = L4 (2 slots), GPU1 = L2 (2 slots) -> 2 procs/GPU (MPS) --------------
echo "===== launch: GPU0={L4 Nf6 || Nf4->Nf2}, GPU1={L2 Nf6 || Nf4->Nf2}  [$(date +%F_%H:%M:%S)] ====="

run_slot () {
  local gpu=$1 L=$2; shift 2
  for nf in "$@"
  do
    echo "### [gpu$gpu L$L] Nf${nf} start [$(date +%F_%H:%M:%S)] ###"
    CUDA_VISIBLE_DEVICES=${gpu} ./hmc_block_L${L}_Nf${nf}.out ${GSQ} ${nf} ${NU0} \
      2>&1 | tee run_L${L}_gsq8_Nf${nf}_nfblock_claude.log
    echo "### [gpu$gpu L$L] Nf${nf} done (status ${PIPESTATUS[0]}) [$(date +%F_%H:%M:%S)] ###"
  done
}

run_slot 0 4 6   &     # GPU0 slot A: L4 Nf6
run_slot 0 4 4 2 &     # GPU0 slot B: L4 Nf4 -> Nf2
run_slot 1 2 6   &     # GPU1 slot A: L2 Nf6
run_slot 1 2 4 2 &     # GPU1 slot B: L2 Nf4 -> Nf2
wait
echo "===== L2/L4 gsq8 extension finished  [$(date +%F_%H:%M:%S)] ====="
