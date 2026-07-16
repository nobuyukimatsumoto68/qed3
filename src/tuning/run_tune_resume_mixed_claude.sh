#!/bin/bash
# _claude: RESUME the massless tuning with the L4 FORCE in MIXED precision (-DMIXED_FORCE, ~1.35x on
# the L4 force solve) + use the FREED gpu0 slot (L1 done) to pack a 2nd L4 there. All (L,gsq) dirs
# AUTO-RESUME from their checkpoints (L2 ~19, L4 ~1). L1 already finished -> skipped.
#
# Packing (2/GPU, all MPS): gpu0 = {L4mix gsq2.0 || L4mix gsq4.0};  gpu1 = {L2 gsq1.0->2.0->3.0 || L4mix gsq6.0}
#
# BEFORE running this, STOP ALL old tuning jobs (this script re-launches L2 + L4, and TWO processes
# writing the SAME output dir would corrupt the checkpoints):
#     pkill -f run_tune_massless_claude.sh   # the parent + its slot subshells
#     pkill -f hmc_tune_L2.out               # old fp64 L2 (this script resumes it)
#     pkill -f hmc_tune_L4.out               # old fp64 L4 (2.0, 6.0)   [mixed binary = hmc_tune_L4_mix.out]
#   (L1 is already done -> not re-run. All dirs AUTO-RESUME from their last checkpoint.)
# Then:  nohup bash run_tune_resume_mixed_claude.sh > run_tune_resume_mixed_claude.log 2>&1 &
#
# CAVEAT: -DMIXED_FORCE (Chunk 3) is UNVALIDATED for the force (the dH~tau^2 check never ran) and this is
# its FIRST compile. The mixed force converges to the SAME TOL_INNER=1e-9, so dH should match fp64 -- but
# SANITY-CHECK the first mixed-L4 dH against the fp64 dH already in run_tune_L4_gsq{2.0,6.0}_claude.log.
# If the mixed dH looks off (floor / blow-up), revert to the fp64 hmc_tune_L4.out.
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/tuning
cd "$SRCDIR" || exit 1
source ../../env.sh

NF=2
NU0=1.0
KMAX=25
KRNG=1

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=hmc_hasenbusch_claude.cu

# ---- ensure ambient MPS daemon up --------------------------------------------------------------
if pgrep -f nvidia-cuda-mps-control >/dev/null
then echo "MPS daemon: already running"
else
  echo "MPS daemon not up -- starting"; nvidia-cuda-mps-control -d
  for i in 1 2 3 4 5; do pgrep -f nvidia-cuda-mps-control >/dev/null && break; sleep 1; done
fi
pgrep -f nvidia-cuda-mps-control >/dev/null || { echo "ERROR: MPS daemon failed"; exit 1; }

# ---- build the MIXED-FORCE L4 binary (L2 fp64 hmc_tune_L2.out already built) --------------------
echo "===== build hmc_tune_L4_mix.out (LREF=4 MIXED_FORCE)  [$(date +%F_%H:%M:%S)] ====="
$NVCC $NVCCFLAGS $INCLUDES -DLREF=4 -DMIXED_FORCE -DKMAX=$KMAX -DKRNG=$KRNG "$SRC" $LDFLAGS \
  -o hmc_tune_L4_mix.out 2>&1 | tee build_tune_L4_mix_claude.log
test -f hmc_tune_L4_mix.out || { echo "BUILD FAILED"; exit 1; }
test -f hmc_tune_L2.out     || { echo "hmc_tune_L2.out missing -- run run_tune_massless_claude.sh build first"; exit 1; }

# ---- launch (auto-resume) ----------------------------------------------------------------------
run_bin () {   # gpu binary L gsq logtag
  local gpu=$1 bin=$2 L=$3 gsq=$4 tag=$5
  echo "### [gpu$gpu $bin L$L gsq$gsq] start [$(date +%F_%H:%M:%S)] ###"
  CUDA_VISIBLE_DEVICES=$gpu ./$bin $gsq $NF $NU0 0.0 0.0 \
    2>&1 | tee run_tune_L${L}_gsq${gsq}_${tag}_claude.log
  echo "### [gpu$gpu $bin L$L gsq$gsq] done (${PIPESTATUS[0]}) [$(date +%F_%H:%M:%S)] ###"
}

echo "===== launch: gpu0={L4mix 2.0 || L4mix 4.0}, gpu1={L2 1.0->2.0->3.0 || L4mix 6.0}  [$(date +%F_%H:%M:%S)] ====="

# gpu0: two MIXED L4 packed (uses the freed L1 slot)
run_bin 0 hmc_tune_L4_mix.out 4 2.0 mix &
run_bin 0 hmc_tune_L4_mix.out 4 4.0 mix &

# gpu1: L2 fp64 queue (resume) || one MIXED L4
( for g in 1.0 2.0 3.0; do run_bin 1 hmc_tune_L2.out 2 "$g" resume; done ) &
run_bin 1 hmc_tune_L4_mix.out 4 6.0 mix &

wait
echo "===== resume-mixed tuning finished  [$(date +%F_%H:%M:%S)] ====="
