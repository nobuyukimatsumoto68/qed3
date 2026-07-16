#!/bin/bash
# _claude: L=1 MASSIVE Hasenbusch ACCEPTANCE TEST -- the "c += rescaled mass" scheme
# (massive_production_impl_plan_claude.md). The auxiliary ladder coefficients are shifted UP by the physical
# rescaled mass resc = m*Abar/abar_s, so the inter-frame gaps (tuned massless) are preserved for any mass ->
# no retuning, no degeneracy even at heavy m (m=1.5 -> frame0 coeff 1.419, frame1 = 1.0+1.419 = 2.419).
# Driver = tuning/hmc_hasenbusch_claude.cu (modified: shift applied after base). gsq = 1.5 (largest L1).
#
#   L1 gsq1.5, physical masses m in {0.1, 0.5, 1.0, 1.5}, KMAX=5 -> 4 trajectories (quick acceptance smoke).
#   (traj count = KMAX-1; bump KMAX for a longer run. From a cold start the first dH are thermalization.)
#
# MPS 2/GPU: gpu0 = {m0.1, m0.5}, gpu1 = {m1.0, m1.5}.
# Read acceptance/dH from the per-traj "# dH : .. is_accept : .. rate : .." lines in the tee'd logs.
#
# Run detached:  nohup bash run_massive_test_L1_claude.sh > run_massive_test_L1_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/tuning
cd "$SRCDIR" || exit 1
source ../../env.sh

NF=2
NU0=1.0
GSQ=1.5
KMAX=5            # traj count = KMAX-1 = 4 (quick acceptance smoke test)
KRNG=1
MASSES="0.1 0.5 1.0 1.5"

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=hmc_hasenbusch_claude.cu

# ---- MPS daemon ---------------------------------------------------------------------------------
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

# ---- build L1 -----------------------------------------------------------------------------------
OUT=hmc_test_L1.out
echo "===== build $OUT (LREF=1 KMAX=$KMAX)  [$(date +%F_%H:%M:%S)] ====="
$NVCC $NVCCFLAGS $INCLUDES -DLREF=1 -DKMAX=$KMAX -DKRNG=$KRNG "$SRC" $LDFLAGS -o "$OUT" 2>&1 | tee build_test_L1_claude.log
test -f "$OUT" || { echo "BUILD FAILED"; exit 1; }

# ---- launch -------------------------------------------------------------------------------------
run_one () {   # gpu mass
  local gpu=$1 mass=$2
  echo "### [gpu$gpu L1 gsq$GSQ m$mass] start [$(date +%F_%H:%M:%S)] ###"
  CUDA_VISIBLE_DEVICES=$gpu ./$OUT $GSQ $NF $NU0 $mass 0.0 \
    2>&1 | tee run_test_L1_gsq${GSQ}_m${mass}_claude.log
  echo "### [gpu$gpu L1 gsq$GSQ m$mass] done (status ${PIPESTATUS[0]}) [$(date +%F_%H:%M:%S)] ###"
}
run_seq () {   # gpu mass... : sequential on one GPU slot
  local gpu=$1; shift
  for m in "$@"; do run_one $gpu "$m"; done
}

echo "===== launch (MPS 2/GPU): gpu0={m0.1,m0.5}, gpu1={m1.0,m1.5}  [$(date +%F_%H:%M:%S)] ====="
run_seq 0 0.1 0.5 &   # gpu0
run_seq 1 1.0 1.5 &   # gpu1
wait
echo "===== L1 massive acceptance test finished  [$(date +%F_%H:%M:%S)] ====="
