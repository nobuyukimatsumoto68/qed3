#!/bin/bash
# _claude: MASSIVE Nf2 overlap production (params_massive_claude.md). Reuses the massless-tuned HMC
# parameters AS-IS (no retuning, per NM) -- frame 0 = physical mass via argv; ladder/steps/tau unchanged
# (includes/hasenbusch_ladder_claude.h) via the "c += rescaled mass" shift (aux coeffs += resc =
# m*Abar/abar_s -> additive inter-frame gaps preserved). Physical masses m in {0.1,0.5,1.0,1.5}, R=1 units.
#
#   L1 gsq1.5 (KMAX 120)   L2 gsq3.0 (KMAX 80)   L4 gsq6.0 (KMAX 60)   -- largest corrected gsq per L
#   3 L x 4 m = 12 ensembles.  The c+=resc shift keeps the split VALID for ANY mass (no rung constraint).
#
# Driver = the Nf-PACKED production driver hmc_hasenbusch_block_claude.cu (-DNF=2 -> NSTACK=1 no-op block,
# bit-identical to serial at Nf2). L4 built with -DMIXED_FORCE (force-only mixed precision, ~1.35x, exact by
# Metropolis). MPS packing (2/GPU): L1+L2 (cheap) each on a slot; L4 (expensive) split across both GPUs.
#   gpu0 = {L1 all m}  || {L4 m=0.1,0.5}      gpu1 = {L2 all m}  || {L4 m=1.0,1.5}
#
# Output dir (auto, AUTO-RESUMES): Nf2_gsq<g>at0.200000nu01.000000mRe<m>mIm0.000000nt128L<L>_hb1.000000/
# per-traj "# dH : .. is_accept : .. rate : .." in the tee'd logs. KMAX = target sample size (incl.
# thermalization); extend by rebuilding with a larger KMAX (same dir auto-resumes).
#
# Run detached:  nohup bash run_massive_claude.sh > run_massive_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/production
cd "$SRCDIR" || exit 1
source ../../env.sh

NF=2
NU0=1.0
KRNG=1                       # FULL checkpointing (keep every ckpoint_rng)
MASSES="0.1 0.5 1.0 1.5"     # physical masses (R=1)
KMAX_L1=120
KMAX_L2=80
KMAX_L4=60

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=hmc_hasenbusch_block_claude.cu

# ---- 1. ensure the ambient MPS daemon is up -----------------------------------------------------
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

# ---- 2. build one binary per L (gsq/mass are runtime args); L4 gets -DMIXED_FORCE -----------------
build () {
  local L=$1 kmax=$2 extra=$3
  local out=hmc_massive_L${L}.out
  echo "===== build $out (LREF=$L NF=$NF KMAX=$kmax $extra)  [$(date +%F_%H:%M:%S)] ====="
  $NVCC $NVCCFLAGS $INCLUDES -DLREF=$L -DNF=$NF -DKMAX=$kmax -DKRNG=$KRNG $extra "$SRC" $LDFLAGS -o "$out" 2>&1 | tee build_massive_L${L}_claude.log
  test -f "$out" || { echo "BUILD FAILED: $out"; exit 1; }
}
build 1 $KMAX_L1 ""
build 2 $KMAX_L2 ""
build 4 $KMAX_L4 "-DMIXED_FORCE"

# ---- 3. launch ----------------------------------------------------------------------------------
run_one () {   # gpu L gsq mass
  local gpu=$1 L=$2 gsq=$3 mass=$4
  echo "### [gpu$gpu L$L gsq$gsq m$mass] start [$(date +%F_%H:%M:%S)] ###"
  CUDA_VISIBLE_DEVICES=$gpu ./hmc_massive_L${L}.out $gsq $NF $NU0 $mass 0.0 \
    2>&1 | tee run_massive_L${L}_gsq${gsq}_m${mass}_claude.log
  echo "### [gpu$gpu L$L gsq$gsq m$mass] done (status ${PIPESTATUS[0]}) [$(date +%F_%H:%M:%S)] ###"
}
run_masses () {   # gpu L gsq mass... : run the listed masses sequentially on one GPU slot
  local gpu=$1 L=$2 gsq=$3; shift 3
  for m in "$@"; do run_one $gpu $L "$gsq" "$m"; done
}

echo "===== launch (MPS 2/GPU): gpu0={L1 || L4 m0.1,0.5}, gpu1={L2 || L4 m1.0,1.5}  [$(date +%F_%H:%M:%S)] ====="

run_masses 0 1 1.5 $MASSES    &   # gpu0 slot A: L1 gsq1.5, all masses
run_masses 0 4 6.0 0.1 0.5    &   # gpu0 slot B: L4 gsq6.0, m=0.1,0.5
run_masses 1 2 3.0 $MASSES    &   # gpu1 slot A: L2 gsq3.0, all masses
run_masses 1 4 6.0 1.0 1.5    &   # gpu1 slot B: L4 gsq6.0, m=1.0,1.5

wait
echo "===== massive runs finished  [$(date +%F_%H:%M:%S)] ====="
