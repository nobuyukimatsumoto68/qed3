#!/usr/bin/env bash
# _claude: build + run the Hasenbusch HMC end-to-end validation (C2c/C3 of
# hasenbusch_massless_impl_plan_claude.md): reversibility + dH~tau^2 scaling + a few trajectories,
# with the REAL MinimumNorm2Hasenbusch integrator + HMCHasenbusch on a massless overlap, ladder
# {0,0.1,0.4}. L2 default (small Nt). Reads back via the .log. No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=test_hasenbusch_hmc_l2_claude.log
GPU=0
export CUDA_VISIBLE_DEVICES=$GPU

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_hasenbusch_hmc_l2_claude.cu

run_phase () {
  local tag="$1"
  local lref="$2"
  local args="$3"
  local bin="test_hb_hmc_${tag}.x"
  echo ""
  echo "==================== BUILD ${tag}  (LREF=${lref}) ===================="
  if ! $NVCC $NVCCFLAGS $INCLUDES -DLREF=${lref} $SRC $LDFLAGS -o "$bin" ; then
    echo "!!! BUILD FAILED: ${tag} -- skipping run"
    return
  fi
  echo "-------------------- RUN ${tag} (GPU ${GPU}) [args: ${args}] --------------------"
  if ./"$bin" ${args} ; then
    echo ">>> ${tag}: exit 0 (PASS)"
  else
    echo ">>> ${tag}: NONZERO exit (some check FAILED -- inspect above)"
  fi
}

{
  echo "######## Hasenbusch HMC validation  $(date)  GPU=${GPU} ########"
  run_phase "L2_gsq8"  2 "8.0 8"    # gsq=8, nsteps=8
  run_phase "L1_gsq8"  1 "8.0 8"    # gsq=8, nsteps=8  (includes section (4) standard-PF comparison)
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
