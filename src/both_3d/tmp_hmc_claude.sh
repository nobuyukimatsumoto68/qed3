#!/usr/bin/env bash
# _claude: build + run the L>2 massive-HMC test (diagonal m_L) for L=2 and L=4.
# Tests per (L x m in {0, 0.1, 0.1i}): (1) reversibility, (2) dH ~ tau^2 scaling, (3) trajectory
# sniff <exp(-dH)>~1. Reads back via test_hmc_diag_mass_claude.log. No rm anywhere (project rule).
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=test_hmc_diag_mass_claude.log
GPU=0                       # MPS up on either GPU; set to 1 for the other card
export CUDA_VISIBLE_DEVICES=$GPU

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_hmc_diag_mass_claude.cu

run_phase () {
  local lref="$1"
  local bin="test_hmc_diag_L${lref}.x"
  echo ""
  echo "==================== BUILD L=${lref} (LREF=${lref}) ===================="
  if ! $NVCC $NVCCFLAGS $INCLUDES -DLREF=${lref} $SRC $LDFLAGS -o "$bin" ; then
    echo "!!! BUILD FAILED: L=${lref} -- skipping run"
    return
  fi
  echo "-------------------- RUN L=${lref} (GPU ${GPU}) --------------------"
  if ./"$bin" ; then
    echo ">>> L=${lref}: exit 0 (ALL PASS)"
  else
    echo ">>> L=${lref}: NONZERO exit (some check FAILED -- inspect above)"
  fi
}

{
  echo "######## massive-HMC diagonal-mass test  $(date)  GPU=${GPU} ########"
  run_phase 2
  run_phase 4
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
