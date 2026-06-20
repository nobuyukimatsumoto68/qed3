#!/usr/bin/env bash
# _claude: build + run the operator-level Ward / current-conservation check for the diagonal mass.
# Per m in {0, 0.1, 0.1i}:  (A) [D_m,Theta]xi = sum_l dtheta_l K^l xi ; (B) sum_z tr(K^{wz} D_m^{-1}) = 0.
# L=1 (analytic reference, deterministic gate) + L=2 (stochastic). Reads back via test_ward_diag_mass_claude.log.
# No rm anywhere (project rule).
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=test_ward_diag_mass_claude.log
GPU=1                       # run on GPU 1 (HMC test is on GPU 0); MPS up
export CUDA_VISIBLE_DEVICES=$GPU

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_ward_diag_mass_claude.cu

run_phase () {
  local lref="$1"
  local bin="test_ward_diag_L${lref}.x"
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
  echo "######## Ward / current-conservation diagonal-mass check  $(date)  GPU=${GPU} ########"
  run_phase 1
  run_phase 2
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
