#!/usr/bin/env bash
# _claude: build + run the MASSIVE Hasenbusch force gate (heatbath identity + force-vs-FD) for the
# "c += rescaled mass" scheme. Tests the L1 ladder {0,1.0} shifted per physical mass m in {0.1,0.5,1.0,1.5}
# (frame0 = phys mass via set_mass; aux = 1.0 + resc, resc = m*Abar/abar_s). Validates the massive force
# path get_force() == +dS/dU (central FD, |grad-fd| < 1e-4). GRAD_L4 = the production block grad. No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=test_hasenbusch_force_massive_l1_claude.log
GPU=0
export CUDA_VISIBLE_DEVICES=$GPU
export OMP_NUM_THREADS=4

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_hasenbusch_force_massive_l1_claude.cu

run_phase () {
  local tag="$1"
  local lref="$2"
  local extra="$3"
  local bin="test_hb_force_massive_${tag}.x"
  echo ""
  echo "==================== BUILD ${tag}  (LREF=${lref} ${extra}) ===================="
  if ! $NVCC $NVCCFLAGS $INCLUDES -DLREF=${lref} ${extra} $SRC $LDFLAGS -o "$bin" ; then
    echo "!!! BUILD FAILED: ${tag} -- skipping run"
    return
  fi
  echo "-------------------- RUN ${tag} (GPU ${GPU}) --------------------"
  if ./"$bin" ; then
    echo ">>> ${tag}: exit 0 (PASS)"
  else
    echo ">>> ${tag}: NONZERO exit (some check FAILED -- inspect above)"
  fi
}

{
  echo "######## MASSIVE Hasenbusch force gate  $(date)  GPU=${GPU} ########"
  run_phase "L1_default" 1 ""            # reference grad
  run_phase "L1_gradl4"  1 "-DGRAD_L4"   # production block grad
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
