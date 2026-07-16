#!/usr/bin/env bash
# _claude: build + run the Hasenbusch C2 stack gate (heatbath identity + force-vs-FD)
# (C2 of hasenbusch_massless_impl_plan_claude.md). Validates the HasenbuschPF manager:
#   (1) heatbath S_i == xi_i^dag xi_i at generation (rel ~1e-7),
#   (2) full summed force get_force() == -dS/dU vs central FD (|grad+fd| ~1e-5),
# for the K=1 {0,0.1} and K=2 {0,0.1,0.4} ladders. Default reference grad required; a -DGRAD_L4
# phase also exercises the production block grad. Reads back via the .log. No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=test_hasenbusch_force_l1_claude.log
GPU=0
export CUDA_VISIBLE_DEVICES=$GPU

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_hasenbusch_force_l1_claude.cu

run_phase () {
  local tag="$1"
  local lref="$2"
  local extra="$3"
  local bin="test_hb_force_${tag}.x"
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
  echo "######## Hasenbusch C2 stack gate  $(date)  GPU=${GPU} ########"
  run_phase "L1_default" 1 ""            # required C2 gate (default reference grad)
  run_phase "L1_gradl4"  1 "-DGRAD_L4"   # also exercise the production block grad
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
