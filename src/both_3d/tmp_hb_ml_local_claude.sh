#!/usr/bin/env bash
# _claude: build + run the multi-timescale MinimumNorm2ML validation (C6c):
#   (1) 2-level regression vs MinimumNorm2Hasenbusch, (2) 3-level reversibility, (3) dH~tau^2.
# Massless, ladder {0,0.1,0.4}, lambda_max frozen. Reads back via the .log. No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=test_hasenbusch_ml_claude.log
GPU=0
export CUDA_VISIBLE_DEVICES=$GPU

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_hasenbusch_ml_claude.cu

run_phase () {
  local tag="$1"
  local lref="$2"
  local bin="test_hb_ml_${tag}.x"
  echo ""
  echo "==================== BUILD ${tag}  (LREF=${lref}) ===================="
  if ! $NVCC $NVCCFLAGS $INCLUDES -DLREF=${lref} $SRC $LDFLAGS -o "$bin" ; then
    echo "!!! BUILD FAILED: ${tag} -- skipping run"
    return
  fi
  echo "-------------------- RUN ${tag} (GPU ${GPU}) --------------------"
  if ./"$bin" 8.0 ; then echo ">>> ${tag}: exit 0 (PASS)"; else echo ">>> ${tag}: NONZERO exit (FAIL)"; fi
}

{
  echo "######## Hasenbusch ML integrator validation  $(date)  GPU=${GPU} ########"
  run_phase "L2"  2
  run_phase "L1"  1
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
