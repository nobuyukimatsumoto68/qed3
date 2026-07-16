#!/usr/bin/env bash
# _claude: build + run the per-term derivative test (localize the L1 dH floor). Compares each
# Hasenbusch force term to the central FD of its scalar with FROZEN vectors (clean FD, no re-solve):
#   Term A per mass {0,0.1,0.4}: grad(D_m,eta) vs -FD[eta^dag DHD_ms eta]
#   Term B: grad_bilinear(phi,eta) vs -FD[2Re<phi|D_ov|eta>]
# |grad+fd| ~ FD floor (<<delta) = exact; ~delta = the culprit. Reads back via the .log. No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=test_hasenbusch_termAB_deriv_claude.log
GPU=1
export CUDA_VISIBLE_DEVICES=$GPU

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_hasenbusch_termAB_deriv_claude.cu

run_phase () {
  local tag="$1"
  local lref="$2"
  local extra="$3"
  local bin="test_hb_termAB_${tag}.x"
  echo ""
  echo "==================== BUILD ${tag}  (LREF=${lref} ${extra}) ===================="
  if ! $NVCC $NVCCFLAGS $INCLUDES -DLREF=${lref} ${extra} $SRC $LDFLAGS -o "$bin" ; then
    echo "!!! BUILD FAILED: ${tag} -- skipping run"
    return
  fi
  echo "-------------------- RUN ${tag} (GPU ${GPU}) --------------------"
  ./"$bin"
}

{
  echo "######## Hasenbusch per-term derivative test  $(date)  GPU=${GPU} ########"
  run_phase "L1_default" 1 ""            # default reference grad
  run_phase "L1_gradl4"  1 "-DGRAD_L4"   # production block grad
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
