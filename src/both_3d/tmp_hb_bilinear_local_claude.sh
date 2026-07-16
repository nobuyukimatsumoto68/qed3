#!/usr/bin/env bash
# _claude: build + run the Hasenbusch external-bra (bilinear) force L=1 FD gate
# (C1 of hasenbusch_massless_impl_plan_claude.md). Validates:
#   (i)  Term B = 2 Re<phi|K|eta> analytic grad vs central-FD of 2Re<phi|D_ov|eta>
#        for bra != ket  (tol ~1e-5), and
#   (ii) phi==eta cross-check: bilinear force == standard massless PF force (tol ~1e-9).
# Default reference grad (NO -DGRAD_L1/2/4) is the required gate; a -DGRAD_L4 phase
# additionally exercises the production block grad. Reads back via the .log. No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=test_hasenbusch_bilinear_l1_claude.log
GPU=0                       # set to 1 to use the other GPU
export CUDA_VISIBLE_DEVICES=$GPU

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_hasenbusch_bilinear_l1_claude.cu

run_phase () {
  local tag="$1"
  local lref="$2"
  local extra="$3"
  local bin="test_hb_bilinear_${tag}.x"
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
  echo "######## Hasenbusch bilinear-force gate  $(date)  GPU=${GPU} ########"
  run_phase "L1_default" 1 ""            # required C1 gate (default reference grad)
  run_phase "L1_gradl4"  1 "-DGRAD_L4"   # also exercise the production block grad
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
