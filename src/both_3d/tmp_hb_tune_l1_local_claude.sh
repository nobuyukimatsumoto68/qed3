#!/usr/bin/env bash
# _claude: build + run the L=1 massless Hasenbusch tuning landscape (C6d).
# Loads the last N_CFG thermalized configs from an existing ensemble and prints, per setup
# (mass ladder + multi-timescale step multipliers), the Osborn cost = <N_CG>/(tau^2 <P_acc>).
# Reads back via the .log. No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=test_hasenbusch_tune_l1_v2_claude.log   # _v2: per-stage step counts (keep the v1 log from the global-MDsteps run)
GPU=0
export CUDA_VISIBLE_DEVICES=$GPU
export OMP_NUM_THREADS=4

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_hasenbusch_tune_l1_claude.cu
BIN=test_hb_tune_l1.x

# scan targets: (gsq Nf) pairs -- edit as desired. First try = one strong coupling, cheap.
SCANS=(
  "16 2"
)
N_CFG=1
N_DRAW=1

{
  echo "######## Hasenbusch L1 massless tuning  $(date)  GPU=${GPU} ########"
  echo "==================== BUILD ===================="
  if ! $NVCC $NVCCFLAGS $INCLUDES -DLREF=1 $SRC $LDFLAGS -o "$BIN" ; then
    echo "!!! BUILD FAILED -- stopping"
    exit 1
  fi
  for s in "${SCANS[@]}"; do
    set -- $s
    gsq="$1"
    Nf="$2"
    echo ""
    echo "-------------------- RUN gsq=${gsq} Nf=${Nf} (GPU ${GPU}) --------------------"
    ./"$BIN" "$gsq" "$Nf" "$N_CFG" "$N_DRAW"
  done
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
