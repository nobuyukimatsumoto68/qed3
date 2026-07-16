#!/usr/bin/env bash
# _claude: L2 setup comparison {0,0.5}{3,3} vs {0,0.5}{2,4} (test_hasenbusch_npole_claude.cu, built with
# -DL2_COMPARE). Run this WHEN A GPU FREES UP: set GPU below to the free device. Fully isolated from the
# other runs -- separate binary (test_hb_npole_L2cmp.x) + separate LOG (test_hasenbusch_npole_L2cmp_claude.log,
# does NOT overwrite the others). N_TRAJ=10 real trajectories -> per-frame force + per-traj dH/P/sec +
# <|dH|>/acceptance/Osborn Cost for each setup. No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=test_hasenbusch_npole_L2cmp_claude.log    # SEPARATE log
GPU=0        # <-- SET to the free GPU (0 or 1) before running
export CUDA_VISIBLE_DEVICES=$GPU
export OMP_NUM_THREADS=4

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -DL2_COMPARE -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_hasenbusch_npole_claude.cu

GSQ=8
NF=2
CONFIG_K=0
N_TRAJ=10
WINFAC=0      # legacy arg (window factor fixed at 2)

L=2
bin="test_hb_npole_L2cmp.x"                    # SEPARATE binary (does not clobber test_hb_npole_L2.x)

{
  echo "######## Hasenbusch L2 comparison {0,0.5}{3,3} vs {0,0.5}{2,4} on GPU${GPU}  $(date)  gsq=${GSQ} Nf=${NF} N_TRAJ=${N_TRAJ} ########"
  echo ""
  echo "==================== BUILD L=${L} (-DL2_COMPARE) ===================="
  if ! $NVCC $NVCCFLAGS $INCLUDES -DLREF=${L} $SRC $LDFLAGS -o "$bin" ; then
    echo "!!! BUILD FAILED -- aborting"
  else
    echo "-------------------- RUN L=${L} (gsq=${GSQ} Nf=${NF}) on GPU${GPU} --------------------"
    ./"$bin" "$GSQ" "$NF" "$CONFIG_K" "$N_TRAJ" "$WINFAC"
  fi
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
