#!/usr/bin/env bash
# _claude: RESTART L4 on GPU1 with the new setup comparison (test_hasenbusch_npole_claude.cu, DEFAULT L4
# branch -- NO -DL4_STEPS224):
#     {0,0.2,0.5}       steps {3,3,3}     (3-stage)
#     {0,0.1,0.3,0.5}   steps {2,2,2,4}   (4-stage, heaviest frame sub-cycled x2)
#     {0,0.2,0.5}       steps {2,2,6}     (3-stage, heaviest frame sub-cycled x3)
# N_TRAJ=4 real trajectories each (fresh momentum) -> per-frame force + per-traj dH/P/sec + <|dH|>/acceptance/
# Osborn Cost. Finalized force (n_f=11 @ 2*lambda_min, action n=31). Isolated: GPU1, SEPARATE binary
# (test_hb_npole_L4v2.x) + SEPARATE log (does NOT overwrite the {2,2,4} gpu1 log or the others). ~1h/traj ->
# ~12h for 3 setups x 4 traj. Make sure GPU1 is free (keep the L2-comparison run on GPU0). No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=test_hasenbusch_npole_L4v2_claude.log     # SEPARATE log
GPU=1
export CUDA_VISIBLE_DEVICES=$GPU
export OMP_NUM_THREADS=4

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_hasenbusch_npole_claude.cu

GSQ=8
NF=2
CONFIG_K=0
N_TRAJ=4
WINFAC=0      # legacy arg (window factor fixed at 2)

L=4
bin="test_hb_npole_L4v2.x"                     # SEPARATE binary

{
  echo "######## Hasenbusch L4 restart {0,0.2,0.5}{3,3,3} vs {0,0.1,0.3,0.5}{2,2,2,4} on GPU${GPU}  $(date)  gsq=${GSQ} Nf=${NF} N_TRAJ=${N_TRAJ} ########"
  echo ""
  echo "==================== BUILD L=${L} ===================="
  if ! $NVCC $NVCCFLAGS $INCLUDES -DLREF=${L} $SRC $LDFLAGS -o "$bin" ; then
    echo "!!! BUILD FAILED -- aborting"
  else
    echo "-------------------- RUN L=${L} (gsq=${GSQ} Nf=${NF}) on GPU${GPU} --------------------"
    ./"$bin" "$GSQ" "$NF" "$CONFIG_K" "$N_TRAJ" "$WINFAC"
  fi
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
