#!/usr/bin/env bash
# _claude: PARALLEL GPU1 run of the L4 setup {0,0.2,0.5} steps {2,2,4} for the TWO-OPERATOR split-pole
# overlap HMC (test_hasenbusch_npole_claude.cu, built with -DL4_STEPS224). Runs alongside the GPU0 run
# WITHOUT touching it: separate GPU (CUDA_VISIBLE_DEVICES=1), separate binary, and a SEPARATE LOG
# (test_hasenbusch_npole_gpu1_claude.log -- does NOT overwrite test_hasenbusch_npole_claude.log). N_TRAJ=10
# real trajectories -> per-frame force + per-traj dH/P/sec + <|dH|>/acceptance/Osborn Cost. No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=test_hasenbusch_npole_gpu1_claude.log     # SEPARATE log (GPU0 run keeps test_hasenbusch_npole_claude.log)
GPU=1
export CUDA_VISIBLE_DEVICES=$GPU
export OMP_NUM_THREADS=4

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -DL4_STEPS224 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_hasenbusch_npole_claude.cu

GSQ=8
NF=2
CONFIG_K=0
N_TRAJ=10
WINFAC=0      # legacy arg (window factor fixed at 2)

LEVELS=(4)    # GPU1: only L4 {0,0.2,0.5}{2,2,4}

{
  echo "######## Hasenbusch L4 {0,0.2,0.5}{2,2,4} on GPU1 (parallel)  $(date)  GPU=${GPU}  gsq=${GSQ} Nf=${NF} N_TRAJ=${N_TRAJ} ########"
  for L in "${LEVELS[@]}"; do
    bin="test_hb_npole_gpu1_L${L}.x"     # SEPARATE binary (does not clobber the GPU0 test_hb_npole_L${L}.x)
    echo ""
    echo "==================== BUILD L=${L} (-DL4_STEPS224) ===================="
    if ! $NVCC $NVCCFLAGS $INCLUDES -DLREF=${L} $SRC $LDFLAGS -o "$bin" ; then
      echo "!!! BUILD FAILED L=${L} -- skipping"
      continue
    fi
    echo "-------------------- RUN L=${L} (gsq=${GSQ} Nf=${NF}) on GPU${GPU} --------------------"
    ./"$bin" "$GSQ" "$NF" "$CONFIG_K" "$N_TRAJ" "$WINFAC"
  done
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
