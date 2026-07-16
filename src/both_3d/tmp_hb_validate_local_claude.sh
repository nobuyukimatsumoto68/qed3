#!/usr/bin/env bash
# _claude: build (L=1,2,4) + run the frozen-production Hasenbusch VALIDATION on a thermalized gsq8 config:
# (1) short real HMC chain -> per-traj dH + acceptance + <P_acc>; (2) dH~tau^2 floor-free check (nsteps 2,4,8).
# Uses the production pipeline (frozen window+ladder+nsteps=2, tau=1.0, 2-level MinimumNorm2Hasenbusch). No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=test_hasenbusch_validate_claude.log
GPU=0
export CUDA_VISIBLE_DEVICES=$GPU
export OMP_NUM_THREADS=4

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_hasenbusch_validate_claude.cu

GSQ=8
NF=2
N_TRAJ=10

LEVELS=(1 2 4)

{
  echo "######## Hasenbusch validation  $(date)  GPU=${GPU}  gsq=${GSQ} Nf=${NF} N_TRAJ=${N_TRAJ} ########"
  for L in "${LEVELS[@]}"; do
    bin="test_hb_validate_L${L}.x"
    echo ""
    echo "==================== BUILD L=${L} ===================="
    if ! $NVCC $NVCCFLAGS $INCLUDES -DLREF=${L} $SRC $LDFLAGS -o "$bin" ; then
      echo "!!! BUILD FAILED L=${L} -- skipping"
      continue
    fi
    echo "-------------------- RUN L=${L} (last thermalized config) --------------------"
    ./"$bin" "$GSQ" "$NF" 0 "$N_TRAJ"
  done
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
