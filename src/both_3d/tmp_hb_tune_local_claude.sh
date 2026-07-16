#!/usr/bin/env bash
# _claude: build (L=1,2,4) + run the Hasenbusch nsteps tuning for the FROZEN ladders (hasenbusch_ladder(L))
# on a gsq8 massless config. Prints, per L: per-frame + gauge FORCE NORMS, then a scan of per-stage MD
# step-count candidates with the Osborn cost = N_CG/(tau^2 P_acc). Reads back via the .log. No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=test_hasenbusch_tune_claude.log
GPU=0
export CUDA_VISIBLE_DEVICES=$GPU
export OMP_NUM_THREADS=4

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_hasenbusch_tune_claude.cu

GSQ=8
NF=2
N_CFG=1     # uses the last config; N_DRAW reserved
N_DRAW=1

LEVELS=(4)   # setup comparison, currently L4

{
  echo "######## Hasenbusch nsteps tuning (frozen ladders)  $(date)  GPU=${GPU}  gsq=${GSQ} Nf=${NF} ########"
  for L in "${LEVELS[@]}"; do
    bin="test_hb_tune_L${L}.x"
    echo ""
    echo "==================== BUILD L=${L} ===================="
    if ! $NVCC $NVCCFLAGS $INCLUDES -DLREF=${L} $SRC $LDFLAGS -o "$bin" ; then
      echo "!!! BUILD FAILED L=${L} -- skipping"
      continue
    fi
    echo "-------------------- RUN L=${L} (gsq=${GSQ} Nf=${NF}) --------------------"
    ./"$bin" "$GSQ" "$NF" "$N_CFG" "$N_DRAW"
  done
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
