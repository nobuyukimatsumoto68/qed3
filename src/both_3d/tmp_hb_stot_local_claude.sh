#!/usr/bin/env bash
# _claude: build (L=1,2,4) + run the s_tot (trajectory length) scan for the frozen Hasenbusch setup on
# the stiffest (min lambda_min) config per L. Reports per s_tot, averaged over 4 refreshed momenta:
# <dH>, smeared plaquette <P>, <|dP|>, decorrelation-per-cost <|dP|>/s_tot and <|dP|>/<N_CG>. No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=test_hasenbusch_stot_claude.log
GPU=0
export CUDA_VISIBLE_DEVICES=$GPU
export OMP_NUM_THREADS=4

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_hasenbusch_stot_claude.cu

GSQ=8
NF=2
N_MOM=4
T_FLOW=1.0

LEVELS=(1 2 4)   # L=1,2,4 order

{
  echo "######## Hasenbusch s_tot scan  $(date)  GPU=${GPU}  gsq=${GSQ} Nf=${NF} N_MOM=${N_MOM} t_flow=${T_FLOW} ########"
  for L in "${LEVELS[@]}"; do
    bin="test_hb_stot_L${L}.x"
    echo ""
    echo "==================== BUILD L=${L} ===================="
    if ! $NVCC $NVCCFLAGS $INCLUDES -DLREF=${L} $SRC $LDFLAGS -o "$bin" ; then
      echo "!!! BUILD FAILED L=${L} -- skipping"
      continue
    fi
    echo "-------------------- RUN L=${L} (min-lambda_min config) --------------------"
    ./"$bin" "$GSQ" "$NF" 0 "$N_MOM" "$T_FLOW"
  done
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
