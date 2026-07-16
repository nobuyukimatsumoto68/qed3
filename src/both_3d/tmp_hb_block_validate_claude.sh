#!/usr/bin/env bash
# _claude: block-vs-serial PARITY gate for the Nf-BLOCK Hasenbusch PF. Builds + runs
# test_hasenbusch_block_validate_claude.cu on GPU1: feeds the SAME phi to NSTACK=Nf/2 serial HasenbuschPF
# stacks (summed) and one block HasenbuschPFBlock<...,NSTACK> on one thermalized config, compares S + force
# to solver tolerance (rel < 1e-6 = PASS). NSTACK=1 (Nf=2) already PASSED; the NONTRIVIAL packing check is
# NF=6 (NSTACK=3), e.g. L2. Set NF/LREF below. Writes NO gauge configs (reads one config read-only). No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=test_hasenbusch_block_validate_claude.log
GPU=1
export CUDA_VISIBLE_DEVICES=$GPU
export OMP_NUM_THREADS=4

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_hasenbusch_block_validate_claude.cu

GSQ=8
NF=6          # NONTRIVIAL packing check: NSTACK = NF/2 = 3. (NF=2 = the NSTACK=1 gate, already PASSED.)
CONFIG_K=0    # 0 = last config of the massless ensemble
LREF=2        # L2 for the Nf6 check (needs an Nf6_gsq8...L2 ensemble config to exist).

{
  echo "######## Nf-BLOCK Hasenbusch PARITY (block vs summed serial, NSTACK=$((NF/2)))  $(date)  GPU=${GPU}  gsq=${GSQ} Nf=${NF} L=${LREF} ########"
  bin="test_hb_block_validate_Nf${NF}_L${LREF}.x"
  echo ""
  echo "==================== BUILD NF=${NF} LREF=${LREF} (NSTACK=$((NF/2))) ===================="
  if ! $NVCC $NVCCFLAGS $INCLUDES -DNF=${NF} -DLREF=${LREF} $SRC $LDFLAGS -o "$bin" ; then
    echo "!!! BUILD FAILED -- aborting"
  else
    echo "-------------------- RUN NF=${NF} L=${LREF} (gsq=${GSQ}) on GPU${GPU} --------------------"
    ./"$bin" "$GSQ" "$NF" "$CONFIG_K"
  fi
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
