#!/usr/bin/env bash
# _claude: COMPILE-ONLY check of the Nf-BLOCK Hasenbusch driver hmc_hasenbusch_block_claude.cu
# (nfblock_hasenbusch_impl_plan_claude.md, chunks 1-4). Instantiates the full block template stack
# (BlockedMat / BlockedForce / HasenbuschPFBlock) at NSTACK = Nf/2 = 1, 2, 3 and L1/L4 -- catches template
# errors. Does NOT run (no ensemble writes). Reads back via the .log. No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=tmp_hb_block_build_claude.log
export OMP_NUM_THREADS=4

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=hmc_hasenbusch_block_claude.cu

# (NF, LREF) combos: NSTACK = NF/2 = 1,2,3 ; L1 (small) and L4 (big lattice).
COMBOS=( "2 1" "4 1" "6 1" "6 4" )

{
  echo "######## Nf-BLOCK Hasenbusch driver COMPILE check  $(date) ########"
  for combo in "${COMBOS[@]}"; do
    set -- $combo
    NFV=$1
    LV=$2
    bin="hmc_hb_block_Nf${NFV}_L${LV}.x"
    echo ""
    echo "==================== BUILD NF=${NFV} LREF=${LV} (NSTACK=$((NFV/2))) ===================="
    if $NVCC $NVCCFLAGS $INCLUDES -DNF=${NFV} -DLREF=${LV} $SRC $LDFLAGS -o "$bin" ; then
      echo ">>> OK: NF=${NFV} LREF=${LV}"
    else
      echo "!!! BUILD FAILED: NF=${NFV} LREF=${LV}"
    fi
  done
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
