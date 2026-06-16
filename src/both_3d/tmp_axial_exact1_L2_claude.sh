#!/usr/bin/env bash
# SINGLE-insertion (--ins 0) EXACT AXIAL at L=2: one tp site + one sp link (NOT the full diagonal sum --
# the L=2 lattice is not vertex/edge-transitive, so this is one representative insertion: shape + sign only,
# matched at a reference dt in the notebook).  Output corr_deter_exact1_axial_L2 (keys h0/t0_b/{tp,sp}/Apm).
#
# COST: the K cache data_free_Kcache_L2/K_ins0.h5 is PRESENT -> cache HIT (no op_K solves).  The one-time
# cost is the G=(1-D_ov) build at L=2 (N=10752 op_oneMinusD applies, ~10-15 min) -> cached to
# data_free_Gcache_L2/G.h5 for future L=2 axial runs.  Plus the dense matmuls + conn_shift2.  GPU1, OMP=4.
#
# RERUN NOTE (no rm in this script): if corr_deter_exact1_axial_L2/corr.0.h5 already exists, delete it
# yourself first.  Read back: jj_exact1_axial_free_L2_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES=1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
LOG=jj_exact1_axial_free_L2_claude.log

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

echo "### compile L=2  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=2 jj_exact_axial_deter_free_claude.cu -o jj_exact_axial_deter_free_L2.o \
  2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}; if [ "$st" -ne 0 ]; then echo "### BUILD FAILED ($st) ###" | tee -a "$LOG"; exit 1; fi

echo "### run --ins 0 L=2 (lattice prop, n-t0=2; K cache HIT, builds G once)  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
./jj_exact_axial_deter_free_L2.o --ins 0 --n-t0 2 2>&1 | tee -a "$LOG"
echo "### done L=2  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
echo "### -> corr_deter_exact1_axial_L2; re-run comp_trio_jj_axial_claude.ipynb for the L=2 exact curve ###" | tee -a "$LOG"
