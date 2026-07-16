#!/usr/bin/env bash
# _claude: build + run the Hasenbusch SETUP COMPARISON for the TWO-OPERATOR split-pole overlap HMC
# (test_hasenbusch_npole_claude.cu), FINALIZED force: ACTION n=31; FORCE n_f=11 @ 2*lambda_min. For L4 it
# compares {0,0.2,0.5}{2,2,2} vs {0,0.5}{2,4} over N_TRAJ real trajectories (fresh momentum each), recording
# per-traj dH/P + <|dH|>/acceptance/Osborn Cost. PRODUCTION numbers: N_TRAJ=10, LEVELS=(4). (N_TRAJ=0 =
# probe-only per-frame cost.) L4 is slow -- ~tens of min/traj. Reads back via .log. No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=test_hasenbusch_npole_claude.log
GPU=0
export CUDA_VISIBLE_DEVICES=$GPU
export OMP_NUM_THREADS=4

NVCC=nvcc
# (add -DIsVerbose2 for per-solve SOLVER/MULTISHIFT #iter; the per-frame CG/sec summary is enough here.)
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_hasenbusch_npole_claude.cu

GSQ=8
NF=2
CONFIG_K=0    # 0 = last config of the massless ensemble
N_TRAJ=10     # PRODUCTION: 10 real trajectories (fresh momentum each) -> dH/P/Cost. 0 = probe-only cost.
WINFAC=0      # legacy arg (window factor is fixed at 2 in the setup comparison).

LEVELS=(1 2 4)   # L1,L2 = header ladder (single chain each, fast); L4 = setup comparison {0,0.2,0.5}{2,2,2} vs {0,0.5}{2,4}.

{
  echo "######## Hasenbusch SETUP COMPARISON (n=31 action / n_f=11 force @ 2*lmin)  $(date)  GPU=${GPU}  gsq=${GSQ} Nf=${NF} N_TRAJ=${N_TRAJ} ########"
  for L in "${LEVELS[@]}"; do
    bin="test_hb_npole_L${L}.x"
    echo ""
    echo "==================== BUILD L=${L} ===================="
    if ! $NVCC $NVCCFLAGS $INCLUDES -DLREF=${L} $SRC $LDFLAGS -o "$bin" ; then
      echo "!!! BUILD FAILED L=${L} -- skipping"
      continue
    fi
    echo "-------------------- RUN L=${L} (gsq=${GSQ} Nf=${NF}) --------------------"
    ./"$bin" "$GSQ" "$NF" "$CONFIG_K" "$N_TRAJ" "$WINFAC"
  done
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
