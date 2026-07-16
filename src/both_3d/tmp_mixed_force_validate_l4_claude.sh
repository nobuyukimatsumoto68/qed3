#!/usr/bin/env bash
# _claude: Chunk-3 validation -- mixed-prec FORCE wired into the overlap force precalc (MIXED_FORCE,
# L4-only). Builds test_hasenbusch_validate_claude.cu TWICE at L4 (fp64 baseline vs -DMIXED_FORCE) and
# runs a real Nf6 gsq8 massless L4 config. COMPARE:
#   (1) short HMC chain: per-traj dH + acceptance + sec  (mixed dH must match fp64 to ~1e-9; sec = speedup)
#   (2) dH ~ tau^2 FLOOR-FREE check (MDsteps 1x/2x/4x, ratio ~4 per doubling): the DECISIVE exactness test
#       -- a step-INDEPENDENT dH floor in the MIXED build would mean the mixed force is not exact.
# Reads back via the .log. No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=tmp_mixed_force_validate_l4_claude.log
export OMP_NUM_THREADS=4
export CUDA_VISIBLE_DEVICES=0

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_hasenbusch_validate_claude.cu

# args: gsq Nf config_k N_TRAJ  (config_k=0 -> last config of Nf6_gsq8..L4)
RUNARGS="8.0 6 0 2"

{
  echo "######## Chunk-3 mixed-FORCE validation (L4, Nf6 gsq8)  $(date) ########"

  echo ""
  echo "==================== BUILD fp64 baseline (LREF=4) ===================="
  if $NVCC $NVCCFLAGS $INCLUDES -DLREF=4 $SRC $LDFLAGS -o test_hb_validate_fp64_L4.x ; then
    echo ">>> BUILD OK (fp64)"
    echo "-------------------- RUN fp64  ($RUNARGS) --------------------"
    ./test_hb_validate_fp64_L4.x $RUNARGS
  else
    echo "!!! BUILD FAILED (fp64)"
  fi

  echo ""
  echo "==================== BUILD MIXED force (LREF=4 -DMIXED_FORCE) ===================="
  if $NVCC $NVCCFLAGS $INCLUDES -DLREF=4 -DMIXED_FORCE $SRC $LDFLAGS -o test_hb_validate_mixed_L4.x ; then
    echo ">>> BUILD OK (mixed)"
    echo "-------------------- RUN mixed ($RUNARGS) --------------------"
    ./test_hb_validate_mixed_L4.x $RUNARGS
  else
    echo "!!! BUILD FAILED (mixed)"
  fi

  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
