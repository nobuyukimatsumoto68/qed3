#!/bin/bash
# Phase 0: build the L=2 and L=4 wmass HMC binaries (prerequisite for benchmarking).
# Targets are built in place in both_3d; the benchmark scripts point straight here
# (/project/qed3/qed3/src/both_3d/hmc_fermilab_wmass_L{2,4}_claude.o).
set -u
SRCDIR=/home/nmatsum/project_qed3/qed3/src/both_3d
LOG=$SRCDIR/build_L2L4_claude.log
cd "$SRCDIR" || exit 1
source /lustre2/qed3/env.sh
{
  echo "===== build hmc_fermilab_wmass_L2_claude.o ====="; date
  make -f Makefile.fnal hmc_fermilab_wmass_L2_claude.o || { echo "L2 BUILD FAILED"; exit 1; }
  echo "===== build hmc_fermilab_wmass_L4_claude.o ====="; date
  make -f Makefile.fnal hmc_fermilab_wmass_L4_claude.o || { echo "L4 BUILD FAILED"; exit 1; }
  echo "===== built binaries ====="
  ls -la hmc_fermilab_wmass_L2_claude.o hmc_fermilab_wmass_L4_claude.o
  date
} 2>&1 | tee "$LOG"
