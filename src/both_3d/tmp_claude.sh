#!/bin/bash
# Handoff: force-rebuild the L=2 and L=4 wmass HMC binaries after the 2026-06-16
# parameter changes (nsteps L2=9/L4=12, kmax L4=200, npole L2=17/L4=13,
# NPARALLEL_DUPDATE=1 -> NSTREAMS=1, new max_sec wall-budget arg).
# Removes the stale .o first so the edits are guaranteed to be picked up.
# Run:  bash tmp_claude.sh     (reads back via rebuild_L2L4_claude.log)
set -u
SRCDIR=/home/nmatsum/project_qed3/qed3/src/both_3d
LOG=$SRCDIR/rebuild_L2L4_claude.log
cd "$SRCDIR" || exit 1
source /lustre2/qed3/env.sh
{
  echo "############ FORCE REBUILD L2/L4 wmass HMC ############"; hostname; date
  echo
  echo "===== locked params in source (sanity) ====="
  for L in L2 L4; do
    f=hmc_fermilab_wmass_${L}_claude.cu
    echo "--- $f ---"
    grep -nE "NPARALLEL_DUPDATE=|NSTREAMS=|Fermion D\(DW|const int kmax=|const int k_ckpoint_rng=|if\(Nf==2\) nsteps =|max_sec = atof" "$f"
  done
  echo
  echo "===== remove stale .o ====="
  rm -fv hmc_fermilab_wmass_L2_claude.o hmc_fermilab_wmass_L4_claude.o
  echo
  echo "===== build L2 ====="; date
  make -f Makefile.fnal hmc_fermilab_wmass_L2_claude.o || { echo "L2 BUILD FAILED"; exit 1; }
  echo
  echo "===== build L4 ====="; date
  make -f Makefile.fnal hmc_fermilab_wmass_L4_claude.o || { echo "L4 BUILD FAILED"; exit 1; }
  echo
  echo "===== built binaries ====="
  ls -la hmc_fermilab_wmass_L2_claude.o hmc_fermilab_wmass_L4_claude.o
  echo "ALL OK"; date
} 2>&1 | tee "$LOG"
