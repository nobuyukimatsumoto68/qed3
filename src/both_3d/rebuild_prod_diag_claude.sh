#!/usr/bin/env bash
# _claude: rebuild the L2/L4 production HMC binaries with the diagonal measure-weighted mass
# (includes/overlap_wmass_claude.h). Run on the FNAL cluster (nvcc/GPU env). Reads back via
# rebuild_prod_diag_claude.log. Removes only the two stale .o.
#
# Run AFTER cancelling the old scalar chains (see the launch sequence), so no afterany
# successor execs this new binary against an old-mass checkpoint dir.
set -u
SRCDIR=/project/qed3/qed3/src/both_3d
LOG=$SRCDIR/rebuild_prod_diag_claude.log
cd "$SRCDIR" || { echo "ERROR: no $SRCDIR"; exit 1; }
source /lustre2/qed3/env.sh
{
  echo "######## rebuild diagonal-mass L2/L4  $(date) ########"; hostname
  echo "== sanity: driver treats mass as PHYSICAL m + operator builds diagonal m_L =="
  grep -nE "physical_m =|mass_coeff = physical_m|m_L = m\*A_y" hmc_fermilab_wmass_L2_claude.cu | head
  grep -nE "Complex mass_coeff|apply_mL\(|apply_mLdag\(" includes/overlap_wmass_claude.h | head -3
  echo "== remove stale .o =="
  rm -fv hmc_fermilab_wmass_L2_claude.o hmc_fermilab_wmass_L4_claude.o
  echo "== build L2 =="; date
  make -f Makefile.fnal hmc_fermilab_wmass_L2_claude.o || { echo "L2 BUILD FAILED"; exit 1; }
  echo "== build L4 =="; date
  make -f Makefile.fnal hmc_fermilab_wmass_L4_claude.o || { echo "L4 BUILD FAILED"; exit 1; }
  echo "== built =="; ls -la hmc_fermilab_wmass_L2_claude.o hmc_fermilab_wmass_L4_claude.o
  echo "######## ALL OK ########"
} 2>&1 | tee "$LOG"
