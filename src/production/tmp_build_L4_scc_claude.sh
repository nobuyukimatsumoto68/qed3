#!/bin/bash -l
# tmp_build_L4_scc_claude.sh  (_claude handoff, 2026-07-16)
# Build-only SMOKE TEST for the SCC L4 setup: compiles ONE binary (Nf2, k400, sm_70) via the production
# wrapper in DRYRUN (no submit) + single-arch/single-Nf mode, to validate that hmc_hasenbusch_block_scc_claude.cu
# (absolute geometry paths) + the SCC module include/lib flags + HighFive all compile cleanly.
# NO qsub, NO submit, NO run (login node has no GPU). Output tee'd to the log below for Claude to read back.
#
# Run on an SCC login node:   bash tmp_build_L4_scc_claude.sh
set -u

cd /projectnb/qfe/nmatsum/qed3/src/production || exit 1
LOG=build_smoketest_L4_scc_claude.log

{
  echo "######## SCC L4 build smoke test  $(date) ########"
  echo "-------- geometry data present? --------"
  ls -la /projectnb/qfe/nmatsum/qed3/geometry/data/*n4*.dat

  echo "-------- build one binary (Nf2 k400 sm_70, DRYRUN=no-submit) --------"
  DRYRUN=1 \
  WHICH=massless \
  NF_LIST=2 GSQ_LIST=2.0 MASS_LIST=0.0 \
  ARCH_LIST=sm_70 SUBMIT_ARCHS=sm_70 \
  bash run_wrapper_L4_scc_claude.sh
  rc=$?

  echo "-------- result --------"
  if [ -f hmc_L4_scc_Nf2_k400_sm_70.out ]
  then
    echo ">>> COMPILE OK: hmc_L4_scc_Nf2_k400_sm_70.out"
    ls -la hmc_L4_scc_Nf2_k400_sm_70.out
  else
    echo ">>> COMPILE FAILED (see build_L4_scc_Nf2_k400_sm_70_claude.log above); wrapper rc=$rc"
  fi
  echo "######## done ########"
} 2>&1 | tee "$LOG"
