#!/usr/bin/env bash
# tmp_glue_at01_gap_claude.sh -- fill the 4 EMPTY L1 at0.1 glue ensembles (NM 2026-09-14).
# These 4 have 1999 gauge ckpts but <20 glue h5 (measurement barely started): Nf4 g0.5, Nf4 g1.5,
# Nf6 g0.5, Nf6 g1.5 (all L1, at0.1).  Dumps BOTH glueball obs to match the fermion at0.1 coverage.
#
# The drivers now derive a_t from the ens-dir, so these are flowed with the correct beta_s(0.1).
# NORMAL prefix (glue_f2_v2_shapes / glue_msm_shapes) -- the flow-smearing a_t does NOT affect the
# extracted masses (spatial per-timeslice smearing only), so no _v2.1 label is needed here.
#
# AT_ONLY=0.100000 L_ONLY=1 processes all at0.1 L1 ensembles; the already-complete ones (Nf2 g0.5/1/1.5,
# Nf4/6 g1.0) are SKIPPED by the per-config "complete" gate, so effectively only the 4 empties are
# measured (~1994 cfg each, both obs).  CPU-only; NO rm, NO process kill.
#
# Run detached:
#   nohup bash tmp_glue_at01_gap_claude.sh > tmp_glue_at01_gap_claude.log 2>&1 &
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1

NWORK="${NWORK:-8}"

echo "################ at0.1 L1 glue GAP dump START $(date) ################"

echo "======== phase 1: F^2 / F^4 (glue_f2_v2_shapes), at0.1 L1 ========"
AT_ONLY=0.100000 L_ONLY=1 NWORK="$NWORK" bash run_glue_f2_v2_sweep_claude.sh > run_glue_f2_v2_sweep_at01gap_claude.log 2>&1
RC1=$?
echo "phase 1 exit=$RC1 $(date)"
if [ "$RC1" -ne 0 ]
then
  echo "ABORT: phase 1 failed rc=$RC1 -- see run_glue_f2_v2_sweep_at01gap_claude.log"
  exit 1
fi

echo "======== phase 2: linear F_12 (glue_msm_shapes), at0.1 L1 ========"
AT_ONLY=0.100000 L_ONLY=1 NWORK="$NWORK" bash run_glue_msm_topup_claude.sh > run_glue_msm_topup_at01gap_claude.log 2>&1
RC2=$?
echo "phase 2 exit=$RC2 $(date)"
if [ "$RC2" -ne 0 ]
then
  echo "ABORT: phase 2 failed rc=$RC2 -- see run_glue_msm_topup_at01gap_claude.log"
  exit 1
fi

echo "################ at0.1 L1 glue GAP dump DONE $(date) ################"
