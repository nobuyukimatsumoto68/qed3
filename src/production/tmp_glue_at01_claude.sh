#!/usr/bin/env bash
# ============================ SHELVED 2026-09-14 (NM) ============================
# NOT NEEDED for the physics: the driver's a_t enters ONLY the spatial per-timeslice Wilson-flow
# SMEARING (beta_s = at/(vol*gsq)).  Smearing changes the interpolator overlap / signal-to-noise but
# NOT the temporal transfer matrix, so the GEVP/plateau MASSES are smearing-level-independent.  The
# existing at=0.2-flowed at0.1 h5 therefore give the SAME masses as this _v2.1 re-dump would.  The
# only genuinely a_t-dependent step for the result is the mass conversion Delta_eff = -log(lambda)/
# (dt*at), which uses the ANALYSIS a_t (CLI AT=0.1) and was already correct.
# => Kept available (driver fix + _v2.1 label are harmless) but DO NOT run unless we specifically
#    want the smearing-self-consistent at0.1 set.  See findings_gluonic_claude.md.
# ================================================================================
# tmp_glue_at01_claude.sh -- at0.1 glueball obs dump handoff (NM 2026-09-14).
# Runs BOTH glueball dumpers restricted to the at0.1 ensembles (AT_ONLY=0.100000; at0.1 exists only
# at L1/L2).  The drivers now derive a_t from the ensemble-dir string (at_from_ensdir), so at0.1
# configs are flowed with the CORRECT beta_s(0.1) -- not the old hardcoded beta_s(0.2).
#
# NEW LABEL: OUT_SUFFIX=_v2.1 tags this a_t-corrected data with a DISTINCT h5 prefix, so it neither
# clobbers nor is skipped-against the existing (old at=0.2-flowed) at0.1 h5.  NO deletion needed.
#   phase 1 : glue_f2_v2_shapes_v2.1.<k>.h5   (F^2 0++ and F^4, a_t-corrected)
#   phase 2 : glue_msm_shapes_v2.1.<k>.h5     (linear F_{12} l=1,2, a_t-corrected)
# The old glue_f2_v2_shapes.<k>.h5 / glue_msm_shapes.<k>.h5 stay in place, untouched.
# Analysis of the corrected set: pass prefix "glue_f2_v2_shapes_v2.1" / "glue_msm_shapes_v2.1".
#
# RESUME-SAFE per label: the "complete" gate keys on the _v2.1 filename, so a re-run only fills its
# own deficit.  CPU-only; does NOT touch the GPUs.  NO rm, NO process kill anywhere.
#
# Run detached:
#   nohup bash tmp_glue_at01_claude.sh > tmp_glue_at01_claude.log 2>&1 &
# then read run_glue_f2_v2_sweep_at01_claude.log and run_glue_msm_topup_at01_claude.log.
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1

NWORK="${NWORK:-8}"
LABEL="${LABEL:-_v2.1}"   # h5-prefix suffix for the a_t-corrected at0.1 set

echo "################ at0.1 glue obs dump START $(date)  label=$LABEL ################"

echo "======== phase 1: F^2 / F^4 (glue_f2_v2_shapes${LABEL}), at0.1 only ========"
AT_ONLY=0.100000 OUT_SUFFIX="$LABEL" NWORK="$NWORK" bash run_glue_f2_v2_sweep_claude.sh > run_glue_f2_v2_sweep_at01_claude.log 2>&1
RC1=$?
echo "phase 1 exit=$RC1 $(date)"
if [ "$RC1" -ne 0 ]
then
  echo "ABORT: phase 1 (f2_v2 sweep at0.1) failed rc=$RC1 -- see run_glue_f2_v2_sweep_at01_claude.log"
  exit 1
fi

echo "======== phase 2: linear F_12 (glue_msm_shapes${LABEL}), at0.1 only ========"
AT_ONLY=0.100000 OUT_SUFFIX="$LABEL" NWORK="$NWORK" bash run_glue_msm_topup_claude.sh > run_glue_msm_topup_at01_claude.log 2>&1
RC2=$?
echo "phase 2 exit=$RC2 $(date)"
if [ "$RC2" -ne 0 ]
then
  echo "ABORT: phase 2 (msm topup at0.1) failed rc=$RC2 -- see run_glue_msm_topup_at01_claude.log"
  exit 1
fi

echo "################ at0.1 glue obs dump DONE $(date) ################"
