#!/usr/bin/env bash
# tmp_claude.sh -- L4 glueball obs top-up handoff (NM 2026-09-14).
# Runs BOTH resumable glueball dumpers restricted to L=4, so the recently-landed L4 configs are
# measured.  Both underlying sweeps build their own per-L binary and skip already-complete configs
# (the driver's "complete" h5 gate), so this is safe to re-run after the still-generating Nf6 gsq6
# stream finishes -- it will only compute the deficit.
#
#   phase 1 : glue_f2_v2_shapes  (F^2 0++ and F^4)      -> data_<ens>/glue_f2_v2_shapes.<k>.h5
#   phase 2 : glue_msm_shapes     (linear F_{12} l=1,2)  -> data_<ens>/glue_msm_shapes.<k>.h5
#
# Current L4 deficit (ckpoints 799, except Nf6 gsq6 = 778 still generating):
#   f2_v2 : Nf2 complete; Nf4 lags 43/54/77 (g2/4/6); Nf6 lags 61/71/138 (g2/4/6).
#   msm   : all lag (599 of 799; Nf6 gsq6 = 445).
#
# CPU-only (gradient flow + host Wilson-loop holonomies); does NOT touch the GPUs, so it is safe
# alongside the fermionic jj jobs.  NO rm, NO process kill anywhere.
#
# Run detached:
#   nohup bash tmp_claude.sh > tmp_glue_L4_topup_claude.log 2>&1 &
# then read back run_glue_f2_v2_sweep_L4_claude.log and run_glue_msm_topup_L4_claude.log.
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1

# one worker per L4 ensemble (9 of them); user-side workload, not the 4-core Claude cap.
NWORK="${NWORK:-8}"

echo "################ L4 glue obs top-up START $(date) ################"

echo "======== phase 1: F^2 / F^4 (glue_f2_v2_shapes), L4 only ========"
L_ONLY=4 NWORK="$NWORK" bash run_glue_f2_v2_sweep_claude.sh > run_glue_f2_v2_sweep_L4_claude.log 2>&1
RC1=$?
echo "phase 1 exit=$RC1 $(date)"
if [ "$RC1" -ne 0 ]
then
  echo "ABORT: phase 1 (f2_v2 sweep) failed rc=$RC1 -- see run_glue_f2_v2_sweep_L4_claude.log"
  exit 1
fi

echo "======== phase 2: linear F_12 (glue_msm_shapes), L4 only ========"
L_ONLY=4 NWORK="$NWORK" bash run_glue_msm_topup_claude.sh > run_glue_msm_topup_L4_claude.log 2>&1
RC2=$?
echo "phase 2 exit=$RC2 $(date)"
if [ "$RC2" -ne 0 ]
then
  echo "ABORT: phase 2 (msm topup) failed rc=$RC2 -- see run_glue_msm_topup_L4_claude.log"
  exit 1
fi

echo "################ L4 glue obs top-up DONE $(date) ################"
