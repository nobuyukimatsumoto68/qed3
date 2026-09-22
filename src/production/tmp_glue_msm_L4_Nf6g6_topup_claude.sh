#!/usr/bin/env bash
# tmp_glue_msm_L4_Nf6g6_topup_claude.sh -- HANDOFF (NM runs this).
#
# GOAL: fill the 21 missing glue_msm_shapes h5 on the ONE incomplete L4 ensemble Nf6 gsq6.0
#       (currently 778/799; the other 8 L4 ensembles are already at the full 799).
#       glue_msm_shapes = a7's linear-F (F l=1 + Fl2) glueball driver -- the ONLY glue obs needed on L4.
#       glue_f2_v2_shapes is DELIBERATELY SKIPPED on L4 (NM: F^2/F^4 0++ stays L1+L2, no clean L4 plateau).
#
# METHOD: re-run the existing complete-gated sweep run_glue_msm_topup_claude.sh scoped to L4 at=0.2.
#       Complete-gate => the 8 already-done L4 ensembles skip in seconds (no compute, no memory);
#       only Nf6 g6.0 computes its 21 missing configs (glue2_msm_shapes_L4_claude.o already built).
#
# RESOURCE SIZING (per NM -- a memory-spiking job recently OOM-killed a peer; do NOT repeat):
#       NWORK=1 (a SINGLE worker -> at most ONE L4 flow process resident at a time) and the sweep's
#       own OMP_NUM_THREADS=1 => 1 CPU thread, minimal bounded memory, no parallel spike. CPU-only,
#       no GPU contention with the running conn on GPU0.
#
# Run:  bash tmp_glue_msm_L4_Nf6g6_topup_claude.sh
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1

LOG=tmp_glue_msm_L4_Nf6g6_topup_claude.log
echo "### L4 glue_msm topup (Nf6 g6.0, 21 cfg) START $(date) ###" | tee "$LOG"
echo "### NWORK=1  L_ONLY=4  AT_ONLY=0.200000  (complete-gated; only Nf6 g6.0 computes) ###" | tee -a "$LOG"

NWORK=1 L_ONLY=4 AT_ONLY=0.200000 bash run_glue_msm_topup_claude.sh >> "$LOG" 2>&1
rc=$?

echo "### glue_msm L4 topup done (status $rc) $(date) ###" | tee -a "$LOG"
if [ "$rc" -ne 0 ]
then
  echo "### NONZERO EXIT -- check $LOG and run_glue_msm_topup_claude.log ###" | tee -a "$LOG"
  exit "$rc"
fi
# quick verify
n=$(find data_Nf6_gsq6.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000 -maxdepth 1 -name 'glue_msm_shapes.*.h5' 2>/dev/null | wc -l)
echo "### Nf6 g6.0 L4 glue_msm_shapes now: $n / 799 ###" | tee -a "$LOG"
