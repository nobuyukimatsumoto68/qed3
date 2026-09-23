#!/bin/bash -l
# run_conn_s2_hit2_scc_claude.sh  (_scc, 2026-08-29, NM)  -- SECOND-HIT variant of run_conn_s2_scc_claude.sh
# =============================================================================================
# BU SCC SGE *BATCH* for ONE stride-2 connected-Ylm unit = one (ensemble, offset) on one GPU, run with
# --nhits 2 --outdir-nhits 1 to ADD the second stochastic hit (h1) on top of the existing h0.
#   - --outdir-nhits 1 writes h1 as corr.<k>.h1.h5 INTO THE EXISTING data_<ESNID>/corr_ylm_conn_t00_nhits1_s1/
#     dir (same dir as h0; NO separate nhits2 dir, NO hardlinks). The driver's per-hit "complete"-gate finds h0
#     already there and SKIPS it (a stat, no solve) -> only h1 is solved. This keeps the LOCAL (barracuda22)
#     pull-back working unchanged (its glob is corr_ylm_conn_t00_nhits1_s1/corr.*.h5, which now also grabs h1).
#   - h1 RNG is DETERMINISTIC, per the EXISTING per-hit seed convention: seed_from_string(esnid_k<k>_h1)
#     (the same scheme as h0's esnid_k<k>_h0, just hit index 1 -> an INDEPENDENT, reproducible draw). Matches h0's
#     determinism and stays poolable with FNAL/LOCAL.
#   - RESUMABLE / idempotent: per-config "complete"-gate -> a wall kill or a re-run skips finished (k,hit) pairs.
# Measurement only (reads ckpoint_lat.*). Parent plan: conn_stride2_three_site_impl_plan_claude.md ;
# SCC plan: scc_conn_s2_L4_impl_plan_claude.md.
#
# Parameters arrive as environment variables (qsub -v), IDENTICAL to run_conn_s2_scc_claude.sh:
#   APP  GSQ NF  NU0  ENSDIR  KMIN STRIDE KMAX
# =============================================================================================

#$ -P qfe
#$ -M nmatsum@bu.edu
#$ -j y
#$ -N conn_s2h2

set -u

echo "=========================================================="
echo "Start date : $(date)"
echo "Job name   : ${JOB_NAME:-?}   Job ID : ${JOB_ID:-?}"
echo "Host name  : ${HOSTNAME:-?}   NSLOTS : ${NSLOTS:-?}"
echo "CUDA_VISIBLE_DEVICES (SGE-assigned) : ${CUDA_VISIBLE_DEVICES:-unset}"
echo "=========================================================="

# ---- environment: cuda/gcc + repo Eigen (env.sh), and HDF5/GSL runtime libs (module) ----
source /projectnb/qfe/nmatsum/qed3/env.sh
module load hdf5/1.10.10 2>/dev/null
module load gsl 2>/dev/null

SRCDIR=/projectnb/qfe/nmatsum/qed3/src/production
cd "$SRCDIR" || { echo "ERROR: cannot cd $SRCDIR"; exit 1; }

NU0=${NU0:-1.0}
: "${APP:?APP (conn binary) must be set via qsub -v}"
: "${GSQ:?GSQ must be set}"
: "${NF:?NF must be set}"
: "${ENSDIR:?ENSDIR must be set}"
: "${KMIN:?KMIN must be set}"
: "${STRIDE:=10}"
: "${KMAX:?KMAX must be set}"

# host-side OpenMP (driver uses a few threads for gauge/sort); 1 job/GPU so give it the job's slots
export OMP_NUM_THREADS="${NSLOTS:-4}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-4}"

test -f "$APP" || { echo "ERROR: binary $APP not found (build via the wrapper first)"; exit 1; }
test -d "$ENSDIR" || { echo "ERROR: ens dir $ENSDIR not found"; exit 1; }

LOG="conn_s2h2_scc_L4_Nf${NF}_g${GSQ}_kmin${KMIN}_$(date +%y%m%d%H%M).log"
{
  echo "### CONN-S2 HIT2 $ENSDIR  Nf${NF} g${GSQ}  --kmin $KMIN --stride $STRIDE --kmax $KMAX  nhits2 outdir-nhits1 t0 0 spin-dilution"
  echo "###   APP=$APP  OMP=$OMP_NUM_THREADS  [$(date +%F_%H:%M:%S)]  (--outdir-nhits 1: h1 written INTO the nhits1 dir; gate skips h0)"
  ./"$APP" --gsq "$GSQ" --Nf "$NF" --nu0 "$NU0" --ens-dir "${ENSDIR}/" \
    --kmin "$KMIN" --stride "$STRIDE" --kmax "$KMAX" --nhits 2 --outdir-nhits 1 --t0 0 --spin-dilution
  echo "### CONN-S2 HIT2 $ENSDIR kmin${KMIN} done (status ${PIPESTATUS[0]:-$?})  [$(date +%F_%H:%M:%S)]"
} 2>&1 | tee "$LOG"

echo "end date : $(date)"
echo "finished"
