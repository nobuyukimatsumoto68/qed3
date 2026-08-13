#!/bin/bash -l
# run_conn_s2_scc_claude.sh  (_scc, 2026-08-10, NM)
# =============================================================================================
# BU SCC SGE *BATCH* script for ONE stride-2 connected-Ylm unit = one (ensemble, offset) on one GPU.
# Submitted by run_wrapper_conn_s2_scc_claude.sh (one qsub per unit; resource requests on the qsub line).
# Measurement only -- reads ckpoint_lat.*, writes data_<ESNID>/corr_ylm_conn_t00_nhits1_s1/corr.<k>.h0.h5.
# RESUMABLE: the driver is "complete"-gated per config, so a wall-time kill / a re-run just skips finished
# configs (a stat, not a solve). Chained links (-hold_jid) walk a unit's whole k-range to completion.
#
# Parameters arrive as environment variables (qsub -v):
#   APP     conn binary (arch-specific: jj_conn_s2_L4_sm_70.o / _sm_80.o)
#   GSQ NF  ensemble id (--gsq/--Nf); NU0 (default 1.0)
#   ENSDIR  sea config dir basename (contains ckpoint_lat.*); driver adds the trailing slash
#   KMIN STRIDE KMAX   config range: --kmin <off+1> --stride 10 --kmax <last+1>
# =============================================================================================

#$ -P qfe
#$ -M nmatsum@bu.edu
#$ -j y
#$ -N conn_s2

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

LOG="conn_s2_scc_L4_Nf${NF}_g${GSQ}_kmin${KMIN}_$(date +%y%m%d%H%M).log"
{
  echo "### CONN-S2 $ENSDIR  Nf${NF} g${GSQ}  --kmin $KMIN --stride $STRIDE --kmax $KMAX  nhits1 t0 0 spin-dilution"
  echo "###   APP=$APP  OMP=$OMP_NUM_THREADS  [$(date +%F_%H:%M:%S)]"
  ./"$APP" --gsq "$GSQ" --Nf "$NF" --nu0 "$NU0" --ens-dir "${ENSDIR}/" \
    --kmin "$KMIN" --stride "$STRIDE" --kmax "$KMAX" --nhits 1 --t0 0 --spin-dilution
  echo "### CONN-S2 $ENSDIR kmin${KMIN} done (status ${PIPESTATUS[0]:-$?})  [$(date +%F_%H:%M:%S)]"
} 2>&1 | tee "$LOG"

echo "end date : $(date)"
echo "finished"
