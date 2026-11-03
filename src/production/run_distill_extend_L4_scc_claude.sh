#!/bin/bash -l
# run_distill_extend_L4_scc_claude.sh -- BU SCC SGE batch script: ONE incremental basis-EXTENSION worker unit
#   (distill_peram_extend_claude.cu, Nv NV_OLD -> NV) = one k-offset of ONE ensemble on ONE GPU.  Copy of
#   run_distill_sym_L4_scc_claude.sh (untouched).  Submitted by run_wrapper_distill_extend_L4_scc_claude.sh.
#   Reads data_<ENSDIR>/distill_Nv<NV_OLD><OLD_SUFFIX>/peram.<k>.h5 (stored, read-only); writes
#   data_<ENSDIR>/distill_Nv<NV><OUT_SUFFIX>/peram.<k>.h5.  RESUMABLE (skip-if-present).  A k without a stored file is
#   skipped LOUDLY.  NO rm.  Never touches the stored dirs.
# qsub -v variables: APP GSQ NF AT ENSDIR KMIN STRIDE KMAX NV NV_OLD [NU0=1.0] [TSRC0=0] [TWIN=32] [NSRC=2]
#                    [OLD_SUFFIX=_sym] [OUT_SUFFIX=_sym_ext]
#$ -P qfe
#$ -M nmatsum@bu.edu
#$ -j y
#$ -N dext
set -u
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name   : ${JOB_NAME:-?}   Job ID : ${JOB_ID:-?}   Host : ${HOSTNAME:-?}   NSLOTS : ${NSLOTS:-?}"
echo "CUDA_VISIBLE_DEVICES (SGE) : ${CUDA_VISIBLE_DEVICES:-unset}"
echo "=========================================================="
source /projectnb/qfe/nmatsum/qed3/env.sh
module load hdf5/1.10.10 2>/dev/null
module load gsl 2>/dev/null
SRCDIR=/projectnb/qfe/nmatsum/qed3/src/production
cd "$SRCDIR" || { echo "ERROR: cannot cd $SRCDIR"; exit 1; }
: "${APP:?APP (distill binary) must be set}"
: "${GSQ:?}"; : "${NF:?}"; : "${AT:?}"; : "${ENSDIR:?}"; : "${KMIN:?}"; : "${STRIDE:?}"; : "${KMAX:?}"; : "${NV:?}"; : "${NV_OLD:?}"
NU0=${NU0:-1.0}
TSRC0=${TSRC0:-0}
TWIN=${TWIN:-32}
NSRC=${NSRC:-2}
OLD_SUFFIX=${OLD_SUFFIX:-_sym}
OUT_SUFFIX=${OUT_SUFFIX:-_sym_ext}
export OMP_NUM_THREADS=${NSLOTS:-4}
echo "### dext $ENSDIR  kmin=$KMIN stride=$STRIDE kmax=$KMAX  Nv=$NV_OLD->$NV nsrc=$NSRC twin=$TWIN old$OLD_SUFFIX out$OUT_SUFFIX  app=$APP  [$(date +%F_%H:%M:%S)] ###"
./"$APP" --gsq "$GSQ" --Nf "$NF" --nu0 "$NU0" --at "$AT" \
  --nv "$NV" --nv-old "$NV_OLD" --old-suffix "$OLD_SUFFIX" \
  --tsrc0 "$TSRC0" --twin "$TWIN" --nsrc "$NSRC" --fused-sources --out-suffix "$OUT_SUFFIX" \
  --ens-dir "$ENSDIR/" --kmin "$KMIN" --stride "$STRIDE" --kmax "$KMAX"
rc=$?
echo "### dext unit kmin=$KMIN exit status $rc  [$(date +%F_%H:%M:%S)] ###"
echo "### perams now in data_${ENSDIR}/distill_Nv${NV}${OUT_SUFFIX}: $(ls data_${ENSDIR}/distill_Nv${NV}${OUT_SUFFIX}/peram.*.h5 2>/dev/null | wc -l) ###"
exit $rc
