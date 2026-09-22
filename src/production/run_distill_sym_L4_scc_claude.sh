#!/bin/bash -l
# run_distill_sym_L4_scc_claude.sh -- BU SCC SGE batch script: ONE symmetrized-basis distillation worker unit =
#   one k-offset of ONE ensemble on ONE GPU.  Submitted by run_wrapper_distill_sym_L4_scc_claude.sh (qsub -v ...).
#   Output: data_<ENSDIR>/distill_Nv<NV>_sym/peram.<k>.h5.  RESUMABLE: the binary skips k whose peram.<k>.h5 already
#   exists, so a wall-time kill / re-run / chained -hold_jid link just continues.  NO rm.  Never touches _v2 dirs.
# qsub -v variables: APP GSQ NF AT ENSDIR KMIN STRIDE KMAX NV [NU0=1.0] [TSRC0=0] [TWIN=32] [NSRC=2] [OUT_SUFFIX=_sym]
#$ -P qfe
#$ -M nmatsum@bu.edu
#$ -j y
#$ -N dsym
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
: "${GSQ:?}"; : "${NF:?}"; : "${AT:?}"; : "${ENSDIR:?}"; : "${KMIN:?}"; : "${STRIDE:?}"; : "${KMAX:?}"; : "${NV:?}"
NU0=${NU0:-1.0}
TSRC0=${TSRC0:-0}
TWIN=${TWIN:-32}
NSRC=${NSRC:-2}
OUT_SUFFIX=${OUT_SUFFIX:-_sym}
export OMP_NUM_THREADS=${NSLOTS:-4}
echo "### dsym $ENSDIR  kmin=$KMIN stride=$STRIDE kmax=$KMAX  Nv=$NV nsrc=$NSRC twin=$TWIN out$OUT_SUFFIX  app=$APP  [$(date +%F_%H:%M:%S)] ###"
./"$APP" --gsq "$GSQ" --Nf "$NF" --nu0 "$NU0" --at "$AT" \
  --nv "$NV" --tsrc0 "$TSRC0" --twin "$TWIN" --nsrc "$NSRC" --fused-sources --out-suffix "$OUT_SUFFIX" \
  --ens-dir "$ENSDIR/" --kmin "$KMIN" --stride "$STRIDE" --kmax "$KMAX"
rc=$?
echo "### dsym unit kmin=$KMIN exit status $rc  [$(date +%F_%H:%M:%S)] ###"
echo "### perams now in data_${ENSDIR}/distill_Nv${NV}${OUT_SUFFIX}: $(ls data_${ENSDIR}/distill_Nv${NV}${OUT_SUFFIX}/peram.*.h5 2>/dev/null | wc -l) ###"
exit $rc
