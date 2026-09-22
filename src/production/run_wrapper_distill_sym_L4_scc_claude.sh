#!/bin/bash -l
# run_wrapper_distill_sym_L4_scc_claude.sh -- BU SCC login-node WRAPPER: symmetrized-basis (BASIS_SYM=1) distillation
#   perambulators for L4 Nf2 gsq2.0 (at0.2, massless).  (1) builds the binary (both arches), (2) splits the stride grid
#   into NUNITS staggered k-offset units (unit u: kmin=first+u*STRIDE, stride=STRIDE*NUNITS -> disjoint k sets, no
#   overwrite race), (3) computes a LIVE --kmax from the current ckpoint_lat max (HMC may still be running),
#   (4) submits N_CHAIN dependent (-hold_jid) links per unit; the binary's skip-if-present makes each link resume.
#   Pins strong-FP64 GPUs (V100 / A100) like the conn/HMC jobs.  Measurement only.  NO rm.  NO kill.
# Usage:  bash run_wrapper_distill_sym_L4_scc_claude.sh             # build + submit
#         DRYRUN=1 bash run_wrapper_distill_sym_L4_scc_claude.sh    # print qsub lines only
#         NOBUILD=1 ...  (binaries present)     ARCH=sm_70|sm_80 (default sm_70=V100)
set -u
SRCDIR=/projectnb/qfe/nmatsum/qed3/src/production
cd "$SRCDIR" || { echo "ERROR: cannot cd $SRCDIR"; exit 1; }
source /projectnb/qfe/nmatsum/qed3/env.sh
module load hdf5/1.10.10 2>/dev/null
module load gsl 2>/dev/null

ENSDIR=${ENSDIR:-Nf2_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000}
NF=2
GSQ=2.0
AT=0.2
NV=${NV:-24}
STRIDE=${STRIDE:-2}        # L4 has ~800 configs -> stride 2 = ~400 perams
NUNITS=${NUNITS:-4}        # parallel GPU units (staggered offsets)
N_CHAIN=${N_CHAIN:-4}      # 12h links per unit; skip-if-present resumes
H_RT=${H_RT:-12:00:00}
PE_OMP=${PE_OMP:-4}
ARCH=${ARCH:-sm_70}
DRYRUN=${DRYRUN:-0}
NOBUILD=${NOBUILD:-0}
case "$ARCH" in
  sm_70) GPUC=7.0; GPUT=V100 ;;
  sm_80) GPUC=8.0; GPUT=A100 ;;
  *) echo "unknown ARCH $ARCH"; exit 1 ;;
esac
APP=distill_peram_mrhs_L4_Nv${NV}_sym_${ARCH}.out

if [ "$NOBUILD" -eq 0 ] || [ ! -f "$APP" ]; then
  NV=$NV bash tmp_build_distill_sym_L4_scc_claude.sh 2>&1 | tee tmp_build_distill_sym_L4_scc_claude.log
  [ -f "$APP" ] || { echo "ERROR: $APP missing after build"; exit 1; }
fi
[ -d "$ENSDIR" ] || { echo "ERROR: no ensemble dir $ENSDIR"; exit 1; }
first=$(ls "$ENSDIR"/ckpoint_lat.* | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n | head -1)
last=$( ls "$ENSDIR"/ckpoint_lat.* | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n | tail -1)
KMAX=$(( last + 1 ))
echo "### $ENSDIR : ckpoints $first..$last -> kmax=$KMAX ; Nv=$NV stride=$STRIDE units=$NUNITS chain=$N_CHAIN arch=$ARCH ($GPUT) ###"
for (( u=0; u<NUNITS; u++ )); do
  kmin=$(( first + u * STRIDE ))
  wstride=$(( STRIDE * NUNITS ))
  name="dsym_L4g${GSQ}_u${u}"
  vars="APP=${APP},GSQ=${GSQ},NF=${NF},AT=${AT},ENSDIR=${ENSDIR},KMIN=${kmin},STRIDE=${wstride},KMAX=${KMAX},NV=${NV}"
  prev=""
  for (( c=0; c<N_CHAIN; c++ )); do
    qsub_cmd=( qsub -terse -N "${name}_c${c}" -l gpus=1 -l "gpu_c=${GPUC}" -l "gpu_type=${GPUT}" -l "h_rt=${H_RT}" -pe omp "${PE_OMP}" )
    [ -n "$prev" ] && qsub_cmd+=( -hold_jid "$prev" )
    qsub_cmd+=( -v "${vars}" run_distill_sym_L4_scc_claude.sh )
    echo "+ ${qsub_cmd[*]}"
    if [ "$DRYRUN" -eq 0 ]; then
      jid=$( "${qsub_cmd[@]}" | tr -d '[:space:]' )
      echo "  -> submitted jid=$jid (unit $u link $c)"
      prev=$jid
    else
      prev="<jid_u${u}_c${c}>"
    fi
  done
done
echo "### submitted $NUNITS units x $N_CHAIN links.  Monitor: qstat -u nmatsum ; ls data_${ENSDIR}/distill_Nv${NV}_sym | wc -l ###"
