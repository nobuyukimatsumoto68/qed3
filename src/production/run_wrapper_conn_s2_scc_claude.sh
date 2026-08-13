#!/bin/bash -l
# run_wrapper_conn_s2_scc_claude.sh  (_scc, 2026-08-10, NM)
# =============================================================================================
# BU SCC login-node WRAPPER for SCC's L4 slice of the stride-2 connected-Ylm completion.
#   (1) BUILDS the L4 conn binary for sm_70 (V100) and sm_80 (A100),
#   (2) reads SCC's L4 units (ensemble, offset) from conn_stride2_assign_claude.txt,
#   (3) computes a LIVE --kmax = (current ckpoint_lat max)+1 per ensemble (configs still growing),
#   (4) round-robins the 12 units across the arches, PINS strong-FP64 GPUs (gpu_type sm_70->V100,
#       sm_80->A100; the conn solve is double-precision -> OFF weak-FP64 L40S, same as the HMC), and
#   (5) submits a dependent CHAIN per unit (run_conn_s2_scc_claude.sh); complete-gating makes each 12h
#       link resume where the last stopped, walking the unit's k-range to completion.
# Measurement only: reads configs, no chain fork hazard with the LIVE L4 HMC on these same 3 ensembles.
# Plan: scc_conn_s2_L4_impl_plan_claude.md ; parent: conn_stride2_three_site_impl_plan_claude.md.
#
# Usage:  bash run_wrapper_conn_s2_scc_claude.sh            # build + submit
#         DRYRUN=1 bash run_wrapper_conn_s2_scc_claude.sh   # print qsub lines, submit nothing
#         NOBUILD=1 bash run_wrapper_conn_s2_scc_claude.sh  # skip build (binaries present)
#         FORCE_BUILD=1 bash ...                            # rebuild even if present
# =============================================================================================
set -u

SRCDIR=/projectnb/qfe/nmatsum/qed3/src/production
cd "$SRCDIR" || { echo "ERROR: cannot cd $SRCDIR"; exit 1; }
source /projectnb/qfe/nmatsum/qed3/env.sh          # cuda/12.8, gcc, $QED3_INC (repo Eigen)
module load hdf5/1.10.10 2>/dev/null
module load gsl 2>/dev/null

ASSIGN=${ASSIGN:-conn_stride2_assign_claude.txt}
SRC_CONN=jj_local_ylm_scalar_conn_stoch_claude.cu
LREF=4                               # SCC does L4 only (n4 geometry complete; n1/n2/n3 gap -> LOCAL/FNAL)

PE_OMP=${PE_OMP:-4}
NPACK=1                              # 1 job/GPU, no MPS
H_RT_FIRST=${H_RT_FIRST:-12:00:00}
H_RT=${H_RT:-12:00:00}
N_CHAIN=${N_CHAIN:-4}                # links/unit (complete-gating: each resumes; over-provision is free)
ARCH_LIST=${ARCH_LIST:-"sm_70 sm_80"}
SUBMIT_ARCHS=${SUBMIT_ARCHS:-"$ARCH_LIST"}
GPUT_SM70=${GPUT_SM70:-V100}
GPUT_SM80=${GPUT_SM80:-A100}
DRYRUN=${DRYRUN:-0}
NOBUILD=${NOBUILD:-0}
FORCE_BUILD=${FORCE_BUILD:-0}

# ---- build recipe (per arch) --------------------------------------------------------------------
NVCC=nvcc
NVCCBASE="-g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ ${QED3_INC} -I/projectnb/qfe/nmatsum/opt/highfive/include -I${SCC_HDF5_INCLUDE} -I${SCC_GSL_INCLUDE}"
LDFLAGS="-L${SCC_HDF5_LIB} -L${SCC_GSL_LIB} -lhdf5 -lgsl -lgslcblas -lm"

binname () { echo "jj_conn_s2_L${LREF}_${1}.out"; }   # $1 = arch

gpuc_of () {
  case "$1" in
    sm_70) echo "7.0" ;; sm_80) echo "8.0" ;; sm_90) echo "9.0" ;; *) echo "7.0" ;;
  esac
}
gput_of () {
  case "$1" in
    sm_70) echo "$GPUT_SM70" ;; sm_80) echo "$GPUT_SM80" ;; *) echo "" ;;
  esac
}

build_all () {
  local arch out
  for arch in $ARCH_LIST
  do
    out=$(binname "$arch")
    if [ "$FORCE_BUILD" -ne 1 ] && [ -f "$out" ]
    then
      echo "  $out present -> skip (FORCE_BUILD=1 to rebuild)"
      continue
    fi
    echo "===== build $out (N_REFINE_CLI=${LREF} arch=$arch) [$(date +%F_%H:%M:%S)] ====="
    $NVCC -arch="$arch" $NVCCBASE -DN_REFINE_CLI="$LREF" $INCLUDES $LDFLAGS \
      "$SRC_CONN" -o "$out" 2>&1 | tee "build_conn_s2_L${LREF}_${arch}_claude.log"
    test -f "$out" || { echo "BUILD FAILED: $out"; exit 1; }
  done
}
if [ "$NOBUILD" -eq 1 ]
then
  echo "NOBUILD=1 -- skipping build"
else
  build_all
fi

# ---- SCC L4 units from the assignment table -----------------------------------------------------
# cols: 1=ensemble 2=L 3=Nf 4=gsq 5=at 6=offset 7=kmin 8=stride 9=kmax 10=cfg 11=site 12=worker_h
[ -f "$ASSIGN" ] || { echo "ERROR: assignment table $ASSIGN not found"; exit 1; }
mapfile -t UNITS < <(awk '$0 !~ /^#/ && toupper($11)=="SCC" && $2==4 {print $1" "$3" "$4" "$7}' "$ASSIGN")
echo "===== SCC L4 units: ${#UNITS[@]} (want 12) ====="
[ "${#UNITS[@]}" -eq 0 ] && { echo "ERROR: no SCC L4 units parsed from $ASSIGN"; exit 1; }

read -r -a ARCHES <<< "$SUBMIT_ARCHS"
n_arch=${#ARCHES[@]}

# live kmax = (max ckpoint_lat)+1 for an ensemble dir (configs still growing)
live_kmax () {   # ensdir
  local k
  k=$(ls "$1"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n | tail -1)
  echo $(( ${k:-0} + 1 ))
}

# cross-run anchoring within the conn_s2 namespace (^cs2) -> re-run EXTENDS each unit's chain (new kmax)
EXIST_JOBS=$(qstat -u "$USER" -r 2>/dev/null | awk '
  /^[0-9]+ /{ jid=$1 }
  /Full jobname:/{ if($3 ~ /^cs2/) print jid"\t"$3 }
')
existing_tail () {   # token
  local a=${1//./\\.}
  printf '%s\n' "$EXIST_JOBS" | awk -F'\t' -v pat="(${a})(__|_c)" '$2 ~ pat {print $1}' | sort -n | tail -1
}

hrt_to_sec () { local h m s; IFS=':' read -r h m s <<< "$1"; echo $(( 10#$h*3600 + 10#$m*60 + 10#$s )); }

submit_unit () {   # arch  "ensdir nf gsq kmin"
  local arch=$1 rec=$2
  local ensdir nf gsq kmin
  read -r ensdir nf gsq kmin <<< "$rec"
  local app gpuc gput kmax name tok prev c hrt jid
  app=$(binname "$arch")
  gpuc=$(gpuc_of "$arch")
  gput=$(gput_of "$arch")
  kmax=$(live_kmax "$ensdir")
  # token is ARCH-INDEPENDENT (no cs2/arch prefix) so it matches the "_Nf<nf>g<gsq>k<kmin>" substring of the
  # job name REGARDLESS of arch -> a re-run (e.g. switching sm_80->L40S) anchors the new chain behind the
  # current one and EXTENDS it, instead of racing a concurrent writer. EXIST_JOBS is already ^cs2-filtered.
  tok="Nf${nf}g${gsq}k${kmin}"
  name="cs2${arch#sm_}_Nf${nf}g${gsq}k${kmin}"
  local vars="APP=${app},GSQ=${gsq},NF=${nf},NU0=1.0,ENSDIR=${ensdir},KMIN=${kmin},STRIDE=10,KMAX=${kmax}"
  prev=$(existing_tail "$tok")
  [ -n "$prev" ] && echo "  [anchor] ${name}: existing chain -> hold on tail jid=$prev"
  for (( c=0; c<N_CHAIN; c++ ))
  do
    if [ "$c" -eq 0 ]; then hrt=$H_RT_FIRST; else hrt=$H_RT; fi
    local qsub_cmd=( qsub -terse -N "${name}_c${c}"
                     -l "gpus=1" -l "gpu_c=${gpuc}" -l "h_rt=${hrt}" -pe omp "${PE_OMP}" )
    [ -n "$gput" ] && qsub_cmd+=( -l "gpu_type=${gput}" )
    [ -n "$prev" ] && qsub_cmd+=( -hold_jid "$prev" )
    qsub_cmd+=( -v "${vars}" run_conn_s2_scc_claude.sh )
    echo "+ ${qsub_cmd[*]}"
    if [ "$DRYRUN" -eq 0 ]
    then
      jid=$( "${qsub_cmd[@]}" | tr -d '[:space:]' )
      echo "  -> submitted jid=$jid (link $c, h_rt=$hrt, kmax=$kmax)"
      prev=$jid
    else
      prev="<jid_c${c}>"
    fi
  done
}

# round-robin the 12 units across the submit arches
echo "===== submit ${#UNITS[@]} units x $N_CHAIN links (round-robin ${SUBMIT_ARCHS}, NPACK=1) ====="
ri=0
for rec in "${UNITS[@]}"
do
  a=${ARCHES[$(( ri % n_arch ))]}
  submit_unit "$a" "$rec"
  ri=$(( ri + 1 ))
done

echo "===== conn-s2 wrapper done (DRYRUN=$DRYRUN) [$(date +%F_%H:%M:%S)] ====="
