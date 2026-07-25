#!/bin/bash -l
# run_wrapper_L4_halfat_scc_claude.sh  (_scc, 2026-07-22, NM)
# =============================================================================================
# BU SCC *WRAPPER* for the HALF-a_t L=4 overlap-fermion QED3 probe. Sibling of
# run_wrapper_L4_scc_claude.sh (the full-a_t, a_t=0.2 production wrapper); this one drives the
# a_t = 0.1 set with N_t FIXED at 128 (so the physical time extent N_t*a_t halves, 25.6 -> 12.8).
#
# Campaign = "middle gsq per L, all Nf=2,4,6, half a_t". SCC's slice = Nf2 ONLY at L4 middle gsq
# (gsq=4.0); Nf4/6 half-a_t are handled at FNAL (mirrors the full-a_t L4 split). So the default
# submission here is a SINGLE ensemble: Nf2, L4, gsq4.0, at=0.1, massless, KMAX=200 (= half the
# full-a_t massless target of 400). Design doc: halfat_L4_scc_impl_plan_claude.md.
#
# HMC parameters are REUSED from the a_t=0.2 tuning AS-IS (NM: "I don't believe they'll be affected
# so much"): ladder {0.0,0.4,1.0} steps {4,4,4} tau=1.0 gauge MG=100, Zolotarev frozen (0.008,5)
# n_action=31 n_force=11, -DMIXED_FORCE. The frozen Zolotarev window is REUSED (a_t shift judged
# negligible); startup should be watched for "eval outside window" warnings (smoke-validate only).
#
# a_t is a COMPILE-TIME knob in the driver: -DAT_VAL=0.1 (default 0.2 preserves the full-a_t builds).
# The driver embeds at into the output dir name (at0.100000...), distinct from the at0.200000 dirs --
# no checkpoint collision, auto-resume intact.
#
# SEPARATION from the full-a_t run (CORRECTNESS): the full-a_t massless run already has an ensemble
# token Nf2g4.0m0.0 (job L4<arch>_Nf2g4...). This wrapper uses a DISTINCT job-name namespace "L4h<arch>"
# and scopes cross-run anchoring to ^L4h, so a half-a_t chain NEVER -hold_jid's onto the full-a_t one.
# Binaries carry an at tag (hmc_L4_scc_Nf<nf>_k<kmax>_at<AT_VAL>_<arch>.out); batch logs carry AT_TAG.
#
# Usage:  bash run_wrapper_L4_halfat_scc_claude.sh            # build + submit
#         DRYRUN=1 bash run_wrapper_L4_halfat_scc_claude.sh   # print qsub commands, do NOT submit
#         NOBUILD=1 bash run_wrapper_L4_halfat_scc_claude.sh  # skip build (binaries already present)
# Refs: M. Hasenbusch hep-lat/0107019; B. Jegerlehner hep-lat/9612014; A. D. Kennedy hep-lat/0402038.
# =============================================================================================
set -u

SRCDIR=/projectnb/qfe/nmatsum/qed3/src/production
cd "$SRCDIR" || { echo "ERROR: cannot cd $SRCDIR"; exit 1; }
source /projectnb/qfe/nmatsum/qed3/env.sh          # unified env: cuda/12.8, gcc/13.2.0, $QED3_INC (repo Eigen)

# ---- temporal spacing knob ----------------------------------------------------------------------
AT_VAL=${AT_VAL:-0.1}                # a_t (half). Passed as -DAT_VAL and embedded in binary/log tags.
AT_TAG="at${AT_VAL}"                 # campaign tag for binary names + batch log names

# ---- ensemble set (EDIT THESE) ------------------------------------------------------------------
# SCC slice = Nf2 only at L4 middle gsq=4.0, massless, half KMAX. (Nf4/6 half-a_t = FNAL.)
WHICH=${WHICH:-massless}             # massless | massive | both

ML_NF=${ML_NF:-"2"}                  # MASSLESS: SCC does Nf2 only
ML_GSQ=${ML_GSQ:-"4.0"}             # middle gsq at L4
ML_MASS=${ML_MASS:-"0.0"}
ML_KMAX=${ML_KMAX:-200}             # half the full-a_t massless target (400/2)

MV_NF=${MV_NF:-"2"}                  # MASSIVE: not part of the half-a_t probe by default (kept for symmetry)
MV_GSQ=${MV_GSQ:-"4.0"}
MV_MASS=${MV_MASS:-"0.2 0.3 0.4"}
MV_KMAX=${MV_KMAX:-30}              # half the full-a_t massive target (60/2)

NU0=${NU0:-1.0}
KRNG=${KRNG:-4}                      # keep every KRNG-th rng checkpoint

# ---- GPU architectures --------------------------------------------------------------------------
ARCH_LIST=${ARCH_LIST:-"sm_70 sm_80"}      # build these (V100=sm_70, A100=sm_80)
SUBMIT_ARCHS=${SUBMIT_ARCHS:-"$ARCH_LIST"} # submit to these node pools (round-robin)
PE_OMP=${PE_OMP:-16}                       # CPU slots per job
NPACK=${NPACK:-1}                          # 1 = one ensemble/GPU, NO MPS (per NM 2026-07-18). NPACK=2 re-enables MPS.

# ---- job chaining -------------------------------------------------------------------------------
H_RT_FIRST=${H_RT_FIRST:-12:00:00}
H_RT=${H_RT:-12:00:00}
N_CHAIN=${N_CHAIN:-4}
BUFFER_SEC=${BUFFER_SEC:-900}

# HH:MM:SS -> seconds
hrt_to_sec () {
  local h m s
  IFS=':' read -r h m s <<< "$1"
  echo $(( 10#$h*3600 + 10#$m*60 + 10#$s ))
}

DRYRUN=${DRYRUN:-0}
NOBUILD=${NOBUILD:-0}

SRC=hmc_hasenbusch_block_scc_claude.cu
NVCC=nvcc
NVCCBASE="-g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ ${QED3_INC}"
LDFLAGS="-lm"

# gpu_c (SGE minimum compute capability) per arch
gpuc_of () {
  case "$1" in
    sm_70) echo "7.0" ;;
    sm_75) echo "7.5" ;;
    sm_80) echo "8.0" ;;
    sm_86) echo "8.6" ;;
    *) echo "7.0" ;;
  esac
}

# L4 MD steps per Nf: Nf2 -> {4,4,4} (default 4); Nf4/6 -> {5,5,5} (-DL4_MDSTEP=5), same as full-a_t.
mdstep_of () {
  if [ "$1" -eq 2 ]; then echo 4; else echo 5; fi
}

# binary name per (Nf, KMAX, arch) -- carries the at tag so half-a_t binaries never overwrite the
# full-a_t ones (hmc_L4_scc_Nf<nf>_k<kmax>_<arch>.out from the sibling wrapper).
binname () { echo "hmc_L4_scc_Nf${1}_k${2}_${AT_VAL}_${3}.out"; }

# ---- selected ensemble SET(S): each = "name|nf_list|gsq_list|mass_list|kmax" --------------------
SETS=()
case "$WHICH" in
  massless) SETS+=( "massless|${ML_NF}|${ML_GSQ}|${ML_MASS}|${ML_KMAX}" ) ;;
  massive)  SETS+=( "massive|${MV_NF}|${MV_GSQ}|${MV_MASS}|${MV_KMAX}" ) ;;
  both)     SETS+=( "massless|${ML_NF}|${ML_GSQ}|${ML_MASS}|${ML_KMAX}" )
            SETS+=( "massive|${MV_NF}|${MV_GSQ}|${MV_MASS}|${MV_KMAX}" ) ;;
  *) echo "ERROR: WHICH must be massless|massive|both"; exit 1 ;;
esac
read -r -a ARCHES <<< "$SUBMIT_ARCHS"
n_arch=${#ARCHES[@]}

# unique (Nf, KMAX) build targets across all selected sets
declare -A NKSEEN
for s in "${SETS[@]}"
do
  IFS='|' read -r sname nfs gsqs masses kmax <<< "$s"
  for nf in $nfs
  do
    NKSEEN["$nf $kmax"]=1
  done
done

# ---- 1. build one binary per (Nf, KMAX) per arch (with -DAT_VAL and per-Nf -DL4_MDSTEP) ----------
build_all () {
  local arch nf kmax out key mds
  for arch in $ARCH_LIST
  do
    for key in "${!NKSEEN[@]}"
    do
      read -r nf kmax <<< "$key"
      out=$(binname "$nf" "$kmax" "$arch")
      mds=$(mdstep_of "$nf")
      echo "===== build $out (LREF=4 NF=$nf KMAX=$kmax AT_VAL=$AT_VAL L4_MDSTEP=$mds arch=$arch -DMIXED_FORCE) [$(date +%F_%H:%M:%S)] ====="
      $NVCC -arch="$arch" $NVCCBASE $INCLUDES \
        -DLREF=4 -DNF="$nf" -DKMAX="$kmax" -DKRNG="$KRNG" -DAT_VAL="$AT_VAL" -DL4_MDSTEP="$mds" -DMIXED_FORCE \
        "$SRC" $LDFLAGS -o "$out" 2>&1 | tee "build_L4halfat_scc_Nf${nf}_k${kmax}_${arch}_claude.log"
      test -f "$out" || { echo "BUILD FAILED: $out"; exit 1; }
    done
  done
}
if [ "$NOBUILD" -eq 1 ]
then
  echo "NOBUILD=1 -- skipping build"
else
  build_all
fi

# ---- dependency safety ACROSS wrapper re-runs -- SCOPED TO THE HALF-a_t NAMESPACE (^L4h) ---------
# Snapshot existing HALF-a_t jobs only (Full jobname starts L4h). Anchoring onto a full-a_t L4<arch>
# chain would be WRONG (different at, different dir) -- the ^L4h filter prevents that cross-campaign match.
EXIST_JOBS=$(qstat -u "$USER" -r 2>/dev/null | awk '
  /^[0-9]+ /{ jid=$1 }
  /Full jobname:/{ if($3 ~ /^L4h/) print jid"\t"$3 }
')
existing_tail () {   # tokA [tokB]
  local a=${1//./\\.} b=${2:-} pat
  pat="(${a})(__|_c)"
  if [ -n "$b" ]
  then
    b=${b//./\\.}
    pat="(${a}|${b})(__|_c)"
  fi
  printf '%s\n' "$EXIST_JOBS" | awk -F'\t' -v pat="$pat" '$2 ~ pat {print $1}' | sort -n | tail -1
}

# ---- pair 2/GPU, submit a DEPENDENT CHAIN per pair ----------------------------------------------
submit_chain () {   # arch  A="Nf|KMAX|gsq|mass"  [B=...]
  local arch=$1 A=$2 B=${3:-}
  local gpuc
  gpuc=$(gpuc_of "$arch")
  local nfA kmaxA gsqA massA appA
  IFS='|' read -r nfA kmaxA gsqA massA <<< "$A"
  appA=$(binname "$nfA" "$kmaxA" "$arch")
  local vars="NU0=${NU0},AT_TAG=${AT_TAG},APPA=${appA},GSQA=${gsqA},NFA=${nfA},MASSA=${massA}"
  local name="L4h${arch#sm_}_Nf${nfA}g${gsqA}m${massA}"
  if [ -n "$B" ]
  then
    local nfB kmaxB gsqB massB appB
    IFS='|' read -r nfB kmaxB gsqB massB <<< "$B"
    appB=$(binname "$nfB" "$kmaxB" "$arch")
    vars="${vars},APPB=${appB},GSQB=${gsqB},NFB=${nfB},MASSB=${massB}"
    name="${name}__Nf${nfB}g${gsqB}m${massB}"
  fi

  local tokA="Nf${nfA}g${gsqA}m${massA}" tokB=""
  [ -n "$B" ] && tokB="Nf${nfB}g${gsqB}m${massB}"
  local prev
  prev=$(existing_tail "$tokA" "$tokB")
  [ -n "$prev" ] && echo "  [anchor] ${name}: existing HALF-a_t chain found -> new links hold on tail jid=$prev (no concurrent writer)"
  local c hrt maxsec jid
  for (( c=0; c<N_CHAIN; c++ ))
  do
    if [ "$c" -eq 0 ]; then hrt=$H_RT_FIRST; else hrt=$H_RT; fi
    maxsec=$(( $(hrt_to_sec "$hrt") - BUFFER_SEC ))
    [ "$maxsec" -lt 60 ] && maxsec=60
    local qsub_cmd=( qsub -terse -N "${name}_c${c}"
                     -l "gpus=1" -l "gpu_c=${gpuc}" -l "h_rt=${hrt}" -pe omp "${PE_OMP}" )
    [ -n "$prev" ] && qsub_cmd+=( -hold_jid "$prev" )
    qsub_cmd+=( -v "${vars},MAX_SEC=${maxsec}" run_L4_scc_claude.sh )
    echo "+ ${qsub_cmd[*]}"
    if [ "$DRYRUN" -eq 0 ]
    then
      jid=$( "${qsub_cmd[@]}" | tr -d '[:space:]' )
      echo "  -> submitted jid=$jid (link $c, h_rt=$hrt, max_sec=$maxsec)"
      prev=$jid
    else
      prev="<jid_c${c}>"
    fi
  done
}

# ---- 3. per SET: round-robin ITS ensembles across arches, pair 2/GPU WITHIN the set, submit chains
for s in "${SETS[@]}"
do
  IFS='|' read -r sname nfs gsqs masses kmax <<< "$s"
  recs=()
  for nf in $nfs
  do
    for gsq in $gsqs
    do
      for m in $masses
      do
        recs+=( "${nf}|${kmax}|${gsq}|${m}" )
      done
    done
  done
  [ "${#recs[@]}" -eq 0 ] && continue
  echo "===== set $sname [half-a_t $AT_TAG] : ${#recs[@]} ensembles (Nf {$nfs} x gsq {$gsqs} x m {$masses}, KMAX=$kmax) ====="

  unset SB
  declare -A SB
  ri=0
  for rec in "${recs[@]}"
  do
    a=${ARCHES[$(( ri % n_arch ))]}
    SB[$a]+="${rec}"$'\n'
    ri=$(( ri + 1 ))
  done

  for arch in "${ARCHES[@]}"
  do
    mapfile -t r < <(printf "%s" "${SB[$arch]:-}" | sed '/^$/d')
    [ "${#r[@]}" -eq 0 ] && continue
    echo "--- $sname / $arch : ${#r[@]} ensembles -> $(( (${#r[@]}+NPACK-1)/NPACK )) chains x $N_CHAIN links (NPACK=$NPACK/GPU) ---"
    i=0
    while [ "$i" -lt "${#r[@]}" ]
    do
      if [ "$NPACK" -ge 2 ]
      then
        submit_chain "$arch" "${r[$i]}" "${r[$((i+1))]:-}"
        i=$(( i + 2 ))
      else
        submit_chain "$arch" "${r[$i]}"
        i=$(( i + 1 ))
      fi
    done
  done
done

echo "===== half-a_t wrapper done ($AT_TAG, DRYRUN=$DRYRUN) [$(date +%F_%H:%M:%S)] ====="
