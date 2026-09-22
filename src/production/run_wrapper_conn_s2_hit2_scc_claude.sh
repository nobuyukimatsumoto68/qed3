#!/bin/bash -l
# run_wrapper_conn_s2_hit2_scc_claude.sh  (_scc, 2026-08-29, NM)  -- SECOND-HIT login wrapper
# =============================================================================================
# Adds the SECOND stochastic hit (h1) to SCC's L4 stride-2 connected-Ylm measurement, across ALL 3 SCC
# ensembles (Nf2 g2.0, Nf6 g4.0, Nf6 g6.0) x 4 offsets = 12 units. Sibling of run_wrapper_conn_s2_scc_claude.sh;
# reuses the SAME binaries (jj_conn_s2_L4_sm_{70,80}.out) and the SAME assignment table.
#
# WHAT IS DIFFERENT vs the nhits1 wrapper:
#   (1) submits run_conn_s2_hit2_scc_claude.sh, which runs --nhits 2 --outdir-nhits 1 -> writes h1 as
#       corr.<k>.h1.h5 INTO the EXISTING corr_ylm_conn_t00_nhits1_s1 dir (SAME dir as h0). The driver's per-hit
#       "complete"-gate finds h0 already there and SKIPS it, so only h1 is solved. NO separate nhits2 dir and NO
#       hardlinks -> the LOCAL (barracuda22) pull-back is UNCHANGED (its glob nhits1_s1/corr.*.h5 now grabs h1 too).
#   (2) job namespace c2h<arch> (distinct from cs2<arch>) so its chains anchor independently.
# h1 RNG = deterministic per-hit seed esnid_k<k>_h1 (the EXISTING convention, same scheme as h0; INDEPENDENT,
# reproducible, poolable with FNAL/LOCAL). Output: data_<ESNID>/corr_ylm_conn_t00_nhits1_s1/corr.<k>.h1.h5.
# NB requires the binary rebuilt with the --outdir-nhits flag (jj_local_ylm_scalar_conn_stoch_claude.cu, 2026-08-29).
#
# Usage:  bash run_wrapper_conn_s2_hit2_scc_claude.sh
#         DRYRUN=1 bash ...                                    # print qsub lines, do nothing
#         GPUT_SM80=L40S ... (default)                         # conn is bandwidth-bound -> L40S is fast + off the HMC V100
# =============================================================================================
set -u

SRCDIR=/projectnb/qfe/nmatsum/qed3/src/production
cd "$SRCDIR" || { echo "ERROR: cannot cd $SRCDIR"; exit 1; }
source /projectnb/qfe/nmatsum/qed3/env.sh
module load hdf5/1.10.10 2>/dev/null
module load gsl 2>/dev/null

ASSIGN=${ASSIGN:-conn_stride2_assign_claude.txt}
LREF=4                               # SCC does L4 only
BATCH=run_conn_s2_hit2_scc_claude.sh

PE_OMP=${PE_OMP:-4}
H_RT_FIRST=${H_RT_FIRST:-12:00:00}
H_RT=${H_RT:-12:00:00}
N_CHAIN=${N_CHAIN:-4}                # links/unit (complete-gating: each resumes; over-provision is free)
ARCH_LIST=${ARCH_LIST:-"sm_80"}      # 2nd-hit default = single arch (sm_80 runs native on L40S/A100)
SUBMIT_ARCHS=${SUBMIT_ARCHS:-"$ARCH_LIST"}
GPUT_SM70=${GPUT_SM70:-V100}
GPUT_SM80=${GPUT_SM80:-L40S}         # conn is bandwidth-bound -> L40S fast + off the live HMC V100 (default L40S)
DRYRUN=${DRYRUN:-0}

binname () { echo "jj_conn_s2_L${LREF}_${1}.out"; }   # $1 = arch  (reuse the nhits1 binaries; no rebuild)

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

# ---- SCC L4 units from the assignment table (same as the nhits1 wrapper) ------------------------
# cols: 1=ensemble 2=L 3=Nf 4=gsq 5=at 6=offset 7=kmin 8=stride 9=kmax 10=cfg 11=site 12=worker_h
[ -f "$ASSIGN" ] || { echo "ERROR: assignment table $ASSIGN not found"; exit 1; }
# 2026-09-02 (NM): SCC owns ALL 9 L4 conn -> h1 production over all 36 units (matches the nhits1 wrapper's
# $2==4 widening). Old SCC-only-12 filter kept below for rollback.
# mapfile -t UNITS < <(awk '$0 !~ /^#/ && toupper($11)=="SCC" && $2==4 {print $1" "$3" "$4" "$7}' "$ASSIGN")
mapfile -t UNITS < <(awk '$0 !~ /^#/ && $2==4 {print $1" "$3" "$4" "$7}' "$ASSIGN")
echo "===== L4 conn h1 units (all 9 ensembles): ${#UNITS[@]} (want 36) ====="
[ "${#UNITS[@]}" -eq 0 ] && { echo "ERROR: no SCC L4 units parsed from $ASSIGN"; exit 1; }

# binaries must already exist AND be rebuilt with the --outdir-nhits flag; do NOT rebuild here.
# GUARD: an OLD binary (no --outdir-nhits) would hit getopt '?' -> PrintHelp -> exit(0), i.e. SILENTLY compute
# nothing while logging "done (status 0)". grep the compiled help string to catch that before submitting.
for a in $SUBMIT_ARCHS
do
  app=$(binname "$a")
  test -f "$app" || { echo "ERROR: $app missing -- build it with run_wrapper_conn_s2_scc_claude.sh first"; exit 1; }
  grep -q 'outdir-nhits' "$app" || { echo "ERROR: $app lacks --outdir-nhits (would silently no-op). REBUILD the conn binary with the updated jj_local_ylm_scalar_conn_stoch_claude.cu (FORCE_BUILD=1) first."; exit 1; }
done

# NOTE: no h0 pre-linking needed -- with --outdir-nhits 1 the driver writes h1 INTO the existing nhits1_s1 dir
# where h0 already lives, so the per-hit "complete"-gate skips h0 automatically. (h1-only solve, same dir.)

read -r -a ARCHES <<< "$SUBMIT_ARCHS"
n_arch=${#ARCHES[@]}

# live kmax = (max ckpoint_lat)+1 for an ensemble dir (configs still growing)
live_kmax () {   # ensdir
  local k
  k=$(ls "$1"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n | tail -1)
  echo $(( ${k:-0} + 1 ))
}

# cross-run anchoring within the c2h namespace (^c2h) -> re-run EXTENDS each unit's chain (new kmax)
EXIST_JOBS=$(qstat -u "$USER" -r 2>/dev/null | awk '
  /^[0-9]+ /{ jid=$1 }
  /Full jobname:/{ if($3 ~ /^c2h/) print jid"\t"$3 }
')
existing_tail () {   # token
  local a=${1//./\\.}
  printf '%s\n' "$EXIST_JOBS" | awk -F'\t' -v pat="(${a})(__|_c)" '$2 ~ pat {print $1}' | sort -n | tail -1
}

submit_unit () {   # arch  "ensdir nf gsq kmin"
  local arch=$1 rec=$2
  local ensdir nf gsq kmin
  read -r ensdir nf gsq kmin <<< "$rec"
  local app gpuc gput kmax name tok prev c hrt jid
  app=$(binname "$arch")
  gpuc=$(gpuc_of "$arch")
  gput=$(gput_of "$arch")
  kmax=$(live_kmax "$ensdir")
  # ARCH-INDEPENDENT anchoring token (same as the cs2 wrapper's 2026-08-10 fix)
  tok="Nf${nf}g${gsq}k${kmin}"
  name="c2h${arch#sm_}_Nf${nf}g${gsq}k${kmin}"
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
    qsub_cmd+=( -v "${vars}" "$BATCH" )
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

echo "===== submit ${#UNITS[@]} units x $N_CHAIN links (round-robin ${SUBMIT_ARCHS}, --nhits 2, h1 only) ====="
ri=0
for rec in "${UNITS[@]}"
do
  a=${ARCHES[$(( ri % n_arch ))]}
  submit_unit "$a" "$rec"
  ri=$(( ri + 1 ))
done

echo "===== conn-s2 HIT2 wrapper done (DRYRUN=$DRYRUN) [$(date +%F_%H:%M:%S)] ====="
