#!/usr/bin/env bash
# run_postproc_local_claude.sh -- standalone DISC POST-PROCESSING driver (CPU only).
#
# Runs jj_disc_postproc over corr_nt0<N>_nhits<H>/ dirs: reads the RAW disc traces J(t) in
# corr.<k>.h5 and writes disc_proced.<k>.h5 (the disc TWO-POINT, option A) into the SAME dir.
# [jj_analysis_global_procedure_claude.md stage 2; jj_disc_postproc_impl_plan_claude.md]
#
# This is stage (2) of run_obs_local_claude.sh, PULLED OUT so it can run independently of the
# (slow, GPU) corr inversions: launch it any time over whatever corr.*.h5 already exist -- it is
# config-granular and idempotent (own "complete" sentinel per disc_proced.<k>.h5).  run_obs_local
# STILL runs this stage inline; this standalone is for re-runs / catching up after more configs land.
#
# Two ways to point it at data:
#   (A) ensemble interface (mirrors run_obs_local naming):
#         bash run_postproc_local_claude.sh --ens-root DIR [--nf N] [--n-t0 N] [--nhits N]
#                                           [VAL_RE=.. VAL_IM=.. as env]
#   (B) explicit corr dirs (ad-hoc; any naming):
#         bash run_postproc_local_claude.sh DIR1/corr_nt02_nhits1/ DIR2/corr_nt02_nhits1/ ...
#
# Mass case via env (same as run_obs_local): VAL_RE/VAL_IM + ENS_NAMES select the esnid suffix.
set -u

# ---------------- CONFIG ----------------
BINDIR="."
POSTPROC_BIN="${POSTPROC_BIN:-jj_disc_postproc_claude.o}"   # disc two-point post-processor (CPU)

NT0="${NT0:-2}"            # source-time origins -> corr_nt0<N> dir
NHITS="${NHITS:-1}"        # -> corr_nt0<N>_nhits<H> dir
VAL_RE="${VAL_RE:-0.0}"; VAL_IM="${VAL_IM:-0.0}"            # valence mass -> esnid suffix
ENS_ROOT="${ENS_ROOT:-}"
NF_FILTER="${NF_FILTER:-}"

# Ensembles (basenames; massless-sea default -- same list as run_obs_local).
ENS_NAMES=(
  "Nf2_gsq8.000000at0.200000nu01.000000nt128L1"
  "Nf4_gsq8.000000at0.200000nu01.000000nt128L1"
  "Nf6_gsq8.000000at0.200000nu01.000000nt128L1"
)

DO_BUILD="${DO_BUILD:-0}"; DRYRUN=0
DIRS=()   # explicit corr dirs (mode B), collected from positional args

# ---------------- CLI ----------------
USAGE="usage: $0 [--ens-root DIR [--nf N]] [--n-t0 N] [--nhits N] [--build] [--dry-run] [CORRDIR ...]
  ensemble mode: --ens-root DIR (+ optional --nf N).  Mass case via env VAL_RE/VAL_IM + ENS_NAMES.
  ad-hoc mode  : pass one or more corr_nt0<N>_nhits<H>/ dirs as positional args (no --ens-root)."
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ens-root)   ENS_ROOT="$2"; shift 2 ;;
    --ens-root=*) ENS_ROOT="${1#*=}"; shift ;;
    --nf)         NF_FILTER="$2"; shift 2 ;;
    --nf=*)       NF_FILTER="${1#*=}"; shift ;;
    -T|--n-t0)    NT0="$2"; shift 2 ;;
    --n-t0=*)     NT0="${1#*=}"; shift ;;
    -H|--nhits)   NHITS="$2"; shift 2 ;;
    --nhits=*)    NHITS="${1#*=}"; shift ;;
    --build)      DO_BUILD=1; shift ;;
    --dry-run)    DRYRUN=1; shift ;;
    -h|--help)    echo "${USAGE}"; exit 0 ;;
    --*)          echo "unknown arg: $1"; echo "${USAGE}"; exit 1 ;;
    *)            DIRS+=("$1"); shift ;;   # positional => explicit corr dir
  esac
done

run() { echo "+ $*"; [[ "${DRYRUN}" == 1 ]] || "$@"; }

if [[ "${DO_BUILD}" == 1 ]]; then
  echo "### building ${POSTPROC_BIN} ###"
  run make "${POSTPROC_BIN}"
fi

if [[ ! -x "${BINDIR}/${POSTPROC_BIN}" && "${DRYRUN}" == 0 ]]; then
  echo "POSTPROC_BIN missing: ${BINDIR}/${POSTPROC_BIN} -- run with --build or make it"; exit 1
fi

# ESNID + corr-dir must mirror the corr binary's naming (<ens_basename>_vmRe%.6fvmIm%.6f).
esnid()    { printf "%s_vmRe%.6fvmIm%.6f" "$1" "${VAL_RE}" "${VAL_IM}"; }
corr_dir() { printf "data_%s/corr_nt0%d_nhits%d/" "$1" "${NT0}" "${NHITS}"; }

# Build the dir list: explicit (mode B) wins; else derive from the ensemble list (mode A).
if [[ ${#DIRS[@]} -eq 0 ]]; then
  : "${ENS_ROOT:?--ens-root <dir> required (or pass explicit corr dirs as positional args)}"
  if [[ -n "${NF_FILTER}" ]]; then
    filtered=(); for nm in "${ENS_NAMES[@]}"; do [[ "${nm}" =~ ^Nf${NF_FILTER}_ ]] && filtered+=("${nm}"); done
    ENS_NAMES=("${filtered[@]}")
    [[ ${#ENS_NAMES[@]} -gt 0 ]] || { echo "no ensemble matches Nf${NF_FILTER}"; exit 1; }
  fi
  for name in "${ENS_NAMES[@]}"; do DIRS+=("$(corr_dir "$(esnid "${name}")")"); done
fi

echo "### disc post-proc | valence mRe=${VAL_RE} mIm=${VAL_IM} | n-t0=${NT0} nhits=${NHITS} | ${#DIRS[@]} dir(s) ###"

for cdir in "${DIRS[@]}"; do
  if [[ ! -d "${cdir}" && "${DRYRUN}" == 0 ]]; then echo "## skip (no dir): ${cdir}"; continue; fi
  echo "## ${cdir}"
  run ./"${POSTPROC_BIN}" "${cdir}"
done
echo "### done ###"
