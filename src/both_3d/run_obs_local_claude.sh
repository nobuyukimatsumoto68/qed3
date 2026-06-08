#!/usr/bin/env bash
# run_obs_local_claude.sh -- LOCAL obs driver, MASSLESS case by default; generalizable to m_F / m_P.
#
# Runs the valence-(0,0) obs pipeline over the massless-SEA ensembles (NO free-field case here):
#   Nf{2,4,6}_gsq8.000000at0.200000nu01.000000nt128L1
#   (NOTE: massless-sea dirs carry NO mRe/mIm suffix; the m_F / m_P sea dirs DO.)
#
# Pipeline per ensemble  [jj_analysis_global_procedure_claude.md stages 1-3]:
#   (1) corr : jj_corr (connected correlators + RAW disc traces)  -> data_<ESNID>/corr_nt0<N>_nhits<H>/
#   (2) disc : jj_disc_postproc (raw J(t) -> disc two-point)       -> .../disc_proced.<k>.h5
#   (3) R    : reweighting_R (Eq 2.5) -- PARITY ONLY (VAL_IM!=0); SKIPPED here (massless => R=1).
#
# INTENT / how to generalize (local sibling of run_obs_fermilab_claude.sh):
#   * Switch mass case by setting VAL_RE/VAL_IM and the ENS_NAMES list:
#       massless: VAL=(0,0),  ENS_NAMES = Nf{2,4,6}_..nt128L1                   (no mRe/mIm suffix)
#       m_F     : VAL=(m,0),  ENS_NAMES = Nf{2,4,6}_..mRe<m>mIm0.000000nt128L1
#       m_P     : VAL=(0,m),  ENS_NAMES = Nf{2,4,6}_..mRe0.000000mIm<m>nt128L1  (R fires automatically)
#   * ENS_ROOT = directory holding the sea-ensemble config dirs (each with ckpoint_lat.*).  REQUIRED.
#   * NHITS: interacting ensembles jackknife over CONFIGS -> 1 hit/config (NHITS=1).  Must match the
#     `nhits` set in the analysis notebook so its loaders find corr_nt0<NT0>_nhits<H>/.
#   * Binaries run from THIS dir (both_3d) so the relative geometry path ../../geometry/data resolves
#     (Fermilab uses the *_fermilab* binaries with an ABSOLUTE geometry path instead).
#
# Usage: bash run_obs_local_claude.sh --ens-root DIR [--gpu N] [--nhits N] [--n-t0 N] [--ninter N]
#                                     [--build] [--dry-run]
set -u

# ---------------- CONFIG ----------------
BINDIR="."
CORR_BIN="${CORR_BIN:-jj_corr_block_t_claude.o}"            # connected+disc; fallback jj_corr_claude.o
POSTPROC_BIN="${POSTPROC_BIN:-jj_disc_postproc_claude.o}"   # disc two-point post-processor (CPU)
R_BIN="${R_BIN:-reweighting_R_claude.o}"                    # parity reweighting (local geometry build)

GSQ=8.0; NU0=1.0; NU1=1.0
NT0="${NT0:-2}"            # source-time origins (--n-t0); names the corr_nt0<N> dir
NINTER="${NINTER:-10}"     # ensemble config stride ckpoint_lat.k for k=0,N,2N,...
NHITS="${NHITS:-1}"        # interacting => jackknife over configs => 1 hit/config (match the notebook)
GPU="${GPU:-0}"

# ---- mass case: MASSLESS (valence 0,0).  See header to switch to m_F / m_P. ----
VAL_RE="${VAL_RE:-0.0}"; VAL_IM="${VAL_IM:-0.0}"

# Where the sea-ensemble config dirs live (each contains ckpoint_lat.*).  REQUIRED (no free case).
ENS_ROOT="${ENS_ROOT:-}"

# Ensembles to run (basenames under ENS_ROOT).  MASSLESS default = the three massless-sea ensembles:
ENS_NAMES=(
  "Nf2_gsq8.000000at0.200000nu01.000000nt128L1"            # massless sea, Nf=2
  "Nf4_gsq8.000000at0.200000nu01.000000nt128L1"            # massless sea, Nf=4
  "Nf6_gsq8.000000at0.200000nu01.000000nt128L1"            # massless sea, Nf=6
)

DO_BUILD="${DO_BUILD:-0}"; DRYRUN=0

# ---------------- CLI ----------------
USAGE="usage: $0 --ens-root DIR [--gpu N] [--nhits N] [--n-t0 N] [--ninter N] [--build] [--dry-run]
  Mass case via env: VAL_RE/VAL_IM + ENS_NAMES (default massless)."
while [[ $# -gt 0 ]]; do
  case "$1" in
    -g|--gpu)     GPU="$2"; shift 2 ;;
    --gpu=*)      GPU="${1#*=}"; shift ;;
    --ens-root)   ENS_ROOT="$2"; shift 2 ;;
    --ens-root=*) ENS_ROOT="${1#*=}"; shift ;;
    -H|--nhits)   NHITS="$2"; shift 2 ;;
    --nhits=*)    NHITS="${1#*=}"; shift ;;
    -T|--n-t0)    NT0="$2"; shift 2 ;;
    --n-t0=*)     NT0="${1#*=}"; shift ;;
    -I|--ninter)  NINTER="$2"; shift 2 ;;
    --ninter=*)   NINTER="${1#*=}"; shift ;;
    --build)      DO_BUILD=1; shift ;;
    --dry-run)    DRYRUN=1; shift ;;
    -h|--help)    echo "${USAGE}"; exit 0 ;;
    *)            echo "unknown arg: $1"; echo "${USAGE}"; exit 1 ;;
  esac
done

: "${ENS_ROOT:?--ens-root <dir> required (directory holding the sea-ensemble config dirs)}"
export CUDA_VISIBLE_DEVICES="${GPU}"
export NVIDIA_VISIBLE_DEVICES="${GPU}"

parity=0; awk "BEGIN{exit !(${VAL_IM}>0)}" && parity=1
echo "### local obs | valence mRe=${VAL_RE} mIm=${VAL_IM} (parity=${parity}) | GPU=${GPU} nhits=${NHITS} n-t0=${NT0} ###"
echo "### ens-root = ${ENS_ROOT} ###"

run() { echo "+ $*"; [[ "${DRYRUN}" == 1 ]] || "$@"; }

if [[ "${DO_BUILD}" == 1 ]]; then
  echo "### building ${CORR_BIN} ${POSTPROC_BIN} ${R_BIN} ###"
  run make "${CORR_BIN}" "${POSTPROC_BIN}" "${R_BIN}"
fi

# ESNID + corr-dir must mirror the binaries' naming (jj_corr: <ens_basename>_vmRe%.6fvmIm%.6f).
esnid()    { printf "%s_vmRe%.6fvmIm%.6f" "$1" "${VAL_RE}" "${VAL_IM}"; }
corr_dir() { printf "data_%s/corr_nt0%d_nhits%d/" "$1" "${NT0}" "${NHITS}"; }

for name in "${ENS_NAMES[@]}"; do
  edir="${ENS_ROOT%/}/${name}/"
  if [[ ! -d "${edir}" && "${DRYRUN}" == 0 ]]; then echo "## skip (no dir): ${edir}"; continue; fi
  [[ "${name}" =~ ^Nf([0-9]+)_ ]] && Nf="${BASH_REMATCH[1]}" || Nf=2
  es="$(esnid "${name}")"; cdir="$(corr_dir "${es}")"
  echo "## ${name}  Nf=${Nf}  ->  ${cdir}"

  # (1) connected + raw disc traces
  if [[ -x "${BINDIR}/${CORR_BIN}" || "${DRYRUN}" == 1 ]]; then
    run ./"${CORR_BIN}" --gsq "${GSQ}" --Nf "${Nf}" --nu0 "${NU0}" --nu1 "${NU1}" \
        --mass-re "${VAL_RE}" --mass-im "${VAL_IM}" --nhits "${NHITS}" --n-t0 "${NT0}" \
        --ninter "${NINTER}" --ens-dir "${edir}"
  else
    echo "   (CORR_BIN missing: ${BINDIR}/${CORR_BIN} -- run with --build or make it)"
  fi

  # (2) disc two-point from the raw J(t)
  if [[ -x "${BINDIR}/${POSTPROC_BIN}" || "${DRYRUN}" == 1 ]]; then
    run ./"${POSTPROC_BIN}" "${cdir}"
  else
    echo "   (POSTPROC_BIN missing: ${BINDIR}/${POSTPROC_BIN})"
  fi

  # (3) reweighting R -- parity valence only (massless/flavor: skipped, R trivial)
  if [[ "${parity}" == 1 ]]; then
    if [[ -x "${BINDIR}/${R_BIN}" || "${DRYRUN}" == 1 ]]; then
      run ./"${R_BIN}" --nu0 "${NU0}" --nu1 "${NU1}" \
          --mass-re "${VAL_RE}" --mass-im "${VAL_IM}" --ninter "${NINTER}" --ens-dir "${edir}"
    else
      echo "   (R_BIN missing: ${BINDIR}/${R_BIN})"
    fi
  fi
done
echo "### done ###"
