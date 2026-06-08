#!/usr/bin/env bash
# run_obs_fermilab_claude.sh -- FERMILAB wrapper: loop over ensembles and compute corr (+ R for parity).
#
# Place a COPY in /lustre2/qed3 and in /lustre2/affine (it auto-detects its own directory as the
# ensemble root), OR run it from anywhere with --ens-root <dir>.  For each ensemble subdir it parses
# Nf / mRe / mIm from the name and runs:
#   - the connected+disc correlators (CORR_BIN, valence = sea mass), and
#   - the reweighting factor R (R_BIN) -- ONLY for PARITY ensembles (mIm != 0).
#
# Ensemble dirs:  Nf<N>_gsq8.000000at0.200000nu01.000000mRe<re>mIm<im>nt128L1/  (each holds ckpoint_lat.*)
#   /lustre2/qed3   = masses {0.1, 0.2};   /lustre2/affine = masses {0.01, 0.05};  Nf in {2,4,6}.
#
# The binaries are run with CWD = BOTH3D so their RELATIVE geometry path ("../../geometry/data")
# resolves (= /project/qed3/qed3/geometry/data).  They therefore WRITE their output under
#   BOTH3D/data_<ESNID>/{corr_nt0<N>_nhits<H>/, R/}     (ESNID = <ensname>_vmRe<re>vmIm<im>)
# i.e. in the repo, not on /lustre2.  (--ens-dir is passed as the absolute /lustre2 path.)
#
# Usage:  bash run_obs_fermilab_claude.sh [options]
#
# ============================================================================================
# CLUSTER AGENT -- this script is the OUTLINE (what to compute per ensemble); refine the DETAILS:
#   [ ] BOTH3D path (repo on the cluster) -- CONFIG below or --both3d.
#   [ ] geometry path: the *_fermilab*.cu binaries hardcode an ABSOLUTE geometry root
#       (reweighting_R_fermilab_claude.cu: /project/qed3/qed3/geometry/...). Fix that line if different.
#   [ ] CORR_BIN: if the corr binary also needs absolute geometry, build a jj_corr_block_t_fermilab
#       variant and set CORR_BIN to it (R_BIN already defaults to the _fermilab build).
#   [ ] GPU selection: under SLURM --gres=gpu:N the allocated GPU is visible as device 0; GPU=0 is fine.
#   [ ] build: pass --build (runs `make` in BOTH3D) or prebuild the binaries.
#   [ ] batch: submit via submit_obs_fermilab_claude.sh (SLURM template).  For a SLURM job ARRAY,
#       use --list to enumerate ensembles and --index $SLURM_ARRAY_TASK_ID to run one per task.
# The corr/R logic, parity-only-R rule, name parsing, and flags below are the stable outline.
# ============================================================================================
set -u

# ---------------- CONFIG (edit if the cluster paths differ) ----------------
BOTH3D="${BOTH3D:-/project/qed3/qed3/src/both_3d}"   # repo to run/update (project_qed3/qed3/src/both_3d)
CORR_BIN="${CORR_BIN:-jj_corr_block_t_claude.o}"     # connected+disc; fallback jj_corr_claude.o.
                                                     #   If corr also hits the geometry-path issue,
                                                     #   make a jj_corr_block_t_fermilab variant (absolute
                                                     #   geometry, like R_BIN) and set CORR_BIN to it.
R_BIN="${R_BIN:-reweighting_R_fermilab_claude.o}"    # reweighting factor R (parity only); FERMILAB build
                                                     #   (absolute geometry paths). Built via pattern rule:
                                                     #   make reweighting_R_fermilab_claude.o
ENS_GLOB="${ENS_GLOB:-Nf*nt128L1}"                   # ensemble subdir pattern under the ensemble root

GSQ=8.0; NU0=1.0; NU1=1.0                            # ensemble labels (match the dir naming)
NHITS="${NHITS:-1}"                                  # stochastic hits per config (production: 1)
NT0="${NT0:-2}"                                      # source-time origins (--n-t0)
NINTER="${NINTER:-10}"                               # config stride: ckpoint_lat.k for k=0,N,2N,...
GPU="${GPU:-0}"                                      # CUDA device (via CUDA_VISIBLE_DEVICES)
DO_BUILD="${DO_BUILD:-0}"                            # 1 => `make` the two binaries in BOTH3D first
DO_CORR=1; DO_R=1                                    # toggles (see --corr-only/--r-only)
DRYRUN=0

# ensemble root = this script's own directory unless overridden by --ens-root
ENS_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------- CLI ----------------
USAGE="usage: $0 [--ens-root DIR] [--both3d DIR] [--gpu N] [--nhits N] [--n-t0 N] [--ninter N]
                 [--build] [--corr-only|--r-only] [--dry-run]
  Defaults: ens-root = this script's dir; both3d = ${BOTH3D}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ens-root)  ENS_ROOT="$2"; shift 2 ;;
    --ens-root=*) ENS_ROOT="${1#*=}"; shift ;;
    --both3d)    BOTH3D="$2"; shift 2 ;;
    --both3d=*)  BOTH3D="${1#*=}"; shift ;;
    -g|--gpu)    GPU="$2"; shift 2 ;;
    --gpu=*)     GPU="${1#*=}"; shift ;;
    -H|--nhits)  NHITS="$2"; shift 2 ;;
    --nhits=*)   NHITS="${1#*=}"; shift ;;
    -T|--n-t0)   NT0="$2"; shift 2 ;;
    --n-t0=*)    NT0="${1#*=}"; shift ;;
    -I|--ninter) NINTER="$2"; shift 2 ;;
    --ninter=*)  NINTER="${1#*=}"; shift ;;
    --build)     DO_BUILD=1; shift ;;
    --corr-only) DO_R=0; shift ;;
    --r-only)    DO_CORR=0; shift ;;
    --dry-run)   DRYRUN=1; shift ;;
    -h|--help)   echo "${USAGE}"; exit 0 ;;
    *)           echo "unknown arg: $1"; echo "${USAGE}"; exit 1 ;;
  esac
done

export CUDA_VISIBLE_DEVICES="${GPU}"
export NVIDIA_VISIBLE_DEVICES="${GPU}"

[[ -d "${BOTH3D}" ]]   || { echo "ERROR: BOTH3D not found: ${BOTH3D}"; exit 1; }
[[ -d "${ENS_ROOT}" ]] || { echo "ERROR: ens-root not found: ${ENS_ROOT}"; exit 1; }
echo "### Fermilab obs wrapper ###"
echo "  ens-root = ${ENS_ROOT}"
echo "  both3d   = ${BOTH3D}   (binaries run here; geometry = ../../geometry/data)"
echo "  GPU=${GPU}  nhits=${NHITS}  n-t0=${NT0}  ninter=${NINTER}  corr=${DO_CORR} R=${DO_R}"

if [[ "${DO_BUILD}" == 1 ]]; then
  echo "### building ${CORR_BIN} ${R_BIN} in ${BOTH3D} ###"
  make -C "${BOTH3D}" "${CORR_BIN}" "${R_BIN}" || { echo "BUILD FAILED"; exit 1; }
fi

run() { echo "+ $*"; [[ "${DRYRUN}" == 1 ]] || "$@"; }

shopt -s nullglob
ens_list=( "${ENS_ROOT}"/${ENS_GLOB}/ )
[[ ${#ens_list[@]} -gt 0 ]] || { echo "no ensembles matching ${ENS_GLOB}/ under ${ENS_ROOT}"; exit 1; }
echo "### ${#ens_list[@]} ensemble(s) under ${ENS_ROOT} ###"

for ens in "${ens_list[@]}"; do
  name="$(basename "${ens}")"
  if [[ ! "${name}" =~ Nf([0-9]+)_.*mRe([0-9.]+)mIm([0-9.]+)nt ]]; then
    echo "## skip (name not parseable): ${name}"; continue
  fi
  Nf="${BASH_REMATCH[1]}"; mRe="${BASH_REMATCH[2]}"; mIm="${BASH_REMATCH[3]}"
  parity=0; awk "BEGIN{exit !(${mIm}>0)}" && parity=1   # parity ensemble iff mIm != 0
  echo "## ${name}  ->  Nf=${Nf} mRe=${mRe} mIm=${mIm}  parity=${parity}"

  # connected + disc correlators (valence mass = the sea mass of this ensemble)
  if [[ "${DO_CORR}" == 1 ]]; then
    if [[ -x "${BOTH3D}/${CORR_BIN}" ]]; then
      run bash -c "cd '${BOTH3D}' && ./'${CORR_BIN}' --gsq ${GSQ} --Nf ${Nf} --nu0 ${NU0} --nu1 ${NU1} \
        --mass-re ${mRe} --mass-im ${mIm} --nhits ${NHITS} --n-t0 ${NT0} --ninter ${NINTER} --ens-dir '${ens}'"
    else
      echo "   (CORR_BIN missing: ${BOTH3D}/${CORR_BIN} -- run with --build or make it)"
    fi
  fi

  # reweighting factor R -- PARITY ensembles only (massless/flavor: R trivial)
  if [[ "${DO_R}" == 1 && "${parity}" == 1 ]]; then
    if [[ -x "${BOTH3D}/${R_BIN}" ]]; then
      run bash -c "cd '${BOTH3D}' && ./'${R_BIN}' --nu0 ${NU0} --nu1 ${NU1} \
        --mass-re ${mRe} --mass-im ${mIm} --ninter ${NINTER} --ens-dir '${ens}'"
    else
      echo "   (R_BIN missing: ${BOTH3D}/${R_BIN} -- run with --build or make it)"
    fi
  fi
done
echo "### done ###"
