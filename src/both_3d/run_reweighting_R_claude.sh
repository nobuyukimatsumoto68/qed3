#!/usr/bin/env bash
# run_reweighting_R_claude.sh -- PRODUCTION driver for the parity reweighting factor R (PDF Eq. 2.5).
#
# Standalone: run SEPARATELY from run_jj_analysis_claude.sh.  Builds reweighting_R_claude.o (if needed),
# then computes R per gauge config into  data_<ESNID>/R/R.<k>.h5  (the SAME obs dir as the correlators,
# same ESNID = <ens|free>_vmRe<re>vmIm<im>).  Stage 3 of jj_analysis_global_procedure_claude.md.
#
# R = conj( prod_i ((1-mP)lam_i + mP)/(lam_i + mP) ),  {lam_i} = eigenvalues of the MASSLESS D_ov.
# R is ONLY needed for PARITY (m_P, imaginary) masses: massless gives R=1 (mP=0) and flavor (real mass)
# is trivial, so the free/sweep modes below run the m_P cases only.  The spectrum is MASS-INDEPENDENT and
# is stored in each R.<k>.h5 (`lam`), so you can also reweight any other m_P in post from `lam`.
#
# Usage: bash run_reweighting_R_claude.sh [--gpu N] [--ninter N] [--build|--no-build] [free|sweep|ens]
#   ens mode (single explicit ensemble): also --ens-dir <path/> [--mass-re X] [--mass-im Y] [--Nf N]
# Heavy: each config = N=3072 overlap applies (dense build) + one 3072^2 geev.

set -u

# ---------------- configuration ----------------
BINDIR="."
SRC="${BINDIR}/reweighting_R_claude.cu"
BIN="${BINDIR}/reweighting_R_claude.o"
GPU="${GPU:-0}"
NINTER="${NINTER:-}"      # ensemble config stride k=0,N,2N,...; empty => program default (10)
DO_BUILD="auto"           # auto (build if missing/stale) | force | skip

# ensemble labels (match run_jj_analysis_claude.sh / the hmc sea-dir naming)
GSQ=8.0; NU0=1.0; NU1="${NU0}"; AT=0.2; NT=128; L=1
NFS=(2 4 6)
MASSES=(0.01 0.05 0.1 0.2)

# ens-mode knobs
ENS_DIR="${ENS_DIR:-}"; ENS_NF="${ENS_NF:-2}"; VAL_RE="${VAL_RE:-0.0}"; VAL_IM="${VAL_IM:-0.0}"

# build flags (mirror the Makefile %.o rule)
NVCC="nvcc"
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

# ---------------- helpers ----------------
build() {
  case "${DO_BUILD}" in
    skip) echo "### build: skipped (--no-build) ###"; return 0 ;;
    auto) if [[ -x "${BIN}" && "${BIN}" -nt "${SRC}" ]]; then
            echo "### build: up to date (${BIN}) ###"; return 0
          fi ;;
  esac
  echo "### building ${BIN} ###"
  # shellcheck disable=SC2086
  ${NVCC} ${NVCCFLAGS} ${INCLUDES} ${LDFLAGS} "${SRC}" -o "${BIN}" \
    || { echo "BUILD FAILED"; exit 1; }
  echo "### build OK ###"
}

# sea config directory (matches hmc_w_mass dir naming, std::to_string => 6 decimals)
ens_dir() {  # args: Nf seaRe seaIm
  printf "Nf%d_gsq%.6fat%.6fnu0%.6fmRe%.6fmIm%.6fnt%dL%d/" \
         "$1" "${GSQ}" "${AT}" "${NU0}" "$2" "$3" "${NT}" "${L}"
}

run_R() {  # args: massRe massIm [ens_dir or "" for free]
  local vr="$1" vi="$2" ed="$3"
  local edflag=(); [[ -n "${ed}" ]]     && edflag=(--ens-dir "${ed}")
  local niflag=(); [[ -n "${NINTER}" ]] && niflag=(--ninter "${NINTER}")
  "${BIN}" --nu0 "${NU0}" --nu1 "${NU1}" --mass-re "${vr}" --mass-im "${vi}" \
           --gpu "${GPU}" "${niflag[@]}" "${edflag[@]}"
}

# ---------------- CLI ----------------
USAGE="usage: $0 [--gpu N] [--ninter N] [--build|--no-build] [free|sweep|ens]   (default: free)
  ens mode: also --ens-dir <path/> [--mass-re X] [--mass-im Y] [--Nf N]"
MODE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -g|--gpu)    GPU="$2"; shift 2 ;;
    --gpu=*)     GPU="${1#*=}"; shift ;;
    -I|--ninter) NINTER="$2"; shift 2 ;;
    --ninter=*)  NINTER="${1#*=}"; shift ;;
    --build)     DO_BUILD="force"; shift ;;
    --no-build)  DO_BUILD="skip";  shift ;;
    --ens-dir)   ENS_DIR="$2"; shift 2 ;;
    --ens-dir=*) ENS_DIR="${1#*=}"; shift ;;
    --Nf)        ENS_NF="$2"; shift 2 ;;
    --Nf=*)      ENS_NF="${1#*=}"; shift ;;
    --mass-re)   VAL_RE="$2"; shift 2 ;;
    --mass-re=*) VAL_RE="${1#*=}"; shift ;;
    --mass-im)   VAL_IM="$2"; shift 2 ;;
    --mass-im=*) VAL_IM="${1#*=}"; shift ;;
    free|sweep|ens) MODE="$1"; shift ;;
    -h|--help)   echo "${USAGE}"; exit 0 ;;
    *)           echo "unknown arg: $1"; echo "${USAGE}"; exit 1 ;;
  esac
done
MODE="${MODE:-free}"

export CUDA_VISIBLE_DEVICES="${GPU}"
export NVIDIA_VISIBLE_DEVICES="${GPU}"
echo "### R driver | GPU=${GPU} | ninter=${NINTER:-<prog default 10>} | mode=${MODE} ###"

build

# ---------------- modes ----------------
if [[ "${MODE}" == "free" ]]; then
  # R is only needed for PARITY (m_P) masses: massless -> R=1 (mP=0), flavor (real mass) -> trivial.
  echo "### free-field R (single U=1 config): parity m_P only ###"
  run_R 0.0 0.1 ""     # parity (imaginary mass) -- the only case where R is needed

elif [[ "${MODE}" == "sweep" ]]; then
  echo "### production sweep: R per config, PARITY (m_P) ensembles only ###"
  for Nf in "${NFS[@]}"; do
    for m in "${MASSES[@]}"; do
      run_R 0.0 "${m}" "$(ens_dir "${Nf}" 0.0 "${m}")"   # parity sea (R needed only here)
    done
  done

elif [[ "${MODE}" == "ens" ]]; then
  : "${ENS_DIR:?ens mode needs --ens-dir <path/>}"
  [[ "${ENS_DIR}" != */ ]] && ENS_DIR="${ENS_DIR}/"
  echo "### single ensemble R: ${ENS_DIR} (mRe=${VAL_RE} mIm=${VAL_IM}) ###"
  run_R "${VAL_RE}" "${VAL_IM}" "${ENS_DIR}"

else
  echo "${USAGE}"; exit 1
fi
