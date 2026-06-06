#!/usr/bin/env bash
# run_jj_analysis_claude.sh -- driver for the conserved-current correlator measurements
# (Chunk 10 of conserved_current_correlators_impl_plan_v3_claude.md).
#
# Holds the per-ensemble case logic (v3 Sec. 1.1): it points each executable at the SEA config
# directory and passes the appropriate VALENCE mass (--mass-re/--mass-im) and --current.
#   massless ensemble : vector(conn+disc) valence 0      ; axial(conn) valence 0
#   m_F ensemble       : vector valence m_F (exact)       ; axial valence 0 (massless formulas)
#   m_P ensemble       : vector valence m_P (exact)       ; axial valence m_P (exact)
# Disconnected is vector-only; axial computes only C_{A+-}.
#
# Usage: bash run_jj_analysis_claude.sh [--gpu N[,N...]] [--nhits N] [--n-t0 N] [free|sweep]
#   --gpu / -g  N[,N...]  -- GPU device(s) to expose (sets CUDA_VISIBLE_DEVICES and
#                            NVIDIA_VISIBLE_DEVICES; e.g. --gpu 0,1).  Also via GPU= env.
#   --nhits / -H  N        -- stochastic hits per program call.  Also via NHITS= env.
#   --n-t0 / -T   N        -- number of source-time origins forwarded to the connected programs;
#                            empty => each program's default (2).  Also via NT0= env.
# Two modes (default "free"):
#   free   -- free-field check: run WITHOUT --ens-dir (U=1, single deterministic config).
#   sweep  -- production sweep over the 27 ensembles.
#
# NOTE: of the connected programs, only jj_conn_tpproj_claude.o is implemented so far; the
# sp/ylm and disc entries below are wired for when those land.  Sized for the cluster (edit
# GPUS / run in parallel as needed) -- not artificially core-capped.

set -u

# ---------------- configuration ----------------
BINDIR="."
GPU="${GPU:-0}"                 # default device(s); override with --gpu/-g or GPU= env (e.g. 0,1)

GSQ=8.0
NU0=1.0
NU1="${NU0}"   # valence asymmetry defaults to sea nu0 (knob retained)
AT=0.2
NT=128
L=1
NHITS="${NHITS:-}"   # hits per program call; per-mode default applied after CLI (free 16, sweep 1)
NT0="${NT0:-}"       # # source-time origins; empty => program default (2); override --n-t0/-T or NT0= env

NFS=(2 4 6)
MASSES=(0.01 0.05 0.1 0.2)

# program binaries (CLI: --gsq --Nf --nu0 --nu1 --mass-re --mass-im --current --ens-dir --nhits)
DISC="${BINDIR}/jj_disc_claude.o"                       # (not yet implemented)
CONN_SP="${BINDIR}/jj_conn_spproj_claude.o"             # vector sp  (implemented)
CONN_SP_AXIAL="${BINDIR}/jj_conn_spproj_axial_claude.o" # axial  sp  (implemented)
CONN_TP="${BINDIR}/jj_conn_tpproj_claude.o"             # vector tp  (implemented)
CONN_TP_AXIAL="${BINDIR}/jj_conn_tpproj_axial_claude.o" # axial  tp  (implemented)
CONN_YLM="${BINDIR}/jj_conn_ylmproj_claude.o"           # vector ylm (implemented; vector only)

CONN_PROGS=("${CONN_TP}" "${CONN_SP}" "${CONN_YLM}")     # vector projections
CONN_AXIAL_PROGS=("${CONN_TP_AXIAL}" "${CONN_SP_AXIAL}")  # axial projections; add ylm axial as implemented

# ---------------- helpers ----------------
# sea config directory matching hmc_w_mass_claude.cu dir3 (std::to_string => 6 decimals)
ens_dir() {  # args: Nf seaRe seaIm
  printf "Nf%d_gsq%.6fat%.6fnu0%.6fmRe%.6fmIm%.6fnt%dL%d/" \
         "$1" "${GSQ}" "${AT}" "${NU0}" "$2" "$3" "${NT}" "${L}"
}

# run all connected vector programs (+ disc) for one (ensemble, valence)
run_vector() {  # args: Nf valRe valIm [ens_dir or "" for free]
  local Nf="$1" vr="$2" vi="$3" ed="$4"
  local edflag=(); [[ -n "${ed}" ]]   && edflag=(--ens-dir "${ed}")
  local ntflag=(); [[ -n "${NT0}" ]] && ntflag=(--n-t0 "${NT0}")   # empty => program default (2)
  # disconnected dump (vector only; single-time traces, origin-agnostic, no --n-t0)
  [[ -x "${DISC}" ]] && "${DISC}" --gsq "${GSQ}" --Nf "${Nf}" --nu0 "${NU0}" --nu1 "${NU1}" \
        --mass-re "${vr}" --mass-im "${vi}" --nhits "${NHITS}" "${edflag[@]}"
  # connected, all projections
  for p in "${CONN_PROGS[@]}"; do
    [[ -x "${p}" ]] && "${p}" --gsq "${GSQ}" --Nf "${Nf}" --nu0 "${NU0}" --nu1 "${NU1}" \
          --mass-re "${vr}" --mass-im "${vi}" --current vector --nhits "${NHITS}" "${ntflag[@]}" "${edflag[@]}"
  done
}

# run all connected axial programs for one (ensemble, valence)   [no disc for axial]
run_axial() {  # args: Nf valRe valIm [ens_dir or "" for free]
  local Nf="$1" vr="$2" vi="$3" ed="$4"
  local edflag=(); [[ -n "${ed}" ]]   && edflag=(--ens-dir "${ed}")
  local ntflag=(); [[ -n "${NT0}" ]] && ntflag=(--n-t0 "${NT0}")   # empty => program default (2)
  for p in "${CONN_AXIAL_PROGS[@]}"; do
    [[ -x "${p}" ]] && "${p}" --gsq "${GSQ}" --Nf "${Nf}" --nu0 "${NU0}" --nu1 "${NU1}" \
          --mass-re "${vr}" --mass-im "${vi}" --current axial --nhits "${NHITS}" "${ntflag[@]}" "${edflag[@]}"
  done
}

# ---------------- CLI ----------------
MODE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -g|--gpu)    GPU="$2"; shift 2 ;;
    --gpu=*)     GPU="${1#*=}"; shift ;;
    -H|--nhits)  NHITS="$2"; shift 2 ;;
    --nhits=*)   NHITS="${1#*=}"; shift ;;
    -T|--n-t0)   NT0="$2"; shift 2 ;;
    --n-t0=*)    NT0="${1#*=}"; shift ;;
    free|sweep)  MODE="$1"; shift ;;
    -h|--help)   echo "usage: $0 [--gpu N[,N...]] [--nhits N] [--n-t0 N] [free|sweep]   (default: free)"; exit 0 ;;
    *)           echo "unknown arg: $1"; echo "usage: $0 [--gpu N[,N...]] [--nhits N] [--n-t0 N] [free|sweep]"; exit 1 ;;
  esac
done
MODE="${MODE:-free}"

export CUDA_VISIBLE_DEVICES="${GPU}"
export NVIDIA_VISIBLE_DEVICES="${GPU}"
echo "### GPU(s): CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES} ###"
echo "### n_t0: ${NT0:-<program default 2>} ###"

# ---------------- modes ----------------
if [[ "${MODE}" == "free" ]]; then
  echo "### free-field (U=1) check ###"
  NHITS="${NHITS:-16}"   # free-field precision via hits (default 16; override --nhits N or NHITS= env)
  # one Nf is enough for the free field (gauge-independent); use Nf=2.
  Nf=2
  # massless valence
  run_vector "${Nf}" 0.0 0.0 ""
  run_axial  "${Nf}" 0.0 0.0 ""
  # a flavor-breaking valence (real) and a parity-breaking valence (imaginary)
  run_vector "${Nf}" 0.1 0.0 "";  run_axial "${Nf}" 0.1 0.0 ""
  run_vector "${Nf}" 0.0 0.1 "";  run_axial "${Nf}" 0.0 0.1 ""

elif [[ "${MODE}" == "sweep" ]]; then
  echo "### production sweep over 27 ensembles ###"
  NHITS="${NHITS:-1}"   # one hit per config by default (override --nhits N or NHITS= env)
  for Nf in "${NFS[@]}"; do
    # massless ensemble (sea = 0): vector + axial, valence 0
    ed=$(ens_dir "${Nf}" 0.0 0.0)
    run_vector "${Nf}" 0.0 0.0 "${ed}"
    run_axial  "${Nf}" 0.0 0.0 "${ed}"
    for m in "${MASSES[@]}"; do
      # flavor-breaking ensemble (sea real m_F): vector valence m_F ; axial valence 0
      ed=$(ens_dir "${Nf}" "${m}" 0.0)
      run_vector "${Nf}" "${m}" 0.0 "${ed}"
      run_axial  "${Nf}" 0.0   0.0 "${ed}"
      # parity-breaking ensemble (sea imag m_P): vector + axial valence m_P
      ed=$(ens_dir "${Nf}" 0.0 "${m}")
      run_vector "${Nf}" 0.0 "${m}" "${ed}"
      run_axial  "${Nf}" 0.0 "${m}" "${ed}"
    done
  done

else
  echo "usage: $0 [free|sweep]"; exit 1
fi
