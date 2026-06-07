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
NHITS="${NHITS:-}"   # hits per program call; per-mode default applied after CLI (free 16, sweep/ens 1)
NT0="${NT0:-}"       # # source-time origins; empty => program default (2); override --n-t0/-T or NT0= env
NINTER="${NINTER:-}" # ensemble config stride k=0,N,2N,...; empty => program default (10); --ninter/-I or NINTER=
COMPONENT="${COMPONENT:-all}"  # which program to run: all (combined corr_) | conn | disc; --all/--conn/--disc or COMPONENT=

# 'ens' mode (single explicit ensemble) knobs -- set via CLI:
ENS_DIR="${ENS_DIR:-}"   # explicit sea config directory (must end with /)
ENS_NF="${ENS_NF:-2}"    # Nf label for the ens-mode run
VAL_RE="${VAL_RE:-0.0}"  # valence mass Re for the ens-mode run
VAL_IM="${VAL_IM:-0.0}"  # valence mass Im for the ens-mode run

NFS=(2 4 6)
MASSES=(0.01 0.05 0.1 0.2)

# program binaries.  Three production programs (one per --all/--conn/--disc):
# C6f-c: deploy the t-BLOCKED sink binaries (2.99x whole-jj; same output dirs => analysis unchanged).
# Old lines kept commented as reference fallback (validated bit-identical/within-1e-10 to the block_t).
# CORR="${BINDIR}/jj_corr_claude.o"             # --all : pre-block_t (jj_corr; mrhs = jj_corr_mrhs_claude.o)
CORR="${BINDIR}/jj_corr_block_t_claude.o"       # --all : COMBINED conn + folded disc, t-blocked sink (corr_ dir; 2.99x)
CONN="${BINDIR}/jj_conn_correlators_claude.o"   # --conn: connected only (vector tp/sp/ylm + axial tp/sp; conn_ dir)
DISC="${BINDIR}/jj_disc_claude.o"               # --disc: disconnected only (vector; cheap, no conn solves; disc_ dir)
                                                #   NOTE: --all (jj_corr_block_t) ALREADY includes disc (folded in,
                                                #   block-t'd); jj_disc is only the disc-ONLY cheap path. jj_disc_block_t
                                                #   exists but is REDUNDANT for --all production (left parked, not used here).
                                                #   (none take --current; axial legs auto-selected by mass)
# Legacy per-projection standalone programs (kept as reference; driven by run_vector/run_axial below):
CONN_SP="${BINDIR}/jj_conn_spproj_claude.o"             # vector sp
CONN_SP_AXIAL="${BINDIR}/jj_conn_spproj_axial_claude.o" # axial  sp
CONN_TP="${BINDIR}/jj_conn_tpproj_claude.o"             # vector tp
CONN_TP_AXIAL="${BINDIR}/jj_conn_tpproj_axial_claude.o" # axial  tp
CONN_YLM="${BINDIR}/jj_conn_ylmproj_claude.o"           # vector ylm (m-summed tower)

CONN_PROGS=("${CONN_TP}" "${CONN_SP}" "${CONN_YLM}")     # vector projections (legacy)
CONN_AXIAL_PROGS=("${CONN_TP_AXIAL}" "${CONN_SP_AXIAL}")  # axial projections (legacy)

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

# PRODUCTION (--all): COMBINED connected + disc for one (ensemble, valence).  ONE valence mass suffices:
# the binary computes vector tp/sp/ylm AND axial tp/sp AND the disc traces, with internal flavor/parity
# logic selecting the axial legs (flavor m_F -> massless D_ov axial; parity m_P -> tilde sink).
run_corr() {  # args: Nf valRe valIm [ens_dir or "" for free]  -> data_<ESNID>/corr_nt0<N>_nhits<H>/
  local Nf="$1" vr="$2" vi="$3" ed="$4"
  local edflag=(); [[ -n "${ed}" ]]     && edflag=(--ens-dir "${ed}")
  local ntflag=(); [[ -n "${NT0}" ]]    && ntflag=(--n-t0 "${NT0}")       # empty => program default (2)
  local niflag=(); [[ -n "${NINTER}" ]] && niflag=(--ninter "${NINTER}")  # empty => program default (10)
  [[ -x "${CORR}" ]] && "${CORR}" --gsq "${GSQ}" --Nf "${Nf}" --nu0 "${NU0}" --nu1 "${NU1}" \
        --mass-re "${vr}" --mass-im "${vi}" --nhits "${NHITS}" "${ntflag[@]}" "${niflag[@]}" "${edflag[@]}"
}

# --conn: CONNECTED-only program (no disc) for one (ensemble, valence).
run_conn() {  # args: Nf valRe valIm [ens_dir or "" for free]  -> data_<ESNID>/conn_nt0<N>_nhits<H>/
  local Nf="$1" vr="$2" vi="$3" ed="$4"
  local edflag=(); [[ -n "${ed}" ]]     && edflag=(--ens-dir "${ed}")
  local ntflag=(); [[ -n "${NT0}" ]]    && ntflag=(--n-t0 "${NT0}")       # empty => program default (2)
  local niflag=(); [[ -n "${NINTER}" ]] && niflag=(--ninter "${NINTER}")  # empty => program default (10)
  [[ -x "${CONN}" ]] && "${CONN}" --gsq "${GSQ}" --Nf "${Nf}" --nu0 "${NU0}" --nu1 "${NU1}" \
        --mass-re "${vr}" --mass-im "${vi}" --nhits "${NHITS}" "${ntflag[@]}" "${niflag[@]}" "${edflag[@]}"
}

# --disc: DISCONNECTED-only dump (vector; origin-agnostic, no --n-t0).  Cheap path for many configs.
run_disc() {  # args: Nf valRe valIm [ens_dir or "" for free]  -> data_<ESNID>/disc_nhits<H>/
  local Nf="$1" vr="$2" vi="$3" ed="$4"
  local edflag=(); [[ -n "${ed}" ]]     && edflag=(--ens-dir "${ed}")
  local niflag=(); [[ -n "${NINTER}" ]] && niflag=(--ninter "${NINTER}")  # empty => program default (10)
  [[ -x "${DISC}" ]] && "${DISC}" --gsq "${GSQ}" --Nf "${Nf}" --nu0 "${NU0}" --nu1 "${NU1}" \
        --mass-re "${vr}" --mass-im "${vi}" --nhits "${NHITS}" "${niflag[@]}" "${edflag[@]}"
}

# Dispatcher honouring COMPONENT {all|conn|disc} (--all default / --conn / --disc).
run_connected() {  # args: Nf valRe valIm [ens_dir or "" for free]
  case "${COMPONENT}" in
    all)  run_corr "$@" ;;   # combined corr_ (efficient: shares phi' + sink across conn & disc)
    conn) run_conn "$@" ;;   # connected-only conn_
    disc) run_disc "$@" ;;   # disconnected-only disc_ (cheap)
  esac
}

# ---------------- CLI ----------------
USAGE="usage: $0 [--gpu N[,N...]] [--nhits N] [--n-t0 N] [--ninter N] [--all|--conn|--disc] [free|sweep|ens]   (default: free, --all)
  ens mode (single explicit ensemble): also pass --ens-dir <path/> [--Nf N] [--mass-re X] [--mass-im Y]
  --all  (default) run the COMBINED jj_corr binary (corr_ dir = connected + folded disc; efficient: shared phi'/sink)
  --conn run ONLY the connected binary jj_conn_correlators (conn_ dir)
  --disc run ONLY the disconnected binary jj_disc (disc_ dir; cheap, no connected solves -- for many configs)"
MODE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -g|--gpu)    GPU="$2"; shift 2 ;;
    --gpu=*)     GPU="${1#*=}"; shift ;;
    -H|--nhits)  NHITS="$2"; shift 2 ;;
    --nhits=*)   NHITS="${1#*=}"; shift ;;
    -T|--n-t0)   NT0="$2"; shift 2 ;;
    --n-t0=*)    NT0="${1#*=}"; shift ;;
    -I|--ninter) NINTER="$2"; shift 2 ;;
    --ninter=*)  NINTER="${1#*=}"; shift ;;
    --all)       COMPONENT="all";  shift ;;
    --conn)      COMPONENT="conn"; shift ;;
    --disc)      COMPONENT="disc"; shift ;;
    --only)      COMPONENT="$2"; shift 2 ;;   # alias: --only all|conn|disc
    --only=*)    COMPONENT="${1#*=}"; shift ;;
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
case "${COMPONENT}" in all|conn|disc) ;; *) echo "component must be all|conn|disc (got '${COMPONENT}')"; exit 1 ;; esac

export CUDA_VISIBLE_DEVICES="${GPU}"
export NVIDIA_VISIBLE_DEVICES="${GPU}"
echo "### GPU(s): CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES} ###"
echo "### n_t0: ${NT0:-<program default 2>}   ninter: ${NINTER:-<program default 10>}   component: ${COMPONENT} ###"

# ---------------- modes ----------------
if [[ "${MODE}" == "free" ]]; then
  echo "### free-field (U=1) check ###"
  NHITS="${NHITS:-16}"   # free-field precision via hits (default 16; override --nhits N or NHITS= env)
  # one Nf is enough for the free field (gauge-independent); use Nf=2.
  Nf=2
  run_connected "${Nf}" 0.0 0.0 ""    # massless valence
  run_connected "${Nf}" 0.1 0.0 ""    # flavor-breaking valence (real)
  run_connected "${Nf}" 0.0 0.1 ""    # parity-breaking valence (imaginary)

elif [[ "${MODE}" == "sweep" ]]; then
  echo "### production sweep over 27 ensembles ###"
  NHITS="${NHITS:-1}"   # one hit per config by default (override --nhits N or NHITS= env)
  for Nf in "${NFS[@]}"; do
    # massless ensemble (sea = 0): valence 0
    run_connected "${Nf}" 0.0 0.0 "$(ens_dir "${Nf}" 0.0 0.0)"
    for m in "${MASSES[@]}"; do
      # flavor-breaking ensemble (sea real m_F): valence m_F (axial auto-massless inside the unified binary)
      run_connected "${Nf}" "${m}" 0.0 "$(ens_dir "${Nf}" "${m}" 0.0)"
      # parity-breaking ensemble (sea imag m_P): valence m_P
      run_connected "${Nf}" 0.0 "${m}" "$(ens_dir "${Nf}" 0.0 "${m}")"
    done
  done

elif [[ "${MODE}" == "ens" ]]; then
  : "${ENS_DIR:?ens mode needs --ens-dir <path/>}"
  [[ "${ENS_DIR}" != */ ]] && ENS_DIR="${ENS_DIR}/"   # binary expects a trailing slash
  NHITS="${NHITS:-1}"
  echo "### single ensemble: ${ENS_DIR}  (Nf=${ENS_NF} valence mRe=${VAL_RE} mIm=${VAL_IM}) ###"
  run_connected "${ENS_NF}" "${VAL_RE}" "${VAL_IM}" "${ENS_DIR}"

else
  echo "${USAGE}"; exit 1
fi
