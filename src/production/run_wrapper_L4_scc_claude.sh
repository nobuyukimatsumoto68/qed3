#!/bin/bash -l
# run_wrapper_L4_scc_claude.sh  (_scc, 2026-07-16, NM)
# =============================================================================================
# BU SCC *WRAPPER* for L=4 overlap-fermion QED3 production. Runs on the LOGIN node:
#   (1) module-loads the SCC toolchain,
#   (2) BUILDS one binary per (Nf, KMAX, GPU-arch)  -- Nf and KMAX are compile-time (-DNF, -DKMAX),
#       arch is -arch=sm_XX; separate executables per architecture (per NM),
#   (3) enumerates the requested ensemble SET(S), ROUND-ROBINs the ensembles across the requested
#       GPU architectures (each ensemble runs once; both V100 and A100 pools get used), PAIRS the
#       ensembles two-per-GPU, and submits a DEPENDENT CHAIN of batch jobs (run_L4_scc_claude.sh) per
#       pair -- each job = 1 GPU running 2 MPS-packed streams. The first chain link is 6h (H_RT_FIRST),
#       the rest H_RT; each link runs to its wall budget, checkpoints, exits, and the next link
#       (-hold_jid) resumes -- walking each ensemble to its KMAX target. The driver's graceful wall-time
#       limiter (max_sec, k_ckpoint=1 lossless checkpointing) is ALREADY in the .cu; the batch script
#       passes max_sec = h_rt - BUFFER_SEC.
#
# Finalized L=4 HMC params (tuning, shipped 2026-07-16; hasenbusch_ladder_claude.h / frozen_window_claude.h):
#   ladder {0.0, 0.4, 1.0}  steps {4,4,4}  tau(s_tot)=0.8  gauge MG=100
#   Zolotarev frozen (0.008, 5.0)  n_action=31  n_force=11  (two-op split-pole force, exact by Metropolis)
#   -DMIXED_FORCE (force-only mixed precision, ~1.35x).
# Driver = hmc_hasenbusch_block_scc_claude.cu (absolute SCC geometry paths).
# Refs: M. Hasenbusch hep-lat/0107019; B. Jegerlehner hep-lat/9612014; A. D. Kennedy hep-lat/0402038.
#
# Usage:  bash run_wrapper_L4_scc_claude.sh            # build + submit
#         DRYRUN=1 bash run_wrapper_L4_scc_claude.sh   # print qsub commands, do NOT submit
#         NOBUILD=1 bash run_wrapper_L4_scc_claude.sh  # skip build (binaries already present)
#         WHICH=massless bash ...                      # only one set (massless|massive|both; default both)
# =============================================================================================
set -u

SRCDIR=/projectnb/qfe/nmatsum/qed3/src/production
cd "$SRCDIR" || { echo "ERROR: cannot cd $SRCDIR"; exit 1; }
source /projectnb/qfe/nmatsum/qed3/env.sh          # unified env: cuda/12.8, gcc/13.2.0, $QED3_INC (repo Eigen)

# ---- ensemble sets (EDIT THESE) -----------------------------------------------------------------
# Each set is (Nf list) x (gsq list) x (mass list), with its own KMAX target (trajectories/ensemble).
# KMAX is COMPILE-TIME and is baked into the binary name, so the two sets get distinct binaries even
# at the same Nf (massless Nf2->k400 vs massive Nf2->k60). Auto-resumes; extend by rebuilding larger.
# default massless-only. NOTE: the massive run at the OLD masses {0.1,0.5,1.0,1.5} finished (k=59) but was at
# the WRONG masses (see MV_MASS below) -- the CORRECT-mass massive {0.1,0.2,0.3,0.4} still needs generating:
# run it explicitly with WHICH=massive. (massless L4 is unaffected -- mRe=0.)
WHICH=${WHICH:-massless}              # massless | massive | both

# FOCUS: Nf=2 only (per NM 2026-07-16 -- Nf4/6 dropped; ~2-3x cheaper force, priority physics).
ML_NF=${ML_NF:-"2"}                  # MASSLESS: params_L1L2_claude.md
ML_GSQ=${ML_GSQ:-"2.0 4.0 6.0"}
ML_MASS=${ML_MASS:-"0.0"}
ML_KMAX=${ML_KMAX:-400}

MV_NF=${MV_NF:-"2"}                  # MASSIVE: largest gsq per L = 6.0
MV_GSQ=${MV_GSQ:-"6.0"}
# physical masses mRe -- MUST match the L1/L2 (affine) runs for a same-physical-mass L=1,2,4 comparison.
# Full matched set = {0.1, 0.2, 0.3, 0.4}. CORRECTED 2026-07-18: the earlier {0.1,0.5,1.0,1.5} (from the stale
# params_massive_claude.md) did NOT match L1/L2 -- the mRe0.5/1.0/1.5 dirs are WRONG-mass leftovers.
# m=0.1 already exists at the correct mass (done, k=59), so this list is only the NEW ones to ADD.
# To regenerate the full set instead, set MV_MASS="0.1 0.2 0.3 0.4".  Run:  WHICH=massive NOBUILD=1 bash ...
MV_MASS=${MV_MASS:-"0.2 0.3 0.4"}
MV_KMAX=${MV_KMAX:-60}

NU0=${NU0:-1.0}
KRNG=${KRNG:-4}                      # keep every KRNG-th rng checkpoint (rng ckpt ~0.5 GB each; 4 -> ~1/4 the disk)

# ---- GPU architectures (EDIT / override) --------------------------------------------------------
# Build for each arch in ARCH_LIST; SUBMIT_ARCHS selects the node pools to actually submit to (ensembles
# are round-robined across them). gpu_c is the SGE minimum-compute-capability request per arch.
ARCH_LIST=${ARCH_LIST:-"sm_70 sm_80"}      # build these (V100=sm_70, A100=sm_80)
SUBMIT_ARCHS=${SUBMIT_ARCHS:-"$ARCH_LIST"} # submit to these node pools (round-robin)
PE_OMP=${PE_OMP:-16}                       # CPU slots per job
# streams packed per GPU: 1 = one ensemble/GPU, NO MPS (per NM 2026-07-18: the per-job MPS daemon caused
# frequent cudaMalloc/init "CUDA error" on shared ece/l40s nodes). Set NPACK=2 to re-enable MPS pairing.
NPACK=${NPACK:-1}

# GPU MODEL PIN (per NM 2026-07-22): this HMC is FP64-heavy (cuDoubleComplex). -l gpu_c is only a MINIMUM
# compute capability, so it lets jobs land on weak-FP64 L40S/A40 (1:64 FP64, e.g. L40S ~1.4 TFLOP/s) -- much
# slower than a strong-FP64 GPU. gpu_type is a requestable resource (==); pin each arch to its native,
# full-rate-FP64 model. Empty value => no pin (old behavior).
GPUT_SM70=${GPUT_SM70:-V100}         # sm_70 native, FP64 1:2 ~7 TFLOP/s (26 nodes)
GPUT_SM80=${GPUT_SM80:-A100}         # sm_80 native, FP64 1:2 ~9.7 TFLOP/s (9 nodes)
GPUT_SM90=${GPUT_SM90:-H200}         # sm_90 native, strongest FP64 ~34 TFLOP/s (8 nodes; Hopper, cuda>=12)

# ---- job chaining (each pair is a CHAIN of dependent jobs; each resumes the previous checkpoint) --
# One SGE job runs until its wall budget, checkpoints (k_ckpoint=1, lossless), exits; the next job in the
# SAME chain (-hold_jid) resumes. This walks each ensemble PAIR (fixed Nf/gsq/mass) up to its KMAX target.
H_RT_FIRST=${H_RT_FIRST:-12:00:00}         # walltime of the FIRST link in every chain (bumped 6h -> 12h, 2026-07-17: run is healthy)
H_RT=${H_RT:-12:00:00}                      # walltime of the subsequent chained links
N_CHAIN=${N_CHAIN:-4}                       # links per chain (over-provision is safe: once KMAX is reached the driver exits immediately)
BUFFER_SEC=${BUFFER_SEC:-900}               # wall budget = h_rt - BUFFER_SEC (covers startup + MPS teardown); driver also self-estimates 1.3x last traj

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
# Only Eigen (repo, via $QED3_INC from env.sh) + the build-dir headers are needed; no HDF5/HighFive/GSL.
INCLUDES="-I./includes/ ${QED3_INC}"
LDFLAGS="-lm"

# gpu_c (SGE minimum compute capability) per arch
gpuc_of () {
  case "$1" in
    sm_70) echo "7.0" ;;
    sm_75) echo "7.5" ;;
    sm_80) echo "8.0" ;;
    sm_86) echo "8.6" ;;
    sm_90) echo "9.0" ;;
    *) echo "7.0" ;;
  esac
}

# strong-FP64 GPU model (gpu_type) per arch -- pins the job to native, full-rate-FP64 hardware and OFF the
# weak-FP64 L40S/A40. Empty for arches without a mapping -> no gpu_type pin.
gput_of () {
  case "$1" in
    sm_70) echo "$GPUT_SM70" ;;
    sm_80) echo "$GPUT_SM80" ;;
    sm_90) echo "$GPUT_SM90" ;;
    *) echo "" ;;
  esac
}

# binary name per (Nf, KMAX, arch): Nf and KMAX are compile-time, arch is the SASS target
binname () { echo "hmc_L4_scc_Nf${1}_k${2}_${3}.out"; }

# L4 MD steps per Nf: Nf2 -> {4,4,4} (default 4); Nf4/6 -> {5,5,5} (-DL4_MDSTEP=5). Matches the FNAL Nf4/6
# builds (needed for the FNAL->SCC L4 config resume: same integrator tuning). Nf is 1:1 with the step count,
# and Nf is already in the binary name, so no extra name tag is needed. resume_fnal_L4_scc_impl_plan_claude.md.
mdstep_of () {
  if [ "$1" -eq 2 ]; then echo 4; else echo 5; fi
}

# ---- selected ensemble SET(S): each = "name|nf_list|gsq_list|mass_list|kmax" ---------------------
# Enumeration + pairing happen PER SET (below), so a massless (k400) stream is never MPS-packed with a
# massive (k60) stream -- like-cost streams share a GPU, no idle slot. Ensemble record = "Nf|KMAX|gsq|mass".
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

# ---- 1. build one binary per (Nf, KMAX) per arch ------------------------------------------------
build_all () {
  local arch nf kmax out key mds
  for arch in $ARCH_LIST
  do
    for key in "${!NKSEEN[@]}"
    do
      read -r nf kmax <<< "$key"
      out=$(binname "$nf" "$kmax" "$arch")
      mds=$(mdstep_of "$nf")
      echo "===== build $out (LREF=4 NF=$nf KMAX=$kmax L4_MDSTEP=$mds arch=$arch -DMIXED_FORCE) [$(date +%F_%H:%M:%S)] ====="
      $NVCC -arch="$arch" $NVCCBASE $INCLUDES \
        -DLREF=4 -DNF="$nf" -DKMAX="$kmax" -DKRNG="$KRNG" -DL4_MDSTEP="$mds" -DMIXED_FORCE \
        "$SRC" $LDFLAGS -o "$out" 2>&1 | tee "build_L4_scc_Nf${nf}_k${kmax}_${arch}_claude.log"
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

# ---- dependency safety ACROSS wrapper re-runs (per NM: "the script should catch all the dependencies") ----
# Snapshot every existing L4 job as "jobid <TAB> full_jobname". Before submitting a new chain, submit_chain
# looks here for any chain that already covers EITHER ensemble in the pair and -hold_jid's the new chain onto
# its tail -- so re-running the wrapper EXTENDS existing chains instead of starting a concurrent c0 that would
# write the same checkpoint dir (the collision that cost allocation on 2026-07-17). Match is by ensemble token
# Nf<nf>g<gsq>m<mass> (arch-independent, so it catches the ensemble even if its arch/partner changed).
EXIST_JOBS=$(qstat -u "$USER" -r 2>/dev/null | awk '
  /^[0-9]+ /{ jid=$1 }
  /Full jobname:/{ if($3 ~ /^L4[0-9]/) print jid"\t"$3 }
')
# echo the max existing job-id whose chain covers tokA (or tokB); empty if none. Trailing (__|_c) anchors the
# token to a chain-name boundary so m0.0 does not match m0.05, etc.
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

# ---- pair 2/GPU, submit a DEPENDENT CHAIN per pair -----------------------------------------------
# N_CHAIN jobs: link 0 (h_rt=H_RT_FIRST=6h), links 1.. (-hold_jid on the previous, h_rt=H_RT). If a chain for
# either ensemble already exists, link 0 -hold_jid's its TAIL (extend, never collide). Pair (Nf,gsq,mass) is
# fixed down the chain, so every link resumes the SAME ensembles' checkpoints.
submit_chain () {   # arch  A="Nf|KMAX|gsq|mass"  [B=...]
  local arch=$1 A=$2 B=${3:-}
  local gpuc gput
  gpuc=$(gpuc_of "$arch")
  gput=$(gput_of "$arch")
  local nfA kmaxA gsqA massA appA
  IFS='|' read -r nfA kmaxA gsqA massA <<< "$A"
  appA=$(binname "$nfA" "$kmaxA" "$arch")
  local vars="NU0=${NU0},APPA=${appA},GSQA=${gsqA},NFA=${nfA},MASSA=${massA}"
  local name="L4${arch#sm_}_Nf${nfA}g${gsqA}m${massA}"
  if [ -n "$B" ]
  then
    local nfB kmaxB gsqB massB appB
    IFS='|' read -r nfB kmaxB gsqB massB <<< "$B"
    appB=$(binname "$nfB" "$kmaxB" "$arch")
    vars="${vars},APPB=${appB},GSQB=${gsqB},NFB=${nfB},MASSB=${massB}"
    name="${name}__Nf${nfB}g${gsqB}m${massB}"
  fi

  # anchor onto any existing chain writing THIS pair's ensembles (extend, don't collide) -> c0 -hold_jid's it
  local tokA="Nf${nfA}g${gsqA}m${massA}" tokB=""
  [ -n "$B" ] && tokB="Nf${nfB}g${gsqB}m${massB}"
  local prev
  prev=$(existing_tail "$tokA" "$tokB")
  [ -n "$prev" ] && echo "  [anchor] ${name}: existing chain found -> new links hold on tail jid=$prev (no concurrent writer)"
  local c hrt maxsec jid
  for (( c=0; c<N_CHAIN; c++ ))
  do
    if [ "$c" -eq 0 ]; then hrt=$H_RT_FIRST; else hrt=$H_RT; fi
    maxsec=$(( $(hrt_to_sec "$hrt") - BUFFER_SEC ))
    [ "$maxsec" -lt 60 ] && maxsec=60
    local qsub_cmd=( qsub -terse -N "${name}_c${c}"
                     -l "gpus=1" -l "gpu_c=${gpuc}" -l "h_rt=${hrt}" -pe omp "${PE_OMP}" )
    [ -n "$gput" ] && qsub_cmd+=( -l "gpu_type=${gput}" )
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

# ---- 3. per SET: round-robin ITS ensembles across arches, pair 2/GPU WITHIN the set, submit chains -
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
  echo "===== set $sname : ${#recs[@]} ensembles (Nf {$nfs} x gsq {$gsqs} x m {$masses}, KMAX=$kmax) ====="

  # round-robin THIS set's ensembles across the submit arches
  unset SB
  declare -A SB
  ri=0
  for rec in "${recs[@]}"
  do
    a=${ARCHES[$(( ri % n_arch ))]}
    SB[$a]+="${rec}"$'\n'
    ri=$(( ri + 1 ))
  done

  # within each arch: pair 2/GPU (both streams from THIS set) and submit a dependent chain per pair
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
        submit_chain "$arch" "${r[$i]}" "${r[$((i+1))]:-}"   # 2 streams/GPU (MPS pair)
        i=$(( i + 2 ))
      else
        submit_chain "$arch" "${r[$i]}"                       # 1 stream/GPU (no MPS)
        i=$(( i + 1 ))
      fi
    done
  done
done

echo "===== wrapper done (DRYRUN=$DRYRUN) [$(date +%F_%H:%M:%S)] ====="
