#!/usr/bin/env bash
# Condensates <sigma_PS>,<sigma_FS> for the HEAVY sea masses, all L=1,2,4 (Nf2 only).
# Uses the DEDICATED driver jj_scalar_condensate_eo_stoch_claude.cu (e/o timeslice + spin dilution, 8 solves/
# config = 2 legs x 4 patterns) -- REPLACES the old run_condensate_heavy_claude.sh, which got the condensate as
# a byproduct of the full disc driver (128 solves/config).  Output -> data_<ens>_vmRe<m>.../corr_condensate_eo_nhits1/
# with h0/condensate/{etadag_xi, xidag_1mDdag_eta}.  Analysis: <sigma_PS>=2Re etadag_xi,
# <sigma_FS>=etadag_xi - xidag_1mDdag_eta.  Per-L binaries (-DN_REFINE_CLI={1,2,4}); 2-wide MPS pool.
# NOTE: the driver's h5 skip-gate keys on xidag_1mDdag_eta, so files from an earlier (old-driver) run auto-
# recompute here (no rm).  The .cu changed -> need_build rebuilds the binaries on this run.
#
# Heavy set (Nf2, gsq8): CLI mass_re = physical m (R=1, same at every L); dir mRe is %f-rounded.
#   phys 0.4228996588195 (dir mRe0.422900), 0.845799317639 (mRe0.845799), 1.2686989764584 (mRe1.268699)
#
# Run detached:
#   nohup bash run_condensate_eo_heavy_claude.sh > run_condensate_eo_heavy_claude.log 2>&1 &
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

GPU="${GPU:-1}"
MAXJOBS="${MAXJOBS:-2}"
KMIN_TRAJ=40          # thermalization: skip configs numbered below this
SAMP=4                # sample every SAMP-th available config (stride = SAMP * base spacing)
NHITS=1

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
SRC=jj_scalar_condensate_eo_stoch_claude.cu

# ---- ensure MPS daemon up (2-packing needs it, else jobs serialize via context switch) ----
if pgrep -f nvidia-cuda-mps-control >/dev/null
then
  echo "### MPS daemon: already running ###"
else
  echo "### MPS daemon not up -- starting nvidia-cuda-mps-control -d ###"
  nvidia-cuda-mps-control -d
  for i in 1 2 3 4 5
  do
    pgrep -f nvidia-cuda-mps-control >/dev/null && break
    sleep 1
  done
fi
pgrep -f nvidia-cuda-mps-control >/dev/null \
  || { echo "### ERROR: MPS daemon failed to start -- aborting ###"; exit 1; }

# ---- build the three per-L binaries (skip if up-to-date; FORCE_BUILD=1 to force) ----
need_build () {   # $1 = output binary
  [ -n "${FORCE_BUILD:-}" ] && return 0
  [ ! -f "$1" ] && return 0
  find . -maxdepth 2 \( -name '*.cu' -o -name '*.h' \) -newer "$1" -print -quit 2>/dev/null | grep -q . && return 0
  return 1
}
for L in 1 2 4
do
  BIN="jj_scalar_condensate_eo_L${L}.o"
  if need_build "$BIN"
  then
    echo "### compile L${L} (-DN_REFINE_CLI=${L}) -> $BIN  [$(date +%F_%H:%M:%S)] ###"
    $NVCC $NVCCBASE -DN_REFINE_CLI=${L} $INCLUDES $LDFLAGS "$SRC" -o "$BIN" \
      || { echo "### L${L} BUILD FAILED ###"; exit 1; }
  else
    echo "### L${L} binary up-to-date, skip (FORCE_BUILD=1 to rebuild) ###"
  fi
done

# ---- ensembles: "<ens-dir>|<exact valence mass>|<L>"  (valence = sea = physical m) ----
ENS_LIST=(
  "Nf2_gsq8.000000at0.200000nu01.000000mRe0.422900mIm0.000000nt128L1|0.4228996588195|1"
  "Nf2_gsq8.000000at0.200000nu01.000000mRe0.845799mIm0.000000nt128L1|0.845799317639|1"
  "Nf2_gsq8.000000at0.200000nu01.000000mRe1.268699mIm0.000000nt128L1|1.2686989764584|1"
  "Nf2_gsq8.000000at0.200000nu01.000000mRe0.422900mIm0.000000nt128L2|0.4228996588195|2"
  "Nf2_gsq8.000000at0.200000nu01.000000mRe0.845799mIm0.000000nt128L2|0.845799317639|2"
  "Nf2_gsq8.000000at0.200000nu01.000000mRe1.268699mIm0.000000nt128L2|1.2686989764584|2"
  "Nf2_gsq8.000000at0.200000nu01.000000mRe0.422900mIm0.000000nt128L4|0.4228996588195|4"
  "Nf2_gsq8.000000at0.200000nu01.000000mRe0.845799mIm0.000000nt128L4|0.845799317639|4"
  "Nf2_gsq8.000000at0.200000nu01.000000mRe1.268699mIm0.000000nt128L4|1.2686989764584|4"
)

# ---- run one condensate job (ens|mass|L): auto-detect config range, then run on GPU.  Backgrounded. ----
run_job () {
  local spec="$1"
  local dir mass L
  IFS='|' read -r dir mass L <<< "$spec"

  local mretag bin LOG
  mretag=$(printf '%s' "$dir" | grep -oE 'mRe[0-9]+\.[0-9]+' | head -1 | sed 's/mRe//')
  bin="jj_scalar_condensate_eo_L${L}.o"
  LOG="cond_eo_L${L}_mRe${mretag}_claude.log"

  # available config numbers (sorted ascending)
  local ks
  mapfile -t ks < <(ls "$dir"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n)
  if [ "${#ks[@]}" -eq 0 ]
  then
    echo "### cond L${L} mRe${mretag} SKIP: no ckpoint_lat in '$dir'  [$(date +%F_%H:%M:%S)] ###" >> "$LOG"
    return 0
  fi

  local first last base kmin stride kmax nconf k
  first="${ks[0]}"
  last="${ks[-1]}"
  if [ "${#ks[@]}" -ge 2 ]; then base=$(( ks[1] - ks[0] )); else base=1; fi
  if [ "$base" -le 0 ]; then base=1; fi
  kmin=""
  for k in "${ks[@]}"
  do
    if [ "$k" -ge "$KMIN_TRAJ" ]; then kmin="$k"; break; fi
  done
  if [ -z "$kmin" ]; then kmin="$first"; fi
  stride=$(( SAMP * base ))
  kmax=$(( last + 1 ))
  nconf=$(( (last - kmin) / stride + 1 ))

  {
    echo "### cond L${L} mRe${mretag}  ens=$dir  mass=$mass  [$(date +%F_%H:%M:%S)] ###"
    echo "###   detect: avail=${#ks[@]} (k=$first..$last base=$base) -> kmin=$kmin stride=$stride kmax=$kmax (~$nconf cfg) ###"
    CUDA_VISIBLE_DEVICES=$GPU ./"$bin" --ens-dir "$dir/" --kmin "$kmin" --stride "$stride" --kmax "$kmax" \
      --nhits "$NHITS" --mass-re "$mass"
    echo "### cond L${L} mRe${mretag} done (status $?)  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

# ---- 2-wide MPS pool over the 9 heavy ensembles ----
echo "### START condensate-eo (9 ensembles, ${MAXJOBS} at a time on GPU ${GPU})  [$(date +%F_%H:%M:%S)] ###"
for spec in "${ENS_LIST[@]}"
do
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]
  do
    wait -n
  done
  echo "### launch: $spec  [$(date +%F_%H:%M:%S)] ###"
  run_job "$spec" &
done
wait
echo "### ALL condensate-eo done  [$(date +%F_%H:%M:%S)] ###"
