#!/usr/bin/env bash
# L=2 local-current Y_lm tower (conn + disc) over the MASSLESS Nf6 + the 24 transferred MASSIVE
# (real m_F) sea ensembles.  Valence mass = SEA mass (unitary point) = the mRe in the dir name.
# Massless Nf6 = no --mass-re.  Parity (m_P) NOT supported (programs assert !parity); all here are mIm=0.
#
# PACKING: at most MAXJOBS=2 programs on GPU 0 at once, co-resident under CUDA MPS.  Each ensemble
# contributes 2 jobs (conn, disc) that run independently (different output subdirs).  A 2-wide pool
# (wait -n) keeps exactly 2 procs busy across the whole job list.
#
# CONFIG SAMPLING is AUTO-DETECTED per ensemble (the transferred chains differ in length/numbering;
# the program BREAKS at the first missing config, so the sampled grid must hit existing files):
#   - read sorted ckpoint_lat.<k> ; base = spacing between the first two ; kmin = first k >= KMIN_TRAJ
#     (thermalization skip) ; stride = SAMP * base ; kmax = last+1.  Mirrors the L2-massless choice
#     (KMIN_TRAJ=40, SAMP=4 => stride 4 when base=1).  Tune KMIN_TRAJ / SAMP below.
#
# Output auto-tags from --ens-dir (+ valence mass) ->
#   data_<ens>_vmRe<m>vmIm0.000000/{corr_ylm_conn_t00_nhits1_s1, corr_ylm_disc_tb2_nhits1}/
# Per-config atomic .h5 + resume-skip (safe to re-run; resumes, skips complete configs).
# MPS must be up:  nvidia-cuda-mps-control -d .   Per-job logs: ylm_{conn,disc}_L2_<tag>_claude.log .
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

GPU=0
MAXJOBS=2
KMIN_TRAJ=40          # thermalization: skip configs numbered below this
SAMP=4                # sample every SAMP-th available config (stride = SAMP * base spacing)
NHITS=1
DISC_TB=2

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -DN_REFINE_CLI=2 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

CONN=jj_local_ylm_conn_stoch_L2.o
DISC=jj_local_ylm_disc_stoch_L2.o

# ---- build once (L2 binaries; -DN_REFINE_CLI=2) ----
echo "### compile L2 (-DN_REFINE_CLI=2) -> $CONN, $DISC  [$(date +%F_%H:%M:%S)] ###"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_local_ylm_conn_stoch_claude.cu -o "$CONN" || { echo "### CONN BUILD FAILED ###"; exit 1; }
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_local_ylm_disc_stoch_claude.cu -o "$DISC" || { echo "### DISC BUILD FAILED ###"; exit 1; }

# ---- ensembles: "<ens-dir>|<exact valence mass>"  (massless = empty mass) ----
# Family B = correctly-rescaled (more configs / more recent).  Dir names are %f-rounded (mRe0.010572)
# but the TRUE sea masses carry full precision -> pass exact below; output still tags vmRe0.010572.
# Family A {mRe0.005338, 0.026688, 0.053375, 0.106750} = WRONGLY RESCALED -> DROPPED.
ENS_LIST=(
  "Nf6_gsq8.000000at0.200000nu01.000000nt128L2|"
  "Nf2_gsq8.000000at0.200000nu01.000000mRe0.010572mIm0.000000nt128L2|0.0105724914705"
  "Nf2_gsq8.000000at0.200000nu01.000000mRe0.052862mIm0.000000nt128L2|0.0528624573524"
  "Nf2_gsq8.000000at0.200000nu01.000000mRe0.105725mIm0.000000nt128L2|0.1057249147049"
  "Nf2_gsq8.000000at0.200000nu01.000000mRe0.211450mIm0.000000nt128L2|0.2114498294097"
  "Nf4_gsq8.000000at0.200000nu01.000000mRe0.010572mIm0.000000nt128L2|0.0105724914705"
  "Nf4_gsq8.000000at0.200000nu01.000000mRe0.052862mIm0.000000nt128L2|0.0528624573524"
  "Nf4_gsq8.000000at0.200000nu01.000000mRe0.105725mIm0.000000nt128L2|0.1057249147049"
  "Nf4_gsq8.000000at0.200000nu01.000000mRe0.211450mIm0.000000nt128L2|0.2114498294097"
  "Nf6_gsq8.000000at0.200000nu01.000000mRe0.010572mIm0.000000nt128L2|0.0105724914705"
  "Nf6_gsq8.000000at0.200000nu01.000000mRe0.052862mIm0.000000nt128L2|0.0528624573524"
  "Nf6_gsq8.000000at0.200000nu01.000000mRe0.105725mIm0.000000nt128L2|0.1057249147049"
  "Nf6_gsq8.000000at0.200000nu01.000000mRe0.211450mIm0.000000nt128L2|0.2114498294097"
)

# ---- run one job (prog|ens): auto-detect config range, then run on GPU0.  Backgrounded. ----
run_job () {
  local spec="$1"
  local prog dir mass
  IFS='|' read -r prog dir mass <<< "$spec"

  local nf mretag tag bin LOG
  nf=$(printf '%s' "$dir" | grep -oE '^Nf[0-9]+')
  mretag=$(printf '%s' "$dir" | grep -oE 'mRe[0-9]+\.[0-9]+' | head -1 | sed 's/mRe//')
  if [ -n "$mretag" ]; then tag="${nf}_mRe${mretag}"; else tag="${nf}_massless"; fi
  if [ "$prog" = conn ]; then bin="$CONN"; else bin="$DISC"; fi
  LOG="ylm_${prog}_L2_${tag}_claude.log"

  local mass_flag=()
  if [ -n "$mass" ]; then mass_flag=(--mass-re "$mass"); fi

  local prog_flags=()
  if [ "$prog" = conn ]; then
    prog_flags=(--t0 0 --spin-dilution)
  else
    prog_flags=(--disc-tblock "$DISC_TB")
  fi

  # available config numbers (sorted ascending)
  local ks
  mapfile -t ks < <(ls "$dir"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n)
  if [ "${#ks[@]}" -eq 0 ]; then
    echo "### $prog $tag SKIP: no ckpoint_lat in '$dir' (not transferred?)  [$(date +%F_%H:%M:%S)] ###" >> "$LOG"
    return 0
  fi

  local first last base kmin stride kmax nconf k
  first="${ks[0]}"
  last="${ks[-1]}"
  if [ "${#ks[@]}" -ge 2 ]; then base=$(( ks[1] - ks[0] )); else base=1; fi
  if [ "$base" -le 0 ]; then base=1; fi
  kmin=""
  for k in "${ks[@]}"; do
    if [ "$k" -ge "$KMIN_TRAJ" ]; then kmin="$k"; break; fi
  done
  if [ -z "$kmin" ]; then kmin="$first"; fi
  stride=$(( SAMP * base ))
  kmax=$(( last + 1 ))
  nconf=$(( (last - kmin) / stride + 1 ))

  {
    echo "### $prog $tag  ens=$dir  [$(date +%F_%H:%M:%S)] ###"
    echo "###   detect: avail=${#ks[@]} (k=$first..$last base=$base) -> kmin=$kmin stride=$stride kmax=$kmax (~$nconf cfg) mass='${mass:-massless}' ###"
    CUDA_VISIBLE_DEVICES=$GPU ./"$bin" --ens-dir "$dir/" --kmin "$kmin" --stride "$stride" --kmax "$kmax" \
      --nhits "$NHITS" "${prog_flags[@]}" "${mass_flag[@]}"
    echo "### $prog $tag done (status $?)  [$(date +%F_%H:%M:%S)] ###"
  } >> "$LOG" 2>&1
}

# ---- build the flat 2-wide job list: conn + disc per ensemble ----
JOBS=()
for ens in "${ENS_LIST[@]}"; do
  JOBS+=("conn|$ens")
  JOBS+=("disc|$ens")
done

echo "### START ylm L2 (massless Nf6 + 12 massive, Family B)  ${#JOBS[@]} jobs, $MAXJOBS at a time on GPU $GPU  [$(date +%F_%H:%M:%S)] ###"
for spec in "${JOBS[@]}"; do
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do
    wait -n
  done
  echo "### launch: $spec  [$(date +%F_%H:%M:%S)] ###"
  run_job "$spec" &
done
wait
echo "### ALL ylm L2 massive+Nf6 done  [$(date +%F_%H:%M:%S)] ###"
