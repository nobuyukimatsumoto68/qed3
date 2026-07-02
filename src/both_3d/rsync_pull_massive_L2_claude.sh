#!/usr/bin/env bash
# Pull the MASSIVE L=2 sea ensembles from the Fermilab HMC machine to local.
# Remote: nmatsum@lq.fnal.gov:/home/nmatsum/affine/<dir>
# Local : $DEST/<dir>   (default = this both_3d dir, where the other ensembles live)
#
# 24 ensembles: Nf{2,4,6} x 8 real masses (mIm=0), nt128 L2.
# NON-DESTRUCTIVE: no --delete (never removes local files); --partial + size/mtime skip make it
# fully RESUMABLE -- safe to re-run, it only fetches new / incomplete configs.
# Run it yourself (needs your ssh auth):  bash rsync_pull_massive_L2_claude.sh
# Output is tee'd to rsync_pull_massive_L2_claude.log .
set -u

REMOTE_HOST="nmatsum@lq.fnal.gov"
REMOTE_BASE="/home/nmatsum/affine"
DEST="/mnt/barracuda22/qed3/qed3/src/both_3d"
LOG="rsync_pull_massive_L2_claude.log"

# ssh transport: keepalives so long transfers survive idle stretches.
SSH_OPTS="ssh -o ServerAliveInterval=30 -o ServerAliveCountMax=10"

# rsync flags: archive, human, partial+progress, compress; NO --delete.
# LAT-ONLY: pull only the gauge configs ckpoint_lat.* ; skip ckpoint_rng.* and everything else.
RSYNC_OPTS=(-avhP --compress --include='ckpoint_lat.*' --include='*/' --exclude='*')

# Family B ONLY (correctly rescaled = more configs / more recent).  Family A
# {mRe0.005338, 0.026688, 0.053375, 0.106750} = WRONGLY RESCALED -> NOT pulled.
DIRS=(
  Nf2_gsq8.000000at0.200000nu01.000000mRe0.010572mIm0.000000nt128L2
  Nf2_gsq8.000000at0.200000nu01.000000mRe0.052862mIm0.000000nt128L2
  Nf2_gsq8.000000at0.200000nu01.000000mRe0.105725mIm0.000000nt128L2
  Nf2_gsq8.000000at0.200000nu01.000000mRe0.211450mIm0.000000nt128L2
  Nf4_gsq8.000000at0.200000nu01.000000mRe0.010572mIm0.000000nt128L2
  Nf4_gsq8.000000at0.200000nu01.000000mRe0.052862mIm0.000000nt128L2
  Nf4_gsq8.000000at0.200000nu01.000000mRe0.105725mIm0.000000nt128L2
  Nf4_gsq8.000000at0.200000nu01.000000mRe0.211450mIm0.000000nt128L2
  Nf6_gsq8.000000at0.200000nu01.000000mRe0.010572mIm0.000000nt128L2
  Nf6_gsq8.000000at0.200000nu01.000000mRe0.052862mIm0.000000nt128L2
  Nf6_gsq8.000000at0.200000nu01.000000mRe0.105725mIm0.000000nt128L2
  Nf6_gsq8.000000at0.200000nu01.000000mRe0.211450mIm0.000000nt128L2
)

echo "### rsync pull START  [$(date +%F_%H:%M:%S)]  ${#DIRS[@]} ensembles  ->  $DEST ###" | tee "$LOG"

n_ok=0
n_fail=0
failed=()
for d in "${DIRS[@]}"; do
  echo "### [$((n_ok+n_fail+1))/${#DIRS[@]}] $d  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
  rsync "${RSYNC_OPTS[@]}" -e "$SSH_OPTS" \
    "${REMOTE_HOST}:${REMOTE_BASE}/${d}/" "${DEST}/${d}/" 2>&1 | tee -a "$LOG"
  st=${PIPESTATUS[0]}
  if [ "$st" -eq 0 ]; then
    n_ok=$((n_ok+1))
    echo "###   OK  $d ###" | tee -a "$LOG"
  else
    n_fail=$((n_fail+1))
    failed+=("$d")
    echo "###   FAIL (rsync exit $st)  $d  -- continuing ###" | tee -a "$LOG"
  fi
done

echo "### rsync pull DONE  [$(date +%F_%H:%M:%S)]  ok=$n_ok  fail=$n_fail ###" | tee -a "$LOG"
if [ "$n_fail" -gt 0 ]; then
  echo "### FAILED dirs (re-run the script to retry; resumable): ###" | tee -a "$LOG"
  for d in "${failed[@]}"; do
    echo "###   $d ###" | tee -a "$LOG"
  done
fi
