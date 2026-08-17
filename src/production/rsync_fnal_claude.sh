#!/usr/bin/env bash
# rsync_fnal_claude.sh -- pull the REDO gauge configs (ckpoint_lat.* ONLY) from FNAL to local.
# Source : nmatsum@lq.fnal.gov:/lustre2/affine/redo/   (L1/L2 all streams + L4 Nf4/6 when launched)
# Dest   : /mnt/barracuda22/qed3/qed3/src/production/<ensemble-basename>/
# Run    : bash rsync_fnal_claude.sh          (repeatable; rsync sends only new k's -- re-run as streams fill)
# Log    : rsync_fnal_claude.log (appended, run delimited by date headers)
# Auth   : FNAL kerberos/ssh (kinit first if your ticket expired).
# SINGLE rsync invocation = one connection/prompt. -z: WAN compression. No -c (mtime+size is enough,
# configs are write-once). NEVER pull ckpoint_rng.* / logs / binaries. No rm anywhere in this script.
set -uo pipefail

DEST=/mnt/barracuda22/qed3/qed3/src/production
SRC='nmatsum@lq.fnal.gov:/lustre2/affine/redo/'
LOG="$DEST"/rsync_fnal_claude.log

{
  echo "##################################################################"
  echo "### rsync FNAL -> local  $(date)"
  echo "##################################################################"
} | tee -a "$LOG"

rsync -avz --prune-empty-dirs \
  --include='Nf*_*nt128L*_hb*/' --include='ckpoint_lat.*' --exclude='*' \
  "$SRC" "$DEST"/ 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]
then
  echo "### ERROR: rsync exited with code $rc ###" | tee -a "$LOG"
  exit "$rc"
fi

{
  echo "### local ckpoint_lat counts per ensemble (for the loc_ncfg blackboard column) ###"
  for d in "$DEST"/Nf*_*nt128L*_hb*
  do
    [ -d "$d" ] || continue
    n=$(ls "$d"/ckpoint_lat.* 2>/dev/null | wc -l)
    echo "$(basename "$d")  $n"
  done
  echo "### done  $(date) ###"
} | tee -a "$LOG"

# ---------------------------------------------------------------------------
# CONN stride-2 datasets (added 2026-08-13).  FNAL runs the L=1 slice of the
# three-site stride-2 completion (conn_s2_fnal_L1_impl_plan_claude.md): the four
# residue classes k = 3,5,7,9 (mod 10), which are DISJOINT from the local
# stride-10 grid (k = 1 mod 10).  The pull is therefore purely ADDITIVE.
#
# SEPARATE SRC: the FNAL jobs cd to /lustre2/affine/redo/conn_s2 (WORKDIR in
# run_conn_s2_fnal_claude.sh:24), so the data_* trees sit one level BELOW the
# config root above.  Pulling from .../redo/ would land them as conn_s2/data_*/
# locally, which no reader would find -- hence the deeper source root here.
#
# --ignore-existing: never overwrite a conn h5 that already exists locally.  The
# residue classes should make collisions impossible, so an overlap would mean a
# partition bug; the post-pull census below reports it rather than silently
# clobbering the production stride-10 data.
#
# This is a SECOND rsync invocation = a second connection/prompt.  For a
# prompt-free repeat run add to ~/.ssh/config:
#   Host lq.fnal.gov
#     ControlMaster auto
#     ControlPath ~/.ssh/cm-%r@%h:%p
#     ControlPersist 8h
# ---------------------------------------------------------------------------
SRC_CONN='nmatsum@lq.fnal.gov:/lustre2/affine/redo/conn_s2/'

{
  echo "### rsync FNAL conn stride-2 -> local  $(date)"
  echo "### src = $SRC_CONN"
} | tee -a "$LOG"

rsync -avz --prune-empty-dirs --ignore-existing \
  --include='data_Nf*/' \
  --include='data_Nf*/corr_ylm_conn_t00_nhits1_s1/' \
  --include='data_Nf*/corr_ylm_conn_t00_nhits1_s1/corr.*.h5' \
  --exclude='*' \
  "$SRC_CONN" "$DEST"/ 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]
then
  echo "### ERROR: conn rsync exited with code $rc ###" | tee -a "$LOG"
  exit "$rc"
fi

{
  echo "### conn h5 census by residue class k mod 10 ###"
  echo "### (k%10==1 = the production stride-10 grid; 3,5,7,9 = the stride-2 completion) ###"
  printf "%-58s %6s %6s %6s %6s %6s %7s\n" "ensemble" "k%10=1" "3" "5" "7" "9" "total"
  for dd in "$DEST"/data_Nf*/corr_ylm_conn_t00_nhits1_s1
  do
    [ -d "$dd" ] || continue
    ens=$(basename "$(dirname "$dd")")
    find "$dd" -maxdepth 1 -name 'corr.*.h5' -printf '%f\n' 2>/dev/null \
      | sed 's/^corr\.\([0-9]*\)\..*/\1/' \
      | awk -v e="$ens" '
          {r = $1 % 10; c[r]++; t++}
          END {printf "%-58s %6d %6d %6d %6d %6d %7d\n", e, c[1], c[3], c[5], c[7], c[9], t}'
  done
  echo "### done  $(date) ###"
} | tee -a "$LOG"
