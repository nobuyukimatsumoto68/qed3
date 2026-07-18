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
