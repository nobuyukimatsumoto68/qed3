#!/usr/bin/env bash
# rsync_scc_claude.sh -- pull the L4 REDO gauge configs (ckpoint_lat.* ONLY) from BU SCC to local.
# Source : nmatsum@scc1.bu.edu:/projectnb/qfe/nmatsum/qed3/src/production/
#          SCC now hosts ALL L4: massless Nf{2,4,6} (Nf4/6 handed off from FNAL 2026-07-22) + Nf2 massive.
# Dest   : /mnt/barracuda22/qed3/qed3/src/production/<ensemble-basename>/
# Run    : bash rsync_scc_claude.sh           (repeatable; rsync sends only new k's -- re-run as streams fill)
# Log    : rsync_scc_claude.log (appended, run delimited by date headers)
# SINGLE rsync invocation covering ALL ensembles = ONE ssh connection = ONE password prompt
# (never a per-directory loop). For fully prompt-free repeat runs, optionally add to ~/.ssh/config:
#   Host scc1.bu.edu
#     ControlMaster auto
#     ControlPath ~/.ssh/cm-%r@%h:%p
#     ControlPersist 8h
# (first run of the day asks once, later runs reuse the connection).
# -z: WAN compression. No -c (mtime+size is enough, configs are write-once).
# NEVER pull ckpoint_rng.* / logs / binaries. No rm anywhere in this script.
set -uo pipefail

DEST=/mnt/barracuda22/qed3/qed3/src/production
SRC='nmatsum@scc1.bu.edu:/projectnb/qfe/nmatsum/qed3/src/production/'
LOG="$DEST"/rsync_scc_claude.log

{
  echo "##################################################################"
  echo "### rsync SCC -> local  $(date)"
  echo "##################################################################"
} | tee -a "$LOG"

rsync -avz --prune-empty-dirs \
  --include='Nf*_*L4_hb*/' --include='ckpoint_lat.*' --exclude='*' \
  "$SRC" "$DEST"/ 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]
then
  echo "### ERROR: rsync exited with code $rc ###" | tee -a "$LOG"
  exit "$rc"
fi

{
  echo "### local ckpoint_lat counts per L4 ensemble (for the loc_ncfg blackboard column) ###"
  for d in "$DEST"/Nf*_*L4_hb*
  do
    [ -d "$d" ] || continue
    n=$(ls "$d"/ckpoint_lat.* 2>/dev/null | wc -l)
    echo "$(basename "$d")  $n"
  done
  echo "### done  $(date) ###"
} | tee -a "$LOG"
