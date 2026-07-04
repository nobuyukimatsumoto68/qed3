#!/usr/bin/env bash
# rsync_pull_configs_claude.sh -- pull gauge CONFIGS (ckpoint_lat.*) for every ensemble listed in
# the canonical config-dir list fermilab_dirs.txt, from the Fermilab cluster down to this machine.
#
#   REMOTE : nmatsum@lq.fnal.gov:<abs_path>           (abs paths come from fermilab_dirs.txt)
#   LOCAL  : $DEST/<basename>/                        (flat by basename; all 102 basenames are unique)
#
# LAT-ONLY by default: pulls ckpoint_lat.* only (what jj compute needs); skips ckpoint_rng.* and
# everything else.  Pass --with-rng to also pull ckpoint_rng.* (needed only to CONTINUE HMC locally).
#
# NON-DESTRUCTIVE: no --delete (never removes local files).  RESUMABLE: --partial + size/mtime skip,
# so re-running only fetches new / incomplete configs (also picks up configs that grew remotely since
# the last pull).  Safe to run repeatedly.
#
# Run it yourself (needs your ssh auth):
#   bash rsync_pull_configs_claude.sh                 # pull lat for ALL dirs in fermilab_dirs.txt
#   bash rsync_pull_configs_claude.sh heavy           # only paths matching 'heavy'  (jj-first priority)
#   bash rsync_pull_configs_claude.sh 'L2'            # only L2 ensembles, etc. (FILTER = grep -E on the path)
#   bash rsync_pull_configs_claude.sh -n              # dry-run (list what would transfer)
#   bash rsync_pull_configs_claude.sh --with-rng      # also pull ckpoint_rng.*
# Output is tee'd to rsync_pull_configs_claude.log .
set -u

REMOTE_HOST="nmatsum@lq.fnal.gov"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DIRS_FILE="${SCRIPT_DIR}/fermilab_dirs.txt"
DEST="${SCRIPT_DIR}"
LOG="${SCRIPT_DIR}/rsync_pull_configs_claude.log"

# ssh transport: keepalives so long transfers survive idle stretches.
SSH_OPTS="ssh -o ServerAliveInterval=30 -o ServerAliveCountMax=10"

DRY=""
WITH_RNG=0
FILTER=""
for arg in "$@"; do
  case "$arg" in
    -n) DRY="--dry-run" ;;
    --with-rng) WITH_RNG=1 ;;
    -*) echo "unknown option: $arg" ; exit 1 ;;
    *) FILTER="$arg" ;;
  esac
done

# rsync include/exclude: archive, human, partial+progress, compress; NO --delete.
RSYNC_OPTS=(-avhP --compress --include='ckpoint_lat.*')
if [ "$WITH_RNG" -eq 1 ]; then
  RSYNC_OPTS+=(--include='ckpoint_rng.*')
fi
RSYNC_OPTS+=(--include='*/' --exclude='*')

if [ ! -f "$DIRS_FILE" ]; then
  echo "ERROR: config-dir list not found: $DIRS_FILE"
  exit 1
fi

# Collect the target remote paths (skip blanks/comments; apply optional FILTER).
paths=()
while IFS= read -r p || [ -n "$p" ]; do
  p="${p#"${p%%[![:space:]]*}"}"
  p="${p%"${p##*[![:space:]]}"}"
  [ -z "$p" ] && continue
  case "$p" in \#*) continue ;; esac
  if [ -n "$FILTER" ]; then
    echo "$p" | grep -qE "$FILTER" || continue
  fi
  paths+=("$p")
done < "$DIRS_FILE"

echo "### rsync config pull START  [$(date +%F_%H:%M:%S)]  ${#paths[@]} ensembles  lat-only=$((1-WITH_RNG))  ->  $DEST ###" | tee "$LOG"

n_ok=0
n_fail=0
failed=()
for remote_path in "${paths[@]}"; do
  base="$(basename "$remote_path")"
  echo "### [$((n_ok+n_fail+1))/${#paths[@]}] $base  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
  # trailing slash on src copies the dir CONTENTS into DEST/base
  rsync "${RSYNC_OPTS[@]}" $DRY -e "$SSH_OPTS" \
    "${REMOTE_HOST}:${remote_path}/" "${DEST}/${base}/" 2>&1 | tee -a "$LOG"
  st=${PIPESTATUS[0]}
  if [ "$st" -eq 0 ]; then
    n_ok=$((n_ok+1))
    echo "###   OK  $base ###" | tee -a "$LOG"
  else
    n_fail=$((n_fail+1))
    failed+=("$base")
    echo "###   FAIL (rsync exit $st)  $base  -- continuing ###" | tee -a "$LOG"
  fi
done

echo "### rsync config pull DONE  [$(date +%F_%H:%M:%S)]  ok=$n_ok  fail=$n_fail ###" | tee -a "$LOG"
if [ "$n_fail" -gt 0 ]; then
  echo "### FAILED (re-run to retry; resumable): ###" | tee -a "$LOG"
  for b in "${failed[@]}"; do
    echo "###   $b ###" | tee -a "$LOG"
  done
fi
