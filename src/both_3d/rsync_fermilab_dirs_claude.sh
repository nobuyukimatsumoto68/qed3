#!/usr/bin/env bash
# rsync_fermilab_dirs_claude.sh -- pull the ensemble directories listed in
# fermilab_dirs.txt from the fnal host (nmatsum@lq.fnal.gov) down to this machine
# for local analysis.
#
# fermilab_dirs.txt holds one ABSOLUTE remote path per line, e.g.
#   /lustre2/qed3/Nf2_gsq8.000000at0.200000nu01.000000mRe0.000000mIm0.100000nt128L1
#   /lustre2/affine/Nf2_gsq8.000000at0.200000nu01.000000mRe0.000000mIm0.010000nt128L1
# Blank lines and lines beginning with '#' are skipped.
#
# Each remote directory is pulled to a LOCAL directory of the same basename under
# DEST (the local both_3d dir by default). The basenames in fermilab_dirs.txt are
# unique across /lustre2/qed3 and /lustre2/affine, so a flat layout has no clashes
# and matches what the analysis notebooks expect for --ens-dir.
#
# Transfers are resumable (--partial), so re-running picks up where it left off and
# skips already-complete files. No --delete is used; nothing local is removed.
#
# Usage:  bash rsync_fermilab_dirs_claude.sh           # pull everything
#         bash rsync_fermilab_dirs_claude.sh -n        # dry-run (list only)
set -u

REMOTE_HOST="nmatsum@lq.fnal.gov"

# Directory holding this script (and fermilab_dirs.txt).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DIRS_FILE="${SCRIPT_DIR}/fermilab_dirs.txt"

# Local destination root for the pulled ensemble directories.
DEST="${SCRIPT_DIR}"

DRY=""
for arg in "$@"; do
  case "$arg" in
    -n) DRY="--dry-run" ;;
    *) echo "unknown argument: $arg" ; exit 1 ;;
  esac
done

if [ ! -f "$DIRS_FILE" ]; then
  echo "ERROR: directory list not found: $DIRS_FILE"
  exit 1
fi

while IFS= read -r remote_path || [ -n "$remote_path" ]; do
  # strip leading/trailing whitespace
  remote_path="${remote_path#"${remote_path%%[![:space:]]*}"}"
  remote_path="${remote_path%"${remote_path##*[![:space:]]}"}"
  # skip blank lines and comments
  if [ -z "$remote_path" ]; then
    continue
  fi
  case "$remote_path" in
    \#*) continue ;;
  esac

  base="$(basename "$remote_path")"
  echo "===== pull ${remote_path} -> ${DEST}/${base}/ ====="
  # trailing slash on the source copies the directory CONTENTS into DEST/base
  rsync -avzuh --progress --partial $DRY \
    "${REMOTE_HOST}:${remote_path}/" \
    "${DEST}/${base}/"
done < "$DIRS_FILE"

echo "DONE -> ${DEST}/<Nf*_...nt128L1>/"
