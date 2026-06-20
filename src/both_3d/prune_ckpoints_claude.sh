#!/bin/bash
# Dry-run lister for stale HMC checkpoints in a given trajectory range.
#
# NOTE: per the user's global rule, this script never executes rm / any delete.
# It only LISTS the matching ckpoint_lat.k / ckpoint_rng.k files (dry run) and
# then prints the exact rm one-liner for YOU to run yourself.
#
# Usage:
#   bash prune_ckpoints_claude.sh DIR KMIN KMAX
#
# Example (drop the post-freeze L2/Nf2 checkpoints k=230..325, keeping k<=229):
#   bash prune_ckpoints_claude.sh Nf2_gsq8.000000at0.200000nu01.000000nt128L2 230 325
set -u

if [ "$#" -ne 3 ]
then
  echo "usage: bash $0 DIR KMIN KMAX"
  exit 1
fi

DIR=$1
KMIN=$2
KMAX=$3

if [ ! -d "$DIR" ]
then
  echo "ERROR: directory not found: $DIR"
  exit 1
fi

echo "=== DRY RUN: checkpoints in $DIR with k in [$KMIN, $KMAX] ==="
count=0
total_kb=0
for k in $(seq "$KMIN" "$KMAX")
do
  for kind in lat rng
  do
    fn="$DIR/ckpoint_${kind}.$k"
    if [ -f "$fn" ]
    then
      sz=$(du -k "$fn" | cut -f1)
      count=$((count + 1))
      total_kb=$((total_kb + sz))
      printf "  %s  (%s KB)\n" "$fn" "$sz"
    fi
  done
done

echo "--- summary ---"
echo "files matched : $count"
echo "total size    : $total_kb KB"
echo
echo "=== to actually delete, run this yourself (NOT run by this script): ==="
echo "for k in \$(seq $KMIN $KMAX); do rm -f \"$DIR/ckpoint_lat.\$k\" \"$DIR/ckpoint_rng.\$k\"; done"
