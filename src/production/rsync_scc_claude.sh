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

# ---------------------------------------------------------------------------
# CONN stride-2 datasets (added 2026-08-13).  SCC runs the L=4 slice of the
# three-site stride-2 completion (scc_conn_s2_L4_impl_plan_claude.md): the four
# residue classes k = 3,5,7,9 (mod 10), DISJOINT from the local stride-10 grid
# (k = 1 mod 10), so the pull is purely ADDITIVE.
#
# Same SRC as the configs: the SCC jobs cd to /projectnb/qfe/nmatsum/qed3/src/
# production (SRCDIR in run_conn_s2_scc_claude.sh:36), which is already the
# transfer root above -- so only the include rules differ.  Kept as its own
# invocation anyway so --ignore-existing applies to the conn h5 ONLY and never
# to the configs.
#
# --ignore-existing: never overwrite a conn h5 that already exists locally.  A
# collision would mean the residue-class partition was violated; the census
# below surfaces it instead of silently clobbering production data.
#
# NOTE: cross-site results agree only to the CG tolerance (1e-8), NOT bit-for-
# bit, since the arch differs -- validate numerically, never by md5.
# ---------------------------------------------------------------------------
{
  echo "### rsync SCC conn stride-2 -> local  $(date)"
} | tee -a "$LOG"

rsync -avz --prune-empty-dirs --ignore-existing \
  --include='data_Nf*/' \
  --include='data_Nf*/corr_ylm_conn_t00_nhits1_s1/' \
  --include='data_Nf*/corr_ylm_conn_t00_nhits1_s1/corr.*.h5' \
  --exclude='*' \
  "$SRC" "$DEST"/ 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]
then
  echo "### ERROR: conn rsync exited with code $rc ###" | tee -a "$LOG"
  exit "$rc"
fi

{
  echo "### conn h5 census by residue class k mod 10 (L4 rows are SCC's slice) ###"
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
