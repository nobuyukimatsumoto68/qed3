#!/usr/bin/env bash
# rsync_scc_claude.sh -- pull from BU SCC to local, ONE password prompt for all sections:
#   [1] L4 REDO gauge configs      (ckpoint_lat.* ONLY)
#   [2] L4 conn stride-2 h5        (corr_ylm_conn_t00_nhits1_s1/corr.*.h5 ; --ignore-existing)
#   [3] L4 distill perambulators   (data_Nf2_gsq2.0 L4 / distill_Nv48_sym/peram.<k>.h5 ; added 2026-09-23)
# Source : nmatsum@scc1.bu.edu:/projectnb/qfe/nmatsum/qed3/src/production/
# Dest   : /mnt/barracuda22/qed3/qed3/src/production/
# Run    : bash rsync_scc_claude.sh           (repeatable; rsync sends only new files -- re-run as SCC fills)
# Log    : rsync_scc_claude.log (appended, run delimited by date headers)
#
# ONE PASSWORD: the script opens a single ssh ControlMaster connection first (that is the only prompt),
# then every rsync section rides on it via -e "ssh -o ControlPath=...".  The master lingers for
# ControlPersist=2h and expires by itself (nothing is killed here); a re-run within 2h asks nothing.
# -z: WAN compression. No -c (mtime+size is enough, all pulled files are write-once).
# NEVER pull ckpoint_rng.* / logs / binaries. No rm anywhere in this script.
set -uo pipefail

DEST=/mnt/barracuda22/qed3/qed3/src/production
HOST=nmatsum@scc1.bu.edu
SRC="$HOST":/projectnb/qfe/nmatsum/qed3/src/production/
LOG="$DEST"/rsync_scc_claude.log
CP="$HOME/.ssh/cm-scc-rsync"
SSH="ssh -o ControlPath=$CP"

{
  echo "##################################################################"
  echo "### rsync SCC -> local  $(date)"
  echo "##################################################################"
} | tee -a "$LOG"

# ---- [0] single ssh master connection = the ONE password prompt ----
if ssh -o ControlPath="$CP" -O check "$HOST" >/dev/null 2>&1
then
  echo "### ssh master already open (reusing, no prompt) ###" | tee -a "$LOG"
else
  echo "### opening ssh master connection to $HOST (enter password ONCE) ###" | tee -a "$LOG"
  ssh -o ControlMaster=yes -o ControlPath="$CP" -o ControlPersist=2h -fN "$HOST"
  rc=$?
  if [ "$rc" -ne 0 ]
  then
    echo "### ERROR: ssh master failed (code $rc) -- aborting ###" | tee -a "$LOG"
    exit "$rc"
  fi
fi

# ---- [1] L4 gauge configs ----
echo "### [1] L4 ckpoint_lat  $(date)" | tee -a "$LOG"
rsync -avz -e "$SSH" --prune-empty-dirs \
  --include='Nf*_*L4_hb*/' --include='ckpoint_lat.*' --exclude='*' \
  "$SRC" "$DEST"/ 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]
then
  echo "### ERROR: config rsync exited with code $rc ###" | tee -a "$LOG"
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
} | tee -a "$LOG"

# ---------------------------------------------------------------------------
# [2] CONN stride-2 datasets (added 2026-08-13).  SCC owns ALL 9 L4 conn (classes k=3,5,7,9 mod 10 plus,
# since 2026-09-23, the k=1 class completion with --nhits 2 -> corr.<k>.h1.h5 in the SAME dir; the
# corr.*.h5 glob grabs both hits).
# --ignore-existing: never overwrite a conn h5 that already exists locally (a collision would mean the
# residue-class partition was violated; the census below surfaces it instead of clobbering data).
# NOTE: cross-site results agree only to the CG tolerance (1e-8), NOT bit-for-bit -- validate numerically.
# ---------------------------------------------------------------------------
echo "### [2] conn stride-2 h5  $(date)" | tee -a "$LOG"
rsync -avz -e "$SSH" --prune-empty-dirs --ignore-existing \
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
} | tee -a "$LOG"

# ---------------------------------------------------------------------------
# [3] L4 DISTILL PERAMBULATORS (added 2026-09-23; blackboard "UPDATE 2026-09-23 (SCC -- L4 DISTILL
# PERAMBULATORS, Nv=48)").  L4 Nf2 g2.0 only, BASIS_SYM, Nv=48, odd k -> ~400 x 183 MB = ~73 GB.
# Lands in the plain data_<ens>/ dir (no _vmRe suffix), same layout as the local L1/L2 distill_Nv24_sym.
# No --ignore-existing here: a re-run refreshes any peram whose size/mtime changed on SCC (rsync writes
# to a temp name and renames, so an interrupted transfer never leaves a truncated peram.<k>.h5).
# Consumers: NVDIR=distill_Nv48_sym LREF=4 ; tau (NOT tau_gw) ; MODE_CONTACT=1.
# ---------------------------------------------------------------------------
PERAM_ENS=data_Nf2_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000
echo "### [3] L4 distill perams  $(date)" | tee -a "$LOG"
rsync -avz -e "$SSH" --prune-empty-dirs \
  --include="$PERAM_ENS/" \
  --include="$PERAM_ENS/distill_Nv48_sym/" \
  --include="$PERAM_ENS/distill_Nv48_sym/peram.*.h5" \
  --exclude='*' \
  "$SRC" "$DEST"/ 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]
then
  echo "### ERROR: peram rsync exited with code $rc ###" | tee -a "$LOG"
  exit "$rc"
fi
{
  pd="$DEST/$PERAM_ENS/distill_Nv48_sym"
  np=$(find "$pd" -maxdepth 1 -name 'peram.*.h5' 2>/dev/null | wc -l)
  sz=$(du -sh "$pd" 2>/dev/null | cut -f1)
  echo "### L4 Nv48 perams local: $np / ~400 (odd k 1..799), $sz ###"
  echo "### done  $(date) ###"
} | tee -a "$LOG"
