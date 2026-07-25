#!/bin/bash -l
# tmp_rsync_fnal_L4_claude.sh  (_claude handoff, 2026-07-22)
# =============================================================================================
# Pull the FNAL L4 Nf4/6 massless gauge configs to SCC and stage them for RESUME on SCC.
# NM runs this on an SCC login node (scc1/scc2); it needs your FNAL ssh/kerberos auth (lq.fnal.gov).
# READ/COPY ONLY: rsync is ADDITIVE (no --delete), and there is NO rm/rmdir/overwrite anywhere.
#
# For each ensemble dir it pulls:
#   * ALL  ckpoint_lat.*        (~0.64 MB each -- the full chain; also usable for measurement later)
#   * the LATEST ckpoint_rng.<k_max>  (~0.5 GB -- REQUIRED to resume; the driver resumes from the highest
#                                      k that has BOTH ckpoint_lat.<k> AND ckpoint_rng.<k>)
# The latest rng k is queried over ssh per ensemble (FNAL thins rng, so rng k_max may trail lat k_max;
# the driver simply resumes from the rng k_max).
#
# After this completes, resume on SCC by BUILDING + submitting the Nf4/6 chains (see the tail echo).
# Precondition (per NM 2026-07-22): FNAL has STOPPED these ensembles (clean handoff, no double-writer).
#
# AUTH: FNAL access is kerberos -- run your usual `kinit` FIRST (e.g. `kinit <you>@FNAL.GOV`) so the ssh/rsync
# below inherit the ticket (GSSAPI). This script does NOT run kinit (it prompts for a password); it only
# WARNS via `klist -s` if no valid ticket is present.
# =============================================================================================
set -u

FNAL_HOST=${FNAL_HOST:-nmatsum@lq.fnal.gov}
FNAL_BASE=${FNAL_BASE:-/lustre2/affine/redo}
DEST=${DEST:-/projectnb/qfe/nmatsum/qed3/src/production}
# rsync flags: -a archive, -v verbose, -z COMPRESS (per NM). Deliberately NO -c (no per-file checksum;
# default size+mtime is enough and far cheaper) and NO --delete (additive pull only, never removes).
RSYNC_FLAGS=${RSYNC_FLAGS:-"-avz"}
RSYNC_DRYRUN=${RSYNC_DRYRUN:-0}          # 1 => add -n (rsync dry-run: list what WOULD transfer, move nothing)
LOG=rsync_fnal_L4_$(date +%y%m%d%H%M)_claude.log

DRYFLAG=""
if [ "$RSYNC_DRYRUN" -ne 0 ]
then
  DRYFLAG="-n"
fi

# kerberos preflight (warn only) -- you run kinit yourself before this script
if ! klist -s 2>/dev/null
then
  echo "WARNING: no valid kerberos ticket (klist -s failed). Run 'kinit <you>@FNAL.GOV' first, else ssh/rsync will fail."
fi

# ---- the 6 FNAL L4 Nf4/6 massless ensembles (blackboard rows p13/14/15) --------------------------
NFS="4 6"
GSQS="2.000000 4.000000 6.000000"
HBTAG="hb0.400000-1.000000"

ensdir () {   # nf gsq  -> the shared ensemble dir basename (identical on FNAL and SCC)
  echo "Nf${1}_gsq${2}at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_${HBTAG}"
}

{
  echo "############ FNAL->SCC L4 config pull  $(date) ############"
  echo "FNAL : ${FNAL_HOST}:${FNAL_BASE}"
  echo "DEST : ${DEST}"
  echo "dry-run (rsync -n) : ${RSYNC_DRYRUN}"
  echo

  for nf in $NFS
  do
    for gsq in $GSQS
    do
      d=$(ensdir "$nf" "$gsq")
      echo "==== $d ===="
      mkdir -p "${DEST}/${d}"

      # latest rng k on FNAL (empty if none)
      kmax_rng=$(ssh "$FNAL_HOST" "ls ${FNAL_BASE}/${d}/ckpoint_rng.* 2>/dev/null | sed 's/.*\.//' | sort -n | tail -1")
      if [ -z "$kmax_rng" ]
      then
        echo "  WARNING: no ckpoint_rng.* found on FNAL for $d -- cannot resume this one; skipping rng (lat still pulled)."
      else
        echo "  latest FNAL rng : ckpoint_rng.${kmax_rng}"
      fi

      # all gauge configs (cheap; full chain)
      echo "  + rsync ckpoint_lat.*"
      rsync $RSYNC_FLAGS $DRYFLAG "${FNAL_HOST}:${FNAL_BASE}/${d}/ckpoint_lat.*" "${DEST}/${d}/"

      # the single latest rng (required to resume)
      if [ -n "$kmax_rng" ]
      then
        echo "  + rsync ckpoint_rng.${kmax_rng}"
        rsync $RSYNC_FLAGS $DRYFLAG "${FNAL_HOST}:${FNAL_BASE}/${d}/ckpoint_rng.${kmax_rng}" "${DEST}/${d}/"
      fi
      echo
    done
  done

  echo "############ pull done ############"
  echo
  echo "NEXT (resume on SCC -- builds Nf4/6 {5,5,5} binaries, then submits chains to KMAX=400):"
  echo "  PE_OMP=4 ML_NF=\"4 6\" ML_GSQ=\"2.0 4.0 6.0\" ML_KMAX=400 WHICH=massless bash run_wrapper_L4_scc_claude.sh"
  echo "  (PE_OMP=4: each single-GPU job eats only 4 CPU slots, so the public GPU per-user slot cap (16) does"
  echo "   not bottleneck concurrency -- the per-queue GPU caps do. Nf2 production default stays PE_OMP=16.)"
} 2>&1 | tee "$LOG"
