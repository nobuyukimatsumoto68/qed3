#!/usr/bin/env bash
# rsync_pull_L4_claude.sh -- pull the L=4 results from the A100 host (fnal) back to this machine.
#   - the small corr.*.h5 outputs (method B and, if run, method D);
#   - the 55 GB propagator prop_deter_L4/Dinv.0.h5 (P = D_ov^{-1} and Dov), so the now-cheap
#     INSERTION-DIAGONAL local/cont contraction can be rerun HERE instead of on the A100.
#
# Pulled (under data_free_vmRe0.000000vmIm0.000000/):
#   corr_deter_local_L4/corr.0.h5       (B = local op + lattice P)
#   corr_deter_local_cont_L4/corr.0.h5  (D = local op + continuum G, if it was run)
#   prop_deter_L4/Dinv.0.h5             (~55 GB; resumable via --partial/--append-verify)
#
# Usage:  bash rsync_pull_L4_claude.sh            # corr + propagator
#         bash rsync_pull_L4_claude.sh -n         # dry-run (list what would transfer)
#         bash rsync_pull_L4_claude.sh --no-prop  # corr only (skip the 55 GB propagator)
set -u

REMOTE_HOST="nmatsum@lq.fnal.gov"
# Remote working dir = src/both_3d/ on the A100 host. EDIT if it differs there.
REMOTE_DIR="~/project_qed3/qed3/src/both_3d"
REMOTE_DATA="data_free_vmRe0.000000vmIm0.000000"

LOCAL_DIR="/mnt/barracuda22/qed3/qed3/src/both_3d"

DRY=""
GET_PROP=1
for arg in "$@"; do
  case "$arg" in
    -n) DRY="--dry-run" ;;
    --no-prop) GET_PROP=0 ;;
  esac
done

# --- small correlator outputs (method B / D) ---
# Per-directory pull: include the dir + its corr.*.h5, exclude everything else.
# A missing remote dir (e.g. method D / cont was never run) just errors out that one pull -- harmless.
for sub in corr_deter_local_L4 corr_deter_local_cont_L4; do
  echo "===== pull ${REMOTE_DATA}/${sub} ====="
  rsync -avh --progress $DRY \
    --include="corr.*.h5" \
    --exclude="*" \
    "${REMOTE_HOST}:${REMOTE_DIR}/${REMOTE_DATA}/${sub}/" \
    "${LOCAL_DIR}/${REMOTE_DATA}/${sub}/"
done

# --- 55 GB propagator (resumable) ---
if [ "$GET_PROP" -eq 1 ]; then
  echo "===== pull ${REMOTE_DATA}/prop_deter_L4/Dinv.0.h5 (~55 GB, resumable) ====="
  mkdir -p "${LOCAL_DIR}/${REMOTE_DATA}/prop_deter_L4"
  rsync -avh --progress --partial --append-verify $DRY \
    "${REMOTE_HOST}:${REMOTE_DIR}/${REMOTE_DATA}/prop_deter_L4/Dinv.0.h5" \
    "${LOCAL_DIR}/${REMOTE_DATA}/prop_deter_L4/"
fi

echo "DONE -> ${LOCAL_DIR}/${REMOTE_DATA}/{corr_deter_local_L4,corr_deter_local_cont_L4,prop_deter_L4}/"
