#!/usr/bin/env bash
# =============================================================================
# _claude -- submit + afterany-chain the FNAL L=1 stride-2 conn job on qed3.
# Chunk 2 of conn_s2_fnal_L1_impl_plan_claude.md. Mirrors run_wrapper_redo_claude.sh.
#   bash run_wrapper_conn_s2_fnal_claude.sh [NCHAIN]      # default NCHAIN=8
# env overrides: ACCT QOS WALL NCHAIN. Resumable via the driver's per-config
# "complete" gate, so each chained job just resumes where the last stopped.
#
# SMOKE (Chunk 3): one ensemble, few configs, qos=test 30min --
#   QOS=test WALL=00:30:00 NCHAIN=1 CONN_ENS_FILTER='Nf2_gsq0.500000at0.200000' \
#     CONN_KMAX_CAP=42 bash run_wrapper_conn_s2_fnal_claude.sh 1
#   (CONN_* are exported to the job via --export below.)
# =============================================================================
set -u

SCRIPTDIR=/project/qed3/qed3/src/production
RUN=$SCRIPTDIR/run_conn_s2_fnal_claude.sh
WORKDIR=/lustre2/affine/redo/conn_s2
ACCT=${ACCT:-qed3.lq2_gpu}
QOS=${QOS:-normal}
WALL=${WALL:-04:00:00}
NCHAIN=${1:-${NCHAIN:-8}}
SLURM_USER=${USER:-nmatsum}

# TWO 1-GPU chains, 2 offsets each (the run script requests --gpus=a100:1 and packs its offsets 2/GPU
# under MPS). Scheduler-friendly (backfills like the HMC 1-GPU jobs) vs one hard-to-place 2-GPU job.
# Each element = "<suffix>:<space-separated offsets>". For a smoke, both (tiny) jobs run.
# NB: do NOT name this `GROUPS` -- that is a bash special variable (the user's unix group IDs).
CHAINS=("a:2 4" "b:6 8")

mkdir -p "$WORKDIR"
CONN_L=${CONN_L:-1}                              # lattice L (1/2/3/4); binary jj_conn_s2_L${CONN_L}_fnal_claude.o must exist
echo "# conn_s2 L${CONN_L} -> qed3 (acct=$ACCT qos=$QOS time=$WALL, 1 GPU/job) x NCHAIN=$NCHAIN per group"
for grp in "${CHAINS[@]}"; do
  sfx=${grp%%:*}; offs=${grp#*:}                # sfx=a ; offs="2 4"
  jobname="conn_s2_L${CONN_L}${sfx}"
  export_env="ALL,CONN_L=${CONN_L},CONN_OFFSETS=${offs// /_}"   # underscore-join: SLURM --export splits on commas
  [ -n "${CONN_ENS_FILTER:-}" ] && export_env="$export_env,CONN_ENS_FILTER=$CONN_ENS_FILTER"
  [ -n "${CONN_KMAX_CAP:-}" ]   && export_env="$export_env,CONN_KMAX_CAP=$CONN_KMAX_CAP"
  anchor=$(squeue -u "$SLURM_USER" -h --name="$jobname" -o "%i" | sort -n | tail -1)
  echo "## group $sfx (offsets $offs) jobname=$jobname anchor=${anchor:-none}"
  for (( c=0; c<NCHAIN; c++ )); do
    dep=""; [ -n "$anchor" ] && dep="--dependency=afterany:${anchor}"
    jid=$(sbatch --parsable --job-name="$jobname" \
        --account="$ACCT" --qos="$QOS" --time="$WALL" ${dep} \
        --export="$export_env" \
        --output="$WORKDIR/slurm_%x_%j.out" \
        "$RUN")
    [ -n "$jid" ] || { echo "ERROR: sbatch failed ($jobname)"; exit 1; }
    echo "  submitted $jobname jid=$jid${anchor:+  after $anchor}"
    anchor=$jid
  done
done
echo "# check: $WORKDIR/conn_s2_L${CONN_L}_off*_jid*_claude.log  and  data_*/corr_ylm_conn_t00_nhits1_s1/"
