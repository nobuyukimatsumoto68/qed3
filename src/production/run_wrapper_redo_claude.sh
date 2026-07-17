#!/bin/bash
# ---------------------------------------------------------------------------
# REDO production wrapper (corrected-gsq, bug-fixed). Submits the L=1,2 streams -- MASSLESS Nf{2,4,6}
# + MASSIVE Nf2 -- as 2-client MPS jobs on affine, TWO streams packed per A100 (run_redo_mps2_claude.sh).
# First try: --time 4h per job. Jobs self-chain via afterany per pair-slot. All streams write into ONE
# OUTDIR (/lustre2/affine/redo); each binary's dir3 encodes Nf/gsq/mRe/L -> no collision, auto-resume.
# Mimics the both_3d fermilab MPS convention (run_massless_mps2_claude.sh + its wrapper).
# Params: /project/qed3/qed3/src/production/params_L1L2_claude.md , params_massive_claude.md .
#
# Usage:  bash run_wrapper_redo_claude.sh [NCHAIN]              (default NCHAIN=1)
#         SLOTS="0 1 2" bash run_wrapper_redo_claude.sh 1        (submit only listed pair-slots)
#   NCHAIN = 4h jobs queued per pair-slot this invocation, each afterany-chained after the last.
# Binaries are built by ~/launch_redo_claude.sh into BINDIR (below) before this runs.
# ---------------------------------------------------------------------------
SLURM_USER=nmatsum
# This wrapper + run script + docs live in the git-synced production tree; RUN DATA (checkpoints, binaries,
# logs, slurm .out) stay in OUTDIR under /lustre2 (off git). RUNSCRIPT = the sbatch, submitted from production.
SCRIPTDIR=/project/qed3/qed3/src/production
OUTDIR=/lustre2/affine/redo
BINDIR=$OUTDIR
RUNSCRIPT=$SCRIPTDIR/run_redo_mps2_claude.sh
ACCT=${ACCT:-affine.lq2_gpu}
QOS=${QOS:-opp}             # affine allows only {opp,test} (default opp); smoke overrides to 'test'
WALL=${WALL:-04:00:00}      # smoke overrides to 00:30:00
WSEC=${WSEC:-14400}         # MUST match WALL (passed to the binary as WALL_SEC); smoke -> 1800
NCHAIN=${1:-1}

# binary resolvers (built by launch_redo_claude.sh)
mlbin() { echo "$BINDIR/hmc_fermilab_redo_massless_L$1_Nf$2_claude.o"; }   # $1=L $2=Nf (massless, mRe=0)
mvbin() { echo "$BINDIR/hmc_fermilab_redo_massive_L$1_claude.o"; }         # $1=L      (massive Nf2)

# ---- stream list: "TYPE L Nf gsq mRe" -- 18 massless (Nf{2,4,6}) then 8 massive (Nf2) -------------
# Ordered so consecutive pairing keeps massless-with-massless and massive-with-massive (never mixed:
# their trajectory targets differ ~20x). 26 streams -> 13 jobs (2 clients each).
# ---- pair table: slot index -> two clients "TYPE L Nf gsq mRe" (CA=client1, CB=client2). EXPLICIT slot
# indices (not consecutive pairing of a flat list) so a COMPLETED slot can be commented out WITHOUT
# renumbering the survivors. Renumbering is unsafe while chains are live: a re-apply under a shifted index
# would anchor-miss and launch a chain-root that runs a stream ALREADY owned by a still-queued job (same
# output dir3) -> two procs writing one checkpoint dir -> corruption. Indices here MUST match the originally
# submitted hmc_redo_p<idx> job names. 26 streams / 13 slots.
declare -A CA CB
# massless L1 gsq{0.5,1.0,1.5} x Nf{2,4,6}
CA[0]="ml 1 2 0.5 0.0";  CB[0]="ml 1 2 1.0 0.0"
CA[1]="ml 1 2 1.5 0.0";  CB[1]="ml 1 4 0.5 0.0"
CA[2]="ml 1 4 1.0 0.0";  CB[2]="ml 1 4 1.5 0.0"
CA[3]="ml 1 6 0.5 0.0";  CB[3]="ml 1 6 1.0 0.0"
CA[4]="ml 1 6 1.5 0.0";  CB[4]="ml 2 2 1.0 0.0"
# massless L2 gsq{1.0,2.0,3.0} x Nf{2,4,6}
CA[5]="ml 2 2 2.0 0.0";  CB[5]="ml 2 2 3.0 0.0"
CA[6]="ml 2 4 1.0 0.0";  CB[6]="ml 2 4 2.0 0.0"
CA[7]="ml 2 4 3.0 0.0";  CB[7]="ml 2 6 1.0 0.0"
CA[8]="ml 2 6 2.0 0.0";  CB[8]="ml 2 6 3.0 0.0"
# massive L1 gsq1.5 Nf2 m{0.1..0.4} -- slots 9,10 COMPLETED 2026-07-16 (all four at 119/120);
# commented out (NOT deleted) so re-apply skips them and indices below stay put:
# CA[9]="mv 1 2 1.5 0.1";  CB[9]="mv 1 2 1.5 0.2"
# CA[10]="mv 1 2 1.5 0.3"; CB[10]="mv 1 2 1.5 0.4"
# massive L2 gsq3.0 Nf2 m{0.1..0.4} -- slots 11,12 COMPLETED 2026-07-16 (all four at 79/80);
# commented out (NOT deleted) so re-apply skips them. ALL massive done -> re-apply now runs the 18 massless only.
# CA[11]="mv 2 2 3.0 0.1"; CB[11]="mv 2 2 3.0 0.2"
# CA[12]="mv 2 2 3.0 0.3"; CB[12]="mv 2 2 3.0 0.4"

binof() {  # $1=TYPE $2=L $3=Nf  -> binary path
  if [ "$1" = ml ]; then mlbin "$2" "$3"; else mvbin "$2"; fi
}

declare -A PREV=()

submit_slot() {  # $1=slot-index (uses CA[idx], CB[idx]); skips undefined (completed/removed) slots
    local idx=$1
    if [ -z "${CA[$idx]:-}" ]; then echo "slot $idx: not defined (completed/removed) -- skip"; return 0; fi
    read -r t1 l1 n1 g1 m1 <<< "${CA[$idx]}"
    read -r t2 l2 n2 g2 m2 <<< "${CB[$idx]}"
    local b1 b2; b1=$(binof "$t1" "$l1" "$n1"); b2=$(binof "$t2" "$l2" "$n2")
    local jobname="hmc_redo_p${idx}"
    local anchor="${PREV[$idx]:-}"
    if [ -z "$anchor" ]; then
        anchor=$(squeue -u "$SLURM_USER" --name="$jobname" -h -o "%i" | sort -n | tail -1)
    fi
    local dep=""; [ -n "$anchor" ] && dep="--dependency=afterany:${anchor}"
    local jid
    jid=$(sbatch --parsable --job-name="$jobname" \
        --account="$ACCT" --qos="$QOS" --time="$WALL" ${dep} \
        --output="$OUTDIR/slurm_%x_%j.out" \
        --export=ALL,BIN1="$b1",GSQ1="$g1",NF1="$n1",MRE1="$m1",BIN2="$b2",GSQ2="$g2",NF2="$n2",MRE2="$m2",OUTDIR="$OUTDIR",WALL_SEC="$WSEC" \
        "$RUNSCRIPT")
    if [ -z "$jid" ]; then echo "ERROR: sbatch failed for $jobname" >&2; return 1; fi
    PREV[$idx]="$jid"
    local after=""; [ -n "$anchor" ] && after=" after $anchor"
    printf 'slot %-2s jid=%-10s [%s L%s Nf%s gsq%s mRe%s] + [%s L%s Nf%s gsq%s mRe%s]%s\n' \
        "$idx" "$jid" "$t1" "$l1" "$n1" "$g1" "$m1" "$t2" "$l2" "$n2" "$g2" "$m2" "$after"
}

DEFAULT_SLOTS=$(printf '%s\n' "${!CA[@]}" | sort -n)   # only DEFINED (non-commented) slots, ascending
for idx in ${SLOTS:-$DEFAULT_SLOTS}; do
    for ((c=0; c<NCHAIN; c++)); do
        submit_slot "$idx"
    done
done
echo "# submitted slots [${SLOTS:-$(echo $DEFAULT_SLOTS | tr '\n' ' ')}] x NCHAIN=$NCHAIN (qos=$QOS, --time=$WALL) into $OUTDIR"
