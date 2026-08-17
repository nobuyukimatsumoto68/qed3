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
ACCT=${ACCT:-qed3.lq2_gpu}  # DEFAULT flipped affine->qed3 2026-08-03 (NM: affine allocation gone -> all L3 on qed3).
                            # qed3-account jobs still write to /lustre2/affine/redo (OUTDIR is on affine Lustre; only
                            # COMPUTE is charged to qed3). Override ACCT=affine.lq2_gpu only if affine ever comes back.
QOS=${QOS:-normal}          # qed3.lq2_gpu allows {normal,opp,test}. Default normal (allocation-backed, NOT preempted).
                            # Override QOS=opp for bulk when qed3 lquota avail nears 0.
WALL=${WALL:-08:00:00}      # opp MaxWall = 8h (used for L1/L2 and L4 alike, per NM 2026-07-17); smoke -> 00:30:00
WSEC=${WSEC:-28800}         # MUST match WALL (passed to the binary as WALL_SEC); smoke -> 1800
NCHAIN=${1:-1}

# binary resolvers (built by launch_redo_claude.sh)
mlbin() { echo "$BINDIR/hmc_fermilab_redo_massless_L$1_Nf$2_claude.o"; }        # $1=L $2=Nf (massless at=0.2)
mlbin01() { echo "$BINDIR/hmc_fermilab_redo_massless_L$1_Nf$2_at0p1_claude.o"; } # $1=L $2=Nf (massless, HALF a_t=0.1)
mvbin() { echo "$BINDIR/hmc_fermilab_redo_massive_L$1_claude.o"; }              # $1=L      (massive Nf2)

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
# massless L1 gsq{0.5,1.0,1.5} x Nf{2,4,6} -- REOPENED 2026-07-20 for the L1/L2 stats bump 2000->4000
# (KMAX 4000, KRNG 20). All 9 L1 were done@1999; they resume from disk and run on to 3999.
# slot 0 CLOSED 2026-07-23: both streams done@3999 (L1 Nf2 g0.5 + g1.0, KMAX 4000). Commented out.
# CA[0]="ml 1 2 0.5 0.0";  CB[0]="ml 1 2 1.0 0.0"
# slots 1,2,3 CLOSED 2026-07-26: all L1 at=0.2 done@3999 (KMAX 4000). Commented out (all 9 L1 now complete).
# CA[1]="ml 1 2 1.5 0.0";  CB[1]="ml 1 4 0.5 0.0"
# CA[2]="ml 1 4 1.0 0.0";  CB[2]="ml 1 4 1.5 0.0"
# CA[3]="ml 1 6 0.5 0.0";  CB[3]="ml 1 6 1.0 0.0"
# slots 4,5,6 CLOSED 2026-07-28: all done@3999 (KMAX 4000). p4=L2 Nf2 g1.0(+L1 Nf6 g1.5); p5=L2 Nf2 g2/g3;
# p6=L2 Nf4 g1/g2. Commented out.
# CA[4]="ml 1 6 1.5 0.0";  CB[4]="ml 2 2 1.0 0.0"
# massless L2 gsq{1.0,2.0,3.0} x Nf{2,4,6} -- slots 5-8, bumped to KMAX 4000 (2026-07-20; were 2000)
# CA[5]="ml 2 2 2.0 0.0";  CB[5]="ml 2 2 3.0 0.0"
# CA[6]="ml 2 4 1.0 0.0";  CB[6]="ml 2 4 2.0 0.0"
# slot 7 CLOSED 2026-07-28: both done@3999 (L2 Nf4 g3.0 + Nf6 g1.0, KMAX 4000). Commented out.
# CA[7]="ml 2 4 3.0 0.0";  CB[7]="ml 2 6 1.0 0.0"
# slot 8 CLOSED 2026-07-30: both done@3999 (L2 Nf6 g2.0 + g3.0, KMAX 4000). Last at=0.2 slot -> at=0.2 campaign COMPLETE. Commented out.
# CA[8]="ml 2 6 2.0 0.0";  CB[8]="ml 2 6 3.0 0.0"
# massive L1 gsq1.5 Nf2 m{0.1..0.4} -- slots 9,10 COMPLETED 2026-07-16 (all four at 119/120);
# commented out (NOT deleted) so re-apply skips them and indices below stay put:
# CA[9]="mv 1 2 1.5 0.1";  CB[9]="mv 1 2 1.5 0.2"
# CA[10]="mv 1 2 1.5 0.3"; CB[10]="mv 1 2 1.5 0.4"
# massive L2 gsq3.0 Nf2 m{0.1..0.4} -- slots 11,12 COMPLETED 2026-07-16 (all four at 79/80);
# commented out (NOT deleted) so re-apply skips them. ALL massive done -> re-apply now runs the 18 massless only.
# CA[11]="mv 2 2 3.0 0.1"; CB[11]="mv 2 2 3.0 0.2"
# CA[12]="mv 2 2 3.0 0.3"; CB[12]="mv 2 2 3.0 0.4"
# L4 massless Nf{4,6} x gsq{2,4,6} -- ADDED 2026-07-17 (NM). Binaries built with steps {5,5,5} (-DL4_MDSTEP=5),
# -DMIXED_FORCE, KMAX 400, KRNG 4 (launch_redo_claude.sh BUILDS). L4 is EXPENSIVE but runs at the same 8h
# default wall as L2 now; a plain re-apply covers them. To run ONLY L4: SLOTS="13 14 15" bash run_wrapper_redo_claude.sh <N>
# (pairs kept ~like-cost: Nf4 g2+g4 / Nf6 g2+g4 / the two g6 together).
# L4 DROPPED 2026-07-20 (NM: pursue L1/L2 stats to 4000 instead of L4). Commented out so re-apply skips
# them; the L4 data (k~30/400) stays on disk. Uncomment + restore the BUILDS lines to revive.
# CA[13]="ml 4 4 2.0 0.0"; CB[13]="ml 4 4 4.0 0.0"
# CA[14]="ml 4 6 2.0 0.0"; CB[14]="ml 4 6 4.0 0.0"
# CA[15]="ml 4 4 6.0 0.0"; CB[15]="ml 4 6 6.0 0.0"
# HALF-a_t (at=0.1, fixed Nt=128) ensembles -- ADDED 2026-07-22 (NM). type "ml01" -> _at0p1 binaries.
# Massless, middle gsq per L (L1 g1.0, L2 g2.0), all Nf2/4/6. KMAX 2000, KRNG 20, SAME HMC params + frozen
# window as at=0.2. MPS pairs kept ~like-cost (cheap L1s together / mid / heavy L2s together). dir3 = at0.100000.
# slot 16 CLOSED 2026-07-26: both at=0.1 L1 done@1999 (Nf2 g1.0 + Nf4 g1.0, KMAX 2000). Commented out.
# CA[16]="ml01 1 2 1.0 0.0"; CB[16]="ml01 1 4 1.0 0.0"   # at0.1: L1 Nf2 g1.0 || L1 Nf4 g1.0
# slot 17 CLOSED 2026-07-30: both done@1999 (at=0.1 L1 Nf6 g1.0 + L2 Nf2 g2.0, KMAX 2000). Commented out.
# CA[17]="ml01 1 6 1.0 0.0"; CB[17]="ml01 2 2 2.0 0.0"   # at0.1: L1 Nf6 g1.0 || L2 Nf2 g2.0
# slot 18 CLOSED 2026-07-30: both done@1999 (at=0.1 L2 Nf4 g2.0 + Nf6 g2.0, KMAX 2000). Last half-a_t slot -> half-a_t COMPLETE. Commented out.
# CA[18]="ml01 2 4 2.0 0.0"; CB[18]="ml01 2 6 2.0 0.0"   # at0.1: L2 Nf4 g2.0 || L2 Nf6 g2.0
# L=3 (at=0.2) ensembles -- ADDED 2026-07-28 (NM). type "ml" L3 -> hmc_fermilab_redo_massless_L3_Nf<Nf>_claude.o
# (3-stage {0,0.4,1.0}/{3,3,3}, MG100, frozen window (0.015,8.0), KMAX 1000 [was 800, bumped 2026-08-14], KRNG 20, -DNO_METROP_UNTIL=2).
# gsq{1.5,3.0,4.5} x Nf{2,4,6} = 9 streams -> 4 like-cost MPS pairs + 1 SOLO (Nf6 g4.5 = most expensive, gets a
# dedicated GPU; CB unset -> run_redo_mps2 runs client 1 only). gsq3.0/Nf2 (CB[19]) already thermalizing (k~4).
# slot 19 REOPENED 2026-08-14: target bumped 800 -> 1000 (KMAX 1000 rebuild). Nf2 g1.5+g3.0 resume 799 -> 999.
CA[19]="ml 3 2 1.5 0.0"; CB[19]="ml 3 2 3.0 0.0"   # L3: Nf2 g1.5 || Nf2 g3.0
# REPACK 2026-08-07 (all jobs pending -> zero-risk moment): 7 live streams -> 3 pairs + 1 solo. Nf2 g4.5 (704,
# finishes soonest) SINGLED OUT as the solo (shortest solo-GPU waste); its old partner Nf4 g1.5 re-paired with
# the old p23 solo Nf6 g4.5. p21/p22 untouched (already balanced pairs). See repack note in the memory snapshot.
# slot 20 REOPENED 2026-08-14: target bumped 800 -> 1000 (KMAX 1000 rebuild). Nf2 g4.5 SOLO resumes 799 -> 999.
CA[20]="ml 3 2 4.5 0.0"                             # L3: Nf2 g4.5 (SOLO -- was paired w/ Nf4 g1.5)
CA[21]="ml 3 4 3.0 0.0"; CB[21]="ml 3 4 4.5 0.0"   # L3: Nf4 g3.0 || Nf4 g4.5
CA[22]="ml 3 6 1.5 0.0"; CB[22]="ml 3 6 3.0 0.0"   # L3: Nf6 g1.5 || Nf6 g3.0
CA[23]="ml 3 4 1.5 0.0"; CB[23]="ml 3 6 4.5 0.0"   # L3: Nf4 g1.5 || Nf6 g4.5 (REPACKED pair -- was Nf6 g4.5 solo)

binof() {  # $1=TYPE $2=L $3=Nf  -> binary path.  ml=massless at0.2 ; ml01=massless HALF a_t=0.1 ; mv=massive
  if [ "$1" = ml ]; then mlbin "$2" "$3"; elif [ "$1" = ml01 ]; then mlbin01 "$2" "$3"; else mvbin "$2"; fi
}

declare -A PREV=()

# g10: encode a gsq or mass as an integer = value x 10 (0.5->5, 1.5->15, 3.0->30). %.0f rounds
# (avoids double-precision truncation, e.g. 0.3*10=2.9999 -> "3"). Used only for the pretty job name.
g10() { awk "BEGIN{printf \"%.0f\", $1*10}"; }

submit_slot() {  # $1=slot-index (uses CA[idx], CB[idx]); skips undefined (completed/removed) slots
    local idx=$1
    if [ -z "${CA[$idx]:-}" ]; then echo "slot $idx: not defined (completed/removed) -- skip"; return 0; fi
    read -r t1 l1 n1 g1 m1 <<< "${CA[$idx]}"
    read -r t2 l2 n2 g2 m2 <<< "${CB[$idx]}"
    local b1 b2; b1=$(binof "$t1" "$l1" "$n1")
    # SOLO slot support: empty CB -> BIN2/GSQ2/NF2/MRE2 all empty (run_redo_mps2 runs client 1 only)
    if [ -n "${CB[$idx]:-}" ]; then b2=$(binof "$t2" "$l2" "$n2"); else b2=""; g2=""; n2=""; m2=""; fi
    # --- pretty job name: p<idx>_[a1]L<l>N<n>g<gsq*10>[m<mRe*10>][_[L..][N..]g<gsq*10>[m..]] ---
    # Single letter per variable (L, N=Nf, g=gsq*10, m=mRe*10); a1 = at=0.1 (ml01), no tag = at=0.2.
    # Stream-2 repeats L/N only when it differs from stream-1 (so a same-Nf pair reads e.g. L3N2g15_g30).
    # p<idx> prefix guarantees uniqueness => afterany chaining (squeue --name below) stays correct.
    local attag=""; [ "$t1" = "ml01" ] && attag="a1"
    local s1="L${l1}N${n1}g$(g10 "$g1")"
    [ "$(awk "BEGIN{print ($m1!=0)}")" = 1 ] && s1="${s1}m$(g10 "$m1")"
    local s2=""
    if [ -n "${CB[$idx]:-}" ]; then
        local p2=""; [ "$l2" != "$l1" ] && p2="${p2}L${l2}"; [ "$n2" != "$n1" ] && p2="${p2}N${n2}"
        s2="_${p2}g$(g10 "$g2")"
        [ "$(awk "BEGIN{print ($m2!=0)}")" = 1 ] && s2="${s2}m$(g10 "$m2")"
    fi
    local jobname="p${idx}_${attag}${s1}${s2}"
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
