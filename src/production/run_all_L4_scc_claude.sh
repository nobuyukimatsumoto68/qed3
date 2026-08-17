#!/bin/bash -l
# run_all_L4_scc_claude.sh  (_scc, 2026-07-22, NM)  -- COMPLETELY UNIFIED L=4 SCC driver
# =============================================================================================
# ONE command builds + submits (or EXTENDS) every L=4 massless stream, each pinned to a chosen GPU:
#   Nf2 = SCC-native production ; Nf4/6 = FNAL handoff (rsync'd configs, resumed here).
# A per-ensemble TABLE (below) maps each (Nf,gsq) to a build arch, and the wrapper maps arch -> gpu_type:
#   sm_70 -> V100 (FP64 1:2), sm_80 -> A100 (FP64 1:2), sm_90 -> H200 (strongest FP64).  [per NM 2026-07-22]
# This HMC is FP64-heavy, so we pin to strong-FP64 GPUs and stay OFF the weak-FP64 L40S/A40 that -l gpu_c
# (a MINIMUM only) would otherwise allow. Rationale + numbers: resume_fnal_L4_scc_impl_plan_claude.md.
#
# CHAINING / RESUME: every ensemble AUTO-RESUMES from its highest checkpoint (Nf2 ~k253, Nf4/6 ~k47-52).
# To CHAIN FURTHER just RE-RUN this script: the wrapper's existing_tail() finds each live chain (by ensemble
# token, arch-independent) and -hold_jid's N_CHAIN new links onto its tail -- so a re-run EXTENDS, never
# starts a concurrent 2nd writer. Each run adds N_CHAIN links (default 4 x 12h). Bump N_CHAIN for more at once.
# Binaries are rebuilt ONLY if missing (or FORCE_BUILD=1); the compile is one-time, scheduling changes are not.
#
# Usage:
#   bash run_all_L4_scc_claude.sh                 # build-if-missing + submit/extend ALL 9 streams
#   DRYRUN=1 bash run_all_L4_scc_claude.sh        # print qsub lines, build nothing, submit nothing
#   N_CHAIN=8 bash run_all_L4_scc_claude.sh       # append 8 links/chain this run (walk further toward KMAX)
#   ONLY="Nf6" bash run_all_L4_scc_claude.sh      # only ensembles whose "Nf gsq arch" row contains Nf6
#   ONLY="sm_90" bash run_all_L4_scc_claude.sh    # only the H200 stream(s)
#   FORCE_BUILD=1 bash run_all_L4_scc_claude.sh   # recompile even if binaries exist (after a source change)
# =============================================================================================
set -u

SRCDIR=/projectnb/qfe/nmatsum/qed3/src/production
cd "$SRCDIR" || { echo "ERROR: cannot cd $SRCDIR"; exit 1; }

# ---- knobs (EDIT) -------------------------------------------------------------------------------
KMAX=${KMAX:-600}                    # trajectory target (all ensembles); bumped 400->600 (2026-08-01). KMAX is
                                     # COMPILE-TIME (-DKMAX, in the binary name k<KMAX>); the output dir does NOT
                                     # encode it, so a k600 binary RESUMES the same dir -> Nf2 (at 399) restarts
                                     # -> 599, Nf4/6 continue past 400. New k600 binaries are built (k400 ones stay).
PE_OMP=${PE_OMP:-4}                  # CPU slots/job (4 -> packs public GPU slots=16/user; jobs are GPU-bound)
N_CHAIN=${N_CHAIN:-4}                # dependent links appended PER ensemble PER run (each ~12h)
DRYRUN=${DRYRUN:-0}
FORCE_BUILD=${FORCE_BUILD:-0}
ONLY=${ONLY:-}                       # optional filter: space/comma list of ensemble tokens (e.g.
                                     # "Nf2g4 Nf6g6") and/or an arch ("sm_70"). Empty = all 9. Each pattern is
                                     # matched as a substring against a row's token "Nf<nf>g<gsq>" + its arch.

WRAP=run_wrapper_L4_scc_claude.sh

# ---- per-ensemble GPU assignment: "Nf  gsq  arch" (arch -> gpu_type in the wrapper) --------------
# sm_70=V100, sm_80=A100, sm_90=H200. Edit an arch here to MOVE that ensemble to a different GPU on the
# next run (qdel its old chain first so it does not just extend on the old GPU).
# 2026-07-22: switched ALL to V100 -- qfe has buy-in priority on ece V100 nodes + broad public V100 access
# (26 nodes), whereas A100/H200 are other groups' buy-in (qfe locked out -> those jobs sat in qw). V100 is
# strong FP64 and needs NO rebuild (sm_70 binaries exist). To retry H200/A100 later, set that row's arch back.
# 2026-08-14: retried A100 on Nf6g6 (the long pole). Measured HMC 3773 vs V100 4543 s/traj = 1.2x faster
# (below the conn benchmark's ~1.5x: HMC's FP64-heavy action/Metropolis solve erodes the gain, and the V100
# baseline is SXM2 not the conn's weak PCIE-16GB). BUT public A100 -pub queues were disabled/full (the 5-6 free
# A100 GPUs sat behind PRIVATE buy-in qfe cannot schedule) -> the chain's 2nd link stalled 3.6h in qw, past the
# ~2.4h/link break-even. Reverted to V100. CONCLUSION: V100 is the DEFAULT. A100 is only worth retrying when the
# -pub A100 queues are ENABLED and actually have a free GPU (qstat -f -l gpu_type=A100 | grep -- '-pub@'); the
# 1.2x is real but not worth the queue risk otherwise. V100 reallocates each chain link in ~2.5 min (measured).
ENSEMBLES=(
  "2 2.0 sm_70"     # Nf2 g2 -> V100   (SCC-native)
  "2 4.0 sm_70"     # Nf2 g4 -> V100   (was A100)
  "2 6.0 sm_70"     # Nf2 g6 -> V100
  "4 2.0 sm_70"     # Nf4 g2 -> V100   (FNAL cont.)
  "4 4.0 sm_70"     # Nf4 g4 -> V100   (was A100)
  "4 6.0 sm_70"     # Nf4 g6 -> V100
  "6 2.0 sm_70"     # Nf6 g2 -> V100   (was A100)   (FNAL cont.)
  "6 4.0 sm_70"     # Nf6 g4 -> V100
  "6 6.0 sm_70"     # Nf6 g6 -> V100   (H200 too contended; A100 1.2x but -pub queues disabled -> see 08-14 note)
)

echo "######## unified L4 SCC driver  [$(date +%F_%H:%M:%S)]  (KMAX=$KMAX N_CHAIN=$N_CHAIN PE_OMP=$PE_OMP DRYRUN=$DRYRUN ONLY='${ONLY}') ########"

for e in "${ENSEMBLES[@]}"
do
  read -r nf gsq arch <<< "$e"
  tok="Nf${nf}g${gsq}"
  if [ -n "$ONLY" ]
  then
    match=0
    for pat in ${ONLY//,/ }
    do
      if [[ "$tok $arch" == *"$pat"* ]]
      then
        match=1
        break
      fi
    done
    if [ "$match" -eq 0 ]
    then
      continue
    fi
  fi

  # build only if the binary is missing (or FORCE_BUILD); DRYRUN never builds
  bin="hmc_L4_scc_Nf${nf}_k${KMAX}_${arch}.out"
  nb=1
  if [ "$DRYRUN" -eq 0 ] && { [ "$FORCE_BUILD" -eq 1 ] || [ ! -f "$bin" ]; }
  then
    nb=0
  fi
  [ "$nb" -eq 0 ] && act="BUILD+submit" || act="submit/extend (binary present)"

  echo
  echo "===== Nf${nf} gsq${gsq} -> ${arch} : ${act} ====="
  PE_OMP="$PE_OMP" N_CHAIN="$N_CHAIN" \
    ARCH_LIST="$arch" SUBMIT_ARCHS="$arch" \
    ML_NF="$nf" ML_GSQ="$gsq" ML_MASS="0.0" ML_KMAX="$KMAX" \
    WHICH=massless NOBUILD="$nb" DRYRUN="$DRYRUN" \
    bash "$WRAP"
done

echo
echo "######## unified L4 SCC driver done  [$(date +%F_%H:%M:%S)] ########"
