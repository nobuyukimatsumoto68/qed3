# Resume FNAL L4 configs on SCC -- impl plan (_claude, 2026-07-22)

**Plan change (NM, 2026-07-22):** the half-$a_t$ study moves to FNAL. Instead, SCC **rsyncs L4 gauge
configs from FNAL and RESUMES HMC generation on them here** (offload FNAL). Supersedes the half-$a_t$ SCC
work (`halfat_L4_scc_impl_plan_claude.md` + `run_wrapper_L4_halfat_scc_claude.sh` are now MOOT -- pending
NM's keep/remove call; the `-DAT_VAL` knob left in the `_scc` driver is harmless (default 0.2) and keeps it
in sync with the FNAL driver).

## What "resume" requires (vs the existing measurement rsync)
The blackboard's FNAL->local rsync pulls **`ckpoint_lat.*` ONLY** (gauge configs, for measurement). To
**continue the Markov chain** the driver resumes from the highest $k$ with BOTH files:
- `ckpoint_lat.<k>` (~0.64 MB) -- gauge config
- `ckpoint_rng.<k>` (~0.5 GB) -- RNG state  <-- **must also be pulled** (measurement rsync skips it)

Only the **latest** $k$'s rng is needed to resume (older rng not required). So: pull ALL `ckpoint_lat.*`
(cheap, and gives the full chain for later measurement) + the **single latest** `ckpoint_rng.<k_max>`.

## Validity of cross-machine resume
- Same L4 geometry on both (`geometry/data/` n4: pts/links/nns/faces/duals + omega_n4/alpha_n4 -- verified
  present on SCC, the Nf2 L4 runs use them). A config is just link angles on that lattice -> portable.
- Same action/integrator params -> the `_scc` driver writes the IDENTICAL output dir name
  (`Nf<nf>_gsq<g>at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000/`), so dropping the
  FNAL checkpoints into that dir on SCC makes the driver auto-resume from $k_max$.
- HMC is Metropolis-exact; a different binary/arch (SCC sm_70/sm_80, mixed-precision force) continuing the
  chain is still correct. The RNG state carries the chain forward.
- **MD steps must match FNAL** to keep the tuning: Nf4/6 -> `-DL4_MDSTEP=5` ({5,5,5}). Note the step count is
  NOT in the dir name (hbtag = ladder coeffs `0.400000-1.000000`, not steps), so a wrong step count would
  still resume the same dir but change acceptance -- so we pass `-DL4_MDSTEP=5` explicitly.

## FORK HAZARD (the key correctness question -- Q2)
If FNAL keeps advancing an ensemble while SCC also resumes the SAME latest config, the chain FORKS: both
share history up to $k_max$, then diverge into two different continuations. That is NOT one chain -- configs
past the fork are correlated near the seam and cannot be naively concatenated. **Clean handoff = FNAL must
STOP that ensemble before SCC resumes it** (or SCC takes ensembles FNAL never launched). Needs NM's word.

## Plan (pending the decisions below)
1. **Build SCC Nf4/6 L4 binaries** from `hmc_hasenbusch_block_scc_claude.cu`:
   `-DLREF=4 -DNF={4,6} -DKMAX=400 -DKRNG=4 -DL4_MDSTEP=5 -DMIXED_FORCE`, both arches (sm_70/sm_80).
   (The full-$a_t` wrapper `run_wrapper_L4_scc_claude.sh` builds per (Nf,KMAX,arch) but does NOT pass
   `-DL4_MDSTEP` -- it must be extended with the per-$N_f$ mdstep logic, same as the now-moot half-$a_t`
   wrapper had, so Nf4/6 build as {5,5,5} not the default {4,4,4}.)
2. **rsync FNAL -> SCC** (NM runs it; he has the FNAL kerberos/ssh auth -- Claude never runs it). Per
   ensemble dir, pull all `ckpoint_lat.*` + the latest `ckpoint_rng.<k_max>` into the matching SCC dir.
   Command drafted in `tmp_rsync_fnal_L4_claude.sh` (NM edits FNAL k_max / dest, then runs). NO rm anywhere.
3. **Resume**: run the (extended) full-$a_t` wrapper with `ML_NF="4 6"` `ML_GSQ="<selected>"` `ML_KMAX=400`
   -- the driver finds $k_max$ (lat+rng) in each dir and continues to 400. `existing_tail()` anchoring +
   `L4<arch>` namespace already handle chaining (these are new Nf4/6 chains, no clash with the Nf2 ones).
4. Update the blackboard L4 rows (Nf4/6 now generated on SCC) + memory.

## Decisions (resolved with NM, 2026-07-22)
1. **All 6** FNAL L4 Nf4/6 massless x gsq{2,4,6} (rows p13/14/15).
2. **Clean handoff** -- FNAL STOPPED these; SCC owns the chains now. No fork hazard.
3. **KMAX=400**, and the rsync pulls the **latest `ckpoint_rng.<k_max>`** per ensemble (+ all `ckpoint_lat.*`).
4. **No MPS** (NPACK=1, one ensemble/GPU) -- already the wrapper default.
5. **Auth = kerberos**: NM runs `kinit` before the rsync; the handoff script only warns (`klist -s`), never
   runs kinit. rsync flags `-avz` (compress, NO `-c` checksum, NO `--delete`).

## Implemented (2026-07-22)
- `run_wrapper_L4_scc_claude.sh`: added `mdstep_of()` + `-DL4_MDSTEP` in `build_all` -> Nf4/6 build as {5,5,5}
  (Nf2 stays {4,4,4}; binaries byte-identical to before). Dry-run verified: `ML_NF="4 6"` yields 6 chains,
  binaries `hmc_L4_scc_Nf{4,6}_k400_<arch>.out`, namespace `L4<arch>_Nf{4,6}g...` (no clash with Nf2 tokens).
- `tmp_rsync_fnal_L4_claude.sh`: NM-run handoff (kinit preflight; per-ensemble pull of all lat + latest rng;
  additive, no rm/--delete). NO submit/build/rm by Claude.

## Run order (NM)
1. `kinit <you>@FNAL.GOV`
2. `RSYNC_DRYRUN=1 bash tmp_rsync_fnal_L4_claude.sh`   (preview) then `bash tmp_rsync_fnal_L4_claude.sh` (pull)
3. `PE_OMP=4 ML_NF="4 6" ML_GSQ="2.0 4.0 6.0" ML_KMAX=400 WHICH=massless bash run_wrapper_L4_scc_claude.sh`
   (builds the Nf4/6 {5,5,5} binaries, then submits 6 chains that auto-resume from the pulled k_max)

## Scheduler / allocation notes (2026-07-22)
- `-P qfe` is set in the batch (`run_L4_scc_claude.sh:22`). qfe has NO dedicated GPU queue -- it is a member
  project on the `ece`/`ece-long` buy-in GPU queues (`engineering` buy-in is CPU-only). ECE per-user RQS caps:
  **gpus=2 in `ece`, gpus=4 across `ece`+`ece-long`** -- the "max" NM hit before.
- The wrapper `qsub` has NO `-q`/`-l buyin`, so once the ECE cap is full, eligible chains SPILL to the public
  GPU queues (a100/v100) automatically. Public caps: shared `a100`-class = gpus=2 AND slots=16 per user;
  `academic-gpu`/`*-pub` = gpus=4 each. Public a100 is acceptable (NM 2026-07-22).
- **GPU vs CPU slots are independent.** `-l gpus=1` = one whole GPU in EXCLUSIVE_PROCESS mode (1 process/GPU;
  MPS was the only way to put 2 streams on 1 GPU -- disabled, NPACK=1). `-pe omp N` = N CPU cores. "More jobs
  on public a100" = more SEPARATE single-GPU jobs (each its own exclusive GPU), NOT GPU sharing. The public
  slots=16/user cap is what limited that -> **PE_OMP=4** for the whole L4 run (Nf2 too, per NM) frees the slot
  budget so the per-queue GPU caps become the only limit. Jobs are GPU-bound, so 4 host threads is fine.

## FP64 / gpu_type pin (2026-07-22) -- IMPORTANT
This HMC is DOUBLE-PRECISION throughout (`cuDoubleComplex`), so GPU **FP64** rate is the bottleneck, NOT FP32.
`-l gpu_c` is only a MINIMUM compute capability, so jobs silently landed on **L40S** (the ece buy-in GPUs,
cc 8.9) which have WEAK FP64 (1:64, ~1.4 TFLOP/s) -> a full trajectory there is far slower than on a
strong-FP64 GPU. The fast "140s on sm80" run was on an actual **A100** (full-rate FP64), not L40S.
- Strong FP64 available on SCC: **A100** (sm_80, ~9.7 TFLOP/s, 9 nodes), **V100** (sm_70, ~7 TFLOP/s, 26
  nodes), **H200** (sm_90, strongest, 8 nodes -- would need a `-arch=sm_90` build). Weak FP64 (AVOID): L40S,
  A40, A6000.
- **Fix (wrapper):** `gput_of()` + `-l gpu_type=${gput}` pins each arch to its native strong-FP64 model:
  **sm_70 -> V100, sm_80 -> A100** (knobs `GPUT_SM70`/`GPUT_SM80`; empty = no pin). Dry-run verified.
- **Tradeoff:** V100/A100 are public-queue only for qfe (ece buy-in is L40S), so the shared-GPU RQS caps
  concurrency tighter (~2 core-shared, up to 4 academic/pub) vs 4 on the L40S buy-in -- but each job is ~5-7x
  faster in FP64, so net throughput still wins. Watch actual concurrency once pinned.
- **Operational:** jobs submitted BEFORE this change carry no gpu_type -> still on L40S. To move them, qdel +
  re-run `run_all_L4_scc_claude.sh` (re-pins, resumes from the same checkpoints). NM qdels; Claude never does.
- **FINAL (2026-07-22): ALL streams -> V100.** A100/H200 free GPUs are in OTHER groups' locked buy-in queues
  (li-rbsp/neuro-autonomy/labci/chapman) -> qfe A100/H200 jobs sat in `qw`. V100 = qfe buy-in priority (ece
  scc-q3x) + 26-node public pool -> 9/9 run concurrently. Driver table is now all `sm_70`; no rebuild. Live:
  Nf2 ~72% (s/traj ~2100-2380, ETA ~3d), Nf4/6 ~17% (~2940-3700 s, ETA ~2wk). Daily check-up via check-threads.

## Unified driver + chaining + H200 (2026-07-22)
- **`run_all_L4_scc_claude.sh` = the single entry point.** A per-ensemble table maps each (Nf,gsq) -> build
  arch -> gpu_type: sm_70=V100, sm_80=A100, sm_90=H200. It builds-if-missing (else NOBUILD) and submits/EXTENDS
  every stream. Current map: Nf2 g2/g6, Nf4 g2/g6, Nf6 g4 -> V100; Nf2 g4, Nf4 g4, Nf6 g2 -> A100; Nf6 g6 -> H200.
- **Chaining = re-run the driver.** `existing_tail()` in `run_wrapper_L4_scc_claude.sh` (qstat -r snapshot,
  token `Nf<nf>g<gsq>m<mass>`, arch-independent) `-hold_jid`s N_CHAIN new links onto each live chain's tail ->
  EXTENDS, never a concurrent 2nd writer. `bash run_all_L4_scc_claude.sh` adds N_CHAIN(=4) links; `N_CHAIN=8 ...`
  for more. To MOVE an ensemble to another GPU: edit its arch in the table AND qdel its old chain first.
- **H200 (sm_90) for Nf6 g6:** one-time `-arch=sm_90` build (cuda/12.8 supports Hopper); submitted via
  `ARCH_LIST=sm_90 SUBMIT_ARCHS=sm_90 ML_NF=6 ML_GSQ=6.0 ...` (namespace L490). Heaviest force -> best FP64
  showcase; compare its per-traj time vs V100/A100 to decide whether to move more ensembles to H200.
- **qdel gotcha:** killing a chain's c0 releases its held c1-c3 (they run). Kill ALL links at once; use the
  one-line self-updating form (`qdel $(qstat -u "$USER" | awk 'NR>2 && $3 ~ /^L4[78h9]/{print $1}')`) so a long
  explicit ID list is not chopped by a terminal line-wrap (that leak happened 2026-07-22).

## Refs
M. Hasenbusch hep-lat/0107019; B. Jegerlehner hep-lat/9612014; A. D. Kennedy hep-lat/0402038.
