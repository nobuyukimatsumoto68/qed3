# Half-$a_t$ L4 ensembles (SCC) -- impl plan

**Goal (NM, 2026-07-22):** add a new set of L=4 ensembles at **half the temporal lattice spacing**
$a_t = 0.1$ (was $0.2$), **with $N_t$ fixed at 128** (so the physical time extent $N_t a_t$ halves,
$25.6 \to 12.8$). All **$N_f = 2, 4, 6$**, but only at the **middle gsq per $L$** (for L4: gsq $= 4.0$).
**Reuse the existing HMC parameters as-is** (NM: "I don't believe they'll be affected so much").
**Target half KMAX** ($400 \to 200$ for massless).

This is the SCC (L4) slice of a broader "middle-gsq, all-$N_f$, half-$a_t$" campaign ("at each $L$"). L1/L2
live at FNAL/affine; SCC does L4 (wrapper is `-DLREF=4`, geometry is L4). See [[scc-l4-production]].

## How $a_t$ enters (delegated) 
- The base driver currently hardcodes `const double at = 0.2;` (`hmc_hasenbusch_block_scc_claude.cu:199`).
- **Another agent is adding a `-DAT_VAL` compile-time knob** to the driver; we adapt to that. `at` is already
  baked into the output dir name (`...at"+std::to_string(at)+"...`, line 253) -> `at0.100000` dirs are distinct
  from the existing `at0.200000` dirs. **No dir collision**, auto-resume still works.

## What changes on our side (wrapper only) -- `run_wrapper_L4_scc_claude.sh`
The wrapper is already fully env-parameterized (`ML_NF`, `ML_GSQ`, `ML_KMAX`). Three surgical edits let it
build+submit the half-$a_t$ set while the default (`AT_VAL=0.2`) reproduces current behavior exactly:

1. **`AT_VAL` knob.** Add `AT_VAL=${AT_VAL:-0.2}`; pass `-DAT_VAL=${AT_VAL}` on the `nvcc` line in `build_all`.
2. **`at` in the binary name.** `binname` -> `hmc_L4_scc_Nf<nf>_k<kmax>_at<AT_VAL>_<arch>.out` so half-$a_t$ binaries
   do NOT overwrite the full-$a_t$ ones. (Output dirs already separate via the dir-name `at` tag.)
3. **Per-$N_f$ MD steps.** In `build_all`, set `mdstep = (nf==2 ? 4 : 5)` and pass `-DL4_MDSTEP=${mdstep}`
   (Nf2 -> {4,4,4}; Nf4/6 -> {5,5,5}, identical to the full-$a_t$ Nf4/6 builds). $N_f$ is already in the
   binary name, so no name collision.

Nothing else changes: tau=1.0, gauge MG100, ladder {0.0,0.4,1.0}, frozen (0.008,5), n_action 31 / n_force 11,
`-DMIXED_FORCE`, both arches, NPACK=1 (no MPS), the `existing_tail()` cross-run anchoring, chaining -- all reused.

## Launch (NM runs; Claude never submits -- see [[never-submit-jobs]])
Massless half-$a_t$, L4 middle gsq, all $N_f$, half KMAX:
```
AT_VAL=0.1 ML_NF="2 4 6" ML_GSQ="4.0" ML_MASS="0.0" ML_KMAX=200 \
  WHICH=massless bash run_wrapper_L4_scc_claude.sh
```
3 ensembles (Nf2/4/6 at gsq4.0), round-robined across sm_70/sm_80, one dependent chain each ($N_{chain}=4$).
Build first (omit `NOBUILD=1`) because these are NEW binaries (different `at`, and Nf4/6 never built on SCC).

## Decisions (resolved with NM, 2026-07-22)
1. **$N_f$ scope = Nf2 ONLY on SCC.** Nf4/6 half-$a_t$ go to FNAL (mirrors the full-$a_t$ L4 split). So the SCC
   slice is a SINGLE ensemble: Nf2, L4, gsq4.0, at=0.1, massless, KMAX=200.
2. **Massless only** (mRe=0).
3. **Half KMAX = 200** (massless 400/2).
4. Nf4/6 MD steps would be {5,5,5} -- but N/A on SCC now (FNAL's). Nf2 keeps {4,4,4} (default). The per-$N_f$
   `-DL4_MDSTEP` build logic is kept in the wrapper anyway (harmless; correct if someone runs Nf4/6 here).
5. **Separate wrapper** `run_wrapper_L4_halfat_scc_claude.sh` (copy of the full-$a_t$ wrapper, adapted). The
   BATCH script `run_L4_scc_claude.sh` is REUSED (one backward-compatible edit -- optional `AT_TAG` in the log
   name). Keeps the two campaigns as separate launchers, one shared batch.

## Namespace separation (CORRECTNESS -- must not cross the full-$a_t$ run)
The full-$a_t$ massless run ALREADY has an ensemble token `Nf2g4.0m0.0` (job `L480_Nf2g4`, running). If the
half-$a_t$ wrapper reused the same job-name namespace, `existing_tail()` would wrongly `-hold_jid` the half-$a_t$
chain onto the full-$a_t$ one. So the half-$a_t$ wrapper:
- Names chains **`L4h<arch>_Nf...`** (prefix `L4h`, vs `L4<arch>` for full-$a_t`) -> distinct scheduler namespace.
- Filters `EXIST_JOBS` to `Full jobname ~ /^L4h/` -> only anchors onto PRIOR half-$a_t$ chains, never full-$a_t`.
- Binary name carries the $a_t$ tag: `hmc_L4_scc_Nf<nf>_k<kmax>_at<AT_VAL>_<arch>.out` (full-$a_t` has no at tag).
- Batch log name gets an optional `_at<AT_VAL>` (via `AT_TAG` -v) so logs are campaign-distinguishable.

The `check-threads` skill monitor (`.claude/skills/check-threads/check_threads_claude.sh`) qstat filter is
widened to also match `L4h<arch>` jobs.

## Launch (NM runs; Claude never submits -- see [[never-submit-jobs]])
```
bash run_wrapper_L4_halfat_scc_claude.sh          # defaults: AT_VAL=0.1 Nf2 gsq4.0 mRe0.0 KMAX200, build+submit
DRYRUN=1 bash run_wrapper_L4_halfat_scc_claude.sh # preview qsub lines
```

## Refs
M. Hasenbusch hep-lat/0107019; B. Jegerlehner hep-lat/9612014; A. D. Kennedy hep-lat/0402038.
