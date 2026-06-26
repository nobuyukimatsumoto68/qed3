# Zolotarev n-bump + reset-window-from-current-config (L4) — impl plan

**Date:** 2026-06-26. **Goal:** stop the recurring L4 force spikes caused by Wilson
zero-crossings (a near-zero eigenvalue of $D_W^\dagger D_W$ escaping the Zolotarev
sign-approximation window). See `[[l2l4-refinement-run]]` (the mRe0.211450 / mRe0.052862
episodes) and `[[mass-measure-factor]]` (sign solve is mass-independent).

Algorithm source: A. D. Kennedy, hep-lat/0402038 (Zolotarev sign-function approximation);
the `Zolotarev` struct in `includes/overlap_wmass_claude.h` cites it (`:55`).

## Diagnosis (confirmed from code)

`OverlapWMass : public Zolotarev`. The overlap sign function $\text{sign}(H_W)$ is a
Zolotarev rational approximation accurate to error $\Delta$ for eigenvalues of
$\sqrt{D_W^\dagger D_W}$ in the window $[k\,\lambda_\text{max},\,\lambda_\text{max}]$, where
$k = $ Zolotarev parameter (constructor `k_`, default 0.01) and $n = $ degree (constructor
`n_`).

- L4 driver builds `Fermion D(DW, mass, 13)` (`hmc_fermilab_wmass_L4_claude.cu:187`):
  $n=13$, $k=0.01$.
- `update(U)` (`overlap_wmass_claude.h:289`) measures $\lambda_\text{min},\lambda_\text{max}$
  each call (`compute_lambda_max`, power + inverse-power iteration) and, IF
  `lambda_min/lambda_max < 0.1 k` AND `is_update`, RE-FITS the window to
  $k \to 0.1\,\lambda_\text{min}/\lambda_\text{max}$, printing `# Smaller Delta Detected`.
- Driver freezes the window: `if(k%20==0) D.is_update = false;` (`:320`) — `is_update`
  starts `true` (`overlap_wmass_claude.h:213`), and once the traj counter passes a multiple
  of 20 it is set false and never re-enabled. So the window adapts over the first ~20 configs
  of a run, then is FROZEN.

**Failure mode:** after the freeze, a later config develops a smaller eigenvalue (a
topology-change attempt / Wilson zero-crossing). The frozen window no longer covers it; at
$n=13$ the re-fit (when it does fire) only reaches $\Delta\sim1.2\times10^{-2}$ (vs design
$8.95\times10^{-5}$) — a sign function that inaccurate gives a wrong/non-reversible overlap
force → the huge dH (1812, 6320, 246, ...). All such trajectories REJECT (ensemble stays
clean) but the stream parks and makes no progress.

## $\Delta$ vs $(n,k)$ — why bump $n$

Empirically (memory fit) $\log_{10}\Delta$ drops ~0.357 per added pole at the L4 window.
A degraded window gave $\Delta\sim1.2\times10^{-2}$; to bring that back to $\sim10^{-4}$:
$$
\Delta n \approx \frac{\log_{10}(1.2\times10^{-2} / 10^{-4})}{0.357} \approx \frac{2.08}{0.357} \approx 6 .
$$
So $n = 13 + 6 \approx 19$ → round to the next ODD value **21** (n must be odd; `size=n/2+1`).
This also restores the pre-2026-06-16 pole count. Cost scales ~linearly in $n$ (inner
multishift solve has `size-1` poles): $n:13\to21$ is $size:7\to11$, i.e. ~**1.6x** solve
cost per $D_\text{ov}$ apply — paid mostly on the expensive L4 Nf6 corner.

## Two-part change (per NM)

1. **Bump $n$** in the L4 driver: `Fermion D(DW, mass, 13)` → `Fermion D(DW, mass, <N_NEW>)`.
2. **Reset the covered eigenvalues from the current config**: on restart, re-derive the
   window from the resume config's measured spectrum (so it covers the present small
   eigenvalues), rather than inheriting a window frozen on early configs. OPEN: exact
   mechanism — see Open questions.

## Marking (per NM) — to add as a dated comment at the driver change site

```
// 2026-06-26: Zolotarev n 13 -> <N_NEW> + reset window from current config, to cure the
// Wilson zero-crossing force spikes (Smaller Delta Detected -> huge REJECTED dH). Change
// takes effect on restart from these checkpoints:
//   Nf2  k=119 (all 4 masses; DONE/at cap)
//   Nf4  mRe0.010572 k=53  mRe0.052862 k=42  mRe0.105725 k=45  mRe0.211450 k=54
//   Nf6  mRe0.010572 k=31  mRe0.052862 k=25  mRe0.105725 k=25  mRe0.211450 k=29
```

## Files
- `hmc_fermilab_wmass_L4_claude.cu` — bump $n$ (`:187`), the reset mechanism (around the
  startup `D.update(U)` `:189` and/or the `k%20` freeze `:320`), and the dated marking.
- (possibly) `includes/overlap_wmass_claude.h` — only if the reset needs a new method
  (e.g. `set_window_from_spectrum()`); preferred to keep the logic in the driver.

## RESOLVED design (NM, 2026-06-26) — FINAL, implemented

NM simplified: **no adaptive tuning at all** — a single FIXED window, set at construction.

1. **$n = 21$** (from 13). ~1.6x inner-solve cost (paid mostly on Nf6).
2. **FIXED window $k_=0.001$** (from 0.01), so eigenvalue dips ~10x deeper are inside the
   window. `Fermion D(DW, mass, 21, 0.001);` (`hmc_fermilab_wmass_L4_claude.cu:198`). $n=21$
   keeps $\Delta$ small at this wider window.
3. **Adaptive re-fit REMOVED from `OverlapWMass::update()`** (`overlap_wmass_claude.h:~293`):
   the `if(lambda_min/lambda_max < 0.1*k){ if(is_update) Zolotarev::update(...) }` block is
   gone (old code left commented). Replaced by a PASSIVE `# WARNING: eval below Zolotarev
   window` print (under `InfoDelta`) that fires if a dip still escapes the fixed $k$ — a spike
   predictor + the signal to lower $k$ further. `is_update` is now vestigial; no startup reset
   block. The periodic `if(k%20==0) is_update=false` (`:~330`) is commented out.
4. Startup print `# Zolotarev window k = <k> (fixed for the run)`.

**Why fixed (not per-config / per-trajectory):** the adaptive path was confusing and a latent
reversibility hole (`update()` runs per MD step → could re-fit $k$ mid-trajectory). A single
fixed window is unambiguous and exactly reversible. Robustness to deep dips comes from the
wider $k=0.001$ + $n=21$, not from adapting.

**SHARED-HEADER NOTE:** `overlap_wmass_claude.h::update()` is used by L2 + ~20 measurement
programs. The removal only takes effect when each is REBUILT; the handoff rebuilds **L4 only**,
so the running L2 binary + measurement binaries are unchanged for now. L2 is behavior-neutral
(its adaptive path never fired). Measurement programs (jj/condensate, default $k=0.01$) inherit
a fixed window at their next build — verify their $k$ is adequate before their next use (the
passive WARNING will flag any escape).

**VALIDATION (NM, on first restarted job):** confirm printed `# delta` < ~1e-4 at
$k=0.001$/$n=21$ (else bump $n$); watch first-traj |dH| + any `# WARNING: eval below ...` on the
previously-parked streams (Nf6 pairA heavy mRe0.211450, Nf4 pairB mRe0.052862). If WARNINGs
recur, lower $k$ (e.g. 1e-4) and rebuild.

**Cost:** $n$=21 + wider window raises L4 cost ~1.7-2x (more poles + the extreme-pole CG is
more ill-conditioned) — ties into the parked mixed-precision-CG idea.

**Restart scope:** kill all L4 (Nf4+Nf6; Nf2 done/dropped) → rebuild L4 → relaunch from current
checkpoints (NF_LIST default `4 6`). Handoff in `/lustre2/qed3/tmp_claude.sh`.
