# Half-$a_t$ ensembles (fixed $N_t$) — impl plan (_claude, 2026-07-22)

## Goal
Add a new set of massless ensembles at **half the temporal spacing**, $a_t = 0.1$ (currently
$a_t = 0.2$), **keeping $N_t = 128$ fixed**. All $N_f = 2,4,6$, but only at the **middle $g^2$ at
each $L$**. This is a temporal-discretization probe at fixed $N_t$.

Middle $g^2$ per $L$ (from the production grid):
- $L=1$: $g^2 \in \{0.5, 1.0, 1.5\}$ → **middle $g^2 = 1.0$**
- $L=2$: $g^2 \in \{1.0, 2.0, 3.0\}$ → **middle $g^2 = 2.0$**
- $L=4$: $g^2 \in \{2.0, 4.0, 6.0\}$ → **middle $g^2 = 4.0$**  (include? — see Open Q1)

Ensemble count: 3 $N_f$ $\times$ {L1,L2} = **6**, or $\times$ {L1,L2,L4} = **9**.

## Physical note (confirm)
$T = N_t\, a_t$ **halves**: $128 \times 0.2 = 25.6 \;\to\; 128 \times 0.1 = 12.8$. Fixing $N_t$ while
halving $a_t$ shrinks the physical temporal extent by 2. Assumed intended ("fixing $N_t$").

## How $a_t$ is set today (investigated)
- `hmc_hasenbusch_block_fermilab_claude.cu:182`: `const double at = 0.2;` — **hardcoded**, NOT a CLI
  arg (CLI = gsq, Nf, nu0, mass_re, mass_im, max_sec) and NOT a `-D` macro. `Nt=128` is `constexpr`
  (`:62`).
- `at` feeds: `WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0)` (`:184`) and `Action SW(gsq, at, base)`
  (`:214`).
- dir3 name embeds `at` via `std::to_string(at)` (`:236`) → $a_t=0.1$ auto-writes `at0.100000...`
  dirs, distinct from `at0.200000` — **no collision** with existing ensembles. Good.
- Admissibility assert (`:183`): $\sqrt{3}\,\overline{\ell}/a_t - 4/\sqrt{3} > 0$. Halving $a_t$
  **increases** the LHS → still satisfied. No issue.

## THE physics risk: Zolotarev frozen window is $a_t$-specific
`frozen_window_claude.h` sets, per $L$, the fixed $(\lambda_\text{min}, \lambda_\text{max})$ that
bracket the $D_W$ singular-value range — **tuned at $a_t = 0.2$** (L1 (0.1,13), L2 (0.06,8),
L4 (0.008,5)). The Wilson-Dirac UV edge carries a temporal hopping $\sim 1/a_t$, so halving $a_t$
pushes $\lambda_\text{max}$ up (roughly $\times 2$ in the worst case); $\lambda_\text{min}$ (IR,
physical) shifts less. **Reusing the $a_t=0.2$ window at $a_t=0.1$ risks "eval above/below Zolotarev
window" → sign-approx degradation → wrong/rejected trajectories.**

Precedent: the shelved L8 study ran $a_t=0.1$, $N_t=128$ (npole 13) and passed startup physics — so
$a_t=0.1$ is feasible; the L1/L2/L4 windows just need (re)setting. Procedure = smoke-thermalize with a
provisionally WIDE window, read the measured $\lambda_\text{min}/\lambda_\text{max}$ + watch the
"eval outside window" warnings, then freeze new per-$(L,a_t)$ values. **NM owns the window-tuning
procedure (he set the $a_t=0.2$ ones) — see Open Q3.**

## Files to modify
1. `hmc_hasenbusch_block_fermilab_claude.cu` — make `at` a compile-time macro:
   `#ifndef AT_VAL` / `#define AT_VAL 0.2` / `const double at = AT_VAL;` (mirrors -DKMAX/-DKRNG/-DNF/
   -DLREF). Default 0.2 keeps every existing build byte-identical; new binaries pass `-DAT_VAL=0.1`.
2. `includes/frozen_window_claude.h` — add $a_t=0.1$ windows. Either an $a_t$-aware overload
   `frozen_window(L, at, lmin, lmax)` or a parallel set guarded on a build flag. VALUES = TBD from the
   smoke (Open Q3).
3. `launch_redo_claude.sh` (BOTH ~/ and production copies) — add BUILDS entries for the half-$a_t$
   binaries, e.g. `hmc_fermilab_redo_massless_L<L>_Nf<Nf>_at0p1_claude.o ... -DAT_VAL=0.1` (+ the L1/L2
   with the new window). New name suffix so they don't clobber the $a_t=0.2$ binaries.
4. `run_wrapper_redo_claude.sh` (production) — add slots for the new ensembles (new indices, MPS-paired
   by like cost). Keep the $a_t=0.2$ slots untouched.
5. `params_L1L2_claude.md` (+ new `params_half_at_claude.md`) — record the $a_t=0.1$ params, window,
   nsteps, KMAX.

## Ordered chunks
- **C1 — code param**: `-DAT_VAL` macro in the driver (+ keep default 0.2). Files: driver.
- **C2 — window**: smoke one $(L,a_t{=}0.1)$ per $L$ with a wide window, read extrema, set frozen
  values. Files: frozen_window_claude.h, a smoke build. (Blocks C3/C4.)
- **C3 — build**: add BUILDS entries, compile the 6 (or 9) binaries. Files: both launchers.
- **C4 — wrapper + launch**: add slots, MPS pairs, `--apply 2` (opp, chain length 2). Files: wrapper.
- **C5 — validate**: startup delta/admiss/alat/acc/|dH|, no "eval outside window". Files: none.

## Resolved (NM 2026-07-22: "same HMC parameters, won't affect much")
- **Same HMC params** as production: KMAX 4000, KRNG 20, current per-$L$ nsteps, $\nu_0=1.0$, massless
  ($m_\text{Re}=0$). → Q2, Q4, Q5 closed.
- **Reuse the existing frozen Zolotarev window** per $L$ (NM judges the $a_t$ shift negligible). → Q3
  closed: `frozen_window_claude.h` is NOT edited. C2 becomes smoke-VALIDATE only (watch for "eval
  outside window" warnings + measured $\lambda$ extrema at startup; revisit window ONLY if they fire).

## Still open
1. **Include $L=4$** (middle $g^2=4.0$) → 9 ensembles, or only $L=1,2$ → 6? ($L=4$ was just dropped
   from the 4000-stat production; this is a separate probe, so it may or may not want $L=4$.)

## Revised chunks (window retune dropped)
- **C1** — `-DAT_VAL` macro in driver (default 0.2 preserves all existing builds). Files: driver.
- **C2** — build the 6/9 half-$a_t$ binaries with `-DAT_VAL=0.1` (distinct `_at0p1` name). Files: both
  launchers (BUILDS).
- **C3** — wrapper: add slots + MPS pairs for the new ensembles; extend binary resolution for the
  `_at0p1` tag. Files: run_wrapper_redo_claude.sh.
- **C4** — smoke 1 slot/$L$ (qos=test): startup delta/admiss/alat/acc/|dH| + NO "eval outside window".
- **C5** — `--apply 2` (opp, chain length 2). Record params in `params_half_at_claude.md`.
