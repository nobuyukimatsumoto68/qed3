# Frozen Zolotarev window for the massless overlap (per L)

**Fixed 2026-07-13 (NM).** These are the production **frozen** $(\lambda_\text{min}, \lambda_\text{max})$ for the
massless overlap operator, one window per lattice refinement $L$. They are held constant (no per-config
recompute) so the sign approximation is fixed, ensemble-covering, exact, and reproducible.

## The numbers ($n = 21$ poles)

| $L$ | $\lambda_\text{min}$ | $\lambda_\text{max}$ | $k = \lambda_\text{min}/\lambda_\text{max}$ | $\Delta$ |
|-----|------|------|------|------|
| 1 | 0.1   | 13 | $7.69\times10^{-3}$ | $2.54\times10^{-7}$ |
| 2 | 0.06  | 8  | $7.50\times10^{-3}$ | $2.72\times10^{-7}$ |
| 4 | 0.008 | 5  | $1.60\times10^{-3}$ | $7.08\times10^{-6}$ |

- $k$ is the Zolotarev window ratio; $\Delta = (\lambda_\text{inv}-1)/(\lambda_\text{inv}+1)$ is the sign-approx
  error, a function of $(n,k)$ only (NOT of the gauge config). Computed by `zolotarev_delta_claude.cpp`.
- All three beat the old fixed window $k=10^{-3}$ ($\Delta = 1.50\times10^{-5}$), because they are narrower.
- $n=21$ is comfortable for all; L1/L2 are overkill ($\Delta\sim10^{-7}$), L4 is the loosest ($7\times10^{-6}$).

## Corrected (weaker) gsq range -- L1 measured 2026-07-15

The gsq8 ensembles were found to be too strong; the corrected massless couplings are L1 {0.5,1.0,1.5},
L2 {1.0,2.0,3.0}, L4 {2.0,4.0,6.0}. **L1 Wilson D_W spectrum measured** on the fresh tuning configs
(24 configs each, `eig_wilson_lowmode_claude.cu`, `run_eig_L1_tune_claude.sh`; free U=0 ref (1.154, 10.82)):

| gsq | $\lambda_\text{min}$ range | $\lambda_\text{max}$ (max) | $k=\lambda_\text{min}/\lambda_\text{max}$ |
|-----|------|------|------|
| 0.5 | [0.775, 0.886] | 11.00 | ~0.070-0.081 |
| 1.0 | [0.577, 0.781] | 11.15 | ~0.052-0.070 |
| 1.5 | [0.383, 0.611] | 11.12 | ~0.034-0.055 |

Spectrum sits NEAR FREE (weak coupling): $\lambda_\text{min}$ HIGH (0.38-0.89, vs gsq8's 0.2-0.3),
$\lambda_\text{max}$ ~11 (below free's tail). **The frozen L1 window (0.1, 13) COVERS safely**
($\lambda_\text{min}\ge0.38\gg0.1$, $\lambda_\text{max}\le11.15<13$) -- zero "eval below window" warnings
across the HMC runs. It is over-conservative on the low edge (0.1 vs actual $\ge0.38$); a tighter window
(~0.3, 11.5) would give $k\approx0.026$ (vs 0.0077) -> fewer poles for the same $\Delta$. **Per NM: window
KEPT AS-IS for now.** L1 HMC dH is O(0.1), floor-free, acceptances ~0.74/0.91/0.88 (gsq 0.5/1.0/1.5).
L2/L4 spectra pending (their tuning runs in progress). See tuning sandbox `qed3/src/tuning/`.

## Where they live / how they are applied

- Single source of truth: `includes/frozen_window_claude.h` -> `frozen_window(L, lmin, lmax)`.
- Applied via `OverlapWMass::set_lambda(lmin, lmax)` (in `includes/overlap_wmass_claude.h`): sets both edges,
  rebuilds the Zolotarev fit on $k = \lambda_\text{min}/\lambda_\text{max}$, and sets `is_lambda_fixed = true`
  so `update()` never recomputes. Call once after construction, before the first `update()`.
- Wired: `test_hasenbusch_tune_l1_claude.cu`. **PENDING: production HMC drivers** (hmc_hasenbusch /
  hmc_fermilab_wmass massless L1/L2/L4).

## Provenance -- Wilson low-mode distribution

Chosen from `eig_wilson_lowmode_claude.cu` scans (massless gsq8 + gsq12, Nf2 + Nf6), which record per config
$\lambda_\text{min}$ = smallest singular value of $D_W$ (= $\sqrt{\min\text{ev}(D_W^\dagger D_W)}$; there is
NO Hermitian $H_W=\gamma_5 D_W$ for the 2-component Wilson operator in this setup), $\lambda_\text{max}$, and
the ratio. Rule: $\lambda_\text{min}$ set below the observed smallest $\lambda_\text{min}$, $\lambda_\text{max}$
above the observed largest, so the window covers the ensemble.

Measured (Nf2 gsq8 / Nf6 gsq8 / Nf6 gsq12), $\lambda_\text{min}$ [min, max] and min ratio:

| $L$ | $\lambda_\text{min}$ range (Nf2 g8) | min ratio (Nf2 g8) | note |
|-----|------|------|------|
| 1 | [0.205, 0.301], $\lambda_\text{max}$ up to 14.6 | $1.5\times10^{-2}$ | far from zero; not stiff |
| 2 | [0.049, 0.108], $\lambda_\text{max}$ up to 7.7  | $7.2\times10^{-3}$ | far from zero |
| 4 | [0.0099, 0.295], $\lambda_\text{max}$ up to 6.5 | $2.3\times10^{-3}$ | **low tail -> the stiff case** |

Only $L=4$ has a genuine near-zero tail (down to $\sim0.01$); $L=1,2$ sit far from zero. Surprise: gsq12
has *higher* $\lambda_\text{min}$ than gsq8 at every $L$ (i.e. less low-mode stiffness at stronger coupling)
-- so the strong-coupling "freezing" is likely topological / additive-mass-renormalization, not $D_W$
near-zero modes.

### Free-field reference (U = 0, gauge-independent, per $L$)

| $L$ | $\lambda_\text{min}^\text{free}$ | $\lambda_\text{max}^\text{free}$ |
|-----|------|------|
| 1 | 1.154 | 10.82 |
| 2 | 0.987 | 6.62 |
| 4 | 0.911 | 4.33 |

Gauge fluctuations push $\lambda_\text{min}$ down hard (free $\sim1 \to$ interacting $\sim0.2$ at L1, $\sim0.01$
tail at L4) while $\lambda_\text{max}$ barely moves (free vs mean within a few percent; tails +10-35%).

## Singular-value distribution plots (SV = $\lambda_\text{min}(D_W)$)

Per-L histograms of the smallest singular value of $D_W$ (fraction of configs; each panel on its own range).
Plot script: `wilson_lowmode_claude.gp` (PNG via `gnuplot -e "TERM='png'; NF=..; GSQ=.." wilson_lowmode_claude.gp`;
PDF default also has the $\lambda_\text{max}$ panels + the history page).

**Nf2, gsq8** (the reliable anchor):

![Wilson low-mode distribution Nf2 gsq8](wilson_lowmode_Nf2_gsq8_lmin.png)

**Nf6, gsq8:**

![Wilson low-mode distribution Nf6 gsq8](wilson_lowmode_Nf6_gsq8_lmin.png)

**Nf6, gsq12** (stronger coupling -- note $\lambda_\text{min}$ sits HIGHER than gsq8):

![Wilson low-mode distribution Nf6 gsq12](wilson_lowmode_Nf6_gsq12_lmin.png)

Read: L1 tight near 0.26, L2 near 0.09, **L4 the only broad low tail** (down to ~0.01) -> the stiff case and
where the smallest ladder rung ($c_1$) matters. ($\lambda_\text{max}$ panels: `wilson_lowmode_Nf*_gsq*_lmax.png`.)

## Coverage caveats (accepted per NM "looks good")

Worst-case measured extrema (Nf2 gsq8) slightly exceed the frozen edges:
- L1: $\lambda_\text{max}=13 <$ observed 14.6.
- L2: $\lambda_\text{min}=0.06 >$ observed 0.0485 (lowest config clipped).
- L4: $\lambda_\text{max}=5 <$ observed 6.5.

If production configs cross an edge the operator prints `# WARNING: eval below Zolotarev window`. Revisit if
that fires often (widen the edge / add poles).

## Hasenbusch ladders + step counts (per L) -- finalized 2026-07-13 (NM)

Production Hasenbusch mass ladders (entries [1..K] = **M_mass coefficients**, `set_mass_coeff`; frame 0 =
PHYSICAL target, `set_mass`) + per-stage MD step counts. Source of truth:
`includes/hasenbusch_ladder_claude.h` (`hasenbusch_ladder(L)` + `hasenbusch_steps(L)`).

Trajectory length $\tau$ is now PER-L (`hasenbusch_tau(L)`): **1.2 for L1/L2, 1.0 for L4** (no bump for L4 --
near-zero-mode stiffness). FINALIZED 2026-07-15 (NM):

| $L$ | stages | ladder | steps | $\tau$ | status |
|-----|--------|--------|-------|-----|--------|
| 1 | 2-stage | $\{0,\ 0.5\}$ | $\{2,2\}$ | 1.2 | **FINAL** |
| 2 | 2-stage | $\{0,\ 0.5\}$ | $\{3,3\}$ | 1.2 | **FINAL** ($\{3,3\}$ beat $\{2,4\}$: Cost 3.32e5 vs 3.66e5) |
| 4 | 4-stage | $\{0,\ 0.16,\ 0.32,\ 0.5\}$ | $\{2,2,2,4\}$ | 1.0 | **FINAL** (heavy frame sub-cycled x2) |

**TWO-OPERATOR SPLIT-POLE is now the standard** (`two_operator_force_impl_plan_claude.md`): ACTION op D uses
$n_\text{act}(L)=25/25/31$ (`hasenbusch_naction`) on the full frozen window (heatbath + accept/reject); FORCE op
Df uses $n_f=11$ (`hasenbusch_nforce`) on the NARROWED window $[2\lambda_\text{min},\lambda_\text{max}]$ (MD force
only). Exact by Metropolis. **ALWAYS the ML integrator `hmc_hasenbusch_ml_claude.h`** (FermionGroupLevel drives
Df; `run()` re-solves the action eta at final $U$ before $h_1$; `GaugeLevel` MG=20). Drivers:
`hmc_hasenbusch_claude.cu` (serial); NEW `hmc_hasenbusch_block_claude.cu` (Nf-BLOCK `-DNF`, all $N_f\ge2$).

Cost/tuning via `test_hasenbusch_npole_claude.cu` (SETUP COMPARISON). Findings: pole count is a WEAK speed knob
(CG iters npole-independent); window $\times2$ bites only at L4 (near-zero mode); the heavy frame carries the
larger force ($\sim18\times$ the massless). Full results `hasenbusch_tuning_results_claude.md`.

PENDING: re-validate the finalized config (L4 is a NEW 4-stage; $\tau$ now per-L) via
`test_hasenbusch_validate_claude.cu`; then production. NEXT FOCUS: CG cost = mixed precision (fp32
reliable-update) + deflation of the near-zero Wilson modes (the L4 0.67-acceptance stiffness).

## Related

- Ladder ($c_1, c_2$) strategy and the tuning: `hasenbusch_massless_impl_plan_claude.md`.
- Scan driver `eig_wilson_lowmode_claude.cu`; plot `wilson_lowmode_claude.gp`; $\Delta$ calc
  `zolotarev_delta_claude.cpp`.
