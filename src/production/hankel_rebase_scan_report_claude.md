# Block-Hankel + rebase knob scan for the $0^{++}$ two-$\sigma$ / two-meson GEVP

Ensemble `Nf2_gsq1.000000at0.200000...nt128L1_hb1.000000`, `NVDIR=distill_Nv24`, 400 configs,
connected vacuum-free basis (`CONN=1`), `SPLIT=1`, jackknife `BINSIZE=10` (autocorrelation).
Reference scales: $m_{PS}\approx0.322$, $2m_{PS}\approx0.644$.

Physics target: $m_0$ (one-meson-rich ground) and $m_1$ (two-meson $\approx 2m_{PS}$) with a **FLAT $m_1$**
plateau (a sloped $m_1$ disqualifies a choice), reached early, at balanced S/N, both states resolvable.

## How this was run (economy)

The expensive per-config correlator store is computed ONCE and cached; every knob combination is then
evaluated in-process (cheap linear algebra). Files (all new, nothing existing was overwritten):

- `hankel_rebase_scan_claude.py` — store builder + evaluation engine. **New code paths** vs the driver:
  - `hankel_off(Cts, offsets)` — generalized block-Hankel with an EXPLICIT offset list,
    $\hat C(t)_{(a,i),(b,j)}=C(t+\text{off}[a]+\text{off}[b])$; supports non-uniform ladders
    (Dt=1,2,3 $\to$ offsets $[0,1,2,3]$; Dt=1..5 $\to$ $[0,1,2,3,4,5]$). $t_\max=$ twin$-2\max(\text{off})$.
  - `staged_project` — **multiple (staged) rebase**: composes per-stage rebase projections
    $V_\text{tot}=V_1V_2\cdots$ (algebraically identical to successive project-then-rebase).
  - `rebased_effmass_fixed` — pads unresolved levels with NaN so rank-deficient jackknife samples
    (e.g. a fine ladder on a rank-2 base) keep a homogeneous array shape.
- `hankel_rebase_runscan_claude.py` — the controlled one-knob-at-a-time scan (prints the tables below).
- `hankel_rebase_plots_claude.py` — writes the shortlist / reject figures to `figs/`.

$m_0,m_1$ below are jackknife window-averages; "slope" is the jackknife least-squares slope of the
effmass over the $m_1$ window $t\in[6,12]$ (a choice TILTS if $|{\rm slope}|>2\sigma$ and $>0.008$).

---

## Shortlist (ranked)

| rank | config (env) | $m_0$ | $m_1$ | $m_1$ onset | quality note |
|---|---|---|---|---|---|
| **1** | `NOOA=1 HANKEL=1 NSH=3 SHIFT=1 REBT=5 NKEEP=2 T0=3` (Dt=1,2) | 0.470(10) | **0.644(3)** | $t\simeq7$ | flat $m_1$ right at $2m_{PS}$, slope $-0.003(0)$; best all-around, range to $t\simeq20$ |
| **2** | `NOOA=1 HANKEL=1 NSH=2 SHIFT=3 REBT=5 NKEEP=2 T0=3` (Dt=3) | 0.479(9) | 0.648(3) | $t\simeq8$ | 2-block, cheapest; flat, slope $-0.004(0)$; slightly tighter $m_0$ at mid-$t$ |
| **3** | `NOOA=1 HANKEL=1 NSH=3 SHIFT=1 REBT=5 NKEEP=2 T0=2` (Dt=1,2) | 0.483(9) | 0.649(3) | $t\simeq6$ | earliest onset, tightest $m_1$ S/N; slope $-0.004(0)$ still flat |
| **4** | `NOOA=1 HANKEL=1 NSH=3 SHIFT=2 REBT=5 NKEEP=2 T0=3` (Dt=2,4) | 0.430(12) | 0.633(3) | $t\simeq4$ | fastest $m_1$ plateau, but $m_1$ pulled slightly UNDER $2m_{PS}$; range only to $t\simeq13$ |
| **5** | `NOOA=0 HANKEL=1 NSH=3 SHIFT=1 REBT=5 NKEEP=3 T0=3` (Dt=1,2, 3-op) | 0.438(22) | 0.637(3)$^\dagger$ | $t\simeq6$ | resolves THREE levels; two-meson is **level 2** ($0.637$, clean flat); use only if the intermediate $O_A$ state (level 1, $\sim0.48$, tilts) is wanted |

$^\dagger$ For PICK 5 the clean two-meson is `m2`, not `m1` — see the 2-op-vs-3-op trend below.

Full env (prepend always): `OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 ENS=... NVDIR=distill_Nv24 CONN=1 BINSIZE=10`.

**Recommendation.** PICK 1 (Dt=1,2, 2-op, $T_0=3$) is the reference: clean flat two-meson exactly at
$2m_{PS}$, widest usable range, robust. PICK 2 (Dt=3) is an equally good cross-check with fewer Hankel
blocks. If early onset / tightest two-meson error matters, PICK 3 ($T_0=2$). Use PICK 4 (Dt=2,4) only when
the earliest possible $m_1$ plateau is the priority and the mild $m_1$ undershoot and short range are
acceptable. Reach for the 3-op PICK 5 only to expose the intermediate $O_A$-rich state; the two-meson
itself is cleaner and correctly labelled in any 2-op pick.

---

## Rejected examples (for contrast)

| config | why rejected |
|---|---|
| `NOOA=1 NSH=2 SHIFT=1 T0=3` (**Dt=1**) | $m_1$ **TILTS**: slope $-0.0097(3)$, $m_1$ still falling (0.77 -> 0.62). One shift under-cleans; contamination survives in state 1. |
| `NOOA=1 NSH=3 SHIFT=4 T0=3` (**Dt=4,8**) | over-cleaned: $m_1$ flat but pulled to **0.612** ($<2m_{PS}$), and $m_0$ has **no plateau** (dips to 0.28 then rises); range only $t\le13$. |
| `NOOA=1 offsets=[0..5]` (**Dt=1..5, 2-op**) | numerically degenerate: a 6-block ladder on a rank-2 base is massively over-complete; the reduced metric is singular and the effmass is garbage (NaN/wild). Fine ladders need a >=3-rank base. |
| `NOOA=0 NKEEP=2` (**3-op, keep only 2**) | mislabels: with $O_A$ in the basis the leading-2 rebased states are {ground, tilted $O_A$-state $\sim0.48$}; the two-meson is pushed to level 2, so `m1` here TILTS ($-0.013$) and is NOT the two-meson. Must keep `NKEEP=3`. |

Figures for PICK 1-5 and the first two rejects are in `figs/hankel_scan_*_claude.png`.

---

## Per-knob TREND (what each parameter does)

**`SHIFT`/`NSH` = the Hankel shape (dominant knob).** Larger time-shifts remove contamination faster
(earlier $m_1$ plateau) but eat $t$-range and pull $m_1$ below $2m_{PS}$; too small under-cleans and $m_1$
tilts.
- Dt=1 $[0,1]$: under-cleaned, $m_1$ TILTS ($-0.0097$), $m_1=0.672$. REJECT.
- Dt=1,2 $[0,1,2]$: sweet spot, $m_1=0.644$ flat, $t_\max=27$.
- Dt=3 $[0,3]$: equivalent sweet spot, $m_1=0.648$ flat, 2 blocks only.
- Dt=2,4 $[0,2,4]$: earliest plateau ($t\sim4$) but $m_1=0.633$ (undershoot), $t_\max=23$.
- Dt=4,8 $[0,4,8]$: over-cleaned, $m_1=0.612$, $m_0$ no plateau, $t_\max=15$. REJECT.
- Fine ladders (Dt=1,2,3 or 1..5) on the **2-op** base add blocks WITHOUT new rank (base rank 2), so they
  only add noise (Dt=1,2,3: $m_0$ error doubles) or go singular (Dt=1..5). Only worth it on a $\ge3$-op base.

**`T0` (GEVP metric).** Weak, monotone: $T_0=2\to3\to4$ lowers $m_0$ ($0.483\to0.470\to0.461$) and $m_1$
($0.649\to0.644\to0.642$, toward $2m_{PS}$) and **flattens** the residual $m_1$ slope
($-0.0040\to-0.0030\to-0.0024$), but pushes the onset later and grows the errors. $T_0=3$ is the balance;
$T_0=2$ if you want the earliest onset and tightest error.

**`NKEEP` / 2-op vs 3-op base.** This changes WHICH state is `m1`.
- 2-op $\{\sigma^2,O_{2\sigma}\}$: `m0`=ground $\approx0.47$, `m1`=two-meson $\approx0.644$ (clean). This is
  the natural basis for the two-meson.
- 3-op (add $O_A$): a new intermediate state $\sim0.48$ appears; with `NKEEP=3` the levels are
  $m_0\approx0.44$, $m_1\approx0.48$ ($O_A$-rich, **tilts/noisy**), $m_2\approx0.637$ (two-meson, clean flat).
  `NKEEP=2` on the 3-op base is a trap — it keeps the ground + the tilted $O_A$ state and drops the
  two-meson.

**`REBT` (rebase time).** Essentially irrelevant. For the **2-op base the $4\times4$ (or $6\times6$) Hankel
is rank 2**, so the leading-2 subspace is fixed regardless of $t_\text{reb}$ (documented, re-confirmed).
For the 3-op base, $t_\text{reb}=4\to9$ moves the tilted $m_1$ only mildly ($0.483\to0.471$) and leaves
$m_0$ unchanged; the two-meson (level 2) is stable throughout.

**Multiple (staged) rebase.** Implemented (`staged_project`, e.g. rebase $6\to$ then $\to2$ states) and
tested on the 3-op fine ladders. **No gain here** — the physical content is only $\sim3$ states, so the
extra Hankel blocks the staging feeds on are largely redundant and the staged $m_1$ comes out noisy
(errors $0.02$–$0.3$) with no better plateau than a single rebase. Reported as a negative result: staging
would only pay off with a genuinely larger operator basis.

## Caveat on $m_0$

In every pick the ground `m0` (green) does not fully plateau — it drifts slowly downward
(slope $\sim-0.01$) before the late-$t$ noise, i.e. it still carries some excited contamination and is the
weaker of the two determinations. The window-averaged $m_0\approx0.47(1)$ is a fair central value but the
two-meson `m1` is the robust, flat quantity. This matches the established baseline
($m_0\approx0.46(1)$, $m_1\approx0.644(3)$).
