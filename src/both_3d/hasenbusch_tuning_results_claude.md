# Hasenbusch split-pole HMC -- logged tuning numbers (for reference)

Recorded 2026-07-14 (NM). Production cost/acceptance numbers from `test_hasenbusch_npole_claude.cu`
(two-operator split-pole overlap HMC). See `two_operator_force_impl_plan_claude.md` for the method and
`frozen_window_claude.md` / `hasenbusch_ladder_claude.h` for the frozen parameters.

## Finalized config (as of this recording)

- **Action op** $D$: $n_\text{act}(L)=25$ (L1/L2), $31$ (L4), full frozen window; heatbath + accept/reject.
- **Force op** $D_f$: $n_f=11$ poles on the NARROWED window $[2\lambda_\text{min},\lambda_\text{max}]$; MD force only.
- **Ladders / steps** (`hasenbusch_ladder`/`_steps`), $\tau=1.2$ (bumped from 1.0, all L, NM 2026-07-15),
  gauge sub-cycle $\times 20$:
  - L1: $\{0,0.5\}$ steps $\{2,2\}$
  - L2: $\{0,0.5\}$ steps $\{3,3\}$
  - L4: $\{0,0.15,0.3,0.5\}$ steps $\{3,3,3,3\}$ (4-stage, FINALIZED 2026-07-15)
- **Delta** (sign-approx error): action $\sim10^{-8}$ (all L); force $\sim2\times10^{-4}$ (L1/L2), $2\times10^{-3}$ (L4).

## L1 and L2 production numbers (gsq8, Nf2)

Measured with action $n=31$ (these runs predate the $n_\text{act}=25$ change for L1/L2; the force $n_f=11$
part -- which dominates -- is unchanged, so Cost is ~unaffected; the "action re-solve" CG would drop slightly
at $n=25$). 10-trajectory chains, fresh momentum + phi per trajectory.

### Per-frame force (one force eval on the loaded config)

| L | frame | $c$ (M_mass coeff) | CG iters | sec | $\lVert F\rVert_2$ | $\lVert F\rVert_\infty$ |
|---|-------|------|------|------|------|------|
| 1 | 0 | 0   | 10446 | 3.82 | 0.550 | 0.057 |
| 1 | 1 | 0.5 | 9952  | 2.36 | 9.899 | 1.138 |
| 2 | 0 | 0   | 15722 | 15.41 | 1.104 | 0.054 |
| 2 | 1 | 0.5 | 15051 | 9.14  | 18.93 | 1.083 |

### 10-trajectory chain (finalized force, config k = 3426 (L1) / 1601 (L2))

| L | ladder / steps | $\langle|dH|\rangle$ | acceptance | $\langle P\rangle$ | force CG/traj | total CG/traj (action re-solve) | sec/traj | Osborn Cost |
|---|------|------|------|------|------|------|------|------|
| 1 | $\{0,0.5\}\{2,2\}$ | 0.0200 | 1.00 | 0.9967 | $1.34\times10^5$ | $1.76\times10^5$ ($4.15\times10^4$) | 48.5 | $1.35\times10^5$ |
| 2 | $\{0,0.5\}\{3,3\}$ | 0.0180 | 1.00 | 0.9949 | $3.31\times10^5$ | $4.14\times10^5$ ($8.35\times10^4$) | 248 | $3.32\times10^5$ |

L1 all 10 traj accepted ($dH\in[-0.039,+0.025]$); L2 all 10 accepted. Both saturate acceptance at 1.0 with
$\langle|dH|\rangle\sim0.02$ -- **large headroom** (over-integrated).

## Key observations

- **Force imbalance (drives the step tuning):** the HEAVY frame carries the larger force -- L1 frame-1
  ($c=0.5$) $\lVert F\rVert_2=9.9$ vs frame-0 (massless) $0.55$ ($\sim18\times$); L2 $18.9$ vs $1.1$. So the
  heavy frame wants finer stepping; the small-force *and expensive* massless frame can be coarse.
- **Cost-optimal step ratio** $n_i\propto(\lVert F_i\rVert^2/c_i)^{1/3}$ gives $n_\text{massless}/n_\text{heavy}\approx0.12$
  at L1 -- i.e. the heavy frame "should" have $\sim8\times$ more steps. Current $\{2,2\}$ is far from that.
- **$\{1,2\}$ was considered (would move the right direction) but NOT adopted** -- NM not comfortable with a
  single outer MD step ($MD_\text{steps}=1$ = one coarse MN2 step over the whole $\tau$; dH-tail risk on
  stiff configs). **Keep L1 at $\{2,2\}$.** (Heavier sub-cycling of the heavy frame, e.g. $\{2,4\}$, stays
  on the table since it keeps $MD_\text{steps}\ge2$.)
- **Two-operator split at L1/L2 is ~a wash** (pole count is npole-independent for CG; window$\times2$ only
  ~6% at L1 since the physical low mode dominates the conditioning there). The split genuinely pays at L4.
- **L4 is the expensive case:** ~1 h/traj (frame-0 CG $\sim8.6\times10^4$, sec $\sim110$). Multiprecision
  ($\sim2\times$, reversibility-clean) is the planned next speed lever.

## L2 $\{3,3\}$ vs $\{2,4\}$ comparison (`tmp_hb_npole_L2cmp_claude.sh`, gsq8 Nf2, config k=1601)

**DONE 2026-07-14 -- $\{3,3\}$ CHOSEN** (better on both $\langle|dH|\rangle$ and Osborn Cost). 10 traj each.

| steps | $\langle|dH|\rangle$ | acceptance | $\langle P\rangle$ | force CG/traj | sec/traj | Osborn Cost |
|-------|------|------|------|------|------|------|
| **$\{3,3\}$ (chosen)** | **0.0180** | 1.00 | 0.9949 | $3.31\times10^5$ | 229 | **$3.32\times10^5$** |
| $\{2,4\}$ | 0.0233 | 1.00 | 0.9874 | $3.62\times10^5$ | 211 | $3.66\times10^5$ |

$\{2,4\}$ is ~8% faster per traj (211 vs 229 s) but does MORE total force work ($3.62$ vs $3.31\times10^5$ CG:
the heavy frame sub-cycled $\times2$ means more evals of the ~equally-expensive heavy solve) and has larger dH,
so its Osborn Cost is ~10% worse. Naively the force imbalance (heavy $\lVert F\rVert$ $\sim17\times$ the massless)
favors sub-cycling the heavy frame more, but dropping the massless frame $3\to2$ costs more dH than the heavy
$3\to4$ saves -- the small-norm massless frame's dH contribution is not negligible. **L2 = $\{0,0.5\}\{3,3\}$.**

## L4 finalized + low-mode observation (2026-07-15)

**L4 = $\{0,0.15,0.3,0.5\}$ steps $\{3,3,3,3\}$, $\tau=1.2$** (4-stage; chosen after the $\{2,2,2,4\}$ /
$\{2,2,6\}$ / $\{3,3,3\}$ comparisons). The L4 chain log is informative: on at least one trajectory the MD
**hit the near-zero Wilson mode** (the L4 $D_W$ tail, $\lambda_\text{min}\sim0.008$), giving a stiff
trajectory with **acceptance $\approx0.67$** on that traj -- the expected L4 behavior (the frozen-window doc
already flagged L4 as the only near-zero-tail case). The 4-stage ladder + $\tau=1.2$ manage it, but L4 stays
the expensive/stiff L.

**Future CG-cost work (L4):** the cost is dominated by the inner $D_\text{ov}$ multishift iteration count
(the near-zero mode). Two reversibility-clean levers to cut it:
- **Mixed precision** (fp32 Krylov + fp64 reliable-update): reach the same tol faster, no reversibility tax
  (Clark et al. arXiv:0911.3191). From-scratch fp32 path (see `two_operator_force_impl_plan_claude.md`).
- **Deflation** of the near-zero Wilson low modes: project them out of the inner solve so the CG condition
  number drops. Directly targets the L4 stiffness that caused the 0.67-acceptance trajectory.

## Pending / running comparisons (results to append when done)
- **L4 restart** -- `tmp_hb_npole_L4v2_claude.sh` (GPU1), 4 traj each:
  $\{0,0.2,0.5\}\{3,3,3\}$, $\{0,0.1,0.3,0.5\}\{2,2,2,4\}$, $\{0,0.2,0.5\}\{2,2,6\}$.
