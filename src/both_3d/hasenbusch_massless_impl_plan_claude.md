# Hasenbusch mass preconditioning for the massless overlap HMC

## HANDOFF STATUS (2026-07-11) -- READ FIRST (for the LOCAL agent)

REMOTE (fnal) agent has finished **C1 code**; validation + C2 onward are open. This doc is the single
source of truth. Chunk order: **C1 (done, unvalidated) -> C2+C2b -> C3 -> C4 -> C5 -> C6(deferred)**.

**What is DONE (C1 code, on disk in `src/both_3d/`):**
- `includes/overlap_wmass_claude.h` -- NEW method `precalc_grad_bilinear_deviceAsyncLaunch_ms(U,
  d_bra, d_ket)` (just after `precalc_grad_deviceAsyncLaunch_ms`, ~line 778). Feeds an EXTERNAL bra
  (bare, no $(1{+}m_L)$ fold) + ket into the validated pole/COO path so `grad_deviceAsyncLaunch(link,
  U, d_ket)` returns Term B $=-2(C/\lambda_\text{max})\Re\langle\phi|K|\eta\rangle$. Nothing else in
  this file changed; existing force paths untouched.
- `test_hasenbusch_bilinear_l1_claude.cu` -- L=1 force-vs-FD gate: random $\phi\neq\eta$, analytic
  (`precalc_grad_bilinear` + `Force::compute`) vs central-diff of $S_B=2\Re\langle\phi|D_\text{ov}|
  \eta\rangle$ ($\phi,\eta$ FROZEN); plus a $\phi{=}\eta$ cross-check vs the standard massless force.
- `tmp_hb_bilinear_claude.sh` -- fnal sbatch handoff (affine, was qos=test -> NM set qos=normal, 20
  min) building via the `Makefile.fnal` nvcc line, tees `test_hasenbusch_bilinear_l1_claude.log`.

**IMMEDIATE NEXT (C1 gate -- do this first, LOCALLY):** build+run
`test_hasenbusch_bilinear_l1_claude.cu` (it already uses the LOCAL geometry paths `../../geometry/`
+ `../../geometry/data/`, default reference grad, no `-DGRAD_L1/2/4`). Expect: `phi==eta` cross-check
`~1e-9` (the real correctness gate on the plumbing) and the bra$\neq$ket rows `|grad-fd| ~1e-5`.
- If the bra$\neq$ket rows match: the overload gives $2\Re\langle\phi|K|\eta\rangle$ for a GENERIC bra
  (the key unknown -- the massive path only ever used bra $=(1{+}m_L)\eta$). Proceed to C2.
- If `grad \approx -fd` or off by a constant: NOT a kernel bug -- it just pins the SIGN/FACTOR the C2
  manager uses to sum Term A + Term B. Record the observed sign/factor and carry it into C2.
- If `phi==eta` cross-check FAILS (`\gg 1e-9`): real plumbing bug in the overload -- debug that first
  (most likely the bra-side `d_Ys[0]=X^\dagger\!\cdot`bra / `d_eta_bra` wiring, or a Z/Y block reuse).

**Then C2 (the substance):** NEW `includes/pseudofermion_hasenbusch_claude.h` -- a `HasenbuschPF`
manager holding the mass ladder $\{0.1,0.4\}$ (configurable vector), one $(\phi_i,\eta_i)$ per frame,
per-frame heatbath/action/force, forces SUMMED into the EXISTING 2-level `MinimumNorm2(Block)` outer
fermion kick (NO integrator rewrite -- Grid-canonical: all frames on the outer timescale, gauge
inner). Per frame: Term A $=$ `precalc_grad(D_i,\eta_i)` + `grad` (validated production path); Term B
$=$ `precalc_grad_bilinear(U,\phi_i,\eta_i)` + `grad`, with the sign/factor from the C1 gate. Wire
through `both/hmc.h` (`H()`, reject-restore loop over frames). See the "Per-frame pseudofermion" and
"C2" sections below for the exact heatbath ($\phi_i=D_{i+1}^{-\dagger}D_i^\dagger\xi_i$, check
$S_i=\xi_i^\dagger\xi_i$ at gen) and validation (reversibility, $dH\sim\tau^2$, split-vs-single-PF).

**Key settled decisions (do NOT re-litigate):** measure-weighted mass per frame (NOT scalar shift);
2-level integrator FIRST (3-level split = deferred C6, light OUTER $\varepsilon{\approx}1/4$ / heavy
inner $/12$ / gauge $/100$); ladder $\{0.1,0.4\}$; NO $|dH|$ guard; massless study only, serial Nf2
first then Nf4/6 block. Sources: Hasenbusch `hep-lat/0107019`; ordering per Torsiello/Fleming/Jin/
Osborn PoS LATTICE2024 052 + Grid `Test_hmc_Mobius2p1f.cc` (fermions outer/coarse, gauge inner/fine).

**Local vs fnal:** the test + C2 code use LOCAL relative geometry paths; build with the local
nvcc/Makefile as for `test_diag_mass_l1_claude.cu`. `tmp_hb_bilinear_claude.sh` is fnal-only (ignore
locally). The overload/manager headers are shared source (no path assumptions inside them).

## Goal / physics

The massless ($m=0$) strong-coupling study (Nf6 at gsq16/L2 and gsq12/L4) keeps hitting
near-zero-mode / topology-barrier **blowups**: a single molecular-dynamics trajectory produces a
wildly non-reversible $|\Delta H|\sim 400$, and because the spurious value is large-*negative* it is
Metropolis-*accepted*, locking a garbage config into the ensemble (the chain then rejects forever).
This is intrinsic to a single-pseudofermion light quark: the fermion force $\sim D^{-1}$ diverges as
$D_\text{ov}$ develops a near-zero eigenvalue.

**Hasenbusch mass preconditioning** (M. Hasenbusch, *Phys. Lett. B* **519** (2001) 177,
`hep-lat/0107019`) splits the light determinant across a ladder of heavier auxiliary masses. Each
factor carries a much smaller, better-conditioned force; the dangerous light contribution is
preconditioned by the next-heavier frame, so the integrator stays stable through zero-mode regions
and the $|\Delta H|$ spikes are strongly suppressed. This is the standard cure and needs no
measure-weighted mass -- the auxiliary shifts are plain scalar masses (NM).

## Algorithm (cite in code + here)

Reference: **M. Hasenbusch, `hep-lat/0107019`** (two-mass preconditioning), generalized to a
multi-mass ladder (cf. Hasenbusch--Jansen `hep-lat/0211042`).

### Operators (measure-weighted mass -- NM 2026-07-11, reuses the validated force)

Target massless overlap $D_0 \equiv D_\text{ov}$. Auxiliary ladder
$$
0 = m_0 < m_1 < \dots < m_K,\qquad D_i \equiv D_\text{ov} + m_{L,i},\quad
m_{L,i}=\text{mass\_coeff}(m_i)\cdot M_\text{mass},
$$
each frame a `OverlapWMass` at physical mass $m_i$ with the SAME measure-weighted diagonal
$m_L=\text{diag}(m\,A_y/\bar a_s)$ as the physical action (NOT a plain scalar shift $+m_i\mathbb{1}$).
Cleaner theoretically AND reuses the already-validated massive operator + force. Only $D_\text{ov}$ is
$U$-dependent ($dm_L/dU=0$), so $dD_i/dU = dD_\text{ov}/dU \equiv K$ (the conserved current) for every
frame. Supersedes the earlier scalar-shift plan (line 17).

### Determinant split (telescoping)

$$
\det(D_0^\dagger D_0)
=\left[\prod_{i=0}^{K-1}\frac{\det(D_i^\dagger D_i)}{\det(D_{i+1}^\dagger D_{i+1})}\right]
 \det(D_K^\dagger D_K).
$$
$K$ auxiliary masses $\Rightarrow$ $K$ ratio pseudofermions $+$ 1 heaviest-frame pseudofermion
$= K+1$ pseudofermion fields per $N_f/2$ stack.

### Per-frame pseudofermion (ket $\eta$, action $S$, heatbath, force)

**Heaviest frame** $i=K$ (standard PF for $D_K$):
$$
S_K=\phi_K^\dagger (D_K^\dagger D_K)^{-1}\phi_K,\qquad
\text{heatbath } \phi_K=D_K^\dagger\,\xi_K,\quad
\eta_K=(D_K^\dagger D_K)^{-1}\phi_K .
$$
Force $=-\,\eta_K^\dagger\big(dD_K^\dagger\,D_K + D_K^\dagger\,dD_K\big)\eta_K$ -- identical shape to the
current PF force, evaluated with operator $D_K$ and vector $\eta_K$ (reuse existing DHD force).

**Ratio frame** $i=0..K-1$ (light $D_i$ preconditioned by heavier $D_{i+1}$):
$$
S_i=\phi_i^\dagger\,D_{i+1}(D_i^\dagger D_i)^{-1}D_{i+1}^\dagger\,\phi_i .
$$
Let $\chi_i \equiv D_{i+1}^\dagger\phi_i$ and $\eta_i\equiv (D_i^\dagger D_i)^{-1}\chi_i$
(one DHD solve). Then $S_i=\chi_i^\dagger\eta_i=\phi_i^\dagger D_{i+1}\eta_i$.
Heatbath: draw $\xi_i$ Gaussian, set
$$
\phi_i = D_{i+1}^{-\dagger} D_i^\dagger\,\xi_i
\;\;\Longrightarrow\;\; S_i=\xi_i^\dagger\xi_i\ \text{ at generation}
$$
(one $D_i^\dagger$ apply + one solve of $D_{i+1}^\dagger$). Force has TWO pieces:
$$
dS_i=\underbrace{-\,\eta_i^\dagger\big(dD_i^\dagger D_i + D_i^\dagger dD_i\big)\eta_i}_{\text{A: reuse DHD force, op }D_i,\text{ vec }\eta_i}
\;+\;\underbrace{2\,\mathrm{Re}\,\phi_i^\dagger\,(dD_\text{ov})\,\eta_i}_{\text{B: NEW bilinear single-}D\text{ force, bra }\phi_i,\text{ ket }\eta_i}.
$$
Term A = the existing overlap DHD force with a different operator/vector. Term B = a new bilinear
overlap-derivative contraction ($dD_\text{ov}$ between two DISTINCT vectors) -- the one genuinely new
kernel (the current force always has bra $=$ ket).

Total fermion force per stack $=\sum_{i=0}^{K}$ (frame forces); the integrator sums them into one
gauge momentum kick.

## Current code (what changes)

- `includes/overlap_wmass_claude.h` -- `OverlapWMass` stores COMPLEX `mass` and applies the
  measure-weighted $m_L$ (`apply_mL`, `mult/adj/DHD` + `_ms`). The pre-diagonal scalar path
  `+mass*v` is present but commented (`:441`,`:476`). Force is `precalc_grad`/`grad` (pole loop) with
  bra $=$ ket $=\eta$. Frozen scalar reference lives in `overlap_wmass_obsolete_claude.h`.
- `includes/pseudofermion_claude.h` -- `PseudoFermion<Fermion>`: ONE `d_phi/d_eta`,
  `gen`=$D^\dagger\xi$, `update_eta`=$(D^\dagger D)^{-1}\phi$, `S`=$\phi^\dagger\eta$,
  `get_force`=`pi.compute(u,eta,D)`. (Nf2 path.)
- `includes/pseudofermion_Nfblocked_claude.h` -- `PseudoFermionBlock<Fermion,NSTACK>`: NSTACK columns,
  block heatbath/solve/S, `BlockedForce` (`-DBLOCK_FORCE`). (Nf4/6 path.)
- `both/hmc.h` -- `HMC` holds a SINGLE `pf`; `H()=pi^2/2 + Sg + pf->S()`; `run()` does `pf->gen`,
  integrate, and on reject restores `pf->d_eta` from `d_eta_saved`.
- `both/integrator.h` -- `integrate(U,pi,Sg,fermion,pf)` (leapfrog; per-step fermion force from `pf`).
- Drivers: `hmc_fermilab_wmass_L{2,4}_massless_claude.cu` (Nf2 serial),
  `hmc_Nfblocked_claude.cu` (Nf4/6 block).

## Files to modify

1. `includes/overlap_wmass_claude.h` (or a new `overlap_shift_claude.h`): a **scalar-shift** operator
   mode -- $D_\text{ov}+m_i\,\mathbb{1}$ with $m_i$ real scalar (reactivate the commented `+mass*v`
   path; scalar force = the obsolete reference). Must expose: `adj`, `mult`, DHD multishift solve, and
   the **new bilinear force** $2\Re\,\phi^\dagger dD_\text{ov}\,\eta$.
2. NEW `includes/pseudofermion_hasenbusch_claude.h` -- a `HasenbuschPF` manager holding the mass
   ladder $\{m_i\}$, one $(\phi_i,\eta_i)$ per frame, per-frame `gen`/`update_eta`/`S`/`get_force`,
   and the summed force. (Serial Nf2 first.)
3. `both/hmc.h`, `both/integrator.h` -- generalize the single `pf` to the frame-summed action/force
   (loop over frames in `H()`, `run()` reject-restore, and the integrator force).
4. Block: extend `pseudofermion_Nfblocked_claude.h` / `overlap_force_Nfblock_claude.h` to the ladder
   (frames $\times$ NSTACK) once the serial version validates.
5. Drivers: pass the ladder $\{m_i\}$ (CLI or compile-time), wire the Hasenbusch PF.

## Multilevel (multi-timescale) integrator + MDstep tuning

**START 2-LEVEL (NM 2026-07-11).** The full 3-level split below is a big rewrite; begin instead with the
**existing 2-level integrator UNCHANGED** -- put ALL 3 Hasenbusch frames on the ONE outer fermion
timescale (their forces summed into the single fermion kick), gauge on the inner timescale
(`nsteps_inner`), exactly like Grid's canonical 2+1f test. This needs **NO integrator rewrite** --
only the "fermion force" fed to `MinimumNorm2(Block)` becomes the frame-summed force. Get the split
correct + validated here FIRST; the 3-level multi-timescale split (next paragraph) is a DEFERRED
optimization (chunk C6) once the 2-level Hasenbusch runs and force norms are measured.

The production integrator `MinimumNorm2(Block)` (`includes/integrator*.h`) is already **2-level**
Omelyan ($\lambda=0.1931833275037836$): fermion force on the coarse OUTER timescale (`nsteps`), gauge
force nested on the fine INNER timescale (`nsteps_inner`), kick pattern
$\lambda,\,(1{-}2\lambda),\,\lambda$ at each level. Hasenbusch pays off ONLY with a matching
multi-timescale split: the (now small) light-frame force is evaluated FEW times on the coarsest
timescale, the cheap large gauge force MANY times on the finest.

### Structure -- recursive nested Omelyan MN2, one level per force group

Generalize the hardcoded 2-level MN2 into a recursive `MinimumNorm2ML` taking an ORDERED list of
(force-group, steps) from OUTERMOST (small force / expensive) to INNERMOST (large force / cheap):
$$
\underbrace{\text{light ratio }D_0/D_1}_{\text{outer, fewest steps}}
\supset \underbrace{\text{ratio }D_1/D_2}_{} \supset
\underbrace{\text{heavy }D_2}_{} \supset
\underbrace{\text{gauge } S_g}_{\text{inner, most steps}} .
$$
Each level runs the MN2 kick pattern with ITS force; its "drift" is a full sub-trajectory of the
next-finer level (instead of `U += tau*pi`). The innermost drift is the actual gauge `U += tau*pi`.
Palindromic/symmetric construction => reversibility preserved (the decisive gate on the rewrite). The
2-level case must reduce byte-for-byte to the current `MinimumNorm2Block` (regression check). Frames
with similar force norms MAY share a timescale (fewer levels); first try groups the two lighter frames
on one outer level -> **3 timescales: {light-pair} superset {heavy 0.4} superset {gauge}**.

**Concrete 3-level structure (NM 2026-07-11 corrected + Torsiello/Fleming/Jin/Osborn PoS LATTICE2024
052 + Grid `tests/hmc/Test_hmc_Mobius2p1f.cc`):** outer -> inner, with NM's seed step sizes
$$
\{\text{light-pair: ratios }D_0/D_1,\,D_1/D_2\}\ (\varepsilon\approx 1/4)
\ \supset\ \{\text{heavy }D_2=D_{0.4}\}\ (\varepsilon\approx 1/12)
\ \supset\ \{\text{gauge}\}\ (\varepsilon\approx 1/100).
$$
**Ordering (light = FEWEST evals / OUTERMOST; gauge = MOST / INNERMOST):** the lightest frame is the
most expensive solve ($N_\text{CG}\sim 1/m_s$, paper Eq. 5-6) and, once preconditioned, has the
SMALLEST + smoothest force -> it takes the COARSEST step ($\varepsilon\approx1/4$, only ~4 light solves
/traj = the cost win). The heavier/cheaper frame runs finer ($\varepsilon\approx1/12$), gauge finest
($\varepsilon\approx1/100$). Matches the paper's `G5F2,3 = AHALAHALAHA` (A x6 > H x3 > L x2) and Grid's
2+1f test (all fermions Level1=outer/coarse/fewest, gauge Level2=inner/fine/most).

In nested-multiplier terms: **light outer $=4$ steps ; heavy mult $\approx 3$ ($\to 12$ total) ; gauge
mult $\approx 8$ ($\to \approx 96$ total)** -- so the light force is the FEWEST-evaluated, gauge the
most. Total gauge link-updates $=\texttt{MDsteps}\times\prod_\ell\text{mult}[\ell]$ is a TARGET set by
the innermost multiplier (inserting a fermion level REDISTRIBUTES, does NOT multiply an old
`nsteps_inner` -- there is no "gauge cost blowup"). All per-level counts are runtime parameters;
tune from measured force norms (`-DInfoForce`).

NOTE (Grid caveat): Grid's canonical 2+1f test lumps ALL pseudofermions on ONE outer level and
separates only the gauge (2-level). NM opts for the finer 3-level split (light/heavy on separate
timescales, per the Torsiello paper's `AHALAHALAHA`). If the lightest ratio still spikes at
$\varepsilon\approx1/4$ (residual zero-mode stiffness), the lever is a CLOSER / EXTRA Hasenbusch mass
(e.g. add one below 0.1) to shrink its force -- NOT more light steps. Best integrator in the paper was
force-gradient (`ABACABA`); nested Omelyan MN2 is the baseline. Paper notes little benefit beyond 2
Hasenbusch masses at their $m=0.04$, but lighter ensembles (ours is massless) may benefit from more ->
$\{0.1,0.4\}$ is a sensible start.

### MDstep tuning (per level) -- measure, balance, iterate

Instrument already present: `-DInfoForce` prints each group's force norm per step
(`dSf/dSg.print2log_norm`). Procedure:
1. Thermalize briefly; with `-DInfoForce` measure the typical force norm of each group
   $\langle|F_k|\rangle$ (gauge, heavy 0.4, light frames) on the actual Nf6 gsq16/gsq12 ensembles.
2. **Balance the steps** so each level's integration error is comparable: finer timescale for larger
   force, $n_k \propto \langle|F_k|\rangle$ (gauge largest -> most steps; light frame smallest ->
   fewest). Fix $t_\text{max}\approx1.0$.
3. **Couple with mass tuning:** the auxiliary masses $\{0.1,0.4\}$ are chosen so ADJACENT frames have
   balanced force norms; if $\langle|F_\text{light}|\rangle \gg \langle|F_\text{mid}|\rangle$, shift a
   mass. Step-tuning and mass-tuning are done together from the measured force spectrum.
4. Iterate: candidate (masses, steps/level) -> short run -> watch per-level $dH$ contribution +
   acceptance + cost -> adjust. Gate every candidate with reversibility + $dH\sim\tau^2$ scaling.

First guess before measurement (tmax=1.0): gauge inner $\approx$100 (unchanged), heavy-0.4 level
$\approx$8, light-pair outer level $\approx$4; refine from the force norms.

## Ordered chunks

- **C1 -- external-bra force overload (NO new kernel; measure-weighted frames).** NM 2026-07-11: use
  the MEASURE-WEIGHTED `OverlapWMass(mass=m_i)` per Hasenbusch frame (cleaner theoretically; force
  already validated), NOT a scalar shift. Via the GW relation $D_\text{ov}^\dagger D_\text{ov}=
  D_\text{ov}+D_\text{ov}^\dagger$ the force is $-2\Re\langle(1{+}M)\eta|K\eta\rangle$ with
  $K=dD_\text{ov}/dU$ -- an $\langle\text{bra}|K|\text{ket}\rangle$ contraction that ALREADY takes an
  independent bra (`d_eta_bra`, `d_Ys[0]=X^\dagger\!\cdot`bra; `grad_diag_mass_force_bug_claude.md`
  Sec 5'). The ratio-frame Term B $=2\Re\langle\phi|K|\eta_i\rangle$ is that same contraction with the
  EXTERNAL bra $\phi$. So C1 = add a `precalc_grad_deviceAsyncLaunch(U, d_bra, d_ket)` overload that
  takes the bra as an argument (instead of computing $(1{+}M)\eta$), reusing the pole/COO path
  unchanged; FD-check $\langle\phi|K|\eta\rangle$ (bra $\neq$ ket) on L=1. Files:
  `includes/overlap_wmass_claude.h` (one overload), test `test_diag_mass_l1_claude.cu` pattern.
- **C2 -- serial HasenbuschPF (Nf2) on the EXISTING 2-level integrator.** $K=1$ (one auxiliary mass)
  first to validate heatbath ($S=\xi^\dagger\xi$ at gen), action, and $dH\sim$O(1) reversibility; then
  the $K=2$ ladder $\{0.1,0.4\}$. All frames' forces are SUMMED into the single outer fermion kick of
  the unchanged `MinimumNorm2` -- NO integrator rewrite. Files: NEW
  `includes/pseudofermion_hasenbusch_claude.h`, `both/hmc.h` (loop frames in `H()`/reject-restore).
- **C6 (DEFERRED) -- 3-level multi-timescale split + MDstep tuning.** Only after the 2-level Hasenbusch
  is validated. Recursive nested Omelyan MN2 (light-pair outer $\varepsilon\approx1/4$ superset heavy
  $\varepsilon\approx1/12$ superset gauge $\varepsilon\approx1/100$); 2-level reduction
  regression-checked; measure per-level force norms (`-DInfoForce`) and balance steps/masses. Files:
  `includes/integrator.h`, `includes/integrator_Nfblocked_claude.h`, `both/hmc.h`, drivers.
- **C3 -- validate Nf2 end-to-end.** Force-vs-FD on the full summed force; a short massless Nf2 run;
  confirm the blowup rate collapses (no accepted large-negative $dH$). Files: driver
  `hmc_fermilab_wmass_L{2,4}_massless_claude.cu`.
- **C4 -- port to the block (Nf4/6).** Ladder $\times$ NSTACK in the block PF/force. Files:
  `includes/pseudofermion_Nfblocked_claude.h`, `includes/overlap_force_Nfblock_claude.h`,
  `hmc_Nfblocked_claude.cu`.
- **C5 -- cluster rollout.** Rebuild the 6 massless binaries with the ladder; relaunch (fresh chains
  from the last CLEAN configs of the massless dirs). Files: `~/launch_massless_claude.sh`,
  `~/launch_block_build_claude.sh`, `run_wrapper_massless_claude.sh`.

## RESOLVED (NM, 2026-07-11)

1. **Mass ladder.** First try $\{m_1,m_2\}=\{0.1,0.4\}$ -> $K=2$ auxiliary masses, heaviest frame
   $m=0.4$, **3 pseudofermions** (frames $m=0,\,0.1,\,0.4$). Originally considered adding $1.0$ (option
   b) but that is overkill for the first try. **Code the ladder as a configurable vector of masses**
   (default $\{0.1,0.4\}$) so a $1.0$ frame can be added later without a rewrite.
2. **Scope + order.** Massless study only. Implement **serial Nf2 first**, then port to the **Nf4/6
   block**; both L2 and L4; same ladder everywhere on the first try. **Mimic the existing
   pseudofermion packing of Nf=2 (serial `PseudoFermion`) and Nf=4 (block `PseudoFermionBlock`)** --
   reuse those packing structures for the frames (block: frames $\times$ NSTACK columns). Massive
   campaigns are out of scope for now.
3. **No |dH| reject-guard.** NM: "it only obscures the situation." Rely on Hasenbusch alone.
