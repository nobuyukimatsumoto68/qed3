# Hasenbusch mass preconditioning for the massless overlap HMC

## *** TOP PRIORITY (2026-07-12): lambda_max NON-EXACTNESS BUG FOUND + PARTIALLY FIXED -- READ FIRST ***

**What happened:** while validating Hasenbusch we found the overlap HMC was NOT the exact integrator
of its (Zolotarev-approximate) action. Symptom: a step-INDEPENDENT, tolerance-INDEPENDENT dH FLOOR
(dH stops shrinking as nsteps->inf). At coarse L1 the floor is ~delta (~1.5e-5); tiny but nonzero at
L2/L4. Hasenbusch only EXPOSED it (preconditioning shrinks dH right onto the floor at L1); it is NOT
Hasenbusch-specific.

**ROOT CAUSE (confirmed):** `OverlapWMass::update(U)` recomputed `lambda_max` (the spectral
normalization) EVERY MD step (`overlap_wmass_claude.h:303`), so D_ov silently changed step-to-step and
the FORCE omits `d lambda_max/dU`. The 2026-06-26 change froze only the Zolotarev WINDOW k, NOT
lambda_max. This is a BASE-OPERATOR bug in `overlap_wmass_claude.h` -> affects ALL production overlap
HMC (standard massless/massive `hmc_fermilab_wmass_*`, block `hmc_Nfblocked`, AND the new Hasenbusch).
Reversibility was always fine (lambda_max(U) deterministic), which is why it hid; only exactness broke.

**How it was proven (diagnostic chain, for the record):** force-vs-FD per term with an eps SWEEP
(1e-1..1e-4) to separate a real residual (flat in eps) from the CG/FD floor (~1/eps); mimic the
VALIDATED `test_diag_mass_l1_claude.cu:329` methodology (real PseudoFermion, fixed phi, eta re-solved,
FD of the ACTION S=phi^dag(D^dagD)^-1 phi, |grad-fd|); then FREEZE lambda_max inside the FD -> residual
collapsed to machine for BOTH the standard massless force AND the ratio-frame Term B. => the "Term B
bra!=ket symmetrization" worry was a RED HERRING; lambda_max was the whole story. Standard massless
force is EXACT once lambda_max is frozen. Test: `test_hasenbusch_termAB_deriv_claude.cu` (+ `.sh`).

**FIX APPLIED (operator side, ALL headers) -- FREE-THEN-FREEZE (NM 2026-07-12):** leave lambda_max FREE
(recompute) for the first few THERMALIZATION trajectories, then FREEZE forever. Each overlap header got:
member `bool is_lambda_fixed=false`, method `void freeze_lambda(){ is_lambda_fixed=true; }`, and
`update()` now does `if(!is_lambda_fixed) compute_lambda_max();` (old unconditional call commented).
Applied to ALL FIVE overlap headers (NM: "free lambda everywhere"):
`both_3d/includes/{overlap_wmass_claude.h, overlap.h, overlap_wmass_obsolete_claude.h}`,
`dual/overlap.h`, `both/overlap.h`.

**PLAN TO FINISH (NM decisions 2026-07-13):**
- **TUNE THE OVERLAP PARAMS ON THE EXISTING (fairly THERMALIZED) ENSEMBLES.** We already have thermalized
  configs, so the overlap parameters (the frozen lambda_max value, npole/window, ladder) are tuned using
  them. The `D.freeze_lambda()` driver calls are INSERTED AS PART OF THIS TUNING PROCESS -- i.e. the
  freeze is set from / on the thermalized configs, not from a cold start. (The operator mechanism is in
  place but inert until the tuning wires the call.)
- (2) RESOLVED: freezing lambda (both lambda_max AND lambda_min) is REQUIRED for reversibility anyway
  (a recomputed-lambda-driven refit is a reversibility hole -- same reason the window k was frozen
  2026-06-26). So freezing lambda_min is CORRECT, and the passive spike-warning going stale is accepted.
- (3) RESOLVED: the margin/drift concern (frozen lambda_max must stay an upper bound) is addressed with
  the EXISTING CONFIGS too -- survey lambda_max across the thermalized ensembles and set the frozen value
  (with margin) so it bounds them.
- (4) SCOPE unchanged: shared operator -> all drivers; existing ensembles were non-exact (negligible L2/L4).
- (5) RE-VALIDATE after the tuning wires freeze_lambda: rerun L1 HMC scaling (`tmp_hb_hmc_local_claude.sh`)
  -- standard PF (sec 4) AND Hasenbusch (sec 2) scale cleanly (no floor); confirm C1/C2/C3 still green;
  then C4 block + C5 rollout.

**Is Hasenbusch itself done+validated?** YES for the CODE and its gates AT PRODUCTION SIZE (L2): C1
bilinear gate PASS; C2 force+heatbath PASS (force-vs-FD 1e-6, heatbath exact); C2c/C3 HMC at L2 =
reversibility 8e-12, clean dH~tau^2 (order_p 2.0 at tmax=0.1), 6/6 accept. The "drift to debugging" was
uncovering the lambda_max base-operator NON-EXACTNESS (a SEPARATE, broader bug that Hasenbusch merely
EXPOSED at coarse L1 -- the standard massless force has it identically). So: Hasenbusch = code-complete +
validated at L2; the ONLY thing between here and "fully validated end-to-end" is finalizing the
lambda_max freeze (tune+wire) and one confirming re-run showing L1 also scales cleanly.

Everything below (Hasenbusch C1/C2/C2c) is DONE + validated EXCEPT it inherits the lambda_max floor
until the fix above is finalized and re-validated.

## HANDOFF STATUS (2026-07-11) -- READ FIRST (for the LOCAL agent)

REMOTE (fnal) agent finished **C1 code**; **C1 GATE VALIDATED LOCALLY 2026-07-11** (see below); C2
onward are open. This doc is the single source of truth. Chunk order: **C1 (DONE+VALIDATED) -> C2+C2b
-> C3 -> C4 -> C5 -> C6(deferred)**.

**C1 GATE RESULT (2026-07-11, local, `test_hasenbusch_bilinear_l1_claude.log`):** PASS in substance.
Ran `tmp_hb_bilinear_local_claude.sh` (L1_default + L1_gradl4). Findings:
- bra != ket rows: analytic `grad = -fd` to the FD truncation floor (`|grad+fd| ~ 5e-6` at `eps=1e-5`)
  on ALL sp/tp links -> the overload returns the FORCE `-2(C/lambda_max) Re<phi|K|eta>` =
  `-d/dU[2Re<phi|D_ov|eta>]`, i.e. the standard force = -dS/dU convention (NOT the +dS/dU derivative).
- `phi==eta` cross-check == production massless force at `1.4e-16` (machine) -> plumbing exactly
  reproduces the validated path for a GENERIC bra. Same result under `-DGRAD_L4` (production block grad).
- The test binary reports "FAIL" only because it compared `|an-fd|` (it expected grad==+fd); the true
  relation is grad==-fd. This is the anticipated sign-pin, not a kernel bug.

**SIGN/FACTOR (pinned by the C2 FD gate, 2026-07-11 -- SUPERSEDES the earlier "added directly"):**
Both grad routines return the ACTION GRADIENT +dS/dU (what the integrator feeds to `pi += -tau*dSf`;
production convention). Term A `precalc_grad(D_i,eta_i)+compute = +dS_A/dU` (the inverse variation, with
the extra minus of dM^{-1}=-M^{-1}dM M^{-1} already baked in). BUT Term B's kernel
`precalc_grad_bilinear(phi_i,eta_i)+compute = -d[2Re<phi|D_ov|eta>] = -dS_B/dU`, whereas the numerator
variation is dS_B/dU = +2Re<phi|dD_ov|eta> = +d[2Re<phi|D_ov|eta>] (NO inverse -> NO extra minus). So
Term B enters `get_force` with the OPPOSITE relative sign and must be **SUBTRACTED**:
  get_force = sum_i [ precalc_grad(D_i,eta_i)+compute ]  -  sum_{ratio i} [ precalc_grad_bilinear(phi_i,eta_i)+compute ].
C1 only checked the kernel vs a non-action scalar (grad_bilinear == -d[scalar]); the ASSEMBLY sign was
always flagged as "fixed in the manager + L=1 FD gate" -- this is that fix. Implemented as `dSf += -1.0*ftmp`.
NB the C2 FD gate compares |grad - fd| (analytic = +dS/dU, FD of the actual action = +dS/dU), with a
TIGHT solve tol (1e-12) so the re-solved near-singular massless inverse does not swamp the FD.

**C2 CODE WRITTEN (2026-07-11, LOCAL) -- awaiting the force gate.** Decisions locked: measure-weighted
mass, ONE shared operator + `set_mass`, ladder `{0,0.1,0.4}`. Files:
- `includes/overlap_wmass_claude.h:290` -- NEW `set_mass(Complex m)` (rewrites mass/mass_coeff,
  is_precalc=false; update(U) is mass-free so no re-update needed). [C2a DONE]
- `includes/pseudofermion_hasenbusch_claude.h` -- NEW `HasenbuschPF<Fermion,Force>` manager: ladder,
  per-frame (phi_i,eta_i,chi_i), shared op + set_mass; `gen` (heatbath phi_i=D_{i+1}^{-dag}D_i^dag xi
  via DHD-solve+mult both at m_{i+1}; heavy phi_K=D_K^dag xi), `update_eta`, `S`=sum chi_i^dag eta_i,
  `get_force`=sum(Term A precalc_grad + Term B precalc_grad_bilinear for ratio frames). Heatbath
  identity S_i=xi_i^dag xi_i verified analytically for ALL frames (ratio + heavy). xi_sqnorm instrumented.
- `includes/hmc_hasenbusch_claude.h` -- NEW `MinimumNorm2Hasenbusch` (copy of MinimumNorm2, DROPS the
  external `fermion->precalc_grad(U,pf->d_eta)` -- get_force is self-contained/multi-frame) + `HMCHasenbusch`
  (copy of HMC2, Hasenbusch heatbath, no post-gen external precalc). Originals untouched.
- `test_hasenbusch_force_l1_claude.cu` + `tmp_hb_force_local_claude.sh` -- C2 gate: heatbath identity
  (rel ~1e-7) + full summed force vs central-FD of S=sum chi_i^dag eta_i (phi frozen, eta re-solved;
  grad==-fd, |grad+fd|~1e-5). Ladders K=1 {0,0.1} + K=2 {0,0.1,0.4}. NEXT = build+run this gate.
**C2 GATE PASS (2026-07-11, test_hasenbusch_force_l1_claude.log).** Heatbath identity
S_i=xi_i^dag xi_i validated (mean over M=64 draws): mean_rel ~2-6e-7 ~= per-draw rms_rel, BOTH below
delta=3.3e-6 -> the residual is the SYSTEMATIC Zolotarev/GW-shortcut floor (does NOT average down, as
predicted), NOT a bug; gate on the delta-level averaged identity (<1e-5). Full summed force vs FD
matches ~1e-6 for K=1 {0,0.1} and K=2 {0,0.1,0.4} (after the Term B SUBTRACT sign fix + tight-tol
de-noised FD comparing |grad-fd|).

**dH FLOOR + GW SQUARE (NM 2026-07-11, DECISION = proceed, floor benign).** The L1 dH has a
step-INDEPENDENT, TOLERANCE-INDEPENDENT floor F ~ 3 delta (~5e-5) that the naive dH(n)/dH(2n) ratio
sees as a "failure" at fine steps; L2 (dH >> floor) scales cleanly (~3.8). The overlap `D^dag D` is
STILL computed via the GW shortcut D^dag D = D + D^dag in `OverlapWMass::DHD_deviceAsyncLaunch(_ms)`
(overlap_wmass_claude.h:558/606; comment :561), reached by every action solve/heatbath/eta-update
(LinOpDHDWrapper -> Op_DHD, both PseudoFermion + HasenbuschPF). An honest adj(mult) square
(`sq_deviceAsyncLaunch`) does NOT exist for OverlapWMass (grep empty; the commented wrapper refs point
to a removed method). BUT the overlap FORCE is ALSO GW-derived (standard force = 2Re<eta|K|eta> =
grad of 2Re<eta|D_ov|eta>, which uses GW), so action + force are GW-CONSISTENT and GW is NOT the floor
source (my earlier guess retracted). Floor origin left open; it is benign (acceptance ~1, reversibility
exact 1e-12, L2 tau^2 clean, force-vs-FD 1e-6). NM: keep the GW square, proceed to production.
Scaling scan set to tmax_scan=0.1 nsteps 4..16 (may be floor-dominated -> read as informational).

**C2c CODE WRITTEN (2026-07-11) -- HMC GATE effectively PASSED (L2 reversibility 8e-12, dH~tau^2 ~3.8,
6/6 accept; L1 reversibility 2e-12, floor-limited scaling = benign).** `test_hasenbusch_hmc_l2_claude.cu` +
`tmp_hb_hmc_local_claude.sh`: end-to-end validation with the REAL `MinimumNorm2Hasenbusch` +
`HMCHasenbusch` (massless, ladder {0,0.1,0.4}, L2 Nt=4 small): (1) reversibility (fwd then -pi back ->
recover U0,pi0 to CG floor, tol 1e-6), (2) dH~tau^2 (dH ratio nsteps->2*nsteps ~4 [3,5] = the decisive
gate that the summed force is the true action gradient), (3) a few HMCHasenbusch trajectories (dH +
acceptance). NEXT = build+run this HMC gate; on pass -> production Nf2 massless driver (clone
hmc_fermilab_wmass_L{2,4}_massless pattern, swap HasenbuschPF/MinimumNorm2Hasenbusch/HMCHasenbusch +
ladder), then the blowup-rate collapse study (C3).
>>>>>>> f09d67fcdaa20fe51ea94687ec214c92fdc2cd89

**Allocation (2026-07-11):** **qed3 is DEAD** (0 hrs -- omitted from `lquota`; normal-qos jobs won't
accrue -> effectively opportunistic-only). **affine ALIVE but tiny: ~2 kcore-hrs GPU** (`lq2_gpu`,
qos `normal,test,opp`), 0 used. Fairshare unchanged (~0.54 both). => the quick C1 gate runs fine on
affine (`~0.005` kcore-hr); the massless PRODUCTION campaign is NOT fundable on fnal until allocation
is restored. Fallback for anything bigger: **`--qos=opp`** (opportunistic, zero-allocation,
preemptible) works on BOTH affine and qed3. (`tmp_hb_bilinear_claude.sh` is affine, qos=normal.)

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

#### Where Term B comes from (derivation -- NM asked 2026-07-11)

The subtlety: in the ratio action the NUMERATOR operator $D_{i+1}$ is ALSO $U$-dependent (it contains
$D_\text{ov}$), not a frozen constant. Write $M_i \equiv D_i^\dagger D_i$ and
$\eta_i \equiv M_i^{-1} D_{i+1}^\dagger \phi_i = M_i^{-1}\chi_i$, so
$$
S_i=\phi_i^\dagger\,D_{i+1}\,M_i^{-1}\,D_{i+1}^\dagger\,\phi_i=\phi_i^\dagger D_{i+1}\eta_i .
$$
$U$ appears in THREE places -- the two explicit $D_{i+1}$ factors and the inverse $M_i^{-1}$:
$$
dS_i=\underbrace{\phi_i^\dagger (dD_{i+1})\,\eta_i+\eta_i^\dagger (dD_{i+1}^\dagger)\,\phi_i}_{\text{explicit numerator}}
\;+\;\underbrace{\phi_i^\dagger D_{i+1}\,(dM_i^{-1})\,D_{i+1}^\dagger\phi_i}_{\text{inverse}} .
$$
The inverse piece uses $dM_i^{-1}=-M_i^{-1}(dM_i)M_i^{-1}$ and $dM_i=dD_i^\dagger D_i+D_i^\dagger dD_i$,
giving $-\eta_i^\dagger(dD_i^\dagger D_i+D_i^\dagger dD_i)\eta_i = $ **Term A** (the ordinary PF-force
shape, bra $=$ ket $=\eta_i$, operator $D_i$). The explicit-numerator piece is
$z+z^*=2\Re\,z$ with $z=\phi_i^\dagger(dD_{i+1})\eta_i$, and since the mass part is $U$-independent
$dD_{i+1}=dD_\text{ov}=K$:
$$
2\Re\,\phi_i^\dagger (dD_\text{ov})\,\eta_i=2\Re\langle\phi_i|K|\eta_i\rangle=\textbf{Term B}.
$$
So Term B is the variation of the numerator $D_{i+1}$, a bilinear with the RAW (undifferentiated)
pseudofermion $\phi_i$ as bra and the solved $\eta_i$ as ket -- distinct from Term A's bra $=$ ket.
Because $dD_{i+1}=K$ is mass-independent, Term B can be evaluated by ANY frame's grad object (C1
overload), and the bra $\phi_i$ enters BARE (no $(1{+}m_L)$ fold).

#### Scalar shift vs measure-weighted, and operator sharing (NM asked 2026-07-11)

Verified: `OverlapWMass::update(U)` (`:290`) = `d_DW.update(U)` + `compute_lambda_max()`, BOTH
mass-INDEPENDENT (the Zolotarev poles/window are fixed at construction; the sign function is applied
per-`mult` in the pole loop, NOT cached in `update`). The mass is only a cheap diagonal add in
`mult`/`DHD` (`apply_mL`) with $m_L=\text{mass\_coeff}\cdot M_\text{mass}$, $M_\text{mass}$ built once.

Consequence: **operator sharing is orthogonal to scalar-vs-weighted.** We can share ONE $D_\text{ov}$
across all frames and vary only the mass with a cheap `set_mass(m)` (updates `mass`/`mass_coeff`,
sets `is_precalc=false`) -- this works for EITHER a scalar shift OR the measure-weighted mass, because
in both cases the mass is a $U$-independent diagonal and `update(U)` recomputes only the shared
mass-free part. The scalar shift does NOT unlock sharing by itself.

Two viable designs:
- **(A) vector of operators** (one `OverlapWMass` per frame): simplest logic, no mass-swap bookkeeping,
  but 3x the force-scratch MEMORY (`d_Zs/d_Ys/d_Zblock/d_Yblock/d_XZpre/...`, each $N(\text{size}-1)$ --
  the real cost, matters at L4) + a duplicated `compute_lambda_max` per frame per MD step (cheap).
- **(B) one operator + `set_mass`**: minimal memory, but frames MUST be processed strictly
  sequentially (the precalc scratch + `is_precalc` are single-buffered) and the mass must be re-set
  around each $D_i$ / $D_{i+1}$ apply in heatbath/action.

The scalar shift's genuine merits (independent of sharing): auxiliary Hasenbusch masses are just
regulators (need NOT be physical -- only the target $D_0=D_\text{ov}$ massless matters), it is
Hasenbusch's textbook form, and `set_mass` becomes a real scalar with no geometry factor. Its force is
STILL the same generic-bra $K$-contraction the C1 overload validated (bra $=(1{+}m_i)\eta$ for Term A,
$\phi_i$ for Term B) -- so force-reuse holds for scalar too. NM leaning: reconsider scalar shift +
one shared operator (design B) for minimal L4 memory. DECISION PENDING.

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

### TUNING STRATEGY (NM 2026-07-13, from collaborator Osborn et al. PoS LATTICE2024 052 "Automated
### tuning for HMC mass ratios", inspirehep 7648e35370ffac62bffd7cc2a5b14461)

DECISIONS: (integrator) nested-Omelyan multilevel NOW (`MinimumNorm2ML`, = their ABABA / G5F2,3
AHALAHALAHA family); force-gradient (their best, ABACABA) LATER. (tuning) adopt their COST FUNCTION as
the objective + MANUAL SCAN (grid/hand) -- NOT the full differentiable-HMC/Adam infra.

Their cost model (Eq. 4-7): $T_\text{eff}=\tau\sqrt{\langle P_a\rangle N}$; solve-cost per update
$C_S=\sum_s N_{CG}(m_s)/N_{CG}(m)\approx\sum_s m/m_s$ (sum over EVERY solve, cost $\propto 1/$mass,
verified 15%); **Cost $=C_S/(\tau^2\langle P_a\rangle)$**, minimize (equiv. maximize
$\langle P_a\rangle\tau^2/C_S$). Their integrators: ABABA (Omelyan), ABACABA (force-gradient C), G5F2
= AGABAGABAGA (gauge-only kicks G to cut fermion evals), G5F2,3 = AHALAHALAHA (light L / heavy H on
separate timescales). FG (ABACABA) best at 1 hMass; setup = staggered, thermalized cfg, 200 traj tune
+ 400 measure, n_s=1. Confirms OUR ordering: massless solve is the MOST expensive -> outer/fewest-eval
timescale.

**OUR cost metric (massless overlap adaptation):** the paper normalizes by $N_{CG}(m_\text{dyn})$ which
is ill-defined at $m_\text{dyn}=0$ and overlap $N_{CG}$ is NOT $\propto 1/m$; so use the ABSOLUTE solve
cost $C_S$ = TOTAL CG iterations per trajectory (sum over heatbath + all MD force `update_eta`/precalc
solves), measured by an iteration counter in the CG solver. Then **Cost $=C_S/(\tau^2\langle P_a\rangle)$**
(wall time is the ultimate cost; CG-iters is the clean proxy). TUNE: scan Hasenbusch masses + per-level
step multipliers to minimize Cost on the thermalized ensembles.

TODO for tuning: (i) CG-iteration counter in `MatPoly::solve`/`solve_multishift`, accumulated per
trajectory; (ii) HMC reports $C_S$, $\langle P_a\rangle$, $\tau$ -> Cost; (iii) manual scan.

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

**PENDING (after nsteps) -- tune $s_\text{tot}$ (trajectory length $\tau$) via the SMEARED PLAQUETTE (NM
2026-07-13).** Once the per-stage step counts are fixed, tune the trajectory length by how far the smeared
plaquette moves per unit cost. Wilson-flow smearing from `includes/flow_claude.h` (`Flow(&SW, t_flow,
nsteps); flow(U)` -- spatial-only per-timeslice flow). Order: window+ladder+nsteps (done) -> $s_\text{tot}$.

DESIGN (driver `test_hasenbusch_stot_claude.cu`, L-generic):
- **Config = min $\lambda_\text{min}$ (stiffest)** per L (Nf2 gsq8): L1 k=2105 ($\lambda_\text{min}$=0.205),
  L2 k=5 (0.0485), **L4 k=93 (0.0099, the stiff one)**. CLI-overridable.
- Frozen window + ladder + base steps `hasenbusch_steps(L)`.
- **$s_\text{tot} \in \{1.0, 1.4, 1.8\}$.** Loop **4 refreshed momenta** (redraw $\pi$ + pseudofermions);
  for each momentum use the **SAME $\pi$/$\phi$ across all $s_\text{tot}$** (fair comparison -- longer
  $s_\text{tot}$ integrates the same start further). **MDsteps FIXED to the tuned `hasenbusch_steps(L)`; ONLY
  trajL $=s_\text{tot}$ varies** -> $eps=s_\text{tot}/n_\text{steps}$ grows, so $N_\text{CG}$ (cost) stays
  ~constant while $dH$ grows. Pick the largest $s_\text{tot}$ with acceptable $dH$/acceptance.
- Observable: flow the evolved config by a fixed $t_\text{flow}$, take the **smeared gauge action** $S_g(U_\text{flow})$
  (= $\sum\theta^2$; the non-compact "smeared plaquette"). Report $\langle\Delta S_g\rangle = \langle|S_g(s_\text{tot})
  - S_g(\text{config})|\rangle$ (decorrelation), $\langle dH\rangle$, and $\langle\Delta S_g\rangle/s_\text{tot}$
  (decorrelation per cost), averaged over the 4 momenta. Pick $s_\text{tot}$ maximizing decorrelation-per-cost.

OPEN (confirm before building): (a) observable = flowed $S_g$ vs mean plaquette angle (non-compact avg $\approx0$,
so $S_g/\theta^2$ is the meaningful one); (b) $t_\text{flow}$ value (e.g. 1.0); (c) fine vs coarse step rounding.

### WINDOW + LADDER FROM THE LOW-MODE DISTRIBUTION (NM 2026-07-13)

Measure the $D_W$ low-mode distribution first (`eig_wilson_lowmode_claude.cu`, per config: $\sigma_\text{min}=\lambda_\text{min}$
= smallest singular value of $D_W$, $\lambda_\text{max}$, ratio). Then set, **beforehand and fixed for the run**:

1. **Zolotarev window** ($D_\text{ov}$ uses both edges; the sign approx lives on $[\lambda_\text{min},\lambda_\text{max}]$):
   $$\lambda_\text{max}^\text{set}=s_\text{hi}\max_\text{cfg}\sigma_\text{max}\ (s_\text{hi}\sim1.1\text{--}1.3),\qquad
     \lambda_\text{min}^\text{set}=s_\text{lo}\min_\text{cfg}\sigma_\text{min}\ (s_\text{lo}\lesssim1,\ \text{a bit BELOW the observed min so every config is covered}),$$
   window $k=\lambda_\text{min}^\text{set}/\lambda_\text{max}^\text{set}$ -- expected to land near the current $k=0.001$; keep
   $n_\text{pole}=21$ ($\delta\sim1.5\times10^{-5}$). If the tail drives $\min_\text{cfg}(\sigma_\text{min}/\lambda_\text{max})$ below $k$, lower
   $k$ / add poles. Fixed constants -> covers the ensemble AND gives HMC exactness (no per-config recompute).

2. **3-stage ladder $\{0, c_1, c_2\}$** (`set_mass_coeff`; $m_L=c\,\text{diag}(A_y/\bar A)$):
   - $c_1$ = the min-tail eigenvalue $\times\,O(1)$ (same low-mode scale as $\lambda_\text{min}^\text{set}$, but $c_1$ is the
     Hasenbusch REGULATOR, a separate number from the window edge). Not something to over-parametrize -- pick from
     the tail, keep the $O(1)$ as the tuning knob.
   - $c_2\sim0.2\text{--}0.5$ (well above all low modes; well-conditioned heavy frame), tuned later by cost.

Do this per $L$ ($L=1,2,4$) from each ensemble's own distribution.

**FROZEN VALUES FIXED (NM 2026-07-13)** -- `includes/frozen_window_claude.h` (single source of truth),
applied via `OverlapWMass::set_lambda(lmin,lmax)` (rebuilds Zolotarev on $k$, freezes both edges). $n=21$:

| $L$ | $\lambda_\text{min}$ | $\lambda_\text{max}$ | $k$ | $\Delta$ |
|---|---|---|---|---|
| 1 | 0.1   | 13 | $7.69\times10^{-3}$ | $2.54\times10^{-7}$ |
| 2 | 0.06  | 8  | $7.50\times10^{-3}$ | $2.72\times10^{-7}$ |
| 4 | 0.008 | 5  | $1.60\times10^{-3}$ | $7.08\times10^{-6}$ |

$\Delta$ from `zolotarev_delta_claude.cpp`; all better than the old fixed $k=0.001$ ($\Delta=1.5\times10^{-5}$).
Coverage caveats (accepted): L1 $\lambda_\text{max}=13<$ obs 14.6; L2 $\lambda_\text{min}=0.06>$ obs 0.0485; L4
$\lambda_\text{max}=5<$ obs 6.5. Wired into `test_hasenbusch_tune_l1_claude.cu`; production drivers PENDING.

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
- **C6 (INTEGRATOR DONE+VALIDATED 2026-07-13; tuning pending) -- 3-level multi-timescale split.**
  C6a: `HasenbuschPF` per-group force (`update_eta_frame(s)` + `get_force_frames(dSf,U,i_lo,i_hi)`).
  C6b: `includes/hmc_hasenbusch_ml_claude.h` = `MinimumNorm2ML` (Grid `Integrator_algorithm.h`
  MinimumNorm2::step mimicked VERBATIM) + polymorphic `MDLevel`/`GaugeLevel`/`FermionGroupLevel` (lambda-
  free) + `HMCHasenbuschML`. C6c PASS (`test_hasenbusch_ml_claude.cu`): 2-level regression ==
  MinimumNorm2Hasenbusch (dH+U ~1e-9), 3-level reversibility ~1e-12, dH~tau^2. C6d (PENDING): cost-function
  instrumentation (CG-iter counter -> Cost=C_S/(tau^2<P_a>)) + manual scan of masses+multipliers; then wire
  `MinimumNorm2ML`+`HMCHasenbuschML` into the production driver. See TUNING STRATEGY section above.
  **C6d IN PROGRESS (2026-07-13): tuning-landscape driver.**
  - CG-iter counter added to `includes/matpoly_claude.h`: file-scope `g_matpoly_cg_iters` +
    `reset_cg_iters()`/`get_cg_iters()`; every `solve`/`solve_multishift`/`solveAsync` adds its Krylov
    count `k` on exit. `unfreeze_lambda()` added to `overlap_wmass_claude.h` (re-freeze per config).
  - Driver `test_hasenbusch_tune_l1_claude.cu` (+ `tmp_hb_tune_l1_local_claude.sh`): L1, Nt=128, massless.
    PROBE (per NM): load the last `N_CFG`(=5) `ckpoint_lat.*` of `Nf<Nf>_gsq<gsq>...nt128L1/`; per config
    freeze `lambda_max`, then `N_DRAW`(=3) fresh momenta+pseudofermions integrated ONE traj from the FIXED
    config (U reset each draw, no accept). Cost = `<N_CG>/(tau^2 <P_acc>)`, `<P_acc>=<min(1,e^{-dH})>`,
    tau=1. Reset the CG counter AFTER heatbath so `C_S` = MD integration cost only.
  - Setups surveyed (mass ladder + per-level `MinimumNorm2ML` step multipliers, outer->inner=gauge):
    `{0,.1,.4}` 3-level per-mass `[0]x1 [1]x2 [2]x3` (MD 6/8) vs grouped `[0-1]x1 [2]x2`; mass variants
    `{0,.05,.2}`, `{0,.2,.6}`; K=1 ladders `{0,.1}`/`{0,.2}` `[0]x1 [1]x3`. Default scan gsq{16,12} Nf2.
    Set tuned values INDEPENDENTLY per L (L1 first, then L2/L4). Landscape run pending (user runs the .sh).
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
