# Two-operator (split-pole) overlap HMC -- impl plan

## Goal / physics

Use **two** overlap operators in the massless Hasenbusch HMC, sharing the same Wilson operator
$D_W$ and the same frozen Zolotarev window $(\lambda_\text{min},\lambda_\text{max})$, differing only in the
**pole count** $n$ of the Zolotarev sign approximation:

- **Action operator** $D$ ($n = 21$, accurate): pseudofermion **heatbath** and **accept/reject**
  (the Hamiltonian $H$).
- **Force operator** $D_f$ ($n_f \in \{5,7,9,11\}$, cruder): the **MD force** only.

### Why this is exact

HMC is exact for **any** reversible, area-preserving MD map, provided the accept/reject uses the true
$H$ (Duane, Kennedy, Pendleton, Roweth, Phys. Lett. B 195 (1987) 216). The leapfrog / Omelyan momentum
kick $p \to p - \tfrac{\epsilon}{2}F(U)$ and drift $U \to U + \epsilon\,p$ are unit-Jacobian shears for
**any** deterministic force $F(U)$ -- $F$ need not be the gradient of the accept/reject action. So
replacing the MD force by $-\,dS_f/dU$ built from the cruder $D_f$ leaves detailed balance intact; the
Metropolis step corrects the integration error.

### Why it helps

- **Cost**: fewer poles $\Rightarrow$ larger smallest Zolotarev node $\Rightarrow$ the inner multishift CG
  (in both the force eta-solve and `precalc_grad`) converges faster, and the per-pole reconstruction /
  force COO loop shrinks $\propto n_f$.
- **Tamer force ("avoid infinite slope")**: the exact sign has infinite slope at zero; the rational
  approx slope $\sim 1/\Delta$. Fewer poles $\Rightarrow$ larger $\Delta$ $\Rightarrow$ the smallest force
  resolvent shift is larger $\Rightarrow$ the near-zero-mode force spikes are bounded. This is the main
  integrator-stability win, especially at L4 (near-zero $D_W$ tail).

### Diagnostic verdict (2026-07-14, `test_hasenbusch_npole` probe, gsq8 Nf2)

Per-frame force CG at L1/L2 for $n_f\in\{9,11,13\}$ showed the **pole count is a weak speed knob**: L1 frame-0
went $10710\to11489$ CG ($+7\%$) from $n_f=9\to13$; frame-1 (heavy $c=0.5$) is the *cheaper* frame, not the
bottleneck. The reason: the multishift solves ALL poles in one shared Krylov pass, so the CG **iteration
count is $n_f$-independent** (set by the inner conditioning); more poles add only per-iteration AXPYs. The
real cost is the inner $D_\text{ov}$-apply iteration count $\times$ lattice size. So the effective levers are
the **force window** and the **force solve tolerance**, below.

## Force solve tolerance (the effective cost lever)

The force has two tolerance-gated CG types; the **inner** one dominates:

1. **Inner multishift** -- every $D_\text{ov}$ application (partial-fraction $\text{sign}(H)v=\sum_m A_m(H^2+\sigma_m)^{-1}Hv$)
   and every `precalc_grad` resolvent is a multishift CG to `Comp::TOL_INNER` $=10^{-9}$. Called $\sim 2N_\text{outer}+O(1)$
   times per force eval. **This is the bottleneck.**
2. **Outer eta CG** -- `Op_DHD_f.solve` on $D_f^\dagger D_f$ to `Comp::TOL_OUTER` $=10^{-8}$: $N_\text{outer}$ iterations.

### Cost model

For a CG on a system with condition number $\kappa$ to relative tolerance $\tau$,
$$
N_\text{CG} \;\sim\; \sqrt{\kappa}\,\log(1/\tau).
$$
The inner shifted system $(D_W^\dagger D_W/\lambda_\text{max}^2+\sigma)$ has smallest shift $\sigma_\text{min}\sim k_f^2$
($k_f=\lambda_\text{min}^f/\lambda_\text{max}$), so $\kappa\sim k_f^{-2}$ and
$$
N_\text{inner} \;\sim\; \frac{1}{k_f}\,\log(1/\tau_\text{inner}).
$$
Two independent multipliers:

- **Window** ($k_f\to 2k$ via $\lambda_\text{min}^f=2\lambda_\text{min}$): $N_\text{inner}\to N_\text{inner}/2$.
- **Tolerance** ($\tau_\text{inner}:10^{-9}\to10^{-5}$): $\log(1/\tau)$ from $20.7\to11.5$, factor $0.55$.

Combined $\approx 0.5\times0.55\approx 0.28$ -- a **$\sim3$--$4\times$** inner-cost cut. The L4 verbose trace
confirms the tolerance part: the seed reached relative $10^{-5}$ by $\sim$iter 310 vs iter 547 for $10^{-9}$
($\sim43\%$ wasted over-convergence).

### Reversibility bounds the tolerance (NM 2026-07-14 correction)

Earlier claim "looser force tol is free because Metropolis corrects" was **too strong**. The correct picture:

- **Area preservation**: holds for ANY deterministic $F$ (leapfrog = composition of unit-Jacobian shears).
- **Reversibility, exact arithmetic**: a *zero-start* CG to a fixed relative tolerance is a deterministic
  function of $U$, and the leapfrog is exactly reversible for any such $F$ -- $F(U')$ enters identically as
  the forward step's 2nd kick and the backward step's 1st kick and cancels, however inaccurate $F$ is. So
  tolerance alone does not break exactness *in principle*.
- **Reversibility, finite precision**: it DOES break, at $O(\tau)$. Forward/backward passes don't retrace to
  machine precision; a looser solve amplifies the force's non-smoothness (CG iter count jumps $\pm1$ between
  nearly-equal $U$, giving $O(\tau)$ jumps in $F$). This violation is NOT corrected by Metropolis (the
  acceptance formula *assumes* reversibility) $\Rightarrow$ it biases the chain.
- **Warm-start caveat**: a chronological / previous-solution initial guess makes $F$ history-dependent and
  breaks reversibility even in exact arithmetic. The force solves MUST zero-start (verify `MatPoly::solve`
  / `solve_multishift` memset $x=0$).

**Consequence**: the force tolerance is bounded by a *measured* reversibility test (integrate forward then
backward, check $|\Delta U|,|\Delta p|,|\Delta H_\text{rev}|$), not by acceptance alone. Stiffer L4 demands a
tighter force tol than L1.

### Multiprecision solves (NM 2026-07-14) -- the reversibility-CLEAN speedup

Instead of *loosening* the tolerance (reversibility-bounded), reach the SAME tight tolerance FASTER with a
mixed-precision CG: fp32 Krylov iterations + periodic fp64 residual correction ("reliable updates" /
defect correction). Ref: Clark, Babich, Barros, Brower, Rebbi, Comput. Phys. Commun. 181 (2010) 1517,
arXiv:0911.3191.

- **Reversibility-clean**: the solution is still fp64-accurate and deterministic (given fixed-order
  reductions -- no non-deterministic atomics), so $F(U)$ is as accurate/reversible as a pure fp64 solve.
  This is the key advantage over tolerance loosening: speed *without* the reversibility tax.
- **Gain**: $\sim2\times$ on TITAN V / A100 (fp32 = 2$\times$ fp64 throughput + half the memory bandwidth on
  the bandwidth-bound sparse matvec). Multiplies with the window narrowing (also clean).
- **Scope caveat**: NO fp32 infrastructure exists -- `CSR`/`LinOp`/`MatPoly` are all `CuC` (fp64). Needs an
  fp32 sparse $D_W$ matvec + fp32 CG kernels + the reliable-update outer loop. And the *bottleneck* is the
  inner MULTISHIFT, whose mixed-precision variant (per-shift reliable updates) is more involved than the
  single-shift outer solve. So this is a dedicated implementation, larger than the window/pole/tol knobs.

### Lever cleanliness

- **Reversibility-clean** (fixed deterministic operators; no tolerance dependence): the narrower force
  **window** $\lambda_\text{min}^f=2\lambda_\text{min}$ (halves inner iters at *fixed* tol) and fewer poles
  $n_f$. Lean on these.
- **Reversibility-bounded**: the solve **tolerance**. Loosen modestly and MEASURE. Match-to-$\Delta(n_f)$
  ($\Delta\approx7\times10^{-4}$ L1/L2, $4\times10^{-3}$ L4 at $n_f=11$) is the *floor* set by the pole
  approximation, but the reversibility test sets the *real* limit -- start conservative at $\tau_f=10^{-6}$
  and loosen only while the reversibility violation stays negligible.

### Implementation

- Add member tolerances to `OverlapWMass` (`tol_inner`, default `Comp::TOL_INNER`), used by `precalc_grad`
  and the $D_\text{ov}$ apply instead of the hard-coded constant; set looser on $D_f$ only (`Df.set_force_tol(...)`).
- The outer eta solve is already an explicit argument (`Comp::TOL_OUTER` in `update_eta_force_frame`) -- pass
  a force-specific $\tau_\text{outer}$ there.
- Suggested start: $\tau_\text{inner}^f=\tau_\text{outer}^f=10^{-4}$ (or scan $\{10^{-3},10^{-4},10^{-5}\}$),
  with $\lambda_\text{min}^f=2\lambda_\text{min}$ and $n_f=11$. Recheck force norms + timing + acceptance.

### The force action (Option B -- force operator carries the WHOLE force)

Per Hasenbusch frame $i$, with the SAME drawn $\phi_i$ (heatbath, action operator), the MD force is the
gradient of a **force action** built entirely from $D_f$:

$$
S_{f,i} = \phi_i^\dagger\, D_{f,i+1}\,\big(D_{f,i}^\dagger D_{f,i}\big)^{-1}\, D_{f,i+1}^\dagger\, \phi_i ,
\qquad
\chi_{f,i} = D_{f,i+1}^\dagger \phi_i ,
\qquad
\eta_{f,i} = \big(D_{f,i}^\dagger D_{f,i}\big)^{-1} \chi_{f,i} .
$$

(heaviest frame $i=K$: $\chi_{f,K}=\phi_K$, $S_{f,K}=\phi_K^\dagger(D_{f,K}^\dagger D_{f,K})^{-1}\phi_K$.)
Term A = `precalc_grad(D_f, eta_f)`, Term B (ratio frames) = bilinear($\phi_i$, $\eta_{f,i}$) -- a byte
mirror of the current force, but every operator apply / solve uses $D_f$. The accept/reject action stays
exactly the current $S = \sum_i \chi_i^\dagger \eta_i$ with $\eta_i,\chi_i$ from the **action** operator $D$.

### Reference

- Exactness of Metropolis-corrected inexact MD: Duane-Kennedy-Pendleton-Roweth 1987 (above).
- Different rational approx in MD vs action, corrected by accept/reject: Clark & Kennedy, "Accelerating
  dynamical-fermion computations ... (RHMC)", Phys. Rev. Lett. 98 (2007) 051601, hep-lat/0608015.

## Files to modify

- `includes/hasenbusch_ladder_claude.h` -- add `hasenbusch_nforce(L)` (per-L force pole count; tunable).
- `includes/pseudofermion_hasenbusch_claude.h` -- add force operator $D_f$, its DHD solver, `d_chi_f`,
  `d_eta_f`, `set_frame_mass_force`, `update_eta_force_frame(s)`; retarget `get_force_frames` onto $D_f$.
- `includes/hmc_hasenbusch_ml_claude.h` -- `FermionGroupLevel` gets a force-op pointer; its `get_force`
  updates $D_f$ and solves the FORCE eta; `HMCHasenbuschML::run` re-solves the ACTION eta at the final
  $U$ before $h_1$ (so accept/reject uses the accurate action, no longer a force side effect).
- `hmc_hasenbusch_claude.cu` (production) -- construct $D_f$ ($n_f$ = `hasenbusch_nforce(L)`, same
  `set_lambda`), pass it to `HasenbuschPF` and `FermionGroupLevel`.
- NEW `test_hasenbusch_npole_claude.cu` + `tmp_hb_npole_local_claude.sh` -- scan $n_f\in\{5,7,9,11,21\}$ at
  fixed ladder/steps: force CG iters (`get_cg_iters`), force wall-time, dH, acceptance, Osborn Cost.

## Ordered chunks

**Chunk 1 -- HasenbuschPF split (Files: pseudofermion_hasenbusch_claude.h).**
Add `Fermion& Df`, `LinOpDHDWrapper M_DHD_f`, `MatPoly Op_DHD_f`, `d_chi_f`/`d_eta_f` buffers,
`set_frame_mass_force(i)`, `update_eta_force_frame/_frames`. Point `get_force_frames` at $D_f$/`d_eta_f`.
Constructor takes `(Fermion& D, Fermion& Df, masses, base)`.

**Chunk 2 -- integrator + HMC (Files: hmc_hasenbusch_ml_claude.h).**
`FermionGroupLevel` holds `fermion` (action, for the Wilson-link update parity) + `fermion_force`; its
`get_force` does `fermion_force->update(U)` + `update_eta_frames_force` + `get_force_frames`. Add the
action-eta re-solve in `run()` right before `h1 = H()` (and keep the reject-path `fermion->update(U)`).

**Chunk 3 -- production driver (Files: hmc_hasenbusch_claude.cu).**
Construct `Fermion Df(DW, mass, nforce, 0.001); Df.set_lambda(frozen_window); Df.update(U);`
`nforce = hasenbusch_nforce(N_REFINE)`. Thread `Df` into `HasenbuschPF` + `FermionGroupLevel`.

**Chunk 4 -- scan harness (Files: test_hasenbusch_npole_claude.cu, tmp_hb_npole_local_claude.sh).**
Loop $n_f\in\{5,7,9,11,21\}$; action op fixed $n=21$; fixed `hasenbusch_ladder/steps`. Report per $n_f$:
force CG iters, force sec, $\langle|dH|\rangle$, $P_\text{acc}$, Linf force norm, Cost $=N_\text{CG}/(\tau^2 P_\text{acc})$.

## Decisions (NM 2026-07-14)

1. **Split depth = Option B** (force op $D_f$ carries the whole force: eta-solve + grad; action op does
   heatbath + accept/reject, action-eta re-solved at final $U$ before $h_1$).
2. **Scope = Hasenbusch only** (`hmc_hasenbusch_claude.cu` + scan). Non-Hasenbusch driver untouched.
3. **Default `hasenbusch_nforce(L) = 11` for all L**, until the $n_f$ scan fixes per-L values.
4. **Two persistent etas in parallel** (NM): `d_eta` (action) and `d_eta_f` (force) are BOTH globally
   allocated per-frame members, symmetric to the two `d_chi`/`chi_f` sides. `chi_f` stays a single scratch
   (`d_chif`) since the force never reads it. `update_eta_force_frames` mirrors `update_eta_frames`.

## Implemented (2026-07-14)

- `hasenbusch_ladder_claude.h`: `hasenbusch_nforce(L)` (all L -> 11). Also L2 fixed to `{0,0.5}`, steps `{3,3}`.
- `pseudofermion_hasenbusch_claude.h`: `Fermion& Df` + `Op_DHD_f`/`M_DHD_f`, `d_eta_f[]` (global),
  `d_chif` scratch, `set_frame_mass_force`, `update_eta_force_frame(s)` (solve `d_eta_f` with `Df`).
  `get_force_frames` reads `d_eta_f` (op = `Df`). `get_force()` is now SELF-CONTAINED
  (`update_eta_force_frames` + `get_force_frames`) so legacy single-shot callers work untouched. A
  delegating 3-arg ctor (`Df = D`) preserves old single-op call sites (force == action == n=21).
- `hmc_hasenbusch_ml_claude.h`: `FermionGroupLevel` drives `fermion_force` (updates `Df`, solves the FORCE
  eta, computes the force). `HMCHasenbuschML::run` re-solves the ACTION eta at final $U$ before $h_1$
  (the MD no longer refreshes `d_eta`).
- `hmc_hasenbusch_claude.cu`: constructs `Df` (`n_f = hasenbusch_nforce(L)`, same `set_lambda`), threads
  it into `HasenbuschPF` + `FermionGroupLevel`.
- Manual-`h1` test drivers (tune, validate part 2, stot, ml) get the same post-`integrate` action re-solve
  (`D.update(U); pf->update_eta();`). `hmc_l2`/`force_l1` use the old integrator / self-contained
  `get_force` and are unchanged.
- NEW `test_hasenbusch_npole_claude.cu` + `tmp_hb_npole_local_claude.sh` -- the $n_f\in\{5,7,9,11,21\}$ scan.

## Build / validate (USER runs)

1. `bash tmp_hb_npole_local_claude.sh` -> `test_hasenbusch_npole_claude.log` (the scan; also flushes any
   compile errors from the shared headers).
2. Re-validate exactness with the split ON: rebuild + `bash tmp_hb_validate_local_claude.sh` (n_f=11
   default) -- dH~tau^2 should stay floor-free and acceptance healthy.
3. Rebuild production `hmc_hasenbusch_claude.cu`.
