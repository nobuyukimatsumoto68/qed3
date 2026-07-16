# Deflation + LMA for qed3 (L4 massless only) -- impl plan

## >>> OUTCOME / SUPERSEDED (2026-07-15) -- read `memory/project_overlap_wall_additive_mass.md` first <<<
The deflation/LMA goal is effectively CLOSED by a physics finding.  The overlap `D_ov = 1 + M_DW/|M_DW|` is a
UNITARY PROJECTOR of `D_W` (keeps the PHASE of `M_DW = D_W - M`, M=1); a near-zero overlap mode needs a real
`Re(D_W) < M=1`.  The FREE physical mode sits at `Re~0` (massless dispersion); the GAUGE field adds an
ADDITIVE MASS (~1.04 at L2, ~0.68 at L4, shrinking with L) that shifts the interacting mode RIGHT.  At gsq=8
the additive mass OVERSHOOTS the wall on coarse L -> interacting mode is RIGHT of M=1 at L1/L2 (overlap
GAPPED, no near-zero modes) and only crosses LEFT (light) at L4.  So **gsq=8 is too strong for the coarse
lattices**, which is the true reason for the "L4-only near-zero tail" -- and why overlap deflation/LMA is
MOOT at L1/L2 (nothing to exploit).  Tool built + validated: **shift-invert Arnoldi** `eig_arnoldi_claude.cu`
(`DW_ARNOLDI`; D_W^{-1} via CG on `LinOpWilsonSq`; CSR-only, L4-feasible; matches dense on free L1).  The
Lanczos IRL adopt/validate saga below (chunks 1a-1e) is the earlier, superseded path (overlap A eigensolver +
Rayleigh-Ritz).  Files: `eig_arnoldi_claude.cu`, `eig_lanczos_claude.cu`, `eig_wmass_val_claude.cu`
(`-DSPECTRUM_DW/DOV`), plots `eig_dw_dov_compare_claude.py` / `eig_spectrum_scatter_claude.py` /
`eig_dw_irl_overlay_claude.py`, handoffs `tmp_eig_dw_*`.

## Goal / scope

Low-mode **deflation** (CG speed) and **low-mode averaging / LMA** (variance) for the L4 massless gsq=8
measurements, where a near-zero-mode tail exists. Requires implementing an **IRL (Implicitly Restarted
Lanczos) eigensolver FROM SCRATCH** in the qed3 CUDA code (we have none). Method template = the Grid disc-LMA
suite at `/mnt/baracuda_14/grid_claude/Grid/examples/` (portable methodology, NOT portable code -- Grid/DWF).

## (a) Viability -- CONFIRMED L4 ONLY (2026-07-15, from `wilson_lowmode_Nf2_gsq8.000000_L{1,2,4}_claude.dat`)

`lambda_min` = smallest SINGULAR value of $D_W$ per config:
| L | min | median | max | near-zero tail |
|---|----:|-------:|----:|---|
| L1 | 0.205 | 0.260 | 0.301 | NONE (strongly gapped) -> not worth |
| L2 | 0.049 | 0.087 | 0.108 | none <0.02, 1/320 <0.05 -> not worth |
| L4 | 0.0099 | 0.074 | 0.295 | **min 0.0099, 9% (15/168) <0.05, 3 <0.02** -> VIABLE |

So: **L4 only.** L1/L2 gapped -> no low modes to exploit (matches NM's intuition; blocking already covers L1/L2).

## Lever payoff (from Grid `disc_mrhs_defl_bench_results_claude.md`, 24^3x48 m=0.01, DENSE spectrum)

- **Deflation = modest speed.** Nev 25/50/100 -> 1.27/1.45/1.58x; per-config net ~1.5x. Diminishing returns on
  a dense low cluster (deflating Nev only lifts the floor to eval[Nev]). Their verdict: "bigger Nev is the
  variance/LMA lever, not a speed win." Our L4 tail is SPARSER (fewer modes) -> a few modes should lift the floor
  more per-mode, but total mode count is small -> a modest, config-tail-targeted speed gain.
- **LMA = the statistics prize.** Exact low-mode $L^\text{low}$ (all source points) + high-mode stochastic with the
  low modes PROJECTED OUT; unbiased to solver tol (Grid DET gate $1.6\times10^{-8}$). Helps ONLY in the near-zero
  regime -- which L4 is. Best target = the stochastic-limited **disc / condensate** at m=0.

## REVISED chunk 1 (2026-07-15): the IRL ALREADY EXISTS -- ADOPT + VALIDATE, do NOT rewrite

**Discovery.** `includes/lanczos_claude.h` (2026-05-22, plan `lanczos_impl_plan_claude.md`) is a COMPLETE
Sorensen IRL + Chebyshev eigensolver, already targeting exactly $A=(D_{ov}+m)^\dagger(D_{ov}+m)$ via
`LinOpDHDWrapper` -> `OverlapWMass::DHD_deviceAsyncLaunch`, and already producing the low eigenVECTORS
(`s.d_evec[j]`, j=0..Nk-1) AND eigenvalues (`s.eval`). Components: `ChebyshevFilter`, `lanczos_step`
(periodic re-orth every 4), `diagonalize_Eigen`, `QR_decomp`, `basisRotate`, power-iter $\lambda_{max}$,
per-mode residual convergence. Driver `saved_scripts_claude/eig_lanczos_claude.cu`. So we do NOT write from
scratch (the earlier "we have none" premise was wrong). Caveat: NO validation log exists and the `.o` (Jun 8)
predates the Jul-13 `overlap_wmass_claude.h` edits -> STALE, must rebuild + validate.

**Validation design (normality-free).** The dense reference `saved_scripts_claude/eig_wmass_claude.cu`
(`cusolverDnXgeev`, same `OverlapWMass`) builds the spectrum of $D_m=D_{ov}+m$ ITSELF (`mult`), whose
eigenvalues are complex and equal the IRL's $\sqrt{\lambda_i}$ in modulus ONLY if $D_m$ is normal. To avoid
that assumption, the validation copy `eig_wmass_val_claude.cu` builds $A=D_m^\dagger D_m$ DENSELY (bind
`DHD_deviceAsyncLaunch` in place of `mult`), geev -> REAL eigenvalues, sort ascending, smallest Nk compared
1:1 to IRL `eval`. Both drivers run FREE-FIELD (cold gauge, as both references already do) at matched
$L$/$N_t$ via `-DLREF -DNT_CLI` and identical $M5{=}{-}1$, $r{=}1$, $T{=}4$, mass. Handoff
`tmp_eig_validate_claude.sh` (USER runs) + compare `compare_eig_validate_claude.py`. Pass = lowest Nk agree
to ~solver tol. THEN parametrize config-read + rerun IRL at L4 for the low-mode density (chunk 2).

**CONVENTION FIX (2026-07-15), settled against `/home/nobu/Downloads/lanczos-2.pdf` Eq (3.3).** First
validations did NOT converge; two of my choices were BACKWARDS. The PDF map
$q(\lambda;\alpha,\beta)=\frac{2\lambda^2-(\alpha^2+\beta^2)}{\alpha^2-\beta^2}$ is monotonically increasing
and maps $[\beta,\alpha]\to[-1,1]$, so **$\alpha>\beta$** with **$\alpha$ = spectral UPPER bound**
($\ge|\lambda_1|$) and **$\beta$ = upper edge of the wanted LOW band** ($\gtrsim|\lambda_{n-m+1}|$). Wanted
small modes ($\lambda<\beta$) map to $q<-1$; the PDF states $T_n(z)>0$ for $|z|>1$ **only when $n$ is EVEN**,
so **degree must be EVEN**. My driver had passed $\alpha$=small/$\beta$=large (OPPOSITE) and forced ODD --
both wrong. The `ChebyshevFilter` in `lanczos_claude.h` evaluates Eq (3.3) faithfully; only its *comment*
names $\alpha/\beta$ backwards. **SCALE (confirmed from the dense geev spectrum, massless free-field L1 Nt32):**
$D_{ov}$ eigenvalues lie on the GW circle $|z-1|=1$, so its singular values are in $[0,2]$ ($D_m$: $[|m|,2+|m|]$)
-- **$\alpha=2+|m|$** (measured max singular $=1.9993$; the earlier "$=1$" was WRONG). $\beta$ is the wanted
low-band edge and is SPECTRUM-DEPENDENT: free-field is GAPPED with NO near-zero modes (smallest singular
$\approx0.98$, 8-fold; next $\approx1.285$), so $\beta\sim1.1$--$1.4$ for validation; the interacting L4
near-zero tail wants a SMALL $\beta$. `Dov.lambda_min` is the WILSON $D_W$ bound ($\approx1.36$ here, a
DIFFERENT scale) -- NOT an overlap edge; `2*Dov.lambda_min`$=2.73>\alpha$ tripped the assert. `eig_lanczos_claude.cu`
now: $\alpha=2+|m|$ (auto), $\beta$ input (default $0.6\alpha$), even degree, assert $\alpha>\beta$. Validation
handoff: $\alpha$ auto, $\beta=1.15$, `DEG=12`. (Grid = PDF up to a shift-by-1 in the Cheby order -- an API
detail our honest $T_{\_deg}$ recurrence does NOT inherit; verified the recurrence yields exactly $T_{\_deg}$.)

**CONVERGENCE FIX (2026-07-15): full re-orthogonalization + Nk=1 for the degenerate free field.** With the
convention/window CORRECT, the IRL Ritz values landed right on the answer (Ritz$[0]=0.9608461$ vs dense
$0.9608406$) but the residual test read $0/8$ forever, drifting to ghosts (0.970, 0.980, 1.010). Cause: the
free-field low level is **8-fold EXACTLY degenerate** (dense: eight modes at $0.960840632881936$), and
(i) `lanczos_claude.h` used PARTIAL re-orth (`orth_period=4`) which loses orthogonality on a degenerate
cluster -> ghosts, and (ii) single-vector Lanczos cannot pull $N_k=8$ copies of ONE eigenvalue from one
start vector. FIX: `orth_period` 4 -> **1 (full re-orth)** in `lanczos_claude.h`; validate with **$N_k=1$**
(the single smallest converges cleanly even when degenerate). Handoff now `NK=1 NM=16`. MULTI-mode
validation needs a NON-degenerate spectrum = an INTERACTING config (which is also the real L4 target), so it
folds into chunk 2 (config-read).
2. **LMA split for condensate/disc** -- exact $L^\text{low}=\sum_i$ (all-source low-mode contraction) + projected
   high stochastic ($u_i = M v_i/\sigma_i$ construction so the split cancels term-by-term). Template:
   `examples/disc_lma_v2_common_claude.h` (`BuildLowModes`/`BuildA2ASet`/`BuildPhysicalA2A`) +
   `disc_lma_estimator_bench_impl_plan_claude.md` (unbiasedness derivation). Start with the CONDENSATE (simplest
   one-point), then the vector/axial disc.
3. **(secondary) deflated CG** for L4 conn `op_Dmsq` -- project the Nev modes out of the initial guess/residual.
   Modest per Grid; do only if the eigensolve is already built and cheap to reuse per config.

## Chunks (revised: IRL ADOPTED, not from scratch)

1. **Validate the adopted IRL** (`lanczos_claude.h`) vs dense geev (`eig_wmass_val_claude.cu`).
   1a. Free-field L1 Nt32, $N_k=1$ (bulletproof; degenerate) -- handoff `tmp_eig_validate_claude.sh`. IN PROGRESS.
   1b. INTERACTING L1 Nt128 (N=3072), multi-mode $N_k=8$ (non-degenerate) via `CONFIG_LAT` env.
       **NO DENSE** (dense geev at N=3072 is too expensive): the IRL is SELF-VALIDATING -- the driver
       recomputes, per returned mode, the relative residual $\lVert Au-\lambda u\rVert/(\lambda\lVert u\rVert)$
       (the DEFINITION of an eigenpair; PASS if all $<10^{-7}$). The cheap free-field N=768 dense already
       anchors that the filter targets the SMALLEST. Handoff `tmp_eig_validate_interacting_claude.sh`.
       **beta PLACEMENT (no dense to read $\mu_{min}$):** beta must sit ABOVE $\mu_{min}$ or the wanted band
       is EMPTY and IRL drifts into the bulk (first try $\beta=1.0$ gave Ritz$\sim3.0$, `evalMaxApprox=3.98`).
       Interacting low modes sit HIGHER than free field ($\mu_{min}=0.96$) since $D_{ov}.\lambda_{min}=0.528$
       (Wilson) is a bigger gap. So place beta GENEROUSLY in the lower half of $[0,2]$: $\beta=1.5$
       ($\mu<2.25$), monotone-amplified so the lowest 8 are the top-8 filtered. Widen if still missed.
       **RESULT (2026-07-15): IRL WORKS.** With $\beta=1.5$ it converged the lowest 4 (stable at $\mu=$
       2.2269, 2.2276, 2.2406, 2.2436; 4/8) -- the 5th sits just above the band edge $\mu=2.25$. This config
       is GAPPED: smallest overlap SINGULAR value $\approx\sqrt{2.227}\approx1.49$, i.e. **NO near-zero modes
       at L1 gsq8** (direct overlap confirmation that L1 is too coarse to deflate). Bumped $\beta\to1.8$
       ($\mu<3.24$) to bracket $\ge8$. Machinery is VALIDATED (stable converged eigenpairs). **L1 CLOSED** --
       too coarse for deflation regardless.
   **OPERATOR-MATCH FIX (2026-07-15) -- the "slow convergence" at L2.** Two things were wrong vs production:
   (i) **FROZEN Zolotarev window.** Production does NOT recompute the sign-function window per config -- it
   FREEZES it (`frozen_window_claude.h` + `OverlapWMass::set_lambda`, the exactness fix): L1 (0.1,13), L2
   (0.06,8), L4 (0.008,5), npole 21. My driver let `Dov.update` recompute `lambda_min/max` per config. FIXED:
   both eig drivers now `frozen_window(L,lmin,lmax); Dov.set_lambda(lmin,lmax); Dov.update(U)` -- exactly like
   `hmc_...:190`. (ii) **$a_t=0.2$** (was $T/N_t=0.03125$): the ensemble + all production (`jj_...:262`,
   `hmc_...:178`) use $a_t=0.2$, and the frozen edges were chosen from the $a_t=0.2$ Wilson scan -- wrong
   $a_t$ pushes $D_W$'s spectrum outside [lmin,lmax] -> sign wrong -> slow/failing inner CG. Both FIXED in
   `eig_lanczos_claude.cu` + `eig_wmass_val_claude.cu`. (Free-field L1 machinery validation was $a_t$/window
   -independent and still stands; only physics runs L2/L4 need the match. Note the Chebyshev filter window
   $\alpha=2,\beta$ is a SEPARATE thing -- the A-eigenvalue range for the eigensolver, not the Zolotarev
   Wilson window.)
   1c. **L2 intermediate** (N=10752, dense ALREADY infeasible -> self-validating IRL is the only option).
       Handoff `tmp_eig_lowmode_L2_claude.sh`. **DEBUG TRAIL (2026-07-15):**
       - First tried `ckpoint_lat.5` -> IRL found only bulk $\mu\approx3.5$, no near-zero mode. TWO causes:
         (i) traj-5 is **NON-thermalized**; (ii) its Wilson $\lambda_{min}=0.0485$ is BELOW the frozen L2
         window edge 0.06, so the Zolotarev sign is mis-fit exactly at the near-zero mode -> the |z|~0.05
         overlap mode is distorted up into the bulk.
       - **Filter is NOT the problem** (`cheby_filter_plot_claude.py` -> `cheby_filter_claude.png`): the
         $\alpha=2,\beta=1,n=12$ filter gives gain **2e5 at $\mu{=}0.05$** and $<1$ in the bulk -- it would
         trivially isolate a 0.05 mode if the operator had one. (Header squares its args: $\alpha=2\Rightarrow
         \alpha^2=4=\max\mathrm{eig}(D_{ov}^2)$; so "$\alpha$ should be 4" = the eigenvalue bound, achieved
         with header-$\alpha=2$.)
       - FIX: **THERMALIZED config `ckpoint_lat.740`** (traj 740/1601; its Wilson $\lambda_{min}=0.0643$ sits
         ABOVE the frozen 0.06 -> near-zero mode COVERED) + **DEG 12->8** (per NM: 2e5 gain is excessive
         dynamic range; deg=8 -> 2.7e3, still $\gg$ suppressed bulk $\le1$). $\alpha$ auto (=2), $\beta=1.0$,
         $N_k=12$. READY to re-run.
       - **INVERSE-ITERATION ground-truth check (2026-07-15):** IRL still not finding small modes -> to settle
         whether $A$ even HAS them, added an env-gated (`INV_ITER`) inverse-iteration block to
         `eig_lanczos_claude.cu`: $v\leftarrow A^{-1}v/\lVert\cdot\rVert$ via the multishift solve on the SAME
         frozen-window $D_{ov}$, Rayleigh quotient $\to$ smallest eig($A$). Handoff
         `tmp_eig_inviter_L2_claude.sh` (config 740, $\lambda_{min}^{W}=0.0643$ inside the frozen 0.06 window).
         $\to\sim0.05$ means the mode is REAL (IRL bug); staying $>1$ means the operator has no small modes.
       - Current IRL params being tuned: $k=4$, $n=8$, $m=16$, $\beta=0.5$.

   1d. **DENSE $D_{ov}$ complex spectrum, thermalized L1** (`ckpoint_lat.2105`, `-DSPECTRUM_DOV` in
       `eig_wmass_val_claude.cu` -> `mult`, complex geev; N=3072 fits in ~1 GB). Plot
       `eig_dov_spectrum_claude.png` (`eig_dov_spectrum_plot_claude.py`). **RESULT: L1 overlap is GAPPED.**
       All 3072 eigenvalues lie EXACTLY on the GW circle $|z-1|=1$ but populate ONLY the RIGHT arc
       (top $z\approx1+i$ -> $z=2$ -> bottom $z\approx1-i$); the LEFT half toward $z=0$ is EMPTY.
       **smallest $|z|=1.466$, ZERO modes with $|z|<0.2$.** Direct visual confirmation that $D_{ov}$ has NO
       near-zero modes at L1 -> nothing to deflate/LMA on the overlap here. The near-zero tail is a WILSON
       phenomenon (empty left arc = no Wilson modes near a zero crossing). **-> the deflation/LMA target is
       likely the WILSON operator $D_W^\dagger D_W$ (inner-solve deflation), NOT the overlap.** DECISION PENDING.
   1e. **KEY PHYSICS (per NM): $D_{ov}$ is a UNITARY PROJECTOR of $D_W$.** $D_{ov}=1+V$ with $V$ the unitary/
       sign part -> the overlap keeps only each $D_W$ eigenvalue's PHASE and discards $|\lambda|$. So a small
       $|D_W|$ singular value (the "0.2 at gsq8") does NOT map to small $|z|$ -- it lands on the GW circle by
       its phase. This RESOLVES the "0.2 vs 1.47" oddity (no bug). Reference dense driver:
       `saved_scripts/eig.cu` (its non-overlap branch geev's the bare Wilson via `DiracPf::mult`). Added
       `-DSPECTRUM_DW` to `eig_wmass_val_claude.cu` (binds the overlap's public `M_DW` CSR) -> side-by-side
       $D_W$ vs $D_{ov}$ plot (`eig_dw_dov_compare_claude.py`, handoff `tmp_eig_dw_dov_L1_claude.sh`,
       gsq8 `ckpoint_lat.2105`). **RESULT (`eig_dw_dov_compare_claude.png`):** $D_W$ eigenvalues form the
       classic Wilson RING (centered $\approx(5.4,0)$, Re$\in[0.34,10.4]$, entirely RIGHT half-plane);
       smallest $|\lambda_{DW}|=0.373$ at the ring's left edge (near origin, on the real axis, Re$>0$).
       $D_{ov}$ = phase-projection -> only the RIGHT arc; the 0.373 mode projects to the TOP ($|z|=1.466$),
       NOT to $z=0$. CONCLUSION: overlap near-zero modes are TOPOLOGICAL (need $D_W$ eigenvalues near a
       zero-crossing / negative-real direction); this config's ring never reaches there -> $D_{ov}$ gapped,
       no overlap LMA to exploit at L1. The small-$|\lambda_{DW}|$ tail is magnitude-based (inner-solve cost),
       a DIFFERENT quantity from overlap low modes.
2. Measure the L4 **low-mode density** (rerun IRL on an L4 near-zero-tail config; $\beta$ small) -> sets Nev.
3. LMA for the **condensate** (exact low + projected high), unbiasedness gate, variance-ratio measurement.
4. LMA for the vector/axial **disc**; optional deflated CG for conn.

## Refs

- Grid suite: `/mnt/baracuda_14/grid_claude/Grid/examples/disc_lma_v2_common_claude.h`,
  `disc_lma_impl_plan_claude.md`, `disc_lma_estimator_bench_impl_plan_claude.md`,
  `disc_mrhs_deflation_impl_plan_claude.md`, `disc_mrhs_defl_bench_results_claude.md`,
  `disc_lma_HANDOFF_claude.md`.
- qed3 seed: `eig_wilson_lowmode_claude.cu` (lambda_min/max scan, `compute_lambda_max` power iteration).
- Algorithm: Sorensen IRL (implicitly restarted Lanczos); Chebyshev filter (even degree, auto-order).

## Status

(a) done -> **viable L4 only**. **IRL ADOPTED** (`includes/lanczos_claude.h`, pre-existing) -- NOT written
from scratch. Convention SETTLED vs `lanczos-2.pdf` Eq (3.3): $\alpha$=spectral bound $2+|m|$, $\beta$=wanted
low-band edge, EVEN degree; full re-orth (`orth_period=1`). Free-field $N_k=1$ check in progress; interacting
multi-mode handoff ready (`tmp_eig_validate_interacting_claude.sh`). Files: `eig_lanczos_claude.cu` (IRL
driver, `CONFIG_LAT` env), `eig_wmass_val_claude.cu` (dense geev of $A$, `CONFIG_LAT`),
`compare_eig_validate_claude.py`. Blocking work (separate) is closed: `blocking_impl_plan_claude.md`.
