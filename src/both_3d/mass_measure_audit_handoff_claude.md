# Measure-weighted mass audit — handoff (2026-06-19)

Code audit of the per-site diagonal mass $m_L=\mathrm{diag}(m\,A_y/\bar a_s)$ rollout
(see `mass_measure_factor_impl_plan_claude.md`). Goal of this note: list where the
conversion is complete, the one stale scalar-mass site still live in production, and
independent numerical checks. The header math itself is reviewed and correct; the issue
is a parallel operator path that was NOT converted.

## Conventions

- $m_L = \texttt{mass\_coeff}\cdot M_\text{mass}$, with
  $M_\text{mass}=\texttt{volume\_matrix}(1)=\mathrm{diag}(A_y/\bar A)$ and
  $\texttt{mass\_coeff}=m_\text{phys}\,\bar A/\bar a_s = m\cdot\texttt{mean\_dual\_area}/\texttt{mean\_ell}$.
- Ctor argument `mass` is now the PHYSICAL $m$ (complex). At $L=1$, `volume_matrix(1)=I`,
  so $m_L=\texttt{mass\_coeff}\cdot\mathbb{1}$ (NOT $m\cdot\mathbb{1}$): the effective
  scalar is $m\cdot\bar A/\bar a_s\approx 0.996\,m$ at $L=1$.

## Converted correctly (verified by inspection)

- `includes/overlap_wmass_claude.h`: ctor `:241-244`; `apply_mL`/`apply_mLdag` `:379-389`;
  `mult`/`adj` `:417,:452`; `_ms` variants `:490,:516`; `DHD`/`DDH` (+`_ms`) `:522-608`;
  force `grad`/`grad_l1/l2/l4` fold $(1+m_L^*)$ via `M_mass` `:744-780,871-906,937-974`.
  Squared-op identity $D_m^\dagger D_m=(1+M^*)D+D^\dagger(1+M)+|M|^2$ holds for complex
  diagonal $M$ with no commutator (checked against the `DHD` implementation).
- Zolotarev shifts are mass-free ($\sigma_m=-k^2/cp[m]$, `:476,:507,:677`) — no measure
  factor in the poles, as intended.
- HMC fully on the diagonal path: `pseudofermion_claude.h:19-59` -> `LinOpDHDWrapper(D)`
  -> `DHD_..._ms` + `adj_..._ms`; force -> `grad_*`. BlockedMat is NOT used by HMC.
- Massless ($m=0$): `mass_coeff=0` -> $m_L=0$, both paths agree. Unaffected.

## PRIMARY ISSUE — stale scalar mass in `BlockedMat`  [FIXED 2026-06-19, local Claude]

**RESOLVED.** `BlockedMat::mult/adj/DDH` now apply the diagonal `m_L = mass_coeff*M_mass`,
broadcast over NSTACK via `mult_coo_block(M_mass)` (the `grad_l2` pattern):
- `mult` `:347-350` `+ mass_coeff (M_mass v)`;  `adj` `:362-365` `+ conj(mass_coeff)(M_mass v)`.
- `DDH` `:367-389` REWRITTEN to the diagonal identity `(1+M^*)(D+M)v + (D+M)^dag(1+M)v -
  (M+M^*+|M|^2)v`, mirroring `OverlapWMass::DHD_deviceAsyncLaunch` -- the `(1+M)` is fed INTO
  `adj` (NOT applied as a scalar after), since `[D^dag,M] != 0` for diagonal M (the audit's
  one-line "fix direction" understated this; mult/adj were a drop-in weight, DDH needed the reorder).
- New pool scratch `d_mt`/`d_mw` (N*NSTACK, pole-block-gated) for the DDH `M v` / reusable temps.
- Regression added to `test_diag_mass_l1_claude.cu`: BlockedMat(NSTACK=1) vs `*_ms` ~1e-13 (checks 1-2).
- Contamination check (timestamps): JJ/condensate binaries are 2026-06-13 (PRE-conversion) -> existing
  massive free data used scalar-everywhere = self-consistent, NOT contaminated. Fix is forward-only.
VERIFIED (2026-06-19 rerun, `tmp_claude.sh` 10 phases): BlockedMat(NSTACK=1) vs MatPoly `_ms` = 0.0
(bit-identical) for m = 0, 0.1, 0.1i at L=1/L=2/L=4; all phases ALL PASS, 0 errors.

Original finding (for reference):
`includes/blocked_mat_claude.h` reimplements the operator and applies `D.mass` as a SCALAR:
- `:342` `+ cplx(D.mass) d_xi` (mult);  `:355` `+ cplx(conj(D.mass)) d_xi` (adj)
- `:364-367` DDH old scalar identity `2.0*D.mass.real() + std::norm(D.mass)`

Used in production JJ for the massive operators:
- `jj_corr_dilute_claude.cu:343-344` `blk_tp_Dm(Dm)`, `blk_sp_Dm(Dm)`, `blk_*_Dtil(Dtil)`
  (identical in `jj_corr_block_t_claude.cu`, `_L2`, `_L4`, `_fermilab`).
- drives the mrhs source-leg solves `solve_sq_from_cpu` (e.g. `:585,632,658,745,765,796,835,857`).

Consequence: in one correlator the single solves use `op_Dm/op_Dmsq` (= `Fermion::*_ms`,
correct $m_L$) while the block source legs use `BlockedMat` (scalar) -> the two legs use
DIFFERENT mass operators. Mismatch is present even at $L=1$ (~0.4%, the $\bar A/\bar a_s$
factor) and becomes site-dependent + large at $L>1$. Bites any nonzero-valence-mass JJ
run going forward (real-$m$ ylm family included).

Fix direction: give `BlockedMat::mult/adj/DDH` the same `M_mass`/`mass_coeff` treatment as
`apply_mL`, broadcasting the per-site diagonal across the NSTACK columns — exactly the
pattern `grad_l2` already uses (`mult_coo_block(M_mass, ...)`).

## SECONDARY (parity / $m_P$ path; currently dropped — flag only)

- Tilde-D mass transformed on the scalar before weighting:
  `jj_corr_dilute_claude.cu:295` `Dtil(DW, valence_mass/(1-valence_mass), 11)` gives
  $\mathrm{diag}([m/(1-m)]w)$; the consistent diagonal tilde is $M_L(1-M_L)^{-1}
  =\mathrm{diag}(mw/(1-mw))$. Agree only to $O(m^2 w)$.
- Output $(1+m_P)$ factors are scalar: `jj_corr_dilute_claude.cu:894,901,921` use scalar
  `valence_mass`; for diagonal mass these should be per-site $(1+m_L)^{\pm1}$.

## Provenance note (not a code bug)

Ctor `mass` now means physical $m$; HMC dir names encode it directly
(`hmc_fermilab_wmass_L2_claude.cu:209`, etc.). Old $L=1$ dirs `mRe0.010000` are no longer
the same operator. Existing $L2/L4$ dirs with odd masses (`mRe0.005338`, `mRe0.002705`, …)
look like the SUPERSEDED `mean_ell`-scaling scheme — confirm run scripts pass new physical
$m$ and old-scheme data is not mixed into new-scheme analysis.

## Independent numerical checks (cheap, decisive)

STATUS (2026-06-19, local Claude): checks 1/2/4 are now IN `test_diag_mass_l1_claude.cu`
(BlockedMat NSTACK=1 vs `*_ms`, L=1/L=2, default + `-DGRAD_L4`); 3/5 remain TODO.

1. **[DONE]** BlockedMat vs MatPoly, NSTACK=1: `blk.{mult,adj,DDH}` vs `Dnew.{mult,adj,DHD}_ms`
   on the same random vector; expect ~1e-13 (pre-fix DIFFERS by the `mass_coeff/mass` factor).
2. **[DONE -- equivalent]** the NSTACK=1 op-equality (1) is the operator half; the solve round-trip
   is implied (both legs now share `m_L`). A dedicated block-solve round-trip can be added if wanted.
3. **[TODO]** adj/mult adjointness at $L=2$ (site-varying weights): $\langle u|D_m v\rangle
   =\langle D_m^\dagger u|v\rangle$ via `mult` vs `adj` (catches wrong/missing conj on `mass_coeff`).
   Cheap to add to the test (two random vectors, one `mult` + one `adj` + two dots).
4. **[DONE]** `test_diag_mass_l1_claude.cu` compares `BlockedMat` against the diagonal production at L=1/L=2.
5. **[TODO]** Ward identity / current conservation on a massive $L=2$ config (needs a JJ measurement run;
   defer to the JJ-rerun stage -- now that BlockedMat is fixed it should be SATISFIED, not broken).

## Follow-on test work (requested 2026-06-19)
- **Adjointness + (optionally) round-trip** added to `test_diag_mass_l1_claude.cu` (checks 3/2-dedicated).
- **Separate L>2 massive HMC test `.cu`** (real + imag m), modeled on the existing `hmc_*` massive
  drivers/scripts -- full dH/reversibility/force-vs-FD on a thermalized-ish config at L=2 (and L=4 if cheap).
- **Driver/CLI** physical-`m` plumbing in `hmc_fermilab_wmass_L{2,4}_claude.cu` (see task doc).
