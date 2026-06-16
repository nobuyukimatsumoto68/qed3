# Stochastic Y_lm tower -- handoff for the post-/compact agent

> **SUPERSEDED (2026-06-13).** The "one-end trick" with the SINGLE COMBINED weight
> $W_\ell(n)=A_n\sum_m Y_{\ell m}(n)$ below is **WRONG**: it collapses the m-sum into an OUTER product
> $(\sum_{m_1}Y)(\sum_{m_2}Y)$ = all $(m_1,m_2)$ pairs (convention (b)), NOT the diagonal $\sum_m Y(n)Y(n')$
> with $1/(2\ell+1)$ that Eq.(4.36) requires. It conflated $\ell$-diagonality with $m$-diagonality.
> Verified: the stochastic exact-K ylm in `jj_corr_block_t` (which uses exactly this weight) sits a flat
> factor ~4 above the deterministic at $\ell=2$ (systematic, not noise = the "ell=2 deviation").
> **The CORRECT estimator needs ONE solve per $(\ell,m)$** (not per $\ell$): see the
> `jj_local_ylm_impl_plan_claude.md` Addendum (2026-06-13). Use that plan, not this file's collapsed-weight
> recipe. (Final scope there: per-m, $\ell\le3$, LOCAL VECTOR only, 3 Pauli, deter+stoch.)


Goal: add the **stochastic Y_lm tower** $g_\ell(t)$ ($\ell=0,1,2$, m-summed) for the **LOCAL** current
(vector **and** axial), validated against the deterministic local ylm.  Read `MEMORY.md` ->
`project_state.md` (top) + `project_axial_trio_claude.md` first.

## Where things stand (don't re-derive)
- The conserved-current correlator pipeline is in `src/both_3d`.  The active stochastic measurement program is
  **`jj_corr_dilute_claude.cu`** (master-field superposition + spin/time dilution; `--spin-dilution`,
  `--time-dilution td`; spin-only = `--spin-dilution --time-dilution 1`).  It computes: exact-K conn tp/sp +
  disc + axial tp/sp, and LOCAL vector s1,s2,s3 + LOCAL AXIAL s1,s2,s3 (`h0/axial/s{c}/Apm`, added 2026-06-13).
- **ylm is STRIPPED in the dilute** (`constexpr int N_ELL = 0`): the m-summed ylm machinery (`W_ell[l][n]`,
  `srcL[l]`, `psi_yl[l]`, `PhiLt[l*Nt+t]`, indexers `IYL`/`IPL`, `Ylm_real`) is **left in place but inert**.
  Per the user's decision the ylm is to live in a **SEPARATE one-end-trick file for the LOCAL currents**, NOT
  re-enabled inside the dilute.  (`jj_stoch_local_impl_plan_claude.md` Future-TODOs: `jj_local_ylm_stoch`.)
- The DETERMINISTIC ylm already exists and is the validation ground truth:
  - vector: `corr_deter_{local,disp,exactsum}_*` key `h0/t0_b/ylm/l{0,1,2}/Vpp` (loc via `build_Sigma_ylm`).
  - axial:  `corr_deter_local_axial_L{1,2,4}` key `h0/t0_b/ylm/l{0,1,2}/Apm` (jj_local_axial_deter).
  - cross-checked analytically in `jj_cft_ylm_check_claude.cc` (Eq. 4.35).

## Physics / estimator (qed3int_v2-13, Eq. 4.35/4.36)
m-summed diagonal-$\ell$ Legendre coefficient of the TEMPORAL current:
$$
g_\ell(t)=\frac{1}{2\ell+1}\sum_{m=-\ell}^{\ell}\operatorname{tr}\!\big[\Sigma_{\ell m}(t_0)\,P\,\Sigma_{\ell m}(t)\,P\big],
\qquad \Sigma_{\ell m}(t)=\sum_n A_n\,Y_{\ell m}(\hat n)\,J^{t}(n,t).
$$
For the **local** current $J^t(n)=\sigma_3(n)$ (bare, NO $\kappa$): the weight is $W_\ell[n]=A_n\sum_m Y_{\ell m}(\hat n)$
with $A_n=$ `dual_areas[n]` (this is the m-summed one-end weight; it works because the diagonal-$\ell$ m-sum
factorizes, $\sum_{m_1,m_2}Y_{\ell m_1}(n)Y_{\ell m_2}(n')=W_\ell(n)W_\ell(n')/A_nA_{n'}$ -> ONE solve per $\ell$).
NB the dilute's `W_ell` carries an extra `/kappa_t` (kappa-in, for the EXACT-K ylm, Eq. 4.36) -- for the BARE
LOCAL ylm DROP the `/kappa_t` (use `dual_areas[n]*sum_m Ylm`), matching `jj_local_deter`'s `build_Sigma_ylm`.

**One-end trick (stochastic):** with the diluted source $\eta$ and $\phi=D_m^{-1}\eta$,
$$
\psi_\ell(t_0)=D_m^{-\dagger}\,\Sigma_\ell(t_0)\,\eta,\qquad
g_\ell(t)=\tfrac{1}{2\ell+1}\,\psi_\ell^\dagger\,[\Sigma_\ell(t)\,\phi]
$$
($\Sigma_\ell$ Hermitian; $\Sigma_\ell(t)\phi$ is the Y_lm-weighted $\sigma_3$ sink -- cheap, no extra solve).
This is the local-vector ylm.  **Local AXIAL ylm** = same but the $t_0$-leg propagator is $D_m^{-\dagger}$ (rip
off $(1-D)$, like the local axial we just added): source $\psi^A_\ell=D_m^{-1}\Sigma_\ell\eta$ (solve via
`op_DmH`, cf. the local-axial block), sink $\Sigma_\ell\phi$ shared -> single complex channel.

## Targets (free-field; converge with L)
- decay rates $\ell{=}1\to2$, $\ell{=}2\to3$; amplitude ratio $g_2e^{3t}/g_1e^{2t}\to 12/5=2.4$;
  $g_0\to0$ (charge conservation; icosahedral $\ell{=}6$ aliasing leaks $\sim1/L^4$); $g_1/G_t\to1/3$.
  (Deterministic loc/disp/exact already show this: rates 1.93->1.98, 2.41->2.87, ratio 5.7->3.1->2.4, g0/g1->0
  over L=1,2,4 in `comp_trio_jj_claude.ipynb`.)

## Recommended implementation
1. **New file `jj_local_ylm_stoch_claude.cu`** (one-end trick; copy the dilute's source/solve scaffolding +
   `valence_claude.h` `time_spin_dilution`/`time_dilution`).  Compute the LOCAL vector ylm `g_l` (Vpp) and the
   LOCAL axial ylm `g_l` (Apm, $D_m^{-\dagger}$ t0-leg).  m-summed weight $W_\ell[n]=$ `dual_areas[n]`$\sum_m$`Ylm_real(l,m,base.sites[n])`.
   Reuse the master-field superposition (origins 0, Nt/2) and the spin/time dilution switches.  Output
   `h0/ylm/l{0,1,2}/Vpp` (+ Vmm conj) and `h0/axial/ylm/l{0,1,2}/Apm`.  Per-hit files, atomic + complete-gated.
   - Alternatively (simpler but kappa-in, exact-K): just set `N_ELL` back in the dilute -- but the user wants
     the LOCAL ylm in a separate file, so prefer the new file.
2. **Run script** + **notebook**: extend `jj_dilute_validate_axial_claude.ipynb` (or a new validate nb) with ylm
   loaders `h0/ylm/l{l}` and `h0/axial/ylm/l{l}`, comparing diluted vs `corr_deter_{local,local_axial}` ylm.
   Honest-normalization RULE (the user is strict): ONE constant per estimator for ALL channels (derive from one
   reference, reuse everywhere); raw ratios; NO np.abs / per-curve sign flips on the plotted G's.  Axial: master-
   field HALF-VOLUME shift only (`mirror`, no T-fold), full [0,Nt), plot Re AND Im.

## Files to read
`jj_local_ylm_impl_plan_claude.md`, `axial_ylm_tower_impl_plan_claude.md`, `jj_stoch_local_impl_plan_claude.md`,
`jj_local_stoch_estimator_design_claude.md`, `jj_cft_ylm_check_claude.cc`; the inert ylm code in
`jj_corr_dilute_claude.cu` (N_ELL/W_ell/srcL/psi_yl/PhiLt/IYL/IPL); `jj_local_deter_claude.cu` `build_Sigma_ylm`
(deterministic local-vector ylm reference) and `jj_local_axial_deter_claude.cu` (axial ylm reference).

## Workflow reminders (this user)
- Write/extend a `*_impl_plan_claude.md` before coding; resolve open questions interactively; chunked approval.
- `_claude` suffix on every file Claude creates/edits; no Unicode in code (LaTeX macros bare); no lambdas (prefer
  file-scope/inline); match the original file layout (2-space, one stmt/line); cite the algorithm source.
- NEVER compile/run unilaterally -- hand off via a `tmp_*_claude.sh` that builds (to a DISTINCT binary name if
  another dilute run is live -> avoid ETXTBSY) + tees to a `*_claude.log`; the user runs it.  NO rm/rmdir in any
  .sh (ask the user to rm complete-gated outputs before a rerun).  One plot per cell in notebooks.
