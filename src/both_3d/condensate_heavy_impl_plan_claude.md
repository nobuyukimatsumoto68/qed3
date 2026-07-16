# Condensate for heavy masses (L=1,2,4) + notebook plots — implementation plan

## Goal
Measure **only the $\sigma_{PS}$ condensate** (no ylm current correlators) for the three
**heavy** sea masses on all available L, and add those points to the condensate `.ipynb` plots.

## Physics / quantity
The disc stochastic driver already computes the condensate as a byproduct of the single solve
$\phi = D_m^{-1}\eta$ (shared with the disc current):
$$
\mathrm{etadag\_xi} = -\sum_d \eta_A^\dagger \phi,\qquad \eta_A = A_n\,\eta,
$$
$$
\mathrm{dens} = \frac{\mathrm{etadag\_xi}}{N_t\,(4\pi)},\qquad
\sigma_{PS} = \mathrm{dens} + \overline{\mathrm{dens}},
$$
contact-subtracted ($+2$, `condensate_contact_massive_claude.md` Sec 10) $\to$ SSB order parameter.
"Only the condensate" = keep the $\phi$ solve + this scalar trace, **drop** the per-site ylm current
accumulation and its (large) h5 output.

## Available heavy ensembles (all Nf2, gsq8, at0.2)
| CLI mass_re | dir mRe | eff L1 | eff L2 ($\times0.506305$) | eff L4 ($\times0.259021$) | ncfg L1/L2/L4 |
|---|---|---|---|---|---|
| 0.4228996588195 | 0.422900 | 0.400 | 0.2141 | 0.1095 | 319 / 319 / 103 |
| 0.845799317639  | 0.845799 | 0.800 | 0.4282 | 0.2191 | 319 / 319 / 109 |
| 1.2686989764584 | 1.268699 | 1.200 | 0.6423 | 0.3286 | 319 / 319 / 120 |
(eff = effective fermion mass = the physically comparable $m_F$; from `heavy_mass_L124_impl_plan_claude.md`.)

## Files to modify / create
- **NEW** `condensate_stoch_claude.cu` — condensate-only driver, derived from
  `jj_local_ylm_disc_stoch_claude.cu` (strip ylm-current accumulation; keep $\phi$ solve +
  `etadag_xi` trace + h5). Writes `h0/condensate/etadag_xi` (same key the nb reads).
- **NEW** `run_condensate_heavy_claude.sh` — handoff build+run over the 9 ensembles (USER runs).
- **EDIT** `jj_ylm_condensate_allNf_claude.ipynb` (L1) — add heavy Nf2 points.
- **EDIT** `jj_ylm_condensate_allNf_L2_claude.ipynb` (L2) — add heavy Nf2 points.
- **NEW** `jj_ylm_condensate_allNf_L4_claude.ipynb` (L4) — first creation (mirrors L1/L2).

## Ordered chunks
1. **Condensate-only driver** (`condensate_stoch_claude.cu`). Files: the new .cu.
2. **Measurement handoff** (`run_condensate_heavy_claude.sh`, tee log). Files: the new .sh. USER runs.
3. **Notebook edits** L1 + L2 (surgical: add a heavy-mass loader cell + a heavy series on the
   existing plot, Nf2-only). Files: the two .ipynb.
4. **L4 notebook** (new) + light L4 if its disc data exists, else heavy-only. Files: L4 .ipynb.

## RESOLVED DECISIONS (2026-07-13)
1. **Measurement** = reuse the existing disc driver `jj_local_ylm_disc_stoch_claude.cu` (full disc
   calc; the condensate is the byproduct we keep). Run **disc only** (no conn). Output drops into the
   standard `data_<ens>_vmRe<m>vmIm0.000000/corr_ylm_disc_tb2_nhits1/` -> nb glob works unchanged.
2. **x-axis = effective $m_F$** (the operator's mass_coeff = physical $\times$ measure$_L$). Family-B
   (light) was rescaled so its effective $m_F$ = 0.2/0.1/0.05/0.01 **at every L** (constant -> lines
   up across L). Heavy was NOT rescaled -> effective differs per L (table). So:
   - Light (all L): $m_F$ = 0.2, 0.1, 0.05, 0.01. (L1 nb already this; L2 nb currently plots input
     mRe 0.211450 -> **FIX to 0.2/0.1/0.05/0.01**.)
   - Heavy: L1 {0.40,0.80,1.20}, L2 {0.2141,0.4282,0.6423}, L4 {0.1095,0.2191,0.3286}.
   - **CONFIRM with user** the exact light-vs-heavy scale (heavy uses eff$_L$ = phys/1.05725 at L1;
     light L1 0.2 is dir-mRe). Present the x-table before editing notebooks.
3. **Plot layout** = same axes (x auto-extends; heavy L1 reaches 1.2). Notebooks: EDIT L1 + L2, CREATE L4.
4. **Contact subtraction** = $+2$, **mass-INDEPENDENT** (confirmed `condensate_contact_massive_claude.md`
   lines 36,85: real-$m_F$ $\sigma_{PS}$ contact $=\tfrac12$ both legs, identical to massless). Valid at heavy.
5. **Heavy series styling** = Nf2-only; distinct marker AND color (e.g. red filled diamond 'D') vs the
   light Nf2 'o'. (color-blind rule.)

## Proposed plot x-values (CONFIRM before nb edits)
| $m_F$ (effective) | L1 | L2 | L4 |
|---|---|---|---|
| light Family-B | 0.2, 0.1, 0.05, 0.01 | 0.2, 0.1, 0.05, 0.01 | 0.2, 0.1, 0.05, 0.01 |
| heavy | 0.40, 0.80, 1.20 | 0.2141, 0.4282, 0.6423 | 0.1095, 0.2191, 0.3286 |
