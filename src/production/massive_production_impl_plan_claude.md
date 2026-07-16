# Massive overlap production -- implementation plan (2026-07-15, NM)

Massive Nf2 overlap-fermion QED3 HMC on S2 x time, extending the corrected-gsq massless production
([[project_tuning_corrected_gsq]]) to physical masses $m \in \{0.1, 0.5, 1.0, 1.5\}$ (R=1 units).
**NM scheme: "c += rescaled mass" -- shift aux ladder coeffs up by resc, preserving the massless-tuned
additive gaps. NO retuning. Lives in `src/production/`.**

## Goal / physics
- Physical mass enters as the MEASURE-WEIGHTED diagonal $m_L = \text{mass\_coeff}\cdot\text{diag}(A_y/\bar A)$,
  $\text{mass\_coeff} = m\,\bar A/\bar a_s = m\cdot\text{mean\_dual\_area}/\text{mean\_ell}$
  (`includes/overlap_wmass_claude.h:160-163, 269, 308-310`). Driver takes `mass_re`/`mass_im`; frame 0 =
  physical mass. Massive FORCE = validated $(1+m_L)\eta$-through-resolvent ([[project_measure_weighted_mass]]).
- Frozen Zolotarev windows are MASS-INDEPENDENT (act on $M_{DW}=D_W-M$, $M=1$; mass added after the sign
  function) -> UNCHANGED (L1 (0.1,13), L2 (0.06,8), L4 (0.008,5)).

## Scope (locked with NM)
- Couplings = LARGEST corrected gsq per L: **L1 gsq1.5, L2 gsq3.0, L4 gsq6.0**.
- Masses **{0.1, 0.5, 1.0, 1.5}**, Nf = 2.  -> 3 L x 4 m = **12 ensembles**.
- Target sample sizes (trajectories/ensemble): **L1 -> 120, L2 -> 80, L4 -> 60**.

## Scheme: "c += rescaled mass" (the driver shift)
After `base`, the driver does `masses[0]=mass; for i>=1: masses[i] += resc`, resc $= m\,\bar A/\bar a_s$
($\bar A/\bar a_s$ = 0.946 / 0.506 / 0.259 for L1/L2/L4). Ladder in coeff space = {resc, c_1+resc, ...}:
frame 1 is ALWAYS heavier than frame 0 by exactly the massless gap $c_1$, for ANY mass. This makes the old
"$\text{mass\_coeff}<c_1$" rung constraint MOOT (the fixed ladder would have broken at L1 m=1.0: frame0 0.946
vs c_1 1.0; m=1.5: 1.419 > 1.0). A heavier mass only better-conditions the frames. **Force-validated at L1**
(masses {0.1,0.5,1.0,1.5}, grad-vs-FD ~1e-8, reference + block grad; `test_hasenbusch_force_massive_l1_claude.cu`).

## Parameters (reused from massless, `includes/hasenbusch_ladder_claude.h`)
L1 {0->m, 1.0} steps {2,2} tau 1.0; L2 {0->m, 1.0} steps {3,3} tau 1.0; L4 {0->m, 0.4, 1.0} steps {4,4,4}
tau 1.2. Frame 0's coefficient (was 0) is now the physical mass_coeff (set by the driver via argv). All else
identical. (NOTE: massless L2 gsq3.0 / L4 gsq6.0 tuning is still settling; the massive runs inherit whatever
is current in the tuned includes -- rebuild if the massless pass finalizes differently.)

## Deliverables (all in src/production/)
1. `run_massive_claude.sh` -- LOCAL MPS-packed launcher; builds the Nf-packed driver per L
   (-DNF=2 -DKMAX=target -DKRNG=1; L4 -DMIXED_FORCE), runs all 12 (L1 gsq1.5, L2 gsq3.0, L4 gsq6.0) x
   {0.1,0.5,1.0,1.5}, packed 2/GPU (gpu0={L1 || L4 m0.1,0.5}, gpu1={L2 || L4 m1.0,1.5}).
2. `params_massive_claude.md` -- FNAL-facing param doc (mass list, c+=resc scheme, targets, invocation,
   output-dir naming).

Driver = `hmc_hasenbusch_block_claude.cu` (Nf2 -> NSTACK=1 no-op; NOW carries the c+=resc shift after `base`).
Output dir (auto): `Nf2_gsq<g>at0.200000nu01.000000mRe<m>mIm0.000000nt128L<L>_hb<aux>/` -- the `mRe<m>` tag
separates masses; the `_hb<aux>` tag now reflects the SHIFTED aux coeffs (c_1+resc, mass-dependent, e.g. L1
m=0.1 -> _hb1.094585). No collision with the massless (mRe0, _hb1.000000) runs. Auto-resumes.

## Status
- Scheme implemented (shift in serial + block drivers) and FORCE-VALIDATED at L1 (all 4 masses).
- 4-traj acceptance smoke (`tuning/run_massive_test_L1_claude.sh`) pending. Run script + docs in src/production/.
- SUPERSEDES the earlier tuning/ copy of this plan (massive now lives in production/).
