# Massive overlap production -- implementation plan (2026-07-15, NM)

Massive Nf2 overlap-fermion QED3 HMC on S2 x time, extending the corrected-gsq massless production
([[project_tuning_corrected_gsq]]) to physical masses $m \in \{0.1, 0.2, 0.3, 0.4\}$ (R=1 units).
**NM decision: NO retuning -- reuse the massless HMC parameters as-is. Lives in `src/production/`.**

## Goal / physics
- Physical mass enters as the MEASURE-WEIGHTED diagonal $m_L = \text{mass\_coeff}\cdot\text{diag}(A_y/\bar A)$,
  $\text{mass\_coeff} = m\,\bar A/\bar a_s = m\cdot\text{mean\_dual\_area}/\text{mean\_ell}$
  (`includes/overlap_wmass_claude.h:160-163, 269, 308-310`). Driver takes `mass_re`/`mass_im`; frame 0 =
  physical mass. Massive FORCE = validated $(1+m_L)\eta$-through-resolvent ([[project_measure_weighted_mass]]).
  NO operator/driver/ladder change.
- Frozen Zolotarev windows are MASS-INDEPENDENT (act on $M_{DW}=D_W-M$, $M=1$; mass added after the sign
  function) -> UNCHANGED (L1 (0.1,13), L2 (0.06,8), L4 (0.008,5)).

## Scope (locked with NM)
- Couplings = LARGEST corrected gsq per L: **L1 gsq1.5, L2 gsq3.0, L4 gsq6.0**.
- Masses **{0.1, 0.2, 0.3, 0.4}**, Nf = 2.  -> 3 L x 4 m = **12 ensembles**.
- Target sample sizes (trajectories/ensemble): **L1 -> 120, L2 -> 80, L4 -> 60**.

## Why no retuning is safe -- mass-coeff vs the Hasenbusch first rung $c_1$
$\bar A = 4\pi/N_\text{sites}$, $\bar a_s = \text{mean\_ell}$: $\bar A/\bar a_s$ = 0.946 (L1) / 0.506 (L2) /
0.259 (L4). mass_coeff $= m\,\bar A/\bar a_s$ (must clear $c_1$ = 1.0 L1/L2, 0.4 L4 for a valid split):

| $m$ | L1 (c1=1.0) | L2 (c1=1.0) | L4 (c1=0.4) |
|---|---|---|---|
| 0.1 | 0.095 | 0.051 | 0.026 |
| 0.2 | 0.189 | 0.101 | 0.052 |
| 0.3 | 0.284 | 0.152 | 0.078 |
| 0.4 | 0.378 | 0.203 | 0.104 |

All clear $c_1$ with margin. Worst split ratio (frame1/frame0) = L1 m=0.4: 1.0/0.378 = **2.6** (healthy).
A light mass gap makes frame 0 BETTER conditioned than massless -> the massless ladder/steps/tau hold.

## Parameters (reused from massless, `includes/hasenbusch_ladder_claude.h`)
L1 {0->m, 1.0} steps {2,2} tau 1.0; L2 {0->m, 1.0} steps {3,3} tau 1.0; L4 {0->m, 0.4, 1.0} steps {4,4,4}
tau 1.2. Frame 0's coefficient (was 0) is now the physical mass_coeff (set by the driver via argv). All else
identical. (NOTE: massless L2 gsq3.0 / L4 gsq6.0 tuning is still settling; the massive runs inherit whatever
is current in the tuned includes -- rebuild if the massless pass finalizes differently.)

## Deliverables (all in src/production/)
1. `run_massive_claude.sh` -- LOCAL MPS-packed launcher; builds the Nf-packed driver per L
   (-DNF=2 -DKMAX=target -DKRNG=1; L4 -DMIXED_FORCE), runs all 12 (L1 gsq1.5, L2 gsq3.0, L4 gsq6.0) x
   {0.1,0.2,0.3,0.4}, packed 2/GPU across gpu0/gpu1.
2. `params_massive_claude.md` -- FNAL-facing param doc (mass list, reused params, targets, invocation,
   output-dir naming).

Driver = `hmc_hasenbusch_block_claude.cu` (Nf2 -> NSTACK=1 no-op). Output dir (auto):
`Nf2_gsq<g>at0.200000nu01.000000mRe<m>mIm0.000000nt128L<L>_hb1.000000/` -- the `mRe<m>` tag separates masses;
no collision with the massless (mRe0) runs. Auto-resumes.

## Status
- Constraint verified (all clear). No code/ladder change. Run script + params md written in src/production/.
- SUPERSEDES the earlier tuning/ copy of this plan (massive now lives in production/).
