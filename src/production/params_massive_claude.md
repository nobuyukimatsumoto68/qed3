# Production parameters -- MASSIVE Nf2 (corrected gsq, largest per L), 2026-07-15

Massive overlap-fermion QED3 HMC on S2 x time. Extends the massless production (`params_L1L2_claude.md`) to
four physical masses. **HMC parameters are REUSED from the massless tuning AS-IS -- no retuning (per NM).**

## Physics / lattice
- Overlap fermion (frozen Zolotarev sign approx), Wilson kernel, M5 = -1.0, Nt = 128, at = 0.2, nu0 = 1.0.
- Nf = 2. Sites = 10 L^2 + 2, links = 30 L^2.
- **Physical mass enters as the MEASURE-WEIGHTED diagonal** $m_L = \text{mass\_coeff}\cdot\text{diag}(A_y/\bar A)$,
  $\text{mass\_coeff} = m\,\bar A/\bar a_s = m\cdot\text{mean\_dual\_area}/\text{mean\_ell}$ (NOT a scalar mass).
  The driver sets frame 0 = physical mass from argv `mass_re`/`mass_im`; the massive force is the validated
  $(1+m_L)\eta$-through-resolvent version. Integrator `MinimumNorm2ML`, gauge substeps MG = 20 (L1/L2), 100 (L4).

## Scope
- Physical masses **m = 0.1, 0.5, 1.0, 1.5** (R=1 units), real (mass_im = 0).
- Couplings = LARGEST corrected gsq per L: **L1 gsq1.5, L2 gsq3.0, L4 gsq6.0**.
- 3 L x 4 m = **12 ensembles**. Target trajectories: **L1 -> 120, L2 -> 80, L4 -> 60**.

## HMC parameters -- "c += rescaled mass" scheme (reused massless gaps; `includes/hasenbusch_ladder_claude.h`)
The AUXILIARY ladder coefficients are shifted UP by the physical rescaled mass
$\text{resc} = m\,\bar A/\bar a_s = m\cdot\text{mean\_dual\_area}/\text{mean\_ell}$ (applied in the driver after
`base`). Frame 0 stays the physical mass (set_mass rescales internally); frame $i\ge1$ coeff $= c_i + \text{resc}$.
| L | ladder (coeff space) | MD steps | tau |
|---|------|------|------|
| 1 | {resc, 1.0 + resc}          | {2, 2} | 1.0 |
| 2 | {resc, 1.0 + resc}          | {3, 3} | 1.0 |
| 4 | {resc, 0.4 + resc, 1.0 + resc} | {4, 4, 4} | 0.8 (MG100) |
The ADDITIVE inter-frame gaps ($c_i$, massless-tuned) are PRESERVED at any mass. Frozen windows / pole counts
UNCHANGED (mass-independent; act on $M_{DW}=D_W-M$): L1 (0.1,13) n25/n11, L2 (0.06,8) n25/n11, L4 (0.008,5) n31/n11.

## Why no retuning + no rung constraint
$\bar A/\bar a_s$ = 0.946 (L1) / 0.506 (L2) / 0.259 (L4). With the shift, frame 1 is ALWAYS heavier than
frame 0 by exactly $c_1$ (the massless gap), so the split is valid for ANY mass -- the old
"$\text{mass\_coeff}<c_1$" constraint is MOOT (it would have broken at L1 m=1.0/1.5). A heavier mass gap only
better-conditions the frames. **Force-validated at L1** (all masses {0.1,0.5,1.0,1.5}, grad-vs-FD ~1e-8,
reference + block grad): `test_hasenbusch_force_massive_l1_claude.cu`.

## Driver + invocation
- Source: **`hmc_hasenbusch_block_claude.cu`** (Nf-packed; Nf2 -> NSTACK=1 no-op). Build per L:
  ```
  nvcc -arch=sm_70 -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp \
    -I./includes/ <eigen/hdf5> -DLREF={1|2|4} -DNF=2 -DKMAX={120|80|60} [-DMIXED_FORCE (L4)] \
    hmc_hasenbusch_block_claude.cu <hdf5/gsl libs> -o hmc_massive_L{1|2|4}.out
  ```
  L4 uses `-DMIXED_FORCE` (force-only mixed precision, ~1.35x, exact by Metropolis).
- Run: `./hmc_massive_L{L}.out <gsq> 2 1.0 <mass_re> 0.0`   e.g. `./hmc_massive_L1.out 1.5 2 1.0 0.2 0.0`
- Output dir (auto, AUTO-RESUMES): `Nf2_gsq<g>at0.200000nu01.000000mRe<m>mIm0.000000nt128L<L>_hb<aux>/`
  where `<aux>` = the SHIFTED aux coeffs (mass-dependent; L1 m=0.1 -> `_hb1.094585`). `mRe<m>` separates masses.
  -- the `mRe<m>` tag separates masses (no collision with the massless mRe0 runs). `-DKMAX` = target size
  (incl. thermalization); extend by rebuilding with larger KMAX (same dir resumes).
- Local packed launcher provided: `run_massive_claude.sh` (MPS 2/GPU, builds all three L + runs the 12).
- Design doc: `massive_production_impl_plan_claude.md`.
