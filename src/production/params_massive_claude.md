# Production parameters -- MASSIVE Nf2 (corrected gsq, largest per L), 2026-07-15

Massive overlap-fermion QED3 HMC on S2 x time. Extends the massless production (`params_L1L2_claude.md`) to
four physical masses. **HMC parameters are REUSED from the massless tuning AS-IS -- no retuning (per NM).**

## Physics / lattice
- Overlap fermion (frozen Zolotarev sign approx), Wilson kernel, M5 = -1.0, Nt = 128, at = 0.2, nu0 = 1.0.
- Nf = 2. Sites = 10 L^2 + 2, links = 30 L^2.
- **Physical mass enters as the MEASURE-WEIGHTED diagonal** $m_L = \text{mass\_coeff}\cdot\text{diag}(A_y/\bar A)$,
  $\text{mass\_coeff} = m\,\bar A/\bar a_s = m\cdot\text{mean\_dual\_area}/\text{mean\_ell}$ (NOT a scalar mass).
  The driver sets frame 0 = physical mass from argv `mass_re`/`mass_im`; the massive force is the validated
  $(1+m_L)\eta$-through-resolvent version. Integrator `MinimumNorm2ML`, MG = 20.

## Scope
- Physical masses **m = 0.1, 0.2, 0.3, 0.4** (R=1 units), real (mass_im = 0).
- Couplings = LARGEST corrected gsq per L: **L1 gsq1.5, L2 gsq3.0, L4 gsq6.0**.
- 3 L x 4 m = **12 ensembles**. Target trajectories: **L1 -> 120, L2 -> 80, L4 -> 60**.

## HMC parameters (reused from massless; `includes/hasenbusch_ladder_claude.h`)
| L | ladder (frame0 = physical mass) | MD steps | tau |
|---|------|------|------|
| 1 | {m, 1.0} | {2, 2} | 1.0 |
| 2 | {m, 1.0} | {3, 3} | 1.0 |
| 4 | {m, 0.4, 1.0} | {4, 4, 4} | 1.2 |
Frame 0's coefficient (was 0 for massless) is now $\text{mass\_coeff} = m\,\bar A/\bar a_s$, set by the driver.
Frozen windows / pole counts UNCHANGED (mass-independent; act on $M_{DW}=D_W-M$): L1 (0.1,13) n25/n11,
L2 (0.06,8) n25/n11, L4 (0.008,5) n31/n11.

## Why no retuning is safe -- mass_coeff clears the first rung $c_1$
$\bar A/\bar a_s$ = 0.946 (L1) / 0.506 (L2) / 0.259 (L4). mass_coeff must clear $c_1$ (1.0 L1/L2, 0.4 L4):

| m | L1 (c1=1.0) | L2 (c1=1.0) | L4 (c1=0.4) |
|---|---|---|---|
| 0.1 | 0.095 | 0.051 | 0.026 |
| 0.2 | 0.189 | 0.101 | 0.052 |
| 0.3 | 0.284 | 0.152 | 0.078 |
| 0.4 | 0.378 | 0.203 | 0.104 |

All clear $c_1$ with margin; worst split ratio (frame1/frame0) = L1 m=0.4: 1.0/0.378 = 2.6 (healthy). A light
mass gap makes frame 0 better conditioned than massless -> the massless ladder holds.

## Driver + invocation
- Source: **`hmc_hasenbusch_block_claude.cu`** (Nf-packed; Nf2 -> NSTACK=1 no-op). Build per L:
  ```
  nvcc -arch=sm_70 -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp \
    -I./includes/ <eigen/hdf5> -DLREF={1|2|4} -DNF=2 -DKMAX={120|80|60} [-DMIXED_FORCE (L4)] \
    hmc_hasenbusch_block_claude.cu <hdf5/gsl libs> -o hmc_massive_L{1|2|4}.out
  ```
  L4 uses `-DMIXED_FORCE` (force-only mixed precision, ~1.35x, exact by Metropolis).
- Run: `./hmc_massive_L{L}.out <gsq> 2 1.0 <mass_re> 0.0`   e.g. `./hmc_massive_L1.out 1.5 2 1.0 0.2 0.0`
- Output dir (auto, AUTO-RESUMES): `Nf2_gsq<g>at0.200000nu01.000000mRe<m>mIm0.000000nt128L<L>_hb1.000000/`
  -- the `mRe<m>` tag separates masses (no collision with the massless mRe0 runs). `-DKMAX` = target size
  (incl. thermalization); extend by rebuilding with larger KMAX (same dir resumes).
- Local packed launcher provided: `run_massive_claude.sh` (MPS 2/GPU, builds all three L + runs the 12).
- Design doc: `massive_production_impl_plan_claude.md`.
