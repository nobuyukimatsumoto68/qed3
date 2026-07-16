# Production parameters -- L=1, 2 (corrected massless gsq range), 2026-07-15

Massless overlap-fermion QED3 HMC on the simplicial sphere S2 x time. The earlier gsq8-era ensembles
were too strong; these are the CORRECTED (weaker) couplings, with HMC params tuned locally (`src/tuning/`).
**L=1 and L=2 are ready to generate. L=4 is still being tuned locally -- do NOT run L4 from these yet.**

## Physics / lattice (all L)
- Discretization: overlap fermion (Zolotarev rational sign approx), Wilson kernel `DiracExt<S2Simp>`, M5 = -1.0.
- Nf: HMC params were tuned at Nf2, but production will also run Nf > 2 (Nf = 4, 6). Use the Nf-packed
  block driver (below) for all Nf; it is bit-identical to the serial driver at Nf2. Massless: physical mass = 0.
- Nt = 128, at = 0.2, nu0 = 1.0. Sites = 10 L^2 + 2, links = 30 L^2.
- Integrator: `MinimumNorm2ML` (multi-timescale Omelyan), gauge sub-steps MG = 20.
- Two-operator split-pole force: ACTION op D (accurate, heatbath + accept/reject) vs FORCE op Df (cruder,
  MD force only) -- exact by Metropolis.

## Couplings to generate
- **L1: gsq = 0.5, 1.0, 1.5**
- **L2: gsq = 1.0, 2.0, 3.0**
  (L4 gsq {2.0, 4.0, 6.0} pending local tuning.)

## Target sample sizes (trajectories, per gsq)
| L | target N |
|---|------|
| 1 | 2000 |
| 2 | 1200 |
| 4 | 400 (pending; not for production yet) |

Set `-DKMAX` at or above the target (default 1200); the driver auto-resumes, so a run can be extended by
rebuilding with a larger KMAX and re-launching against the same output dir.

## Hasenbusch ladder / steps / tau  (`includes/hasenbusch_ladder_claude.h`)
| L | ladder (M_mass coeffs, frame0=physical=0) | MD steps | tau (s_tot) |
|---|------|------|------|
| 1 | {0, 1.0}      | {2, 2} | **1.0** |
| 2 | {0, 1.0}      | {3, 3} | **1.0** |

## Zolotarev windows / pole counts  (`includes/frozen_window_claude.h`, `hasenbusch_ladder_claude.h`) -- KEPT AS-IS
| L | frozen window (lmin, lmax) | n_action | n_force | force window [2*lmin, lmax] |
|---|------|------|------|------|
| 1 | (0.1, 13)  | 25 | 11 | [0.2, 13] |
| 2 | (0.06, 8)  | 25 | 11 | [0.12, 8] |

The frozen windows are CONSERVATIVE for these weak couplings (measured L1 spectrum below sits well inside
(0.1, 13)); zero "eval below Zolotarev window" warnings in the local runs. Left unchanged per NM.

## L=1 Wilson D_W spectrum (MEASURED on the local tuning configs; free U=0 ref (1.154, 10.82))
| gsq | lambda_min range | lambda_max (max) | k = lmin/lmax |
|-----|------|------|------|
| 0.5 | [0.775, 0.886] | 11.00 | ~0.070-0.081 |
| 1.0 | [0.577, 0.781] | 11.15 | ~0.052-0.070 |
| 1.5 | [0.383, 0.611] | 11.12 | ~0.034-0.055 |
(L2 spectrum not yet measured.)

## Driver + invocation  (Nf-PACKED block driver -- use this for ALL Nf)
- Source: **`hmc_hasenbusch_block_claude.cu`** (this dir). The Nf-BLOCK (mrhs) driver: the Nf/2 = NSTACK
  pseudofermion flavors are packed into ONE block per frame and driven through `BlockedMat`/`BlockedForce`.
  Unified for all Nf >= 2: **Nf=2 => NSTACK=1 (a no-op block, bit-identical to the serial driver);
  Nf>=4 => genuine mrhs packing** (~2.5x force at L4/Nf6, parity-validated NSTACK 1 & 3). Nf is COMPILE-TIME
  (`-DNF=<Nf>`, and NSTACK = Nf/2 is a template arg) -- the runtime `Nf` arg MUST match `-DNF` (asserted).
  Build per (L, Nf):
  ```
  nvcc -arch=sm_70 -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp \
    -I./includes/ <eigen/hdf5 includes> -DLREF={1|2} -DNF=<Nf> [-DKMAX=<max traj>] [-DKRNG=<rng-thin stride>] \
    hmc_hasenbusch_block_claude.cu <hdf5/gsl libs> -o hmc_L{1|2}_Nf<Nf>.out
  ```
  - `-DNF` = number of flavors (must equal the runtime `Nf` arg; even, >= 2). `-DKMAX` bounds trajectory
    count (default 1200; set >= the target sample size, e.g. -DKMAX=2000 for L1). `-DKRNG` = keep-every-Nth
    rng (default 5; =1 keeps all).
- Run (massless): `./hmc_L{1|2}_Nf<Nf>.out <gsq> <Nf> <nu0> <mass_re=0> <mass_im=0> [max_sec]`
  e.g. `./hmc_L1_Nf2.out 0.5 2 1.0 0.0 0.0`  (Nf must match -DNF)
- Output dir (auto-created, AUTO-RESUMES): `Nf<Nf>_gsq<gsq>at0.200000nu01.000000mRe0.000000mIm0.000000nt128L<L>_hb1.000000/`
  writes `ckpoint_lat.<k>` (+ `ckpoint_rng.<k>` thinned by KRNG). Identical naming to the serial driver, so
  a block-driver run resumes any existing dir.
- Geometry: reads `../../geometry/data/` (relative; must exist). Needs CUDA 12.6+, GCC 13.2+.
- (The serial one-flavor `hmc_hasenbusch_claude.cu` is kept in this dir as reference; the block driver at
  Nf=2 is equivalent, so prefer the block driver everywhere.)

## Local tuning status (context)
- L1 acceptance (tau=1.2): 0.74 / 0.91 / 0.88 (gsq 0.5/1.0/1.5); dH O(0.1), floor-free. tau lowered 1.2->1.0
  to be safe; re-check acceptance at tau=1.0.
- L2 gsq1.0 (tau=1.2): dH O(0.01-0.07), acceptance ~0.95+. tau lowered to 1.0.
- Nf-block driver `hmc_hasenbusch_block_claude.cu` (`-DNF=<Nf>`) is the PRODUCTION driver for all Nf (mrhs
  force packing; ~2.5x force at L4/Nf6). At Nf2 it is a no-op NSTACK=1 block (bit-identical to serial).
- `-DMIXED_FORCE` is a LOCAL L4-only perf option (mixed-precision force); NOT for L1/L2 production.
