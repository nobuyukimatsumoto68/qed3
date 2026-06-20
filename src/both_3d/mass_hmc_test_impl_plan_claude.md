# L>2 massive-HMC test (diagonal m_L) — implementation plan

Standalone HMC-integration test for the measure-weighted diagonal mass
$m_L=\text{diag}(m\,A_y/\bar a_s)$ at $L=2,4$ with **real and imaginary** $m$. Validates the
FULL HMC stack (momentum refresh, pseudo-fermion heatbath, integrator, $\Delta H$, Metropolis)
with the converted operator/force, not just a single operator apply (that is
`test_diag_mass_l1_claude.cu`). New file: `test_hmc_diag_mass_claude.cu`.

Source / conventions: `mass_measure_factor_impl_plan_claude.md`, audit
`mass_measure_audit_handoff_claude.md`. HMC machinery: `includes/hmc.h` (`HMC2`),
`includes/integrator.h` (`MinimumNorm2`), `pseudofermion_claude.h`. Mirrors the driver setup in
`hmc_fermilab_wmass_L2_claude.cu` and the commented force/Hamiltonian checks in
`hmc_fermilab_claude.cu:357-487`.

## HMC2 interface (from includes/hmc.h:28-105) — what the test drives

- public members `U` (Gauge), `pi` (Force/momentum), `pfs` (vec<shared_ptr<PseudoFermion>>),
  `Sg` (Action), `fermion` (OverlapWMass D), `integrator`.
- `double H()` = `0.5*pi.squared_norm() + Sg(U) + sum pf->S()`  (total Hamiltonian).
- `void integrate()` -> `integrator->integrate(U, pi, Sg, fermion, pfs)`.
- `run(r, dH, is_accept, no_reject=false)`: refresh `pi.gaussian(rng)`, heatbath each pf
  (`pf->gen(rng)` + `fermion->precalc_grad`), `h0=H()`, `integrate()`, `h1=H()`, `dH=h1-h0`,
  Metropolis. `no_reject=true` forces acceptance (useful for thermalize / scans).
- momentum negate: `pi = (-1.0)*pi` (friend `operator*(double, Force)`).

## Tests (each per (L in {2,4}) x (m in {0.1, 0.1 i}); m=0 massless baseline optional)

1. **Reversibility (deterministic, decisive).** Refresh `pi` + heatbath pfs (fix `phi`), save
   `U0`, `pi0`; `integrate()`; `pi = -pi`; `integrate()`; compare `U` to `U0`
   (max over links `|sp - sp0|`, `|tp - tp0|`) and `pi` to `-pi0`. Expect ~solver-tol-limited
   (`~1e-7..1e-9`, NOT machine, since the force carries CG residual). A wrong diagonal force
   (e.g. dropped `M_mass` weight) breaks reversibility grossly.
2. **$\Delta H$ scaling (deterministic).** Same fixed (`pi`,`phi`); restore `U=U0` and run
   `integrate()` with `nsteps` and `2*nsteps` (rebuild `MinimumNorm2`); report `|dH|`. MN2 is a
   symmetric 2nd-order integrator => `dH ~ O(tau^2)`, so doubling `nsteps` should cut `|dH|` ~4x.
   Print the ratio. (Diagonal-mass force enters here through the fermion `dSf`.)
3. **Trajectory sniff (statistical, light).** `hmc.run()` for `K~8` trajectories; report per-traj
   `dH`, acceptance, and `<exp(-dH)>` (should be ~1, unbiased Metropolis). Not a tight gate; a
   broken force shows up as runaway `dH` / 0% acceptance.

## Knobs — RESOLVED (NM, 2026-06-19)

- `N_REFINE` via `-DLREF=2|4`.
- **`Nt` = 16** (cheap; `m_L` is t-uniform so integrator behavior is Nt-independent in character).
- **mass set = `{ (0,0), (0.1,0), (0,0.1) }`** (massless baseline + real + imaginary). Machinery
  test, not production physics.
- **integrator = `MinimumNorm2(tmax=1.9, nsteps, nsteps_inner=100)`**; scaling check uses
  `nsteps` and `2*nsteps` (L2: 9/18, L4: 12/24).
- **gauge = gaussian width 0.3, NO thermalization** (each mass restarts from the same `Uinit`).
- **npole = 17 (L2) / 13 (L4)** (match the drivers). `Nf=2` (n_pf=1).
- **force path = `#define GRAD_L4`** (the production force is exercised inside the HMC loop).

## Build / run

- Local: `nvcc` flags as `tmp_claude.sh` (`-arch=sm_70 ... -std=c++20 ...`), `-DLREF=2|4`.
- Add phases to a runner (extend `tmp_claude.sh` or a new `tmp_hmc_claude.sh`): L2/L4 x build+run,
  tee to `test_hmc_diag_mass_claude.log`. Claude does NOT compile/run — handoff to NM.

## Ordered chunks — ALL WRITTEN (2026-06-19); PENDING build+run proof

1. **[DONE]** Scaffold: setup + per-(L,m) loop. *`test_hmc_diag_mass_claude.cu`.*
2. **[DONE]** Test 1 reversibility (forward + flipped-pi backward; max|U-U0|, |pi+pi0| < 1e-5).
3. **[DONE, refined]** Test 2 dH-scaling (ratio > 2.5, |dH| < 5).
   FINDING (2026-06-19 run, Nt=16/tmax=1.9): massive cases gave ratio ~1.3 (m=0.1) / ~2.4 (m=0.1i)
   -> FAIL, massless 3.90 -> PASS. This is the integrator NOT yet in the asymptotic tau^2 regime for
   the massive force (larger higher-order terms), NOT a force bug: L=1 force-vs-FD ~1e-8 for ALL m,
   reversibility ~1e-10, acceptance fine. FIX (NM-confirmed): SMALL traj length `tmax_sc=0.5`,`ns_sc=4`
   for the scaling test only -> clean ~4x for every m. Also Nt=4 and NPARALLEL_DUPDATE=1 (omp/streams=1).
4. **[DONE]** Test 3 trajectory sniff (K=8 `hmc.run`; <exp(-dH)> in [0.5,2], accept > 0).
5. **[DONE]** Runner `tmp_hmc_claude.sh` (build+run L=2, L=4; tee `test_hmc_diag_mass_claude.log`).

PENDING: NM runs `bash tmp_hmc_claude.sh`; Claude reads the log back. Expected: all PASS;
reversibility ~1e-7..1e-9, dH ratio ~4, <exp(-dH)>~1. (Tolerances are deliberately loose;
a force bug would show as gross reversibility failure / runaway dH.)
