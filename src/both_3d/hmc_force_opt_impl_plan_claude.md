# HMC force optimization -- impl plan (2026-06-06)

## Algorithm sources (cite in code)
- HMC: Duane, Kennedy, Pendleton, Roweth, Phys. Lett. B195 (1987) 216.
- Overlap force = derivative of the Zolotarev rational sign approximation; the inner shifted solves are a
  multi-shift CG: **B. Jegerlehner, hep-lat/9612014** (already used for the action/force precompute).
- Multi-timescale (Sexton-Weingarten) integration: **Sexton & Weingarten, Nucl. Phys. B380 (1992) 665**
  -- ALREADY IMPLEMENTED in `integrator.h` (the `ExplicitLeapfrogML` / ML variants); see L3.

## Goal / state
HMC trajectory (n=21, Nf2 gsq4 free field, first traj `# HMC`): pole-loop 1146.7 s; action-solve
multishift 913.1 s (**1.26x**, validated dH ~1e-8). Force-multishift = WASH (the force is dominated by
`grad`, NOT by `precalc_grad`'s pole solves). So the lever is the FORCE's `grad`, not more solve
multishift. Action-solve multishift + the NaN freeze-fix are DONE; this plan targets `grad`.

## Anatomy: `OverlapWMass::grad_deviceAsyncLaunch` (overlap_wmass_claude.h:542) -- the force bottleneck
Called PER GAUGE LINK (`n_links*Nt` links x ~`2*nsteps+1` force evals/traj). Per call:
- builds a single-link COO `coo` for `dD_W/dU(link)` (`:545-547`, host `d_coo_format` + `do_it`);
- loops poles `m=1..size-1` over `nstreams=4` streams (`:553`), each pole:
  - `X.on_gpuAsync(d_XZs[m], d_Zs[m])`  = `X Z_m`   (`:558`)  -- **LINK-INDEPENDENT** (`X=M_DW/lmax`)
  - `coo.Async(d_XYs[m], d_Ys[m])`       = `coo Y_m`  (`:559`)  -- link-dep
  - dot `<coo Y_m, X Z_m>` (`:562`)
  - `X.on_gpuAsync(d_XYs[m], d_Ys[m])`  = `X Y_m`   (`:566`)  -- **LINK-INDEPENDENT**
  - `coo.Async(d_XZs[m], d_Zs[m])`       = `coo Z_m`  (`:567`)  -- link-dep
  - dot `<X Y_m, coo Z_m>` (`:569`)
- 4 matvecs + 2 dots per pole; the 2 `X` applies are link-independent => recomputed
  `n_links*Nt*nsteps...` times REDUNDANTLY. `precalc_grad` (`:484` pole-loop / `:519` _ms) currently
  computes only `Z_m=R_m eta` (`d_Zs[m]`) and `Y_m=R_m X^dag eta` (`d_Ys[m]`), NOT `X Z_m`/`X Y_m`.
The matvecs are also occupancy-starved single-vector (N=3072 -> 13 of 80 SMs); `nstreams=4` only partly
hides it. The per-pole `cudaStreamSynchronize` after each dot serializes.

## L1 -- HOIST the link-independent `X Z_m`/`X Y_m` out of `grad` [DONE+VALIDATED 2026-06-06]
**DONE.** `d_XZpre`/`d_XYpre` buffers + precompute in BOTH precalc_grad variants; `grad_deviceAsyncLaunch_l1`
(drops the 2 X applies, reads the precomputed vectors); `grad_deviceAsyncLaunch` dispatches to `_l1` under
`#ifdef GRAD_L1` (`-DGRAD_L1`, like FORCE_MULTISHIFT; gauge_ext.h call site untouched). All in
overlap_wmass_claude.h, originals kept. Test `test_grad_l1_claude.cu` (per-link grad vs grad_l1, 336
links, mass=0.1): **max|ref-L1| = 0.000e+00 (BIT-IDENTICAL)**, grad sweep **1.651x** (0.455 -> 0.276 s).
Byte-identical force => HMC dH reproduces EXACTLY (no traj run needed for correctness). DEPLOY = define
GRAD_L1 in hmc_claude.cu; optional end-to-end `# HMC` traj time vs the 913 s baseline.

## L1 (original plan text below) -- HOIST the link-independent `X Z_m`/`X Y_m` out of `grad`
`X Z_m` and `X Y_m` depend ONLY on pole `m` (not the link). Compute them ONCE per force eval in
`precalc_grad` (npole vectors each, into NEW persistent buffers `d_XZpre[m]`/`d_XYpre[m]`); `grad` then
drops the two `X.on_gpu` calls and reads `d_XZpre[m]`/`d_XYpre[m]` in the dots. Removes ~half of `grad`'s
matvecs (the occupancy-starved ones); only the per-link `coo Z_m`/`coo Y_m` matvecs remain. Force is
byte-identical (same arithmetic, fewer recomputes) -> validate dH match (~1e-8) like the force-ms bench.

### L1 chunks (side-by-side; keep originals as the dH reference, per project convention)
- **L1a -- buffers + precompute.** Add `std::vector<CuC*> d_XZpre, d_XYpre` (size `size`), alloc in the
  ctor next to `d_Zs`/`d_Ys`/`d_XZs`/`d_XYs`, free in the dtor. In BOTH `precalc_grad` variants
  (pole-loop `:484` + `_ms` `:519`), after `d_Zs[m]`/`d_Ys[m]` are ready, compute
  `d_XZpre[m] = (1/lmax) M_DW d_Zs[m]` and `d_XYpre[m] = (1/lmax) M_DW d_Ys[m]` (one `X.on_gpu` each,
  m=1..size-1; can batch as a block matvec -- see L2). ADDITIVE: the old `grad` ignores them, so this
  alone changes nothing.  Files: includes/overlap_wmass_claude.h
- **L1b -- new `grad` variant.** Add `grad_deviceAsyncLaunch_l1` = copy of `grad_deviceAsyncLaunch` with
  the two `X.on_gpuAsync` lines REMOVED and the dots reading `d_XZpre[m]`/`d_XYpre[m]`
  (`coo Z_m`/`coo Y_m` still computed into `d_XZs[m]`/`d_XYs[m]` scratch). Original `grad` UNTOUCHED.
  Files: includes/overlap_wmass_claude.h
- **L1c -- call-site toggle + validate.** Switch `pseudofermion_claude.h get_force` (or the `grad`
  call site) to `grad_l1` behind a toggle (e.g. `#define GRAD_L1` in the .cu, like FORCE_MULTISHIFT).
  Benchmark a trajectory `# HMC` time + confirm dH matches the original to ~1e-8. Files:
  includes/pseudofermion_claude.h, hmc_claude.cu (toggle + bench).

## L1 DEPLOY (2026-06-06): GRAD_L1 enabled in ALL 6 OverlapWMass hmc files
`#define GRAD_L1` before the overlap_wmass include (byte-identical force => safe on-by-default):
hmc_claude.cu, hmc_w_mass_claude.cu, hmc_w_mass_check_claude.cu, AND (user: include fermilab)
hmc_fermilab_claude.cu, hmc_fermilab_L2_claude.cu, hmc_fermilab_wmass_claude.cu. NOT hmc_debug_claude.cu
(legacy overlap.h -- no GRAD_L1 dispatch there). NOTE fermilab files are make-excluded (*_fermilab*) so
build them manually. pf-block (L3) = WON'T DO per user 2026-06-06.

## PF-BLOCK (user 2026-06-06) -- batch the n_pf = Nf/2 pseudofermions (Nf>=4 only)
`hmc_claude.cu:242` `for f<Nf/2` => n_pf = Nf/2 (Nf2->1 NO benefit, Nf4->2, Nf6->3). The n_pf
pseudofermions are INDEPENDENT (own action solve, precalc_grad, grad; force = sum_pf). Feasible to batch
as nstack=n_pf RHS, but TWO different axes:
- the INVERSIONS (action op_Dmsq solve + precalc_grad's 2 multishift solves) batch over pf via BlockedMat
  (nstack=n_pf) -- the D1/mrhs lever; nstack=2-3 small => modest (~1.04x like jj source solves).
- `grad`'s big occupancy axis is the POLES (npole~10 = L2), NOT pf (2-3). pf could be an EXTRA block axis
  (npole x n_pf) but L2 dominates.
=> pf-block is COMPLEMENTARY (targets inversions) + composes with L1/L2 (target grad). Needs restructuring
the HMC `for(auto pf : pfs)` loop. LOWER priority than L1/L2 (small nstack). DEFER unless Nf>=4 runs need it.

## L2 -- BLOCK `grad`/precompute over poles (compounds L1; fills the SMs)
After L1, `grad`'s remaining per-link work is `npole` matvecs of the SAME `coo` (`coo Z_m`, `coo Y_m`)
+ `npole` dots -- a block over poles: one block COO matvec (npole RHS, `mult_block`-style over the
single-link CSR) + a fused block-dot, instead of `npole` `nstreams`-streamed occupancy-starved kernels.
Likewise L1a's `X Z_m`/`X Y_m` precompute is `npole` RHS of `M_DW` -> one `mult_block`. Reuses the
BlockedMat block kernels (`mult_block`/`axpy`/block-dot). Replaces the `nstreams=4` pole streaming.
Files: includes/overlap_wmass_claude.h (+ maybe a single-link block-COO matvec helper).

## L3 -- multi-timescale integration (ALREADY IN integrator.h; ENABLE/TUNE, not build)
`hmc_claude.cu:298` currently uses `MinimumNorm2(tmax, nsteps, 100)` (nsteps 8-10). The ML integrator
(`ExplicitLeapfrogML`, commented at `:451,481,511`) already implements Sexton-Weingarten nesting (cheap
gauge force on fine steps, expensive fermion force on coarse). Lever = put the fermion force on a coarser
timescale so it is evaluated FEWER times per trajectory. This is a TUNING/enable task (pick the level
split + step counts, retune acceptance), independent of L1/L2. Biggest algorithmic potential. Files:
hmc_claude.cu (integrator choice + step tuning). Validate via acceptance + dH.

## Validation bar
L1/L2: force byte-identical-ish (same math, fewer recomputes / reorder) -> dH matches the original
trajectory to ~1e-8 (the force-ms benchmark standard); per-traj `# HMC` time the speedup. L3: acceptance
rate + dH distribution unchanged within statistics.

## Open questions
- `grad` is `const`; the L1 precompute writes `d_XZpre`/`d_XYpre` (logically-mutable scratch, like
  `d_Zs`) -- keep the same const posture (write through pre-allocated pointers).
- nstreams vs block: L2 replaces the 4-stream pole loop with one block kernel -- confirm no other caller
  relies on the streamed `grad`.
- Memory: +`2*(size-1)` N-vectors (`d_XZpre`/`d_XYpre`) ~ 2*10*49 KB ~ 1 MB -- negligible.
- Per-link COO build (`:545-547`) stays per-link (link-specific); cache only if it profiles significant
  after L1/L2 (secondary).
- Does `pf` loop (multiple pseudofermions) add a batching axis? precalc_grad/grad are per-pf; block over
  pf is a further (separate) lever -- defer.
