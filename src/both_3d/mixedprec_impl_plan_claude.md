# Mixed-precision inner solver -- implementation plan

**Goal.** Cut the inner overlap CG cost (the $D_W^\dagger D_W$ multishift that dominates the L4
force) by running the Krylov iteration in **fp32** while preserving the **fp64** final tolerance via
periodic true-residual **reliable updates**. This is a *reversibility-clean* speed lever: the
converged solution is fp64-accurate to the requested tol (only the search direction is low precision),
so it does NOT loosen the force tolerance -- unlike simply raising `TOL_INNER`, which would corrupt
MD reversibility. (Deflation of the near-zero Wilson modes is handled in a separate thread, on the
disconnected-diagram valence solves where it pays more.)

## CONCLUSION (2026-07-15) -- NET WIN AT L4 (~1.35x); net loss at L1/L2. Corrected on REAL configs.

Mixed-precision inner multishift solver investigated end-to-end (Chunks 0-2 + per-shift cleanup) and
VALIDATED to full 1e-9 on REAL thermalized configs. **Nets ~1.35x at L4 (the expensive production case);
net loss at L1/L2.** Literature-backed (two QUDA papers below). NOTE: an EARLIER "shelve everything"
verdict was WRONG -- an artifact of testing on a tiny well-conditioned gaussian config (nt16); real
configs (nt128 + genuine near-zero Wilson mode) flip L4.

- **Chunk 0 GATE**: isolated fp32 vs fp64 $D_W^\dagger D_W$ matvec = 1.4x/1.9x/2.2x (L1/L2/L4).
- **Chunk 1 (single-shift) PASSED**: reliable-update-to-fp64 works; best `tol_f=1e-3`. It's the INFERIOR
  defect-correction (Alg 1 of 0911.3191, discards Krylov each restart), not the reliable-update Alg 2;
  and it does NOT re-project the search direction `p` at reliable updates (0911.3191 fn 9 / PoS2022 3.1).
- **Chunk 2 (multishift)**: fp32 matvec + periodic fp64 reliable update of the SEED + fp64 crossover tail
  at `tol_f`, then PER-SHIFT fp64 CLEANUP to 1e-9. ROOT ISSUE (confirmed verbatim, PoS LATTICE2022 338
  Sec 3.2): the fp32 matvec injects a NON-PROPORTIONAL per-pole error $\epsilon_m\sim10^{-7}$ that breaks
  the multishift co-linearity $r_m=\zeta_m r_\text{seed}$, so the poles stall at the fp32 floor even while
  the seed converges. Removing $\epsilon_m$ needs a per-shift fp64 cleanup (different RHS per pole -> NO
  shared Krylov; == warm-started per-pole fp64 CG, the two are identical). `MixedSolver::solve_multishift`
  (`cleanup=true`, `tol_f`, `freq`).
- **REAL-CONFIG net speedup to 1e-9** (Nf6 gsq8 massless, nt128, `test_multishift_mixed_claude.cu` +
  `tmp_multishift_mixed_realcfg_claude.sh`, per-shift cleanup ON, best `tol_f=1e-6`):
  | L | N | ref iters | fp32/fp64/cleanup | net speedup |
  |---|---|---|---|---|
  | 1 | 3072 | ~97 | 67/31/72 | **0.59x (loss)** |
  | 2 | 10752 | ~583 | ~400/~180/~195 | **0.81x (loss)** |
  | 4 | 41472 | ~874-1031 | ~700/~260/~160 | **1.33-1.35x (WIN)** |
  Floor-only (poles at ~1e-7, NO cleanup) was L4 ~1.5x; the per-shift cleanup (~160 fp64 iters = npole x
  ~16, warm from ~1e-7) trims it to ~1.35x net. Cleanup is a SMALL fraction of the big stiff L4 solve
  (~1000 iters) but a LARGE fraction of small L1/L2 -> only L4 nets positive.
- **Why L4 wins / L1-L2 lose**: L4 real configs are big (N=41472) AND stiff (near-zero Wilson mode ->
  ~1000 iters), so the solve is matvec-bound and the fp32 matvec surfaces (QUDA's 1.5-2x regime). L1/L2
  are small/latency-bound; the un-shareable cleanup dominates.

## L4 tuning DONE (2026-07-15): nested-fp32 cleanup + tol_f. TUNED = tol_f=1e-6, n_nested=1.
Made each per-pole cleanup MIXED: `MixedSolver::cleanup_pole_mixed` = n_nested fp32 stages (each targeting
the ABSOLUTE crit_d, rel=sqrt(crit_d/||r_m||^2)~1e-2) + fp64 finish. 2D scan `test_multishift_mixed_claude.cu`
(tol_f{1e-4,1e-5,1e-6} x n_nested{0,1,2}) on 3 real L4 configs (`tmp_multishift_tune_l4_claude.sh` ->
`tmp_multishift_tune_l4_v2_claude.log`):
- **BUG FOUND+FIXED**: a FIXED relative tol_f_clean=1e-4 over-solved the cleanup (needs ~2 orders, 1e-4 forces
  ~4) -> clean_lo EXPLODED to 1417 (vs fp64's 181) -> net 0.73x. Fix = fp32 stage targets ABSOLUTE crit_d
  (per-pole rel = sqrt(crit_d/||r_m||^2)); tol_f_clean=1e-4 is now just a floor. After fix clean_lo=181 ==
  the fp64 cleanup count exactly (no over-solve).
- **RESULT (within-run, noise-free nn=0 vs nn=1)**: n_nested=1 beats n_nested=0 (pure fp64 cleanup) by
  **~+7%** at every config/tol_f (e.g. tol_f=1e-6: 1.10x -> 1.17-1.18x). **n_nested=2 == n_nested=1** (one
  fp32 stage already reaches crit_d). **tol_f=1e-6 > 1e-5 > 1e-4** (tighter crossover = more fp32 pass).
- **TUNED L4 config = tol_f=1e-6, n_nested=1, tol_f_clean=1e-4 floor, freq=10** (now the `solve_multishift`
  defaults). Net ~1.4-1.45x at L4 (pure-fp64 ~1.35x + the ~7% fp32-cleanup gain; absolute ms are GPU-clock
  noisy run-to-run, so quote the within-run deltas). All PASS ~9.9e-10.

## NEXT: Chunk 3 -- wire the tuned mixed multishift into the L4 FORCE op Df
L4-ONLY (L1/L2 stay pure fp64 -- net loss there). Route Df's inner multishift (n_f=11 poles, window
[2lmin,lmax]) through `MixedSolver::solve_multishift` (tuned defaults) + `refresh()` per MD step (recast the
fp32 CSR from Df's fp64 D_W, mirroring Grid `precisionChange(Umu_f,Umu)`); guard behind a flag; verify HMC
dH/reversibility UNCHANGED + measure force-solve speedup. Optional later: reliable-update Alg 2 + p-reprojection.

**References** (both from NM): M.A. Clark et al., *"Solving Lattice QCD systems ... mixed precision
solvers on GPUs"*, arXiv:0911.3191 (single-shift reliable-update; Alg 1 defect-correction vs Alg 2
reliable-update; footnote 9 gradient reprojection). M.A. Clark et al., *"Maximizing the Bang Per Bit"*,
PoS LATTICE2022 (2023) 338 (multishift co-linearity breakdown Sec 3.2; per-shift refinement 1.5-2x;
int20/int30 sloppy formats).

Code artifacts kept (do NOT delete; reusable if volumes grow): `includes/sparse_matrix_mixed_claude.h`
(fp32 CSR matvec + casts), `includes/matpoly_mixed_claude.h` (`MixedSolver`: single-shift + multishift
mixed CG), tests `test_{matvec,solve,multishift}_mixed_claude.cu`, handoffs `tmp_*_mixed_*_claude.sh`.

---

**Algorithm source (mandatory citation).** Mixed-precision reliable-update CG:
M.A. Clark, R. Babich, K. Barros, R.C. Brower, C. Rebbi, *"Solving Lattice QCD systems of equations
using mixed precision solvers on GPUs"*, Comput. Phys. Commun. 181 (2010) 1517, **arXiv:0911.3191**.
Reliable-update idea = G.L.G. Sleijpen, H.A. van der Vorst, *"Reliable updated residuals in hybrid
Bi-CG methods"*, Computing 56 (1996) 141.

**Grid reference implementations (pinned; NM pointer to `/mnt/baracuda_14/grid_claude/Grid`).**
Two DISTINCT structures -- mirror the matching one per chunk:
- **Single-shift (Chunk 1)** = `Grid/algorithms/iterative/ConjugateGradientMixedPrec.h`
  (`MixedPrecisionConjugateGradient<FieldD,FieldF>`, C. Kelly). *Restarted defect correction*: OUTER
  loop recomputes the fp64 residual $r=b-A\,x$; if not converged, cast $r\to$fp32, run a WHOLE fp32 CG
  $A\,\delta=r$ to an (adaptive) inner tol, cast $\delta\to$fp64 and accumulate $x{+}{=}\delta$;
  finish with one fp64 "patch-up" CG from $x$ to hit the exact tol. Inner CG lives entirely in fp32
  (minimal precision changes -- matches NM's "lower-prec vector objects" note); only the residual +
  accumulate are fp64.
- **Multishift (Chunk 2)** = `Grid/algorithms/iterative/ConjugateGradientMultiShiftMixedPrec.h`
  (`ConjugateGradientMultiShiftMixedPrec`, CK 2020). Different: ONE multishift pass, matvec + residual
  in fp32, **search directions AND shifted solutions in fp64**; every `ReliableUpdateFreq` iterations
  the residual is recomputed in fp64 (reliable update); final regular fp64 multishift cleanup if
  needed. (Multishift can't cleanly restart -- shifted Krylov spaces are shared -- so it's an in-pass
  reliable update, not a per-restart defect correction.) `ShiftedLinop` adds the shift to the base
  linop for the fp64 cleanup.
The disc thread's `disc_multipleGamma_binary_mixedprec_claude.cc` +
`disc_mixedprec_impl_plan_claude.md` (same Grid dir) already use the single-shift class for the
Mobius valence solves -- prior art for the idiom.

**HMC wiring reference (Chunk 3; NM pointer):**
`/mnt/baracuda_14/grid_claude/Grid_sdm_build/src/gauge_gen_Nc4/dweofa_mobius_HSDM_hasenbusch_claude.cc`
(Hasenbusch-preconditioned EOFA Mobius HMC). Pattern to mirror:
- The mixed-prec CG is wrapped as an `OperatorFunction` (`MixedPrecisionConjugateGradientOperatorFunction`,
  lines 109-226) whose `operator()` **refreshes the fp32 operator from the current fp64 gauge every
  solve** (`precisionChange(Umu_f, Umu_d)` per MD force eval, lines 194-198) then runs the mixed CG.
  -> our analog: call `CSRf::cast_from(Df.M_DW.val)` / `cast_from(Df.M_DWH.val)` right after
  `Df.update(U)` each MD step, before the mixed solve.
- **Separate `ActionStoppingCondition` (tight, heatbath/accept-reject) vs `DerivativeStoppingCondition`
  (force) tolerances**, same mixed CG class for both -- matches our two-operator split (action op D vs
  force op Df) and the reversibility-bounded force tol.
- NOTE: that reference is EOFA (single-shift, no rational) so it uses the single-shift mixed CG
  directly. Our overlap force is the Zolotarev sign function = genuine multishift, so Chunk 3 wires the
  **multishift** mixed CG (Chunk 2) with this per-solve fp32-refresh + separate-tol wrapper.

## Current code (what mixed precision has to replace)

All in NESTED `/mnt/barracuda22/qed3/qed3/src/both_3d/`. Everything is `CuC = cuDoubleComplex` (fp64).

- Inner solve = `MatPoly::solve_multishift<N>` (`includes/matpoly_claude.h:389`). Seed matrix
  $A=(1/\lambda_\text{max}^2)\,D_W^\dagger D_W$ (`MatPoly Aseed` built in
  `overlap_wmass_claude.h:542`), shifts $\sigma_m=-k^2/c_m$ ( `sigma[m-1]=-k*k/cp[m]` ), single RHS,
  all $n_\text{pole}$ poles in one Krylov pass. Seed = smallest shift drives convergence.
- The expensive per-iteration op is the seed matvec `on_gpu<N>(d_q,d_p0)`
  (`matpoly_claude.h:449`) = **two sparse CSR applies** `mult<CuC,N>` of `M_DW` then `M_DWH`
  (`sparse_matrix_claude.h:25`). This is memory-bandwidth-bound in the fp64 `val` array -> fp32
  storage/compute is the ~2x lever.
- Single-shift CG `MatPoly::solve<N>` (`matpoly_claude.h:298`) -- the outer/valence solves.
- Cost accounting: `g_matpoly_cg_iters` (`matpoly_claude.h:10`); Osborn cost
  $C_S/(\tau^2\langle P_a\rangle)$.
- `LinOp` interface (`sparse_matrix_claude.h:9`): `operator()(CuC* res, const CuC* v)` +
  `Async(...)`. CSR val/cols/rows all fp64.

## Mixed-precision algorithm (equations)

Solve the SPD system $A\,x=b$ with $A=D_W^\dagger D_W/\lambda_\text{max}^2+\sigma$ (single shift shown;
multishift below). Keep the accumulated solution $x$ and reference residual $r$ in fp64; run the CG
correction in fp32.

Outer defect-correction loop (target relative tol $\tau$, i.e. stop when
$\Vert r\Vert \le \tau\Vert b\Vert$):

$$
r \;=\; b - A\,x \quad(\text{fp64 matvec, once per restart})
$$

Inner fp32 CG solves the correction system $A\,\delta = r$ starting from $\delta=0$, in single
precision, until the fp32 residual has dropped by a **reliable-update factor** $\rho$ (e.g.
$\rho=10^{-1}\!\dots\!10^{-2}$, or the fp32 floor $\sim\sqrt{\varepsilon_{32}}\approx3\times10^{-4}$),
then:

$$
x \;\leftarrow\; x + \delta_{\text{(fp32)}}\qquad
r \;\leftarrow\; b - A\,x \;\;(\text{fp64})
$$

Repeat. Because $r$ is always the *true* fp64 residual, the loop converges to full fp64 $\tau$; the
bulk matvecs (the many fp32 CG iterations) run at ~2x bandwidth, with only $O(\log_{1/\rho}(1/\tau))$
extra fp64 matvecs for the restarts.

**Multishift version.** Seed $=$ smallest shift. Run one fp32 multishift pass on the current fp64 seed
residual $r$ (shift set $\{\hat\sigma_m=\sigma_m-\sigma_\text{seed}\}$ unchanged), producing fp32
corrections $\delta_m$ for every pole. Reliable update: $x_m\leftarrow x_m+\delta_m$ for all $m$;
recompute the fp64 seed residual $r=b-(A+\sigma_\text{seed})x_\text{seed}$; restart. One fp32
multishift pass per restart corrects all poles simultaneously (same shifts) -- so the multishift
speedup is preserved. Stop when the seed's fp64 residual meets $\tau$ (seed slowest => all poles
converged, same criterion as the current solver).

## fp32 infrastructure needed (there is none today)

1. **fp32 CSR operator** `CSRf` mirroring `CSR<N>`: `val` in `cuFloatComplex`, same `cols/rows`.
   Produced by a device cast kernel from the existing fp64 `d_val` whenever $D_W$ is rebuilt (cheap,
   one pass over the $\sim n_\text{nnz}$ entries). Needed for BOTH `M_DW` and `M_DWH`.
   -- In the FORCE, $D_W$ changes every MD step, so the fp32 mirror is re-cast every step (negligible
   vs the solve).
2. **fp32 vector kernels**: `Taxpy` / dot in `cuFloatComplex` (the fp32 CG inner recurrences). The
   fp64 axpy/dot (`Taxpy<CuC>`, `cublasZdotc`) stay for the fp64 reliable-update residual + accumulate.
3. **fp32 multishift kernels**: fp32 versions of `multishift_x_update` / `multishift_p_update`
   (`multishift_kernels_claude.h`) for the shifted recurrences during the fp32 pass. (The zeta/alpha
   scalar recurrences can stay in fp64 on the host -- they are $O(n_\text{pole})$, not the cost.)

Memory: one extra fp32 copy of $D_W$'s `val` ($\approx$ half the fp64 `val` bytes) per operator.
Non-issue.

## Decisions (NM, 2026-07-15)

1. **Ordering = single-shift beachhead FIRST** (prove fp32 infra + reliable-update logic + measure
   raw matvec win), THEN multishift, THEN wire the force. **Deflation is a separate thread** (valence /
   disconnected-diagram solves, where it pays more).
2. **Micro-benchmark first** (Chunk 0 gate): confirm the raw fp32-vs-fp64 $D_W^\dagger D_W$ matvec win
   per L before building the reliable-update machinery.
3. **Header style (NM):** copy the relevant headers with a `_mixed` suffix and templatize on the
   INNER precision, rather than editing the fp64 headers in place (preserve-original convention).
4. Reliable-update factor $\rho$: default $10^{-2}$ (or fp32 floor $\sqrt{\varepsilon_{32}}\approx
   3\times10^{-4}$); tune after the micro-benchmark.

## Ordered chunks

### GATE RESULT (2026-07-15) -- PASSED, proceed
Raw fp32-vs-fp64 $D_W^\dagger D_W$ matvec (`test_matvec_mixed_claude.cu`, Nt=16, gaussian U w=0.3),
fp32 rel err $\sim7\times10^{-8}$ (= fp32 epsilon, matvec accurate) at every L:
| L | N | nnz | fp64 us | fp32 us | speedup |
|---|---|---|---|---|---|
| 1 | 384 | 6144 | 22.95 | 16.31 | **1.41x** |
| 2 | 1344 | 23424 | 25.56 | 13.76 | **1.86x** |
| 4 | 5184 | 92544 | 27.69 | 12.44 | **2.22x** |
L1 latency-bound (fp64 matvec ~flat in N = launch-dominated at tiny N) -> only 1.41x; the win grows to
**2.22x at L4 (the target)**. Worthwhile. Chunk 1 next.

- **Chunk 0 (infra + GATE):** fp32 CSR matvec + fp64<->fp32 cast kernels in
  `includes/sparse_matrix_mixed_claude.h` (new). Micro-bench `test_matvec_mixed_claude.cu` (new):
  build the real $D_W$ at $L$, cast to fp32, time fp32 vs fp64 $D_W^\dagger D_W$ matvec + report
  rel err ($\sim10^{-6..7}$) and speedup per L. Handoff `tmp_matvec_mixed_build_claude.sh` (user runs
  L1/L2/L4 on a GPU). **GATE: proceed to Chunk 1 only if the win is worthwhile.**
  Access: `Dnew.d_DW.len` = nnz; `Dnew.M_DW.{val,cols,rows}` (fp64 D_W CSR),
  `Dnew.M_DWH.{val,cols,rows}` (D_W^dag = `d_valH`/`d_colsT`/`d_rowsT`).
- **Chunk 1 (single-shift mixed CG) -- DONE + VALIDATED (2026-07-15).** `MixedSolver<N>` in
  `includes/matpoly_mixed_claude.h` (restarted defect correction; fp32 `cg_fp32` + fp64 patch `cg_fp64`
  + `applyA_fp32/_fp64` + `refresh()`); fp32 kernels `Taxpy_f`/`axpby_shift_f/_d` in
  `sparse_matrix_mixed_claude.h`. Validate `test_solve_mixed_claude.cu` (seed/mid/large Zolotarev pole
  vs `MatPoly::solve`). ALL PASS L1/L2/L4: mixed TRUE fp64 residual ~1e-9 (= tol; reliable update
  reaches full precision), sol agrees with fp64 ref ~1e-9.
  - KEY FIX: the fp64 patch-up must target the ABSOLUTE residual `stop = tol^2||b||^2` (reduce r by ~1
    order), NOT relative tol -- a relative tol re-solved the whole system (patch_d 63 -> 1, over-solve
    true-resid 1e-16 -> 1e-9). `cg_fp64(..., abs_crit)`.
  - Speedup (ref/mix, shared GPU1 -> absolute ms noisy across L): L1 ~1.2x, L2 ~2.1x, L4 ~1.1x seed.
    Single-shift RESTART inflates iters (`inner_f ~ 1.24x` the fp64 count -> Krylov superlinearity lost
    per restart), so the fp32 matvec 2x nets only ~1.1-2x. Large-shift poles <1x (finish in ~8 iters,
    fixed fp32 overhead dominates) -- irrelevant: in multishift the SEED drives the count. Expected;
    the real win is Chunk 2 (shared matvec across all poles, in-pass reliable update = no restart loss).
- **Chunk 2 (multishift mixed):** `solve_multishift_mixed<N>` (restarted, fp64 defect correction) in
  `matpoly_mixed_claude.h` + fp32 multishift kernels. Validate vs `solve_multishift<N>`.
  `test_multishift_mixed_claude.cu` (new).
- **Chunk 3 (wire force):** route the FORCE operator Df's inner solve through the mixed multishift
  (guarded; fp64 path untouched). Verify HMC dH / reversibility unchanged; measure force-solve
  speedup at L4. Files: Df construction in `hmc_hasenbusch_claude.cu` + the force headers.

## Chunk 1b -- tol_f tuning (single-shift) -- DONE 2026-07-15
`MixedSolver::solve(d_x,d_b,sigma, tol_f, tol_d, ...)` = explicit TWO tolerances (dropped the adaptive
inner_tol): fp32 stages each solve the fp64-residual system to RELATIVE `tol_f`; fp64 finish to ABSOLUTE
`stop_d = tol_d^2||b||^2`. Scan `test_solve_mixed_claude.cu` (seed pole, tol_f {1e-2..1e-7}, tol_d=1e-9),
ALL PASS L1/L2/L4:
- **tol_f = 1e-3 is the OPTIMUM at every L** (min total fp32 iters, best speedup): 3 stages, fp32it
  131/96/64 vs fp64-ref 121/90/61, speedup ~1.4x/2.4x/1.3x (noisy -- GPU1 contention).
- **fp64it = 0 for ALL tol_f**: the fp32 stages, recomputing the true fp64 residual between them, reach
  tol_d ENTIRELY in fp32 -> the fp64 finish is pure safety (never fires). Reliable-update-between-stages
  IS the full-precision mechanism.
- tol_f optimum = stages-vs-depth balance: 1e-2 -> 5 stages (restart penalty); 1e-4..1e-7 -> tighter
  per-stage overshoots. 1e-3 = 3 stages x ~3 orders = the 9 orders 1->1e-9.
- Even optimal, fp32it slightly EXCEEDS fp64-ref (~1.1x restart penalty -> Krylov superlinearity lost
  per restart). Single-shift ceiling ~1.3-2.4x. Chunk 2 (continuous multishift pass, no restart) removes it.
- LESSON for Chunk 2: fp32 single-pass floor ~1e-6; the fp32->fp64 crossover tol_f should sit ABOVE it
  (~1e-4..1e-5) since the continuous pass has no restarts to punch below the floor -- the fp64 TAIL
  (shared-Krylov continuation) does the sub-floor orders.

## Related
[[project_two_operator_splitpole]], [[project_frozen_window]]. Prior cost findings:
`hasenbusch_tuning_results_claude.md` (pole count weak; window bites only at L4; heavy frame ~18x).
