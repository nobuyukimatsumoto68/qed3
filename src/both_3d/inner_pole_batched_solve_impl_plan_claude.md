# Inner pole-solve batching -- impl plan (multi-shift / block CG over the Zolotarev poles)

**Reference (algorithm):** B. Jegerlehner, "Krylov space solvers for shifted linear systems",
hep-lat/9612014 -- the multi-shift CG used throughout this plan.

## BUGFIX -- freeze-converged (HMC NaN, 2026-06-05)

Running multi-shift inside HMC (massless `OverlapWMass`, `hmc_claude.cu`/`hmc_claude_debug.cu`,
free field) NaN'd at the first solve. **Root cause:** the `zeta` recurrence hits `0/0 = NaN` for
*converged* shifts -- the fast LARGEST-shift pole converges first, `zeta_m -> 0`, and
`zeta_new = (zeta*zeta_old*al_old)/denom` underflows to 0/0 while the slow seed keeps iterating; the
NaN `Z_m` poisons the reduction `sum A[m] Z_m` -> overlap output -> outer `dot2self` NaN. Intermittent;
the pole-loop is immune (no shared `zeta`). Not the pole count (n=11 and n=21 both), not the class
swap, not tolerance, not the lambda-escape print (benign, gives correct 1.154/10.82).

**Fix (in `matpoly_claude.h::solve_multishift`, SHARED header -> protects jj too):** FREEZE-CONVERGED.
`std::vector<char> frozen(npole)`; when `zeta[m]^2 * mu < mu_crit` the shift's residual is below tol ->
`frozen[m]=1`; frozen shifts skip the recurrence (`alm=0`, `zeta_new=0`) and the p-update (`betm=0`),
so `x_m` is kept at ~`TOL_INNER` and `zeta_m` never reaches the 0/0 underflow. The seed (`zeta=1`) is
never frozen early. Confirmed: no nan/assert/breakdown over 316 outer / 196 inner solves.

Diagnostic infra added (DEBUG): the `gam` breakdown guard now `assert(false)` (hard fail) + an
`IsVerbose2` per-shift true-residual/`||Z_m||` dump (at breakdown and end-of-loop). NOTE the
`assert(false)` is in the shared header (affects production) -- decide keep-loud-fail vs revert to
`break`. Side effect: fast poles now stop at ~1e-9 (not over-converged ~1e-11), so re-run
`test_overlap_multishift`/`test_mrhs` to confirm they still PASS (< 1e-6).

## Status

- M1 (device-scalar CG in `matpoly_claude.h`) was **reverted** (git restore; `gpu_header_claude.h`
  deleted) after benchmarking showed it cannot help the *outer* overlap solve (each outer iteration
  costs ~50 ms = ~20 inner pole solves, so removing 2 host syncs/outer-iter is negligible; the
  K-cadence overshoot made it ~24% slower). Conclusion: the real target is the **inner pole solves**.
- `test_mrhs_claude.cu` is kept as the BASELINE harness (reference `op_Dmsq` solve residual +
  end-to-end timing); the batched path will be benchmarked against it.

## The structure we are optimizing (key finding)

`OverlapWMass::mult_deviceAsyncLaunch` / `adj_deviceAsyncLaunch` (`overlap_wmass_claude.h:310,338`)
evaluate the Zolotarev sign function as a sum over `size-1` (~20) poles. Each pole `m` currently runs
its OWN CG via `MatPoly::solveAsync` (`overlap_wmass_claude.h:321,358`), parallelized over
`nstreams=4` OpenMP/CUDA streams, then combined with weights `A[m]` via `Taxpy_gen`.

Crucially, all poles solve the SAME matrix and SAME right-hand side, differing ONLY by a scalar shift:

$$
\left( A + \sigma_m \right) Z_m = \xi,
\qquad
A = \frac{1}{\lambda_\text{max}^2}\, D_W^\dagger D_W \ (\text{Hermitian PD}),
\qquad
\sigma_m = -\frac{k^2}{c'_m} > 0
$$

(`c'_m = cp[m] < 0`, so `\sigma_m > 0`; `m = 1 .. size-1`). The result is
$Z = \xi + \sum_m A[m]\, Z_m$, then $D_\text{ov} = \tfrac{1}{\lambda_\text{max}} D_W Z \cdot C + \dots$.

This is a **multi-shift** system (one matrix, one RHS, many scalar shifts) -- NOT a generic multi-RHS
problem. That distinction drives the algorithm choice below.

## Two algorithm options ("arrays of scalars indexed by thread id" fits both)

### Option A -- multi-shift CG (shared Krylov space; 1 matvec/iter for ALL poles)

Standard shifted/multi-mass CG (Jegerlehner, hep-lat/9612014). Run ONE unshifted-seed CG; the
expensive matvec $A p$ is done ONCE per iteration and shared by every pole. Per-pole solutions are
carried by cheap scalar recurrences. Pick the seed = smallest shift $\sigma_s$ (slowest system);
relative shifts $\hat\sigma_m = \sigma_m - \sigma_s \ge 0$.

Seed CG (on $A + \sigma_s$), scalars $\alpha^{(j)}, \beta^{(j)}$ as usual. For each shift $m$ carry
$\zeta_m^{(j)}, \zeta_m^{(j-1)}, p_m, x_m$:

$$
\zeta_m^{(j+1)} =
\frac{\zeta_m^{(j)}\,\zeta_m^{(j-1)}\,\alpha^{(j-1)}}
     {\alpha^{(j)}\beta^{(j-1)}\!\left(\zeta_m^{(j-1)}-\zeta_m^{(j)}\right)
      + \zeta_m^{(j-1)}\alpha^{(j-1)}\!\left(1-\hat\sigma_m\,\alpha^{(j)}\right)}
$$

$$
\alpha_m^{(j)} = \alpha^{(j)}\frac{\zeta_m^{(j+1)}}{\zeta_m^{(j)}},
\qquad
\beta_m^{(j)} = \beta^{(j)}\!\left(\frac{\zeta_m^{(j+1)}}{\zeta_m^{(j)}}\right)^{2}
$$

$$
x_m \mathrel{+}= \alpha_m^{(j)} p_m,
\qquad
p_m \leftarrow \zeta_m^{(j+1)} r + \beta_m^{(j)} p_m
$$

- **The "arrays of scalars indexed by thread id":** $\zeta_m, \zeta_m^\text{old}, \alpha_m, \beta_m,
  \hat\sigma_m$ are length-`(size-1)` device arrays. The $x_m$/$p_m$ updates run on a grid over
  $N \times (\text{size}-1)$; thread $(i,m)$ reads its shift's scalars by index $m$.
- **Cost:** ~1 matvec/iter total (vs ~20 today), + (size-1) cheap axpys. Matvec is the dominant cost
  => potentially the largest win. Seed (slowest pole) sets the iteration count; converged poles keep
  getting cheap scalar updates (no early-freeze needed in v1).
- **Risk:** the $\zeta$ recurrence and seed choice must be exactly right; all $\sigma_m>0$ keeps every
  system SPD (safe). Validate pole-by-pole vs the current independent solves.

### Option B -- batched block CG (independent per-pole CG, widened kernels; ~20 matvecs/iter)

Keep `size-1` independent CG solves but fuse them: store $x_m, p_m, r_m, q_m$ as $N\times(\text{size}-1)$
column-major blocks; the matvec/axpy/dot become block kernels on a grid over $N\times(\text{size}-1)$;
per-pole scalars $\mu_m,\gamma_m,\alpha_m,\beta_m$ are length-`(size-1)` device arrays indexed by the
thread's pole id. Raises occupancy (12 -> ~12*(size-1) blocks fills the SMs) and removes the
4-stream OpenMP orchestration, but does NOT cut the matvec count (each pole still has its own $p_m$).

- Simpler/more mechanical than A; independent convergence per pole.
- Win = occupancy + launch amortization only (the original perf diagnosis), not matvec reduction.

**Recommendation: Option A.** Same RHS + same matrix = textbook multi-shift; sharing the matvec is the
dominant saving and is exactly what this structure is for. B is the fallback if A's recurrence proves
fiddly or unstable, and B's block kernels are reusable by A's $x_m/p_m$ updates anyway.

## Files (both options)

- `includes/overlap_wmass_claude.h` -- replace the `for(m) solveAsync` loops in `mult_...` / `adj_...`
  with one batched/multi-shift solve; drop the per-pole OpenMP-over-streams orchestration.
- `includes/matpoly_claude.h` -- new `solve_multishift<N>` (A) or `solve_block<N>` (B); block/array
  scalar bookkeeping. Keep the existing single-RHS `solve`/`solveAsync` as reference.
- `includes/sparse_matrix.h` / `gpu_header.h` -- block matvec (option B) and the
  `x_m/p_m`-update kernels with per-pole scalar-array indexing (both options); `dev_div`-style tiny
  kernels for the seed scalars if kept on device.
- `test_mrhs_claude.cu` -- add the batched path next to the reference; diff `op_Dmsq` end-to-end.

## Chunks (proposed; gate = compile/run by user after each)

1. **C1 (math, no GPU):** implement the multi-shift scalar recurrence on the host (or a tiny test) and
   verify against the current independent pole solves on a small random SPD matrix. [Option A]
   **DONE / PASS** (`test_multishift_claude.cpp`): 20 shifts, n=64, all poles match independent CG to
   rel ~1e-10, residuals ~tol. Matvec count 374 (independent sum) -> 34 (multishift = slowest/seed
   pole) => ~11x matvec reduction. Seed = smallest shift confirmed. Sign of the zeta recurrence
   shift term is `(1 + \hat\sigma_m al_j)` (relative shifts >= 0); the seed needs no special case.
2. **C2:** `solve_multishift<N>` in `matpoly_claude.h` (device): seed CG + $N\times$npole block
   $x_m/p_m$ updates with the scalar arrays. Validate block $Z_m$ vs independent `solveAsync` per pole.
   **IMPLEMENTED (pending user compile/run).** Files: new `includes/multishift_kernels_claude.h`
   (`multishift_x_update`, `multishift_p_update` -- block kernels over $N\times$npole, column-major
   $[m N + i]$, per-pole double scalar arrays indexed by `m = gid/N`); `matpoly_claude.h` adds
   `#include` + `solve_multishift<N>(d_X, d_b, sigma_host, npole, tol, maxiter)` after the foreground
   `solve` (existing `solve`/`solveAsync` untouched). Seed scalars on HOST; the $\zeta_m$ recurrence
   runs on host each iter; coeff arrays copied H2D with **synchronous** `cudaMemcpy` (host buffers are
   overwritten every iter -- async would race; documented in an ASYNC NOTE). Fully synchronous on the
   default stream (foreground `on_gpu` matvec + blocking host-pointer dots). Convergence = seed
   (smallest shift) residual `< tol^2 ||b||^2`, matching the foreground `solve` convention (X_m updated
   before the break). Validation driver `test_multishift_gpu_claude.cu`: free-field setup, seed
   $A=(1/\lambda_\text{max}^2)D_W^\dagger D_W$, real poles $\sigma_m=-k^2/c'_m$ ($m=1..\text{size}-1$);
   diffs each block column vs a per-pole foreground `solve` + checks block residual
   $\Vert(A+\sigma_m)X_m-\xi\Vert/\Vert\xi\Vert$. NOT yet wired into overlap (that is C3).
   **DONE / PASS:** 10 poles (size=11), $\sigma\in[8.4\text{e-}6, 2.76]$, all rel `< 5e-11`,
   block resid `< 1e-9` (seed pole 0 largest at `9.3e-10`, as expected). Timing block one-pass
   `0.051s` vs serial 10-pole loop `0.366s` (~7.2x vs serial; production overlaps 4 streams).
3. **C3:** wire it into `mult_deviceAsyncLaunch` / `adj_deviceAsyncLaunch`; validate $D_\text{ov}$
   mult/adj vs current; re-run an existing overlap check (`check_conserved_current_claude.cu`).
   **IMPLEMENTED (pending user compile/run), SIDE-BY-SIDE (user-chosen).** `overlap_wmass_claude.h`:
   added members `d_Zblock`/`d_Yblock` (N*(size-1) blocks; alloc in ctor, free in dtor) and NEW
   methods `mult_deviceAsyncLaunch_ms` / `adj_deviceAsyncLaunch_ms` (single `solve_multishift` pass +
   block reduction, weight `A[m]` <-> column `m-1`, shift `-k^2/cp[m]`; RHS = `d_xi` for mult,
   `d_Ys[0]` for adj). Originals UNTOUCHED -> production still uses the pole loop until sign-off; then
   C3b = swap the `std::bind` call sites (and jj/HMC) to the `_ms` methods. Driver
   `test_overlap_multishift_claude.cu` diffs `mult` vs `mult_ms` and `adj` vs `adj_ms` (expect rel
   `< 1e-6`, both solve to TOL_INNER) + reports loop/ms speedup.
   **DONE / PASS:** D_ov rel `5.8e-12` (4.9x faster), D_ov^dag rel `1.8e-11` (4.4x faster) -- vs the
   REAL 4-stream pole loop (TITAN V), so the ~4.4-4.9x is the realistic per-apply win.
4. **C4:** benchmark `op_Dmsq` end-to-end in `test_mrhs_claude.cu` (expect matvec-count drop -> speedup).
   **IMPLEMENTED (pending user compile/run), side-by-side A/B.** `overlap_wmass_claude.h`: added
   `DHD_deviceAsyncLaunch_ms` / `DDH_deviceAsyncLaunch_ms` (route DHD's (D+m)/(D+m)^dag applies through
   `mult_ms`/`adj_ms`; originals untouched). `test_mrhs_claude.cu`: added `op_Dmsq_ms` (binds
   `DDH_deviceAsyncLaunch_ms`); solves the same system via `op_Dmsq` (ref) and `op_Dmsq_ms`, checks
   resid(ref)/resid(ms) ~tol and rel|x_ms-x_ref| < 1e-6, and benchmarks both (`--nrhs --reps`,
   default reps=10). PASS = resid(ms) < 10*tol and rel < 1e-6. Production call sites still on the
   originals (swap = C4b, after sign-off).
   **DONE / PASS:** resid(ms) `5.876e-6` == resid(ref) `5.876e-6`, rel(ms-ref) `1.6e-11`; end-to-end
   `op_Dmsq` solve `1256 ms` -> `288 ms` = **4.36x** (TITAN V, free field, nrhs=1, reps=10).
   **C4b DONE (Option 2, explicit call-site swap, pending user compile):** kept overlap originals +
   `_ms` intact (live A/B); switched PRODUCTION callers to `_ms`. Old lines left commented beside the
   new (per workflow). Edits: new `includes/sparse_matrix_claude.h` (LinOpDHDWrapper -> `DHD_..._ms`) +
   `includes/pseudofermion_claude.h` (heatbath -> `adj_..._ms`; action CG via the `_claude` wrapper);
   `jj_conn_tpproj_claude.cu` 6 binds -> `_ms`; `hmc_w_mass_claude.cu` + `hmc_w_mass_check_claude.cu`
   swap their `sparse_matrix.h`/`pseudofermion.h` includes -> the `_claude` copies. NOT swapped:
   `hmc_claude`/`hmc_debug`/`hmc_fermilab`/`hmc_fermilab_L2` use legacy `overlap.h` (no `_ms`).
   DEFERRED: `hmc_fermilab_wmass_claude.cu` (uses overlap_wmass but `*_fermilab*` excluded from make).
   HMC FORCE still pole-loop (`precalc_grad_deviceAsyncLaunch` `:468-477`, TWO same-RHS multi-shift
   solves RHS `d_eta`/`d_Ys[0]`) -> separate `solve_multishift` chunk (C5/force).

## Future work -- M2/C6: multi-RHS layered on the multi-shift seed (jj spatial-site loop)

Multi-shift (C1-C4) is ORTHOGONAL to the original multi-RHS (M2) idea: it cut the matvec COUNT
(~20 poles -> 1 shared Krylov), but the one matvec it keeps -- the seed $D_W^\dagger D_W$ apply on
`d_p0` in `solve_multishift` -- is still SINGLE-RHS (~$N/256$ = 12 blocks at $N{=}3072$, vs 80/108
SMs), so it is still occupancy-starved and is now the dominant cost. Multi-RHS composes
multiplicatively: batch several outer RHS so the seed matvec becomes an $N\times$nrhs SpMM that fills
the SMs.

Prime target = the JJ connected estimator's **spatial-site sink loop** (`jj_conn_tpproj_claude.cu`):
the forward leg $\phi'=D_m^{-1}\eta$ is a single shared solve, but the sink leg is
$\psi_n = D_m^{-\dagger} K^\dagger(n,0)\eta$ for $n=0..n_\text{sites}-1$ -- $n_\text{sites}$ solves, SAME
`op_Dmsq`, RHS per site. Batch them into ONE multi-shift-multi-RHS solve (nrhs = $n_\text{sites}$).
Scales the right way: $n_\text{sites}=10L^2+2$ = 12 at $L{=}1$, 42 at $L{=}2$ -- bigger batch exactly
when the lattice is more expensive.

Combined object: `nrhs` seeds $\times$ `npole` shifts; per-`(rhs,pole)` scalar arrays; seed matvec
over $N\times$nstack (column-major); per-seed convergence (iterate-to-slowest in v1, freeze-converged
deferred). `solve_multishift` threads an `nstack` arg; `multishift_x_update`/`_p_update` already block
over $N\times$npole -> extend to $N\times$nstack$\times$npole. DO AFTER C4b/C5; Nsight the seed matvec
first to size the headroom (gain stacks on the 4.36x but is bounded by that matvec's occupancy slack).
**FULL DESIGN (canonical): `mrhs_multishift_impl_plan_claude.md` -- STATUS: PARKED.** Decided not to
build: ~1.1-1.5x expected gain vs a large bug-prone `_block` refactor; multi-shift's 4.36x already
stands. Revisit only if L=2/A100 widens the occupancy gap or a profile shows big headroom. Locked
decisions retained (specialized two-CSR seed, shared worst-column stop, `nstack` runtime, adaptive grid).

## Future work -- conserved-current kernel (GATED on D_ov verification)

The same Zolotarev inner-pole structure (matrix $A=\frac{1}{\lambda_\text{max}^2}D_W^\dagger D_W$,
shift $\sigma_m=-k^2/c'_m$) recurs inside `ConservedCurrent::apply_k` / `apply_k_dag`
(`conserved_current_claude.h`). FOUR inner-solve loops; THREE are clean same-RHS multi-shift,
directly replaceable by `solve_multishift` exactly as in the overlap mult/adj:

- `apply_k` Step 1 (`:128-135`) -- RHS `d_xi`. **same-RHS multi-shift** => `solve_multishift`.
- `apply_k_dag` Step 1 (`:217-223`) -- RHS `d_tmp1` (m-independent). **multi-shift**.
- `apply_k_dag` Step 3 (`:245-249`) -- RHS `d_Zs[0]` (m-independent). **multi-shift**.
- `apply_k` Term B (`:167-189`) -- RHS $u_m = X^\dagger W Z_m - W^\dagger X Z_m$ is **m-DEPENDENT**
  (depends on $Z_m=R_m d_\xi$), shift also varies => genuinely independent (shift,RHS) pairs, NOT
  single-Krylov multi-shift. Candidate for batched-independent block CG (Option B) for occupancy
  only; leave as-is for now.

**C5 IMPLEMENTED (side-by-side, pending user compile).** `conserved_current_claude.h`: added member
`d_Zblock` (N*(size-1); alloc ctor / free dtor) and side-by-side methods `apply_k_ms` /
`apply_k_dag_ms` (originals UNTOUCHED). Each same-RHS loop -> ONE `solve_multishift` over the seed
`{&D.M_DW,&D.M_DWH}`/lambda_max^2 with shifts `-k^2/cp[m]`:
- `apply_k_ms` Step 1 (RHS `d_xi`) -> multishift into `d_Zblock`, COPY OUT to `d_Zs[1..]` so the
  reduction + Term B are byte-identical to `apply_k`. Term B (m-dependent RHS) unchanged.
- `apply_k_dag_ms` Step 1 (RHS `X^dag xi` in `d_tmp1`) -> multishift, copy out to `d_Zs[1..]` (=Y_m).
  Term A^dag (RHS `w0=d_Zs[0]`) -> multishift, reduce `d_result += C*A[m]*d_Zblock[m-1]` directly
  (no copy-out). Term B^dag (m-dependent RHS r_m) unchanged. `d_Zblock` reused across the two
  multishift passes (Y_m already copied to `d_Zs[]` before the 2nd pass).
Driver `test_conserved_multishift_claude.cu`: temporal link el={0,0}, free field, diffs `apply_k` vs
`apply_k_ms` and `apply_k_dag` vs `apply_k_dag_ms` (expect rel < 1e-6) + loop/ms timing.
**C5 VALIDATED (PASS):** test_conserved_multishift_claude.cu -- K rel `1.3e-10` (1.34x), K^dag rel
`5.2e-10` (2.19x) vs the pole-loop originals (partial speedup: Term B/B^dag + COO matvecs aren't
multi-shift). **C5b DONE:** `operator()` dispatch (`:83`) switched to `apply_k*_ms` (old lines
commented) -> jj/meson/disc/sp/axial `op_K` now multi-shift; `check_conserved_current*_claude.cu`
call `apply_k`/`apply_k_dag` DIRECTLY (not via `operator()`) so they still validate the ORIGINAL
reference -- unaffected. Originals kept (not retired).
NOTE (line refs above are pre-C5-edit; apply_k_ms@:297, apply_k_dag_ms@:358 after the edit).

## Decisions (locked)

1. **Algorithm: Option A -- multi-shift CG.** [CONFIRMED by user]
2. **Seed = smallest-shift pole** $\sigma_s=\min_m\sigma_m$; relative shifts $\hat\sigma_m=\sigma_m-\sigma_s\ge0$.
   Guarantees all other poles have converged once the seed meets tolerance. [default]
3. **Layout:** column-major block $x_m[m\cdot N + i]$, $m=0..\text{npole}-1$ (npole = size-1). [default]
4. **Seed scalars on host.** One seed CG; its 2 syncs/iter are negligible against the now-shared
   matvec -- do NOT reintroduce M1 device-scalar machinery for the seed. The per-pole
   $\zeta_m/\alpha_m/\beta_m$ live in device arrays (updated by a tiny npole-thread kernel). [default]
5. **Drop the 4-stream OpenMP pole parallelism** -- one wide solve replaces it. [default]
6. **Convergence:** iterate to the seed (slowest pole) tolerance `TOL_INNER`; no per-pole early-freeze
   in v1. [default]

All defaults above stand unless the user redirects.
