# Multi-RHS on the multi-shift seed (M2/C6) -- impl plan

**Reference (algorithm):** B. Jegerlehner, "Krylov space solvers for shifted linear systems",
hep-lat/9612014 -- the underlying multi-shift CG this layers `nstack` RHS onto.

> **STATUS: PARKED (2026-06-05).** Decided NOT to build now: the expected gain is only ~1.1-1.5x
> (occupancy lever on one matvec), and the full-stack `_block` refactor is a large, bug-prone surface
> -- debugging cost outweighs the return. The multi-shift win (4.36x, C1-C5) already stands. Design
> kept on file. **Revisit only if** (a) a bigger lattice ($L{=}2$, n_sites=42) or A100 widens the
> seed-matvec occupancy gap, or (b) a profile (the Nsight command below) shows large unused headroom.
> If revisited, the locked decisions below (specialized two-CSR seed, shared worst-column stop,
> `nstack` runtime, adaptive grid) still hold.

## Goal

Raise GPU occupancy by batching `nstack` outer right-hand sides ("stack") so the inner multi-shift
**seed matvec** ($A p_0$, $A=\frac{1}{\lambda_\text{max}^2}D_W^\dagger D_W$) becomes an
$N\times$nstack SpMM that fills the SMs. Multi-shift (C1-C5) cut the matvec COUNT (~20 poles -> 1);
this cuts the matvec UNDER-UTILIZATION (~12 blocks of 80/108 SMs -> ~12*nstack blocks). The two
compose. **Realistic target: ~1.5x** (occupancy lever only; the seed matvec's slack bounds it).

Prime consumer NOW: the JJ connected estimator's spatial-site **sink loop** -- `nstack = n_sites`
(= 10L^2+2 = 12 at $L{=}1$, 42 at $L{=}2$). The forward leg $\phi'$ is a single shared solve.

## Decisions (from discussion)

- **`nstack` is a RUNTIME variable, never hardcoded.** It must vary with $L$ (n_sites grows) and
  serve other stacks later. Thread it through as a function arg; recompute every grid from it.
- **Scope = jj for now**, but keep `nstack` fully general (do NOT bake in 12). The `_block` operators
  are written once and are stack-agnostic.
- **Reusable for HMC later (out of scope now):** the same `_block` solve can stack several
  pseudofermions (each flavor-pair's solve is an independent RHS) -> HMC could call the block path
  with `nstack = #pseudofermions`. We do NOT touch HMC now, but nothing in the design precludes it.
- **Adaptive launch sizing (replaces the hardcoded constants):** the global `NThreadsPerBlock`/
  `NBlocks` are hand-tuned for single-RHS `N` and are NOT optimal here. For every `_block` kernel,
  recompute the grid from the actual work size (`N*nstack` or `N*nstack*npole`) and make
  `NThreadsPerBlock` a tunable (sweep 128/256/512/1024; retune for A100). Do NOT reuse the global
  `NBlocks` macro for block launches.
- It is **mostly copy&paste**: each `_block` op mirrors its `_ms` sibling with `nstack` threaded in
  + the grid recomputed; keep the single-RHS `_ms`/originals intact as the reference.

## Key structural fact (sizes the whole effort)

The seed matvec has only ONE inner RHS per overlap apply -- the outer CG's current search
direction. So `nstack` must flow from the OUTER CG down through every layer:

$$
\text{outer block CG (nstack)} \to \text{block } DDH_{ms} \to \text{block } mult_{ms}/adj_{ms}
\to \text{block solve\_multishift (nstack seeds} \times \text{npole shifts)} \to \text{block seed matvec } (N\times\text{nstack}).
$$

This is a FULL-STACK change. It is NOT achievable by only touching `solve_multishift` (there is no
`nstack` at the inner level unless the outer CG carries it).

## Layout (decision: column-major)

- Seed-level blocks (`p0`, `r`, `q`, the nstack outer/seed vectors): `[c*N + i]`, c in [0,nstack).
- Per-(rhs,pole) blocks (multi-shift `X_{c,m}`, `p_{c,m}`): `[(c*npole + m)*N + i]` -- column (c,m)
  contiguous, nstack*npole columns.
- Scalars: per-rhs `mu[c]/gam[c]/al[c]/bet[c]` (length nstack); per-(rhs,pole)
  `zeta[c*npole+m]/zeta_old/alm/betm` (length nstack*npole). Indexed by the thread's (c,m).

## Block kernels (the occupancy lever -- grow the grid, don't loop in a fixed grid)

- **block CSR matvec** `mult_block<T,N>(res, v, csr, cols, rows, nstack)`: thread (i,c) ->
  `res[c*N+i] = sum_k csr[k]*v[c*N+cols[k]]`; CSR structure read once, reused across columns.
  Grid sized `N*nstack`. (Seed $A$ = two CSR applies M_DW then M_DWH -> two block matvecs.)
- **block Taxpy / Taxpy_gen**: elementwise over the relevant block length (N*nstack for seed,
  N*nstack*npole for the per-pole updates); per-column coeff array OR single broadcast coeff.
- **block dot**: column-wise `X^H Y -> nstack` (or nstack*npole) values. Start with strided/offset
  `cublasZdotc` per column (columns contiguous => cheap); custom reduction later if profiled.
- **multishift_x_update/_p_update**: extend the existing N*npole kernels to N*nstack*npole; the
  per-pole scalar arrays become per-(rhs,pole), indexed by `gid -> (c,m) = ((gid/N)/npole? )`.
  (Define index math from the chosen column order.)

## Chunks (each = compile/run gate; bottom-up, validate vs nstack separate single-RHS solves)

- **C6a** block leaf kernels: `mult_block` (CSR), block `Taxpy`/`Taxpy_gen`, block dot. New header
  (do NOT mutate shared `sparse_matrix.h`/`gpu_header.h` -- copy or add `_claude` kernels). Unit-test
  each vs the single-RHS kernel, column by column.
- **C6b** `MatPoly::solve_multishift_block<N>(d_X, d_B, sigma, npole, nstack, tol)` + block `on_gpu`
  for the seed matvec. Validate the block `Z_{c,m}` vs `nstack` separate `solve_multishift` calls
  (identical to tol). STANDALONE-VALIDATABLE -- this is the foundational, self-contained piece.
- **C6c** block `mult_ms`/`adj_ms`/`DDH_ms` in overlap (inherit block-ness from `solve_multishift_block`
  + block `A[m]` combine). Validate `D_m^{-1}`/`D_m^{-dag}` block vs single.
- **C6d** block outer CG: `MatPoly::solve` over `nstack` RHS (per-column scalars, iterate-to-slowest).
  Validate `op_Dmsq` block solve vs nstack single solves.
- **C6e** batch the jj sink loop (`nstack = n_sites`); free-field re-run -> Vpp/Vmm unchanged + faster.

## Resolved (from discussion)

- **Scope = jj sink loop, build bottom-up** (C6a..C6e), each chunk a compile gate. Worth it;
  largely copy&paste off the `_ms` siblings. HMC NOT touched (future stack-of-pseudofermions reuse).
- **`nstack` = runtime arg** (not template); grids recomputed from `N*nstack(*npole)` each launch.
- **Adaptive launch sizing** for all `_block` kernels (own grid; tunable `NThreadsPerBlock`); do not
  reuse the single-RHS `NBlocks` macro.
- **Seed matvec = specialized two-CSR.** `solve_multishift_block` takes `M_DW`, `M_DWH` (+ lambda_max)
  directly and launches a block CSR matvec kernel TWICE (D_W then D_W^dag), reading the CSRs'
  `d_val`/`d_cols`/`d_rows`. NO edit to shared `sparse_matrix.h`; block kernels live in
  `multishift_kernels_claude.h`. (Seed is always D_W^dag D_W in our code, so no generality lost.)
- **Block outer CG stop = SHARED stop on the worst column:** iterate until
  `max_c ( mu_c / ||b_c||^2 ) < tol^2` (RELATIVE -- handles the small-norm sink RHS without the
  machine-precision floor issue). All columns stay live; freeze-converged deferred to a v2.

- **PROFILE-FIRST GATE (decided):** before writing ANY C6 code, Nsight-Compute the C5 seed `mult`
  kernel to confirm the occupancy/BW headroom justifies the ~1.5x. Build C6a.. only if it does.

## Open questions (deferred / refine later)

1. **Interleaved layout** `[i*nstack+c]` (true SpMM coalescing) -- deferred; column-major first,
   switch only if Nsight shows the block-matvec B-access dominates.
2. **Block dot:** start with `nstack` strided/offset `cublasZdotc` (columns contiguous) vs a custom
   `nstack`-accumulator reduction -- start simple, optimize if profiled.

## Parked alternative -- CUDA-stream parallelism (instead of `_block` kernels)

Considered + PARKED (same modest-gain trade). Within one multi-shift solve there is nothing to
overlap (CG is sequential; the block `x_update`/`p_update` already span ~240 blocks). Across
independent solves there are two options:
- **`mult_ms` || `adj_ms` inside `DHD_ms`** (the two legs are independent) -- cheapest, 2-way; the
  seed matvecs overlap (~12 -> ~24 resident blocks).
- **the `nstack` sink solves on streams** -- stream version of mrhs; same ~1.1-1.5x ceiling.

BLOCKER (why parked): `solve_multishift` is intentionally foreground/default-stream/synchronous, and
its seed matvec calls `on_gpu` which does a per-matvec `cudaMalloc`/`free` = a DEVICE-WIDE sync that
serializes streams. Any stream scheme first needs a stream-clean async `solve_multishift`:
preallocate (kill the per-matvec malloc), async memcpy + stream-bound cublas handle, and PER-STREAM
PINNED coeff buffers for `zeta`/`alm`/`betm` -- i.e. exactly the async-race surface the C2 ASYNC NOTE
warns about. Real race-debugging for a gain in the same band as mrhs => not worth it now. If ever
revisited, `mult_ms || adj_ms` is the cheapest first experiment.

## Notes

- Keep the single-RHS `solve_multishift` (C2) intact as the reference; `solve_multishift_block`
  with `nstack=1` must reproduce it bit-for-bit (first validation).
- A100 later: wider GPU => the occupancy lever matters MORE; keep `nstack`/grid parameterized.
