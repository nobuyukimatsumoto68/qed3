# M1 implementation -- device-resident CG scalars (impl plan)

This is the working plan for **Milestone 1 only** of
`matpoly_solver_optimization_impl_plan_claude.md`. M2 (multi-RHS block solves) is
explicitly out of scope and is NOT touched here.

## Goal / physics

Remove the per-iteration host<->GPU round-trip in `MatPoly::solve` /
`MatPoly::solveAsync` (`includes/matpoly_claude.h`). The CG scalars
$\mu, \gamma, \alpha, \beta, \mu_\text{old}, \Vert b\Vert^2$ stay on the device;
the scalar math runs in tiny 1-thread device kernels; the axpy reads a device
coefficient; convergence is tested only every $K$ iterations. **No change to
numerical results** -- validated by diffing the new device-scalar solve against
the retained host-sync reference solve.

## Files to modify

- `includes/gpu_header_claude.h` (NEW): small device kernels for the scalar math
  (`dev_div`, `dev_neg`, `dev_copy_scalar`). Included from `matpoly_claude.h`.
  We do NOT mutate the shared `includes/gpu_header.h`; the existing pointer-overload
  `Taxpy<T,N>(T* d_res, const T* d_a, ...)` in `gpu_header.h` already serves as the
  device-scalar axpy required by the plan, so no new axpy kernel is needed.
- `includes/matpoly_claude.h`: add device scalars to `DeviceMemorySetN`; add new
  device-scalar CG paths `solve_dev` / `solveAsync_dev`; keep the existing
  host-sync `solve` / `solveAsync` UNCHANGED as the reference/fallback. Move the
  `dot2self` imag/real sanity assert out of the hot path.

## Design decisions

1. **Pointer mode.** Rather than a dedicated handle, wrap each device-target
   `cublasZdotc` with
   `cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE)` ... dot ...
   `cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST)`. The set-mode call is
   a cheap host-side handle-attribute write (no device sync) and guarantees every
   OTHER cublas call sharing this handle (each stream has its own handle, see
   `overlap_wmass_claude.h:178`) is unaffected. The host-sync reference path keeps
   `CUBLAS_POINTER_MODE_HOST` and is untouched.

2. **Device scalars** live in `DeviceMemorySetN` (one set per stream), so the
   Async path can use them with no per-call alloc. `d_mu, d_gam, d_mu_old,
   d_al, d_bet` are `CuC`; `cublasZdotc` writes a `CuC`. `mu` is mathematically
   real (norm) but stored as `CuC` for uniformity; the host-side convergence read
   takes `real(.)`. The foreground `solve_dev` uses the same per-stream set when
   `istream>=0`, else a local cudaMalloc.

3. **Scalar kernels** (`gpu_header_claude.h`):
   - `dev_div(CuC* d_out, const CuC* d_a, const CuC* d_b)`: `*d_out = *d_a / *d_b`.
     Used for $\alpha=\mu/\gamma$ and $\beta=\mu/\mu_\text{old}$.
   - `dev_neg(CuC* d_out, const CuC* d_a)`: `*d_out = -*d_a`. Used to form $-\alpha$
     for the residual update without a host round-trip.
   - `dev_copy_scalar(CuC* d_out, const CuC* d_a)`: `*d_out = *d_a`. Used for
     $\mu_\text{old} \leftarrow \mu$.

4. **Loop structure.** All updates device-side, NO host branch inside the loop.
   `mu_crit` is computed once at start (one-time `b_norm_sq` host read + zero-RHS
   guard stay outside the loop). Every `K` iters: `cudaMemcpy(d_mu -> host)`,
   test `real(mu) < mu_crit` (and NaN), break. Up to `K-1` extra iters past
   convergence accepted.

5. **K cadence** exposed as a constant `MATPOLY_CG_CHECK_EVERY` (default 8),
   `#ifndef`-guarded so a `.cu` can override before include.

6. **Assert moved out of hot path.** The device-scalar path does NOT run the
   imag/real sanity assert per iteration. It is available via `dot2self` (the
   reference path still uses it). Under `#ifdef MATPOLY_CG_DEBUG_ASSERT` the
   device path performs the check on the convergence-cadence reads only.

## Chunks

- Chunk 1: `gpu_header_claude.h` with the 3 scalar kernels.
- Chunk 2: `DeviceMemorySetN` device scalars (alloc/free) + include hook.
- Chunk 3: `solve_dev` (foreground) device-scalar CG.
- Chunk 4: `solveAsync_dev` (stream) device-scalar CG.
- Chunk 5: validation notes.

## Validation (user runs)

Diff `solve_dev` vs `solve` (and `solveAsync_dev` vs `solveAsync`) on `op_Dmsq`
and a CSR LinOp: identical result to `tol`, same iter count +/- K granularity.
The reference paths are retained so the diff is possible.
