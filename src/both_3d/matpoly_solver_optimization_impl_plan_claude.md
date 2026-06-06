# MatPoly solver optimization -- impl plan (device-scalar CG + multi-RHS)

> **SUPERSEDED (2026-06-05).** This plan's M1 (device-scalar CG) was implemented then REVERTED
> (it can't help the inner-pole-dominated outer solve). The work pivoted to INNER-POLE MULTI-SHIFT
> CG, now done + validated (4.36x end-to-end). **Canonical plan:
> `inner_pole_batched_solve_impl_plan_claude.md`.** M2 (multi-RHS) lives on there as C6 (multi-RHS
> layered on the multi-shift seed, jj spatial-site sink loop). Kept only for history.

## Goal

Speed up the overlap CG solves (the dominant cost of the JJ / HMC code) by removing the two
distinct sources of GPU under-utilization at the project's small problem size. **No change to
numerical results** -- every milestone is validated by diffing against the current single-RHS,
host-sync solve to solver tolerance.

Two independent, compounding levers (sequenced M1 then M2):

- **M1 -- device-resident CG scalars** (kill the per-iteration host<->GPU sync latency). Lower
  risk, does not touch the operator/data layout. **Do first.**
- **M2 -- multi-RHS block solves + grid growth** (raise occupancy: more resident blocks -> more
  SMs busy -> latency hidden). Higher effort; touches the kernel layout end to end.

## Why (diagnosis, from the perf discussion)

- `N = Comp::N = NS*N_SITES*Nt = 2*12*128 = 3072` at `N_REFINE=1`. A column is ~49 KB; a block of
  `nrhs` columns is ~`nrhs`*49 KB. Both fit in L2 (4.5 MB TITAN V, 40 MB A100) => the matvec/axpy
  are **L2-resident, latency/occupancy/host-sync bound, NOT DRAM-bandwidth bound**.
- `nvidia-smi` util ~60% (= fraction of time a kernel is resident, NOT bandwidth). Streams help
  only ~1.05--1.1x (they overlap different Zolotarev pole solves, but each solve still stalls on
  its own per-iteration host dot).
- Single-RHS launch = ~`N/256 = 12` blocks => ~12 of 80 SMs (TITAN V) / ~12 of 108 (A100) =>
  too few resident warps to hide latency. **A100 is wider => underutilization is worse there =>
  both levers matter MORE on A100.** The host-sync latency (host/PCIe-bound) does not shrink on a
  faster GPU (Amdahl), so M1 is especially important for A100.
- NOTE (user): `NBlocks`/`NThreadsPerBlock` are hand-tuned global constants sized for single-RHS
  `N`; they are NOT optimal and may be freely recomputed/retuned.

## Current CG mechanism (what we are changing) -- `includes/matpoly_claude.h`

`MatPoly::solve<N>(CuC* d_x, const CuC* d_b, tol, maxiter)` is foreground CG. The reductions write
to **host** scalars => `CUBLAS_POINTER_MODE_HOST` => each call blocks until the GPU drains and the
16-byte scalar is copied to host:
```
void dot2self(double* result, const CuC* x, ...) { CuC dummy;                 // HOST var
    cublasZdotc(handle, N, x,1, x,1, &dummy); ... *result = real(dummy); }     // blocks (host ptr)
// in solve():  CuC gam; dot<N>(&gam,d_p,d_q);   const CuC al = mu/gam;        // CPU math
//              Taxpy<CuC,N>(d_x, al, d_p, d_x);  ... dot2self<N>(&mu,d_r);     // 2nd host sync
//              if(mu<mu_crit) break;  const CuC bet = mu/mu_old;              // CPU branch
```
So **2 host round-trips per iteration** + host-side `al`/`bet`/`assert`/convergence branch =>
CPU<->GPU ping-pong, GPU idle during each. `solveAsync` is the same with explicit
`cudaStreamSynchronize` after each `dot*Async`.

Leaf kernels (the things that get widened in M2): `mult<T,N>` (CSR matvec, `sparse_matrix.h`),
`COO::operator()` (matvec, `sparse_matrix.h`), `Taxpy<T,N>` / `Taxpy_gen<...>` (axpy), the cublas
`Zdotc`. The overlap `OverlapWMass::mult/adj/DHD/DDH` (`overlap_wmass_claude.h`) internally run
`~size-1 (=20)` inner pole solves via `MatPoly::solveAsync` + an `A[m]`-weighted `Taxpy_gen` sum;
`ConservedCurrent::apply_k` (`conserved_current_claude.h`) likewise nests pole solves + COO matvecs.
The operator GRAPH (LinOp composition, MatPoly terms, pole sums) is UNCHANGED by both milestones;
we only widen the leaves and the CG scalar bookkeeping.

---

## Milestone 1 -- device-resident CG scalars (host-sync removal)

**Idea:** keep `mu, gam, al, bet, mu_old, b_norm_sq` on the **device**; do the scalar math in tiny
1-thread device kernels; let the axpy read a **device** coefficient; test convergence only **every
K iterations** (copy `mu` to host every K; `K ~ 8--16`). This removes the per-iteration round-trip;
the CG loop becomes device-resident.

**Changes (`includes/matpoly_claude.h`, + a few kernels):**
1. `cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE)` for the solver handle; `cublasZdotc`
   then writes to a device pointer. Add device scalars to `DeviceMemorySetN` (or local cudaMalloc):
   `d_mu, d_gam, d_mu_old, d_bnorm, d_al, d_bet` (each 1x `CuC`/`double`).
2. New kernels (1-thread or tiny): `dev_div(d_out, d_a, d_b)` for `al=mu/gam`, `bet=mu/mu_old`;
   reuse for both. A `Taxpy_devscalar<T,N>(d_y, const CuC* d_al, d_x, d_y)` that reads the device
   coefficient (the current `Taxpy` takes a host `CuC`).
3. Restructure the CG loop: all updates device-side; **no host branch inside the loop**. Every `K`
   iters: `cudaMemcpy(d_mu->host)`, test `mu < mu_crit`, break. (Accept up to `K-1` extra iters
   past convergence; cheap.) Compute `mu_crit` once at start (the one-time `b_norm_sq` host read +
   zero-RHS guard stays -- it's outside the loop).
4. Move the `dot2self` imag/real sanity `assert` out of the hot path (do it once, or under a debug
   macro) -- it currently forces a host read every call.
5. Keep the existing host-sync `solve` as a **reference/fallback path** (rename or `#ifdef`) for
   validation; do not delete until M1 is signed off.

**Files:** `includes/matpoly_claude.h` (solve, solveAsync, dot/dot2self, `DeviceMemorySetN`),
`includes/gpu_header.h` or wherever `Taxpy`/`Taxpy_gen` live (new device-scalar axpy + `dev_div`).

**Validation:** single-RHS solve of `op_Dmsq` (and a CSR) vs the reference host-sync solve --
identical result to `tol`, same iteration count (+/- the K-granularity), wall-time drops. Profile
with Nsight Systems (timeline should densify: gaps between kernels shrink).

**Risks:** off-by-K convergence overshoot; the moved assert; pointer-mode affecting other cublas
calls sharing the handle (use a dedicated handle / restore mode). Mitigated by the reference path.

---

## Milestone 2 -- multi-RHS block solves (occupancy)

**Idea:** solve `A X = B` for `X, B in C^{N x nrhs}` in ONE block solve, and **grow the launch grid
to `N*nrhs`** so the small problem finally fills the SMs. `nrhs` is a **runtime** parameter (not a
template) to avoid template explosion; `NBlocks` is recomputed from `N*nrhs` at each block-kernel
launch (stop reusing the global constant; sweep `NThreadsPerBlock`).

**Layout (decision):** start **column-major** block `X[c*N + i]` (`nrhs` contiguous `N`-vectors):
- columns stay contiguous => per-column cublas uses offset `c*N` stride 1 (no strided penalty),
  host `from_cpu`/`on_gpu` copies stay simple, existing single-vector paths still "see" a vector.
- matvec coalescing is no worse than single-RHS; CSR structure (`rows`,`cols`) read once and reused
  across columns (L2).
- **Interleaved** `X[i*nrhs + c]` (RHS fastest) is the higher-perf alternative (coalesced
  `B[cols[k], :]` chunk per nonzero, true SpMM) but forces rewriting every cublas call + host
  buffers. **Deferred**: switch only if Nsight Compute shows the matvec B-access dominates.

**Grid (the occupancy lever):** do NOT keep `<<<NBlocks, NThreadsPerBlock>>>` (fixed for `N`) and
loop `nrhs` inside a thread -- that gives only launch amortization + L2 reuse, NOT more SMs. Instead
size the grid for `N*nrhs` work: elementwise kernels called with `N' = N*nrhs` (grid auto-grows,
`12 -> 12*nrhs` blocks); matvec as a **2D grid over (rows x nrhs)**; dot as `O(nrhs)` more blocks.
`12 -> ~144` blocks (nrhs=12) oversubscribes 80 SMs (TITAN V) / fills 108 (A100).

**Block kernels (leaves; graph unchanged):**
- block CSR/COO matvec: thread `(i,c)` -> `res[c*N+i] = sum_k val[k]*v[c*N+cols[k]]`; structure read
  once, reused across columns.
- block `Taxpy`/`Taxpy_gen`: elementwise over `N*nrhs`; TWO coeff modes -- (a) per-column coeff
  `coeff[c]` (CG `al`/`bet`, length-`nrhs` device array, composes with M1), (b) single coeff
  broadcast across columns (`MatPoly::on_gpu` term coeffs).
- block dot: column-wise `X^H Y -> nrhs` values. Start with `nrhs` strided/offset `cublasZdotc`
  (columns contiguous -> cheap); optimize to a custom `nrhs`-accumulator reduction later. With M1,
  results go to a length-`nrhs` device array.
- `MatPoly::on_gpu` / `solve`: thread `nrhs` through; scalars become length-`nrhs` device arrays;
  convergence = iterate until the **slowest** column meets its `mu_crit` (v1 keeps all columns live;
  "freeze converged columns" deferred -- it complicates the matvec by shrinking the active set).
- `OverlapWMass::mult/adj/DDH` + `ConservedCurrent::apply_k`: inherit block-ness from the inner
  pole `MatPoly::solve` + block `A[m]` combine -- once `solve`/`on_gpu` are block-aware these come
  almost for free (only the `Taxpy_gen` combines become block).
- `DeviceMemorySetN`: `N*nrhs` work buffers + length-`nrhs` scalar arrays.
- Host: a block `FermionVector` (`N*nrhs`) or `from_cpu`/`on_gpu` overloads taking `nrhs`.

**Consumer:** `jj_conn_tpproj_claude.cu` -- batch the per-site sink solves. Per hit there are
`n_sites = 12` sink solves `psi_n = D_m^{-dag} K^dag(n,0) eta`, all using the SAME `op_Dmsq` with
different RHS `K^dag(n,0) eta` -> a natural `nrhs = 12` block. (The forward `phi'` is a single RHS;
the `(--)` channel mirrors with `op_tilDmsq`.)

**Files:** `includes/matpoly_claude.h`, `includes/sparse_matrix.h` (CSR/COO kernels),
`includes/overlap_wmass_claude.h` (mult/adj/DDH + inner solves), `includes/conserved_current_claude.h`
(apply_k), `includes/gpu_header.h` (kernels, grid macros), `includes/valence_claude.h`
(block FermionVector), `jj_conn_tpproj_claude.cu` (batch sink solves).

**Validation (bottom-up):**
1. block matvec/axpy/dot vs single-RHS, column by column.
2. `op_Dmsq` block solve (random `N x nrhs` B) vs `nrhs` single-RHS solves -- identical to `tol`.
3. through overlap: `D_m^{-1}`/`D_m^{-dag}` block vs single.
4. batch the 12 sink solves in `jj_conn_tpproj`; free-limit re-run -> `Vpp`/`Vmm` unchanged
   (and faster). Diff the `.h5` against a pre-optimization run.

**Risks:** layout/stride bugs; per-column convergence bookkeeping; the deep `apply_k` stack ->
validate incrementally bottom-up; keep single-RHS path available.

---

## Sequencing & checkpoints (each = a compile/run gate for the user)

1. **M1** device-scalar CG (single-RHS), `matpoly` only. Validate + measure. [low risk]
2. **M2a** block matvec/axpy/dot on one LinOp (CSR / `op_Dmsq`). Validate block vs single.
3. **M2b** block through `MatPoly::solve`/`on_gpu`. Validate `op_Dmsq` block solve.
4. **M2c** block through `OverlapWMass` mult/adj. Validate `D_m^{-1}`/`D_m^{-dag}` block.
5. **M2d** batch sink solves in `jj_conn_tpproj`; free-limit re-validate (`Vpp`/`Vmm` identical).
6. (optional) interleaved layout if Nsight says matvec B-access dominates.

Profile before/after each: **Nsight Systems** (timeline gap density) + **Nsight Compute** (matvec/
axpy achieved-vs-peak; occupancy/active-warps).

## Open questions / decisions to confirm before coding

- `nrhs`: tie to the 12 sink solves (or pad to a multiple of warp/occupancy sweet spot)? Runtime.
- Layout: column-major first (decided); interleaved later if profiled.
- M1 convergence cadence `K` (8--16) and whether to expose it as a knob.
- M2 per-column convergence: iterate-to-slowest (v1) vs freeze-converged (deferred).
- Keep the host-sync / single-RHS reference paths during validation (yes) -- remove only after M2d.
- Dedicated cublas handle for `CUBLAS_POINTER_MODE_DEVICE` vs restoring mode on the shared handle.
- `NThreadsPerBlock` retune (sweep 128/256/512) once the grid is `N*nrhs`.

## Notes / pending

- Hardware: target TITAN V now, **A100 later** -- both levers scale up on A100 (wider GPU, same
  host-sync latency). Keep `nrhs`/grid sizing parameterized so A100 can be retuned without code
  changes.
- PENDING (user question, to answer after this draft): do `K` (conserved-current kernel) and
  `D_ov` commute? -- not part of this plan; addressed separately.
