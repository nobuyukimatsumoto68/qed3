# Nf-block HMC force (Phase 2) -- impl plan

**Reference (algorithm):** B. Jegerlehner, hep-lat/9612014 (multishift). Force math = the deployed
serial GRAD_L4 path (`overlap_wmass_claude.h`: `precalc_grad_deviceAsyncLaunch_ms` :719 +
`grad_deviceAsyncLaunch_l4` :940), + the diagonal-mass fix (grad_diag_mass_force_bug_claude.md Sec 5').
Phase 1 (action blocking) is in `hmc_pf_block_parallelize_impl_plan_claude.md`.

## Context / motivation

Phase-1 (action-only) block gave just **~1.1x at L1/Nf6** because it packs only the ~21 action solves
per trajectory, NOT the ~63 FORCE evals (Nf/2 per kick x ~21 kicks). The force is the heavier fermion
cost: per flavor it is **2-3 multishift solves** (`precalc_grad`: Z_m=R_m eta, Y_m=R_m X^dag eta, +
massive W_m) plus a per-link `grad` loop (block matvec over npole + npole dots) over ALL links. Blocking
the force over the Nf/2 flavors is the real Nf=6 lever.

## What blocks, and why it is mostly reuse

The total fermion force per link is `sum_f grad_f(link)`. Both stages of the per-flavor force carry a
clean NSTACK (=Nf/2) axis:

1. **precalc (the expensive part):** `Z_m=R_m eta_f`, `Y_m=R_m (X^dag eta_f)` are same-matrix multishift
   solves differing only by RHS -> **one `BlockedMat::solve_multishift_block` over the NSTACK eta block**
   emits `d_Zg`/`d_Yg` in the `[flavor*npole+pole]` layout the grad kernels already want (NO per-pole
   gather, unlike serial). X-precompute (`X Z_m`, `X Y_m`) = one block `M_DW`-apply over NSTACK*npole cols.
2. **per-link grad:** `mult_coo_block`/`link_matvec_block`/`block_dot`
   (`multishift_block_kernels_claude.h:158,177,209`) all take a runtime `ncol` -> pass **NSTACK*npole**
   instead of `npole`; the `m=0` combine `eta + sum_m A[m] Z_m` is exactly `block_reduce_poles`
   (already used by `BlockedMat::mult`). Host reduction sums `A[m](-Re a_{fm} - Re b_{fm})` over (f,m)
   plus the per-flavor m0 term, then `*= -2 C/lambda_max`.

Massive (Family-B): the bra `(1+m_L)eta` carried through the resolvent -> one extra block multishift
`W_m=R_m X^dag(M_mass eta)`, `hatY = Y + mass_coeff*W`, `eta_bra=(1+m_L)eta` (block versions of the
Sec 5' steps). `M_mass`, `mass_coeff`, `C`, `A[m]`, `cp`, `lambda_max`, `DW.d_coo_format` are all
reachable via the Op (public `struct OverlapWMass`), exactly as `BlockedMat` reaches them.

## Design

Scope clarifications (per discussion 2026-07-09): only the **pseudo-fermion sector** is enlarged
(phi/eta + the force intermediates Z/Y). The **gauge momentum `pi` stays SINGLE** -- all flavor forces
sum into it (`grad_block` returns `sum_f grad_f(link)` per link). The serial `DHD_deviceAsyncLaunch_ms`
/ `precalc_grad_deviceAsyncLaunch_ms` / `grad_l4` are **NOT obsolete**: they remain the validation
reference, the Nf=2 / non-block HMC path, and (DHD_ms) the measurement/valence operator. The block
path is additive machinery beside OverlapWMass, used only by the Nf=4/6 block HMC driver. `BlockedForce`
**reuses `BlockedMat::solve_multishift_block`** (same seed A=(1/lmax^2)M_DWH M_DW, same shifts
-k^2/cp[m]) as the block form of `Aseed.solve_multishift`, and **borrows the action pool's scratch**
(force runs after update_eta, never concurrent) so it owns only the force OUTPUT blocks (~6-7*N*NSTACK
*npole). The block solve emits the `[flavor,pole]` layout directly -> NO per-pole gather.

New helper **`BlockedForce<Idx N, int NSTACK, class Op>`** (new header
`includes/overlap_force_Nfblock_claude.h`), parallel to `BlockedMat` (OverlapWMass + the serial GRAD_L4
stay UNTOUCHED as the reference). Owns the force pool and reuses `BlockedMat::solve_multishift_block`
(pass a `BlockedMat&` or a shared `BlockMemPool`) for the Z/Y/W solves. Methods:
- `precalc_grad_block(U, d_eta_blk)` -> fills `d_Zg,d_Yg,d_XZg,d_XYg` (each N*NSTACK*npole),
  `d_eta_bra_blk` (N*NSTACK); massive-aware.
- `grad_block(link, U)` -> per-link total force scalar (summed over flavors): 2 block matvecs
  (NSTACK*npole) + 2 block_dots + block_reduce_poles m0 term + host reduce.
- `compute(Force& dSf, U)` -> loop links, `dSf.sp/tp = grad_block(...)`; dSf is the TOTAL force.

**`PseudoFermionBlock`** (Phase 1) gains `compute_force(Force& dSf, U)` that owns/drives a
`BlockedForce`. Integrator kick becomes:
`bpf.update_eta(); bpf.compute_force(dSf,U); pi += <coeff>*tau*dSf;` -- ONE force per kick (vs NSTACK).

New integrator `MinimumNorm2BlockF` (side-by-side with the action-only `MinimumNorm2Block`, per the
"preserve variant" convention); driver selects it via `-DBLOCK_FORCE`.

## Memory (the recurring check)

Force pool ~ `6*N*NSTACK*npole` (Zg,Yg,XZg,XYg,CY,CZ) + multishift scratch `3*N*NSTACK*npole` ~
**9*N*NSTACK*npole*16 B**. Worst target L4/Nf6/npole20: N*NSTACK*npole=2.49M -> ~360 MB force +
~145 MB action pool ~ **0.5 GB** (< 12 GB TITAN V). L1 trivial. Startup print + assert vs
`cudaMemGetInfo`, as Phase 1. Option: share the action `BlockMemPool` pole-blocks with the force
multishift scratch (they never run concurrently) to shave ~1/3 -- deferred unless needed.

## Correctness / validation

- **Force-vs-serial at fixed (U, eta):** `BlockedForce::compute` total force == `sum_f grad_l4(eta_f)`
  to ~1e-8 (grad_l4's own block_dot tolerance). Deterministic test driver.
- **Full dH:** block-force trajectory dH == Phase-1 (action-block, serial force) dH to ~1e-6 from an
  identical cold start (same harness pattern as `run_Nfblock_validate_claude.sh`, +`-DBLOCK_FORCE`).
- **Timing (the point):** L4/Nf6 per-traj wall-time: serial vs action-only vs action+force block.

## Chunks

1. DONE (built OK): `includes/overlap_force_Nfblock_claude.h` (`BlockedForce` + `scale_uniform_blk`) +
   `PseudoFermionBlock::compute_force` + `MinimumNorm2BlockF` + driver `-DBLOCK_FORCE` (guarded so
   action-only builds don't allocate the force pool) + memory print + build combos.
2+3 (MERGED): `run_Nfblock_forceval_claude.sh` (USER runs, GPU1) -- builds `-DBLOCK_FORCE` (block force)
   vs default block (action-only, serial force); both share the validated block ACTION, so the FIRST-traj
   dH difference ISOLATES the force -> expect ~1e-6. Same run reports the force-block per-traj speedup at
   **L4/Nf6 (real target, 1 traj)** + L1/Nf6 (few traj). Fresh sibling scratch dirs (no auto-resume,
   production dirs untouched); `nvcc -t 4`.
4. After validation: wire the production `-DBLOCK_FORCE` build targets + update the FNAL blackboard memo.

## Result so far (2026-07-09)

- **Correctness PASS:** block-force vs action-only (serial force) first-traj `|dH|` = **1.5e-9 (L1)**,
  **1.16e-10 (L4)**. Force validated.
- **Timing L4/Nf6 = 2.51x** (serial-force 2754.5 -> block-force 1096.9 s/traj) -- the real Nf=6 win.
  Trusted: long runs (46/18 min) average the `G`-job contention; force-dominated at L4 so 3-flavor
  blocking -> ~3x on the force. **L1 nominal 1.76x = contention artifact** (real ~1x; force is a
  minority at L1 -- action-only was only 1.1x). Caveat: single-traj on a shared (G-only) GPU; a 2-3
  traj avg would pin the exact figure. block+MPS compose (MPS still worth it at L4).
- **Memory:** `npole=10` (n=21); force pool ~17 MiB (L1) / 231 MiB (L4). Non-issue.

## Resolved decisions (2026-07-09)

1. **Housing:** external `BlockedForce<N,NSTACK,Op>` (new header). Rationale beyond style: the block
   scratch is `N*NSTACK*npole`, so whoever owns it must be templated on `NSTACK`. Putting it in
   `OverlapWMass` would force `OverlapWMass<...,NSTACK>` and ripple the template through EVERY consumer
   (jj measurements, prop/eig/glue, ...). The external helper confines compile-time `NSTACK` to the HMC
   block path (same reason `BlockedMat` is external for the action). `OverlapWMass` = operator provider,
   untouched; serial `grad_l4` stays the reference. NOTE `NSTACK` could be made runtime only by
   rewriting `BlockedMat`'s kernels (loses stride specialization) -- out of scope.
2. **Massive INCLUDED now** (Family-B diagonal mass): block `(1+m_L)eta` bra + `W_m=R_m X^dag(M_mass eta)`
   extra block multishift (block form of the Sec 5' fix).
3. **Integrator:** side-by-side `MinimumNorm2BlockF` + driver `-DBLOCK_FORCE` (preserve the action-only
   `MinimumNorm2Block` for A/B).
