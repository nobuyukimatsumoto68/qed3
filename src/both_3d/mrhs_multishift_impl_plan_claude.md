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

## Profile-first gate (2026-06-06 -- reconsidering mrhs; motivated by jj cost on TITAN V, L=1)

Decided to run the PROFILE-FIRST GATE (never executed) before committing to the `_block` refactor.
Free-field proxy chosen (the per-iteration `mult`:vector-op ratio is ~config-independent; only the
iteration count + fixed-setup fraction differ). NOTE (user): on THERMALIZED configs the multi-shift
gain is at most $\sim1.5\times$ (vs $4.36\times$ free-field), so any mrhs gain stacks on a smaller
base -- but mrhs attacks a different axis (the surviving seed-matvec occupancy), and the
ill-conditioned thermalized seed runs MORE iterations, invoking that matvec more often => the gate,
not the ms number, settles it.

**Seed matvec kernel = `mult<CuC,N>` (`gpu_header.h:101`)**, one thread/row, `<<<N/256, 256>>>` =
12 blocks at $L{=}1$ ($N{=}3072$) on an 80-SM (TITAN V) / 108-SM (A100) GPU => the latency-starved
target mrhs would widen to $N\times$nstack.

**Confound found: `MatPoly::on_gpu` (`matpoly_claude.h:108`) does per-matvec 2x `cudaMalloc` +
`cudaMemset` + D2D `cudaMemcpy`s + 2x `cudaFree`** (device-wide syncs). This API/memcpy overhead may
dominate the tiny 12-block `mult` wall time and is fixable by HOISTING the scratch -- a clean win for
EVERY solve, no async-race, no block refactor. Must be separated from the occupancy question.

Gate commands (user runs; free-field `test_overlap_multishift` binary):
```
nsys profile -o ms_prof --stats=true ./<binary>
nsys stats --report cuda_gpu_kern_sum ms_prof.nsys-rep   # mult vs Taxpy vs multishift_* vs cublas dot
nsys stats --report cuda_api_sum     ms_prof.nsys-rep   # cudaMalloc/cudaFree/cudaMemcpy/cudaMemset
ncu --kernel-name-base demangled -k "regex:^mult<" -c 20 -s 200 \
    --section Occupancy --section SpeedOfLight ./<binary>
```
Decision rule: (a) `cudaMalloc/Free/Memcpy` large -> do the prealloc fix first, re-profile (best ROI);
(b) `mult` large slice AND low achieved occupancy + low SM/DRAM throughput -> mrhs justified, ceiling
$\approx$ (`mult` share) $\times$ (matvec speedup, saturating ~nstack 6-8); (c) vector ops dominate +
well-occupied -> skip mrhs.

**`on_gpu_pre` DRAFTED side-by-side (2026-06-06, `matpoly_claude.h`):** new method `on_gpu_pre(d_v,
d_v0, d_tmp, d_Mv0)` (caller-supplied scratch; body byte-identical to `on_gpu` minus malloc/free);
original `on_gpu` UNTOUCHED. `solve_multishift` allocates `d_tmp_mv`/`d_Mv0_mv` once and frees them;
the hot seed matvec keeps `on_gpu` ACTIVE with the `on_gpu_pre` line COMMENTED beneath = one-toggle
A/B (default unchanged so the gate measures the CURRENT overhead). Verbose-only seed matvecs left on
`on_gpu`.

### Gate results (2026-06-06, TITAN V, free-field, `test_overlap_multishift`)

`nsys cuda_gpu_kern_sum`: seed `mult<double2,3072>` = **54.0%** of GPU kernel time, 6824 inst,
**avg 9228 ns for a 12-block launch** (saturated SpMV this size should be ~1-2 us => badly
occupancy-bound). `Taxpy` 19.3%, cublas `dot_kernel` 12.7%, `multishift_x/p_update` 0.4% each.
Per-multishift-iteration accounting (matvec 2x`mult`~18 us vs ~4 `Taxpy` + 2 dots + 2 updates ~16 us)
reproduces 54% => the matvec is ~54% of the solve, and that ratio is CONFIG-INDEPENDENT (holds on
thermalized configs, where only the iteration count grows). **mrhs solve-level ceiling at $L{=}1$ is
$\sim1.7$-$1.9\times$, NOT the earlier 1.1-1.5x** -- the matvec slice is bigger than assumed.

`nsys cuda_api_sum`: `cudaMalloc` 22% (skewed by one 161 ms startup alloc; median 2.8 us),
`cudaMemcpyAsync` 20%, `cudaMemcpy` 8.7%, `cudaFree` 3.8%. `on_gpu`'s per-matvec 2x`cudaMalloc` +
2x`cudaFree` are device-wide syncs that bubble between matvecs => `on_gpu_pre` removes them
(zero-risk ~1.2-1.3x, no block refactor).

`mult_ms` vs pole-loop still PASS (rel ~1e-16); D_ov 4-6x, D_ov^dag ~6x.

**A/B done (2026-06-06):** `on_gpu_pre` removed only ~370 `cudaMalloc`/`cudaFree` (the seed matvecs);
both API totals are dominated by a single ~120-160 ms one-time startup alloc, so the per-matvec malloc
cost was ~1-2 ms and the wall-clock `mult_ms` change (0.01202 -> ~0.0106-0.0111 s) is noise =>
**prealloc is a FREE-FIELD WASH.** It keeps the inherent D2D copies (can't remove). REVERTED to OFF
(original `on_gpu` default; `on_gpu_pre` + the `d_tmp_mv`/`d_Mv0_mv` alloc kept side-by-side for a
later thermalized A/B, where more iterations may expose the sync bubbles). Refined `cuda_gpu_kern_sum`
read: ~87% of GPU time is small occupancy-starved kernels (`mult` 54% @ 9.2 us/12-block, `Taxpy` 19%
@ 1.7 us/12-block, `dot` 13%) -- ALL batch well over nstack => mrhs ceiling for the jj sink loop is
plausibly **~2-4x per-RHS** (until BW saturates ~nstack 6-8 on 80 SMs), higher than the parked
1.1-1.5x.

**ncu DONE (2026-06-06, `sudo ncu`, `mult<double2,3072>` grid 13x256):** Compute(SM) Throughput
**0.79%**, Memory/DRAM Throughput **13.06%**, Achieved Occupancy **12.09%** (theoretical 100%),
Duration 12.3 us. => Occupancy-starved + ~7x DRAM-BW headroom: 13 blocks light ~13/80 SMs. Batching
`nstack` RHS fills the SMs, raising DRAM toward saturation; matvec per-RHS cost drops ~5-6x until BW
saturates (~nstack 7-8; n_sites=12 at L=1 more than fills it). Amdahl on the solve (87% occupancy-
bound) => realistic **~2-4x per-RHS for the jj sink loop**.

**GATE PASSED -> PROCEED with the C6a-e `_block` build** (bottom-up, each chunk a compile gate;
specialized two-CSR seed; shared worst-column stop; originals + `_ms` kept side-by-side). Prealloc
(`on_gpu_pre`) stays OFF/side-by-side (free-field wash).

## C6f -- speed up the sink K applies (DETAILED PLAN, 2026-06-06; the ~90%-of-jj lever)

Reference (algorithm): conserved current K eq. (3.34) `conserved_current_claude.h`; inner multishift
= B. Jegerlehner, hep-lat/9612014. Goal: cut the jj sink passes (`Nt` x `n_sites`/`n_links` x 4
channels ~ 10.7k `K(.,t)phi'` applies ~ 90% of the 1056s/hit). The current sink loop calls
`op_K.from_cpu(kphi, phi')` per `(n,t)` -- i.e. a full `apply_k` each time.

### Anatomy of `apply_k(d_result, d_xi, U, el)` (conserved_current_claude.h:132) -- what depends on what
The link element `el=(t,n)` (the current insertion) enters ONLY through the COO `W^{wz}` (`build_W`,
`:157`). The inner POLE SOLVES split into two sets:
- **Step 1 (`:138-150`): `Z_m = R_m xi` (m=1..size-1), `Z_0 = xi + sum A[m] Z_m`.** `R_m = (X^dag X
  - k^2/cp[m])^{-1}`, `X=D_W/lambda_max`. **xi-ONLY (el-INDEPENDENT).** Also `X Z_m` (`:185`) is
  xi-only.
- **Term A (`:168`): `C W Z_0`** -- el-dependent (W), reuses Z_0. Cheap (1 COO matvec).
- **Term B (`:177-198`, loop m): `-C sum_m A[m] X R_m u_m`, `u_m = X^dag W Z_m - W^dag X Z_m`.**
  `u_m` is **el-DEPENDENT** (W), so the `Op.solve` at `:193` (`R_m u_m`) is a SECOND pole-solve set
  that is **el-dependent => NOT cacheable**, plus the COO matvecs `W Z_m`/`W^dag X Z_m` (el-dependent).
=> `apply_k` ~ `2*(size-1)` pole solves (Step1 + TermB) + ~`(size-1)` COO + Wilson matvecs.

### THE APPROACH (user 2026-06-06): BLOCK THE `t` LOOP (`nstack = Nt`), per fixed site/link `n`
For a fixed `n`, the `Nt` applies `K(n,t)phi'` (t=0..Nt-1) share `xi=phi'` but differ in `el=(t,n)`
(different `W^{wz}` per `t`). Process all `Nt` at once. Why `t` not `n`: (a) `nstack=Nt`=128 is a big
SM-filling batch (vs `n_sites`=12); (b) Step 1 is xi-only so it's computed ONCE per `t`-block = once
per `n` => `n_sites` Step-1 sets total, and since `n_sites`(12) < `Nt`(128) that's FEWER redundant
Step-1 solves than n-blocking would give. So `t`-blocking captures the Step-1 sharing AND the big
batch in one restructure.

A `block_t` sink apply for fixed `n` (xi=phi', the `Nt` el's = `(0..Nt-1, n)`):
- **Step 1 once:** `Z_m = R_m phi'` (xi-only) -- single (size-1) pole-solve set for the whole `t`-block.
- **Term A / Term B block over `t`:** build the `Nt` `W^{wz}(t,n)` (COO per t); Term B's per-`el` inner
  solve `R_m u_m(t)` (`u_m` el-dependent) batches over `t` => one block multishift `nstack=Nt` via
  `BlockedMat::solve_multishift_block` (assemble the `Nt` `u_m` block, one solve, scatter). The
  el-dependent COO matvecs `W(t) Z_m` / `W^dag(t) X Z_m` are a stack over `t` (block-COO / loop).
- Output: `kphi(n,t)` for all `t` (an `N*Nt` block).

### Restructure of the jj sink pass: outer-`t` inner-`n`  ->  outer-`n`, block over `t`
- The accumulators that currently sum over `n` at fixed `t` must hold all `t`: ylm `PhiL[l]` ->
  `PhiL[l][t]` (Nt vectors/l, ~`Nt*N_ELL*N` ~ 19 MB at L1, held across the `n` loop); disc
  `Jtp(t)`/`Jsp(t)` accumulate `+= w[n] eta.dag(kphi(n,t))` over the `t`-block per `n`; connected
  `Ctp[b][dt] += w[n] psi[n].dag(kphi(n,t))`, `dt=(t-t0_b)%Nt`, looped over the `t`-block per `n`.
- Per channel/dag the Step-1 vector is `Z_m=R_m xi` (K) or `Y_m=R_m X^dag xi` (K^dag); `xi` = the
  channel's sink vector (`phi'`, `phimm`, `chi`, ...). Disc traces ride `kphi` unchanged (bit-identical).
- OPTIONAL further: since `xi=phi'` is the SAME across all `n` too, Step 1 could be cached across the
  `n` loop (`n_sites` sets -> 1). Defer; the `t`-block already cuts Step-1 by `Nt`x.

### Chunks (each a compile/run gate; bit-identical to current jj_corr_mrhs = the bar, via h5diff)
- **C6f-a (profile/size):** time one `apply_k` -- Step1 vs TermB pole solves vs COO/Wilson matvecs. Sets
  the ceiling (Step1 ~ half the pole solves; Term-B block over `Nt` is the occupancy lever on the rest).
- **C6f-b (`block_t` apply on ConservedCurrent):** new `apply_k_block_t<NSTACK=Nt>(d_kphi_block, d_xi,
  U, n, dag)` -- Step1 once + Term A/B over the `Nt` el's (Term-B inner solve via
  `BlockedMat::solve_multishift_block`, `nstack=Nt`). Side-by-side; original `apply_k` untouched.
  Unit-test vs the per-`t` `apply_k` (bit-identical).
- **C6f-c (jj sink rewrite):** restructure the 4 channel sink passes (vec tp/sp ++/--, axial tp/sp) to
  outer-`n` + `apply_k_block_t` + per-`t` accumulators. Validate free-field h5diff vs jj_corr_mrhs.
- **C6f-d (optional):** cache Step 1 across `n`; block-COO matvecs -- only if profiled to matter.

### Open questions
- Abstraction: `block_t` bypasses `op_K.from_cpu` -> jj calls `ConservedCurrent::apply_k_block_t`
  directly (it holds `kop`). Confirm no other sink consumer needs op_K.
- `set_temporal/set_spatial` (`:69-80`) only configure `cfg_el`/gauge for `build_W` => Step 1 truly
  el-independent. Verify before relying on Step-1 sharing.
- sp (link) channels: identical with `build_W_wz`, block over `t` (`nstack=Nt`, same as tp).
- Memory: the `N*Nt` `kphi` block + the per-`t` `PhiL[l][t]` buffers + one BlockedMat (`nstack=Nt`)
  scratch; check the L1 budget (Nt=128 block buffers are the new sizable allocation).
- Term B's COO matvecs (`W(t) Z_m`) over the `t`-block: stack of distinct COOs -> block-COO matvec or a
  `t`-loop of single COO matvecs (cheap relative to the batched pole solve); decide in C6f-b.

### REFINEMENT (2026-06-06, from reading apply_k_ms) -- Term B needs a SINGLE-SHIFT block CG
Cost re-check: `apply_k` Step 1 = ONE multishift pass (~1 Krylov solve, all poles share the seed
Krylov sequence). Term B = `npole` INDEPENDENT single-shift CG solves (`conserved_current_claude.h:354`
`Op.solve` per `m`). So Term B ~ `npole` x Step 1 => **Term B dominates `apply_k`**; the t-block lever
is batching TERM B over `t`, and Step-1 sharing is a minor (~`1/(npole+1)`) add-on.

Why NOT `solve_multishift_block` for Term B: Term B's RHS
$$u_m(t) = X^\dagger W(t) Z_m - W^\dagger(t)\, X Z_m$$
depends on BOTH pole `m` and time `t`. For a fixed `m`, the `Nt` RHS $\{u_m(t)\}_t$ share the single
shift $\sigma_m=-k^2/c_p[m]$ => a SINGLE-SHIFT block CG (block over `t`, loop `m`). `solve_multishift_block`
shares ONE RHS across ALL poles -- the wrong structure here. So C6f-b needs a NEW primitive:
- **`BlockedMat::solve_shift_block(d_X, d_B, double sigma, tol)`** -- solves $(A+\sigma)X_c=b_c$ for the
  `NSTACK` columns, $A=(1/\lambda_{max}^2)M_{DW}^\dagger M_{DW}$ (same seed matvec as
  `solve_multishift_block`: apply `M_DW` then `M_DWH`), real coeffs, per-column freeze-converged,
  shared worst-column stop. This is the real-coeff block CG removed in A1, re-added parameterized by an
  explicit single `sigma`. nstack=1 => bit-identical to the single `MatPoly::solve` (the validation bar).

Term B per fixed site/link `n`, over the `Nt` el's `(t,n)`:
- precompute (t-indep) $Z_m$ (Step 1, single-vector `solve_multishift`, once per `n`) and $XZ_m=X Z_m$;
- per `m`: build the `Nt`-wide RHS block $B_m[t]=u_m(t)=X^\dagger(W(t)Z_m)-W^\dagger(t)\,XZ_m$ (Term-B COO
  matvecs over `t`), `solve_shift_block(sigma_m)` -> $V_m[t]=R_m u_m(t)$, then accumulate
  `result[t] -= C A[m] (X V_m[t])`;
- Term A: `result[t] = C W(t) Z_0` (cheap `Nt` COO matvecs).
K^dag (`apply_k_dag_block_t`) is the operator-adjoint mirror: Step-1 seed $Y_m=R_m X^\dagger\xi$, Term
A^dag block over `t` (multishift, RHS `w0(t)=W^dag(t) xi` -- t-DEPENDENT here, so its inner solve also
batches over `t` via `solve_shift_block` per `m`), Term B^dag block over `t`.

### Multishift vs block are ORTHOGONAL -- `solve_shift_block` washes out NOTHING
A natural worry: by making Term B a SINGLE-shift block solve, did we throw away the multishift gain
(`npole` Krylov sequences collapsed to 1)? No. The two levers attach to two DIFFERENT pieces of
`apply_k`, because they exploit two different kinds of RHS sharing:

| Piece of `apply_k` | RHS structure | Applicable lever | how C6f handles it |
|---|---|---|---|
| **Step 1** $Z_m=R_m\xi$ (`:309`) | ONE RHS $\xi$, shared across ALL poles | **multishift** (poles share one Krylov sequence; `npole` seqs -> 1) | kept as a single-vector `MatPoly::solve_multishift`; just computed once per `n` instead of once per `(n,t)` -- strictly better |
| **Term B** $R_m u_m(t)$ (`:354`) | RHS varies per pole `m` AND per `t` | **block** (independent RHS batched for occupancy) | `solve_shift_block`, batched over the `Nt`-wide `t` block, looped over `m` |

Crucial fact about the CURRENT `apply_k_ms`: **Term B already does `npole` INDEPENDENT single-shift
solves** (the `for m` loop of `Op.solve` at `:351-354`) -- it was NEVER multishift, because $u_m$ is
pole-dependent so there is no shared RHS to multishift over. So `solve_shift_block` does not remove a
multishift gain that Term B never had; it ADDS an occupancy (block) gain on top of Term B's irreducible
`npole` solves, while leaving Step 1's multishift fully intact.

- **multishift** = poles share a Krylov sequence (needs ONE RHS for all shifts) -> Step 1.
- **block** = independent RHS batched to fill the SMs (any RHS set at a fixed shift) -> Term B (and
  Step 1 too, if ever blocked).

You would only wash out the multishift gain if you replaced Step 1's multishift with `npole`
single-shift block solves -- which we explicitly do NOT do. **C6f-b1 invariant:** Step 1 in
`apply_k_block_t` stays a single-vector `MatPoly::solve_multishift` (RHS $=\phi'$), NOT
`solve_shift_block`, so the `npole`->1 Krylov collapse is preserved there.

### REFINED chunk order (supersedes the C6f-b bullet above)
- **C6f-b0 [DONE+VALIDATED 2026-06-06]:** `BlockedMat::solve_shift_block` (`blocked_mat_claude.h:170`)
  single-shift real-coeff block CG over `A+sigma`, per-column freeze + worst-column stop. Test
  `test_shift_block_claude.cu` (seed+mid poles, NSTACK=1/4/Nt): ALL PASS, max|block-single| 6e-12 (seed)
  / 1e-14 (mid); speedup **2.44x** (seed) / **2.47x** (mid) at NSTACK=Nt=128, ~1.8-2x at NSTACK=4,
  ~1x at NSTACK=1 (overhead). Ran on GPU 1 via `tmp_claude2.sh`. This is the occupancy gain for Term B's
  inner solve; composes with Step-1 multishift.
- **C6f-b0.5 [DONE 2026-06-06] BlockMemPool + PoolGuard:** extracted the device scratch out of
  `BlockedMat` into `BlockMemPool<N,NSTACK>` (RAII owner; `with_pole_blocks` gates the three big
  `N*NSTACK*npole` arrays = ~189 MB at L1, Nt=128 [NOT 1.9 GB -- earlier 10x arithmetic error]; the
  K^dag block path Term A^dag needs them via solve_multishift_block, so the wrapper uses a FULL pool).
  `BlockedMat` now holds ALIAS
  pointers + a `mempool` ptr: 2-arg ctor `(Op&, BlockMemPool&)` SHARES an external pool, 1-arg ctor
  `(Op&)` owns a full pool (backward compatible -- existing call sites unchanged). Overwrite hazard from
  two clients sharing a pool is caught by `PoolGuard` (RAII lease; owner-tracked depth so self-nesting
  `solve_sq->DDH->mult->solve_multishift_block` does not false-trip; `assert`-only => zero release cost,
  bit-identity untouched) at the top of every public method. Re-validated via C6a-d + b0 (bit-identical).
- **C6f-b1 [K DONE+VALIDATED 2026-06-06]:** `ConservedCurrentBlockT<OverlapOp,Gauge,NSTACK>` wrapper in
  NEW header `includes/conserved_current_block_claude.h` (NOT in conserved_current_claude.h, so that
  header's other consumers are unaffected; it references BlockedMat + block kernels -> include AFTER
  blocked_mat + conserved_current). Owns a LEAN `BlockMemPool<N,NSTACK>(npole,false)` + `BlockedMat` +
  `KBlockScratch` (d_B/d_V/d_MV = N*NSTACK, d_XZ = N). `apply_k_block_t<FixedLink>` (FixedLink=Idx temporal
  / BaseLink spatial): Step1 ONCE via single-vector `MatPoly::solve_multishift` (the b1 INVARIANT); Term A
  loop t; Term B loop m (XZ_m=X Z_m once; assemble Nt-wide u_m(t); `blk.solve_shift_block(sigma_m)`;
  `result -= C A[m] X V_m`). Nt single-link W(t)/W^dag(t) COOs precomputed once via COUNT-constructed
  `std::vector<COO<N>>(NSTACK)` (no realloc => no COO copy/move => no COO edit needed; user OK'd a move
  ctor but count-ctor sidesteps it). Test `test_apply_k_block_t_claude.cu` vs per-t apply_k_ms: TEMPORAL
  2.7e-15 / SPATIAL 3.4e-16 (machine-exact; only Term B's solve_shift_block differs ~1e-11), **3.08-3.10x
  @ NSTACK=Nt=128** (beats b0's 2.44x -- also shares Step1 once vs Nt times). K^dag
  (`apply_k_dag_block_t`) DONE (NOT asymmetric with K: Term A^dag carries S=I+sum A_m R_m on the
  t-DEPENDENT w0(t)=W^dag(t) xi => per-t MULTISHIFT => batches via solve_multishift_block + block_reduce_
  poles, hence FULL pool; Term B^dag = solve_shift_block like K). Test (K+K^dag, temporal+spatial) ALL
  PASS: K 2.7e-15/3.4e-16 @ 3.00x/3.05x; K^dag 1.0e-15/1.9e-16 @ 2.83x/2.87x. C6f-b1 COMPLETE.
- **C6f-c [DETAILED PLAN below]:** jj sink rewrite (outer-`n`, block-over-`t`, per-`t` accumulators) in a
  copy `jj_corr_block_t_claude.cu` + free-field h5diff (WITH TOLERANCE) vs `jj_corr_mrhs`.
- **C6f-d (optional):** cache Step 1 across `n`; block-COO matvecs -- only if profiled to matter.

### C6f-c -- jj sink rewrite (DETAILED, 2026-06-06)
New file `jj_corr_block_t_claude.cu` = copy of `jj_corr_mrhs_claude.cu` (KEEPS the mrhs SOURCE-solve
batching `blk_*.solve_sq_from_cpu`; adds the SINK block-t). One wrapper `ConservedCurrentBlockT<Fermion,
Gauge,Comp::Nt> kblock(kop)` (full pool ~189 MB + KBlockScratch ~19 MB), created once; the existing
`kop`/`op_K` stay for the source `K^dag(.,t0)` applies (single, NOT block).

The SEVEN sink passes (each currently `for t: for n/a: op_K.from_cpu(kphi,sinkvec)`):
| # | pass | sink vec | dag | accumulators |
|---|------|----------|-----|--------------|
| 1 | (++) tp+ylm | phi'  | F | Jtp[t], Ctp[b][dt], PhiL[l]->Jyl/Gyl |
| 2 | (--) tp+ylm | phimm | T | Ctp[b][dt], PhiL[l]->Gyl |
| 3 | (++) sp     | phi'  | F | Jsp[t], Csp[b][dt] |
| 4 | (--) sp     | phimm | T | Csp[b][dt] |
| 5 | disc(--) tilde | tilphi | F | JtpT[t],JylT[l][t],JspT[t] |
| 6 | axial tp    | chi   | T | Atp[b][dt] |
| 7 | axial sp    | chi   | T | Asp[b][dt] |

Each pass -> `for n/a: kblock.apply_k[/dag]_block_t(d_kphi_block, sinkvec, U, fixed); D2H to host kblk
(N*Nt); for t: <existing per-t accumulation on column t>`. The accumulation (host dots
`FermionVector.dag`, axpy) is UNCHANGED -- copy column t into `kphi.field` then run the exact existing
lines -> the ONLY numerical difference from `jj_corr_mrhs` is `apply_k_block_t` vs `op_K` (~1e-15, the
Term-B `solve_shift_block`). dag=F -> `apply_k_block_t`; dag=T -> `apply_k_dag_block_t`.

YLM restructure (passes 1,2 only): `PhiL[l]` (summed over `n` at fixed `t`) -> `PhiL[l][t]` (std::array of
N_ELL std::vector<FermionVector>(Nt) or N_ELL*Nt FermionVectors), accumulated `PhiL[l][t] += W_ell[l][n]
kphi_col(t)` across the `n` loop; the Jyl/Gyl reductions move to a post-`n` `for t` loop. ~N_ELL*Nt*N =
19 MB.

VALIDATION: h5diff is NOT bit-exact here (block K differs ~1e-15) -> use a tolerance:
`h5diff -d 1e-10` (absolute) or `-p 1e-12` (relative) vs `jj_corr_mrhs` free-field, 1 hit. Expect all
datasets within tol.

CHUNKS: (c1) pass 1 [(++) tp+ylm, the PhiL[l][t] pattern-setter] -> h5diff (only tp/ylm Vpp + disc
tp/ylm differ ~1e-15, rest identical). (c2) passes 3,6,7 [sp + axial, simple Ctp-style]. (c3) passes
2,4,5 [parity (--) + disc tilde]. Then full free-field h5diff. Then repoint nothing yet (B deploy uses
jj_corr_mrhs; jj_corr_block_t is the C6f product, benchmarked separately for the whole-jj speedup).

## POST-C6f-c ROADMAP (2026-06-06) -- two fast products now exist; what next

State: mrhs source batching LANDED (`jj_corr_mrhs`, bit-identical). C6f sink t-block: b0/b0.5/b1 DONE+
validated (`solve_shift_block` 2.4x; `BlockMemPool`+`PoolGuard`; `apply_k[_dag]_block_t` 2.8-3.1x).
C6f-c (`jj_corr_block_t`, all 7 sink passes) **PASS 2026-06-06**: h5diff -d 1e-10 within tol (==
jj_corr_mrhs) free-field; **whole-jj 1054.2 -> 352.1 s/hit = 2.99x** (tp+ylm sink 138->48, sp 490->118,
axial tp 161->53, axial sp 397->133). The headline C6f win; `jj_corr_block_t` = production jj binary.
`jj_disc_block_t` (4 sink loops) IMPLEMENTED, validating (tmp_claude3, GPU 0).
B1a (`jj_corr_mrhs` vs `jj_corr` thermalized gate) **PASS 2026-06-06**: h5diff BIT-IDENTICAL on
`Nf2_gsq8...pole11/ckpoint_lat.224`. Thermalized hit times jj_corr 2863 s vs jj_corr_mrhs 2764 s =
**1.036x** (~3.5%) -- confirms K-APPLY-BOUND on a real config (source batching ~3% whole-jj; sink block-t
= the lever). mrhs deploy gate GREEN.

Two fast products: (1) `jj_corr_mrhs` = source-solve batching only, BIT-IDENTICAL to `jj_corr`. (2)
`jj_corr_block_t` = source + SINK block-t, ~1e-15 from `jj_corr_mrhs`. block_t is the headline jj win
(sink K applies were ~90%; block-t ~3x on them => projected ~2-2.5x whole-jj -- the C6f-c log gives the
FIRST measured number).

### Prioritized next steps (pick order)
**P0 -- block-t PROPAGATION SCOPE [CORRECTED 2026-06-06]: ONLY `jj_corr` (= `jj_corr_block_t`, DONE, 2.99x).**
`jj_corr` is the COMBINED program -- it ALREADY computes disc (the J(t) traces ride the connected (++)
sink passes, written under h5 disc/{tp,sp,ylm}/J + Jtil; block-t carries this through). So the --all
production (`jj_corr_block_t`) gives block-t disc for free. The standalone `jj_disc` (--disc) is only the
disc-ONLY cheap many-config path; its block-t copy `jj_disc_block_t_claude.cu` was written + builds clean
but is REDUNDANT for --all => PARKED, not deployed (run script `DISC` left = `jj_disc_claude.o`). The
old/standalone conn programs (`jj_conn_correlators`, `jj_conn_tpproj/spproj/axial/ylm`) are NOT production
-> untouched. Run script repointed: `CORR=jj_corr_block_t_claude.o` (run_jj_analysis_claude.sh:53).

[SUPERSEDED -- jj_disc not needed] earlier scope had also said "add the block-t sink to `jj_disc_claude.cu`
(disc-only: J(t)=sum_a w_a eta^dag K(a,t) phi', VECTOR only, tp/sp/ylm, NO connected source legs => much
simpler than jj_corr -- one sink vector phi', the same per-(n,t) K applies, accumulate eta^dag kphi).
Same block-t pattern: outer-n/a `apply_k_block_t(d_kphi_block, phi', U, fixed)` + D2H + per-t
`J*[t] += w eta^dag kcol(t)`. Gate after the C6f-c speedup is known.

### jj_disc block-t PRE-SCOPE (2026-06-06, read of jj_disc_claude.cu)
SIMPLER than jj_corr. FOUR sink loops, ALL already outer-n/a inner-t, ALL dag=false (forward K) =>
ONLY `apply_k_block_t` needed (no K^dag). Loops: (1) vec tp `for n: for t: T=eta^dag K(n,t)phi';
Jtp[t]+=w_tp[n]T; Jyl[l][t]+=W_ell[l][n]T` (335-341); (2) vec sp `for a: for t: Jsp[t]+=w_sp[il]
eta^dag K(lk,t)phi'` (342-350); (3) parity tp `K(n,t)tilphi; TtilNT=kphi.dag(eta); JtpT/JylT` (362-369);
(4) parity sp `JspT[t]+=w_sp[il] (K(lk,t)tilphi).dag(eta)` (370-378). KEY: disc ylm is a SCALAR sum
`Jyl[l][t]+=W_ell[l][n]Tnt` (Tnt scalar) => NO PhiLt[l][t] restructure (jj_corr needed it only because
its ylm pairs with a source leg). Trivial conversion: swap inner `for t: kop.set;op_K.from_cpu` for
`apply_k_block_t(d_kphi_block,d_sinkvec,U,fixed); D2H; for t: kcol(t); <same accumulation>`. sinkvec =
phi' (vec) / tilphi (parity).
PLUMBING: jj_disc lacks blocked_mat_claude.h -> add it + conserved_current_block_claude.h after
overlap_wmass (multishift_block_kernels already via matpoly_claude.h). Add ConservedCurrentBlockT kblock
(kop) + d_sinkvec/d_kphi_block/kblk + free. Output dir unchanged (disc_nhits<H>/disc.<k>.h5).
FILE: make `jj_disc_block_t_claude.cu` (copy; keeps jj_disc as h5diff reference, matches jj_corr_block_t)
[RECOMMENDED] vs edit in place [user to confirm]. VALIDATE: h5diff -d 1e-10 jj_disc vs jj_disc_block_t
free-field (vector + parity via --mass-im).

**P1 -- LAND + DEPLOY block_t (the headline win).** After C6f-c h5diff passes:
  - P1a. Read the C6f-c log's REF-vs-BLOCK_T wall times => the measured whole-jj speedup (the deploy
    justification). If ~2x+, block_t becomes the production jj binary.
  - P1b. THERMALIZED gate for block_t: h5diff -d 1e-10 `jj_corr_block_t` vs `jj_corr_mrhs` on one real
    config (mirror B1a's single-symlink trick). Confirms block-t on a non-trivial U.
  - P1c. Repoint `run_jj_analysis_claude.sh` / `run_connected` to `jj_corr_block_t_claude.o` (keep
    `jj_corr` / `jj_corr_mrhs` as reference fallbacks). Multi-GPU split already supported (--gpu 0,1).
  - P1d. Run the ensembles (Nf{2,4,6} gsq8, the mass matrix) with block_t.

**P2 -- C6f-d block-t refinements (ONLY if P1a's whole-jj speedup disappoints, e.g. <2x).** After
block-t, re-profile what now dominates a hit (the sink K is cut ~3x; the next bottleneck shifts). Likely
candidates, in order:
  - Step-1 cache across `n`: `Z_m = R_m phi'` is xi-only and `phi'` is shared across ALL `n` in a pass,
    so Step 1 (a single-vector multishift) could be computed ONCE per pass instead of once per `n`
    (`n_sites`/`n_links` -> 1). Needs the wrapper to expose a "set source xi" that precomputes Z_m/Y_m,
    then `apply_k_block_t` skips Step 1. Saves `~n-1` multishift passes per pass.
  - Host-side accumulation: after block-t the per-(n,t) host dots (`FermionVector.dag`) + PhiLt axpy may
    be a non-trivial fraction (they run on host after a D2H). If profiled significant, move the dots to
    device (cublasZdotc over the N*Nt block vs eta/psi on device) -- the D2H then drops too.
  - COO build cost: `Nt` single-link `W(t)/W^dag(t)` COOs rebuilt per `n` (`do_it` host loops + device
    allocs). Precompute/cache per pass or a custom single-link matvec kernel.

**P3 -- HMC (the other big consumer), D1 + C1.** Independent of jj. HMC action solve already multishift
(1.26x). D1 = batch a few pseudofermion SETS as `nstack` RHS via BlockedMat in the action/force
`op_Dmsq` solve. C1 = fused block-dot kernel (replace the `2*NSTACK` host `cublasZdotc`/iter) -- LOW
value for jj (K-apply bound) but the solve IS the HMC cost, so do C1 WITH D1.

**Parked:** B2 (L=2 memory -- per-mass-case instances / shared scratch when L=2 lands); D3 (standalone jj
programs rewire -- `jj_corr`/`jj_corr_block_t` is the production path); C2 (launch-config retune on A100);
the `COO` move-ctor hardening (count-ctor sidesteps it).

## NEXT STEPS (roadmap, 2026-06-06) -- after C7 BlockedMat validation

NOTE: class is now `BlockedMat<Idx N, int NSTACK, class Op>` (renamed from MatPolyBlocked -- it's a
blocked-matrix engine, not a MatPoly). Outer squared-op solve = `solve_sq` / `solve_sq_from_cpu`
(host-block). jj instances `blk_{tp,sp}_{D,Dm,Dtil}` (width-first; tp=N_SITES, sp=N_LINKS).

**C7/C6e VALIDATED 2026-06-06: jj h5diff = IDENTICAL** (jj_corr_mrhs BlockedMat output == single-RHS
jj_corr, free-field 1 hit) + C6a-d bit-identical. The mrhs migration is complete & correct end-to-end.

**A. IMMEDIATE:**
- A1. **DONE + BUILD-CHECK GREEN (2026-06-06):** C6 tests + jj_corr + jj_corr_mrhs + hmc all compile
  after the deletion; C6a-d still bit-identical. REVERTED the now-dead old block code
  (superseded by BlockedMat): removed `MatPoly::solve_multishift_block` + `solve_block_cg`, the
  `DeviceMemorySetN` block buffers + `allocate_block`/`deallocate_block`, the `d_MemorySetParent` global
  (matpoly_claude.h); removed `OverlapWMass::{mult,adj,DDH}_deviceAsyncLaunch_ms_block`
  (overlap_wmass_claude.h). KEPT: the free block kernels include in matpoly (BlockedMat uses them);
  single `_ms`/`solve`/`solveAsync`. (The `on_gpu_pre` prealloc-gate leftover left in place -- unrelated,
  separate cleanup.) Build-check = tmp_claude.sh builds C6 tests + jj_corr + jj_corr_mrhs + hmc + re-runs
  C6 (bit-identical unchanged, BlockedMat untouched).

**B. DEPLOY jj_corr_mrhs to production: [CHOSEN PRIORITY after A1, 2026-06-06]**
  Concrete order:
  - B1a. THERMALIZED h5diff: pick one ensemble config (via `--ens-dir <ENS> --ninter <k>`), run
    `jj_corr_claude.o` and `jj_corr_mrhs_claude.o` on the SAME config, h5diff (bit-identical expected;
    free-field already covers the algebra, thermalized confirms on a real gauge background). OPEN: which
    ensemble/config to use (the n=21 sea dirs vs a `_pole11`; jj builds n=21 operators).
  - B1b. Repoint the run script (`run_jj_analysis_claude.sh` / `run_connected`) to `jj_corr_mrhs_claude.o`
    (or make it the default `--all` binary), keep `jj_corr_claude.o` as the reference fallback.
  - B1c. Run the ensembles with the faster binary.
  - B2 (below) folds in for L=2.
- B1. Validate on a THERMALIZED config (not just free field): run one ensemble config through jj_corr
  vs jj_corr_mrhs, h5diff (bit-identical expected). Then point the run script(s) (run_jj_analysis /
  run_connected) at jj_corr_mrhs (or make it the default), run the ensembles with the faster binary.
- B2. L=2 MEMORY: **DEFERRED (not for a while, per user 2026-06-06).** When L=2 comes: 6 BlockedMat
  instances x exactly-sized scratch grow (n_links=120, larger N) -- check it fits (TITAN V 12GB /
  A100). Lever: create only the instances used per mass-case (massless -> Dm only; parity -> Dm+Dtil;
  flavor -> Dm+D => 2-4 not 6), or share scratch. Not on the current path (L=1 fits comfortably).

**C. OPTIMIZATION -- jj is K-APPLY-BOUND (KEY FINDING from the 1056s free-field ref run):**
  Per-hit breakdown: source solves (`op_Dmsq`, batched by BlockedMat) ~55s (~5%); the SINK
  `K(n,t)phi'` applies (~10.7k applies = `Nt` x sites/links x 4 channels, NOT batched) ~950s (~90%).
  => the inversions we batched are NOT the jj bottleneck; the sink `K` applies are. So the
  source-solve mrhs is only ~3% off the whole jj, and the fused-dot (C1) is <1% -- the real jj lever
  is batching the SINK applies (C6f, promoted to the top here).
- **C6f = THE jj optimization lever [PROMOTED from D2, 2026-06-06].** Batch the sink `K` applies: at
  fixed `t` the `n_sites`/`n_links` `K(.,t)phi'` applies are `apply_k`s with the SAME Zolotarev seed,
  DIFFERENT RHS -- exactly the block-multishift structure BlockedMat already has. Block `apply_k`
  (`nstack = n_sites/n_links`) over the sink pass (`conserved_current_claude.h::apply_k`). ~2x on the
  sink => ~45% off the whole jj (vs ~3% from the source-solve mrhs). The block inner-solve reuses
  BlockedMat's `solve_multishift_block`; the m-dependent Term B / COO matvecs get block kernels. HIGH
  VALUE -- the headline jj optimization.
- C1. FUSED BLOCK-DOT: replace the `2*NSTACK` host-blocking `cublasZdotc`/iter (in `solve_sq` +
  `solve_multishift_block`) with ONE per-column-dot kernel (device output, 1 memcpy) -> removes the
  host serialization. **LOW value for jj** (source solves ~5% => <1% whole-jj), but RELEVANT for HMC
  (D1) where the `op_Dmsq` solve IS the dominant cost -- so do C1 WITH/FOR the HMC batching, not jj.
- C2. Launch-config tuning (threads/block 128/256/512, retune A100): minor.

**D. EXTENSIONS:**
- D1. HMC pseudofermion batching (WANTED): reuse BlockedMat to batch a few pseudofermion SETS as
  `nstack` RHS in the HMC action/force solve (each an independent RHS over op_Dmsq). Wiring in
  hmc_claude.cu / hmc_w_mass. Operator-generic, so it's a wiring task.
- D2. C6f block `apply_k` -- **PROMOTED to the C section (the jj K-apply lever, ~45% off whole jj).**
  No longer "deferred": the sink `K` applies ARE the jj bottleneck (~90% of wall time), so batching
  them is the high-value jj optimization. See C above.
- D3. The standalone jj programs (jj_conn_tpproj/spproj/axial/ylm, jj_disc) still use per-site single
  solves; could get the same BlockedMat rewire, but jj_corr (unified) is the production path -- leave
  unless wanted.

## C7 -- consolidate into `BlockedMat<Idx N, int NSTACK, class Op>` (REFACTOR, 2026-06-06, user)

All mrhs block functionality migrates into ONE class `MatPolyBlocked<N,NSTACK,Op>`
(`includes/matpoly_blocked_claude.h`). It HOLDS a `const Op& D` (the overlap operator), OWNS the
block scratch (ctor allocates, dtor frees -- RAII; sized from N, NSTACK, npole), and contains the
inner sign-solve + the overlap block applies + the outer block CG:
- `solve_multishift_block(dX, dB, tol)` -- seed `(1/lambda_max^2) M_DWH M_DW`, shifts `-k^2/cp[m]`
  computed from D; (was `MatPoly::solve_multishift_block`).
- `mult`/`adj`/`DDH(dres, dxi)` -- the overlap block applies (were `OverlapWMass::*_ms_block`), using
  `D.M_DW/M_DWH/A/C/mass/lambda_max/k/cp/size` + `this->solve_multishift_block`.
- `solve(dX, dB, tol)` -- outer per-column block CG over `this->DDH` (was `MatPoly::solve_block_cg`;
  per-column complex al + per-column freeze + shared worst-column stop).
Method names avoid the free-kernel names (kernels stay `mult_block`/`block_reduce_poles`/`axpy_*`).
**Reverts:** the block buffers come OUT of `DeviceMemorySetN`, `solve_multishift_block`/`solve_block_cg`
come OUT of `MatPoly`, the `*_ms_block` come OUT of `OverlapWMass`, and `d_MemorySetParent` is removed
(all back to pre-C6 state). Op interface required: `M_DW,M_DWH (CSR<N>), lambda_max, k, size, cp[], A[],
C, mass`. WHY a class (not local-in-function): the apply runs once/outer-iter (tens-hundreds of times
/solve) and the ctor allocates ~tens of MB scratch -> a per-call object = per-call device malloc (the
thing we eliminated). So created ONCE. jj: one instance per (operator, width) actually used --
`MatPolyBlocked<N,N_SITES,Fermion> blk_Dm_tp(Dm, npole)`, ... -- call `blk.solve(dX,dB,tol)`. Memory:
each instance owns exactly-sized scratch; at most ~4 instances per run (Dm + one of {Dtil,D}) x {tp,sp}.
**Sequencing: refactor first, validate once** (C6 tests + jj h5diff).
**DONE / VALIDATED (header `includes/matpoly_blocked_claude.h`; C6b/c/d repointed to the class):**
C6a `max|diff|=0`; C6b/c/d bit-identical (`max|block-single|=0`, mult/adj/DDH=0) at nstack=1/4/12;
C6d ~2.45x @ nstack=12. The class consolidation preserved all numerics. PENDING: (1) REVERT the now-
dead old block code (`MatPoly::solve_multishift_block`/`solve_block_cg`, `DeviceMemorySetN` block
buffers + `allocate_block`/`deallocate_block`, `OverlapWMass::*_ms_block`, the `d_MemorySetParent`
global); (2) rewire `jj_corr_mrhs_claude.cu` to `MatPolyBlocked` (drop the named ab_* lambdas +
`d_Bblk`/`d_Xblk`; one `MatPolyBlocked` per (operator,width) used; `blk.solve` per loop); (3) jj
h5diff free-field vs `jj_corr_claude.o`.

## C6 memory & template model (REVISED 2026-06-06, from discussion -- SUPERSEDES the runtime-nstack
## / internal-malloc shape of the first C6a/C6b draft) [SUPERSEDED by C7 consolidation above]

Two decisions reshape how `nstack` and the block scratch are handled (the first C6a/C6b draft used a
runtime `nstack` arg + internal `cudaMalloc`; both are replaced):

1. **`NSTACK` is COMPILE-TIME, native to the OPERATION (not the operator type).** It is a
   `Comp::NSTACK` constexpr defined per `.cu` (value chosen from L / hardware / jj-vs-hmc), threaded
   as a TEMPLATE parameter into the block kernels, `solve_multishift_block`, the block overlap
   methods, and their callers -- exactly as `N` is threaded. `OverlapWMass` is NOT class-templated on
   `NSTACK` (the only reason to was that the single-RHS `d_Zblock` is an Overlap member; decision 2
   moves that scratch out, so the operator type stays single). `npole = size-1 = int(n/2)` stays
   RUNTIME (set at operator construction; passed as a runtime arg where kernels need it). The C6a
   leaf kernels (`mult_block`, `multishift_x/p_update_block`) are templated `<Idx N, int NSTACK>`;
   `axpy_col`/`axpby_col` keep a runtime `ncol` (they serve both the N*NSTACK and N*NSTACK*npole
   lengths). [DONE: kernels templatized.]

2. **The wide block scratch is THREAD-OWNED in the per-stream `d_MemorySets` pool, preallocated once
   -- and REUSED by both the inner block solve AND the conserved-current `K` block apply.** Instead
   of `OverlapWMass` members or per-call `cudaMalloc`, `DeviceMemorySetN` (matpoly_claude.h:7) gains
   block buffers, allocated by a new `allocate_block(int nstack_cap, int npole)` (separate from the
   default `allocate()`, so `.cu` that never batch are untouched).
   **PARENT-STREAM SET (revised per user, 2026-06-06):** the block scratch lives in a dedicated global
   `DeviceMemorySetN d_MemorySetParent` (the default/foreground stream), NOT in the child
   `d_MemorySets[i]` (those stay the threaded per-stream scratch). `allocate_block` is called once on
   the parent; `solve_multishift_block` / block overlap methods / `solve_block_cg` all use
   `d_MemorySetParent` (no `istream` indexing); tests call `d_MemorySetParent.allocate_block(...)`.
   Buffers:
   ```
   d_p0_blk, d_q_blk, d_r_blk, d_seedtmp_blk : N*nstack_cap          (seed/outer)
   d_pm_blk, d_Zblk_blk, d_Yblk_blk          : N*nstack_cap*npole    (per-(rhs,pole))
   d_mp_blk, d_mm_blk                        : N*nstack_cap          (C6c DDH: (D+m)v, (D+m)^dag v)
   d_ocg_p_blk, d_ocg_q_blk, d_ocg_r_blk     : N*nstack_cap          (C6d outer-CG p/q/r)
   d_alm_blk, d_zeta_blk, d_betm_blk         : nstack_cap*npole      (double)
   d_sc_blk, d_sc2_blk                       : nstack_cap            (double, seed-level scalars)
   d_alc_blk                                 : nstack_cap            (CuC, C6d per-column al)
   ```
   Sized from a `Comp::` constexpr PRODUCT of the other `Comp::` constants, e.g. (values chosen per
   `.cu` by Claude):
   ```
   constexpr int N_ZOLO = 21;            // Zolotarev order (matches the Fermion ctor arg)
   constexpr int NPOLE  = N_ZOLO/2;      // = size-1 = 10
   constexpr int NSTACK = N_SITES;       // jj sink-loop batch at this L (=12 at L=1, 42 at L=2)
   // block scratch length multiplier = NSTACK*NPOLE   (the "product Comp:: constexpr")
   ```
   A solve with `NSTACK <= nstack_cap` uses the FRONT N*NSTACK(*npole) portion of the buffers, so one
   `allocate_block(NSTACK_MAX, npole)` (NSTACK_MAX = the largest batch the binary uses) serves all
   smaller batches -- this is how the C6b test exercises NSTACK=1/4/12 against a single allocation.

3. **Launch config stays 256 threads/block for now** (correctness first); `NBlocks`/threads-per-block
   tuning for the wide kernels is delegated to Claude as a later knob (the global `NBlocks` macro is
   untouched; block grids are caller-computed `(work + 256 - 1)/256`).

Consequence: `solve_multishift_block<N,NSTACK>(d_X, d_B, M0, M1, coeff, sigma, npole, tol)` -- drops
the runtime `nstack` arg, uses `d_MemorySets[is]` scratch (no internal malloc), writes the
caller-owned `d_X`.

**REWORK DONE (pending user compile/run):** C6a kernels templatized `<Idx N, int NSTACK>`
(`mult_block`, `multishift_x/p_update_block`; `axpy_col`/`axpby_col` keep runtime `ncol`).
`DeviceMemorySetN` gained the block buffers + `allocate_block(nstack_cap,npole)`/`deallocate_block`
(matpoly_claude.h). `solve_multishift_block<N,NSTACK>` reworked: foreground slot `is=istream>=0?:0`,
`assert(block_set && cap_nstack>=NSTACK && cap_npole>=npole)`, all scratch from `d_MemorySets[is]`,
no internal malloc/free. C6a test calls templatized kernels (NSTACK=3). C6b test: `allocate_block
(Comp::N_SITES,npole)` once, templated lambda `run.operator()<NSTACK>()` for NSTACK=1/4/12 (front
portion of the cap-sized buffers), `deallocate_block` at end. Re-validate: both must still be
bit-identical (`max|diff|=0`).

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
- **Reusable for HMC later (WANTED -- flagged 2026-06-06):** the same `_block` solve should be used
  in the HMC when there are a few pseudofermion SETS -- each pseudofermion's action/force solve is an
  independent RHS over the SAME `op_Dmsq`, so batch them with `nstack = #pseudofermion-sets` via
  `solve_block_cg` (action-solve heatbath/CG) exactly as the jj sink loop does. The block machinery
  (parent-set scratch, `DDH_..._ms_block`, `solve_block_cg`) is operator-generic, so the HMC reuse is
  a wiring task once jj (C6e) is validated. We do NOT touch HMC now, but this is a target application,
  not just a "possible" one.
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
  **IMPLEMENTED (pending user compile/run).** New header `includes/multishift_block_kernels_claude.h`:
  `mult_block<T,N>` (CSR matvec, grid N*nstack, column-major [c*N+i], shared CSR), `axpy_col<N>`
  (per-column real-scalar AXPY a[c]*x+y, covers q+=sig0_c p0 / r-=al_c q / p0=r+bet_c p0),
  `multishift_x_update_block<N>` / `multishift_p_update_block<N>` (extend the single-RHS N*npole
  updates to N*nstack*npole; scalar indexed by col=gid/N=c*npole+m; seed residual r broadcast over
  poles via c=col/npole). nstack RUNTIME arg; grids caller-sized (no NBlocks macro). Block dot
  DEFERRED to C6b (host strided `cublasZdotc` per column). Driver
  `test_multishift_block_kernels_claude.cu` (self-contained random CSR + random vectors, N=64
  nstack=3 npole=5): each block kernel must reproduce its single-RHS sibling on the contiguous
  per-RHS sub-block BIT-IDENTICALLY -> PASS = max|diff| < 1e-12. Originals/`_ms` untouched.
- **C6b** `MatPoly::solve_multishift_block<N>(d_X, d_B, sigma, npole, nstack, tol)` + block `on_gpu`
  for the seed matvec. Validate the block `Z_{c,m}` vs `nstack` separate `solve_multishift` calls
  (identical to tol). STANDALONE-VALIDATABLE -- this is the foundational, self-contained piece.
  **IMPLEMENTED (pending user compile/run).** `matpoly_claude.h` adds `solve_multishift_block<N>`
  after `solve_multishift` (originals untouched): nstack INDEPENDENT multishift solves batched into
  the C6a wide kernels. SEED = specialized two-CSR `coeff*M1*M0` (takes `const CSR<N>& M0,M1` +
  coeff DIRECTLY; mult_block twice + `axpby_col` for `coeff*(M1 M0 p0)+sig0 p0`; NO sparse_matrix.h
  edit). Per-column seed scalars (mu/gam/al/bet) + per-(c,m) zeta host recurrence. freeze-converged
  per (c,m) AND per COLUMN (`colfrozen`): the column freeze is REQUIRED (not the deferred opt) --
  an early-converged column kept live drives p0_c->0 and hits gam_c->0/0 (the same breakdown the
  per-pole freeze prevents); a frozen column does no updates so X_{c,m} is byte-frozen at its
  single-RHS break iter => block == per-column `solve_multishift` BIT-FOR-BIT. Shared worst-column
  stop `max_c mu_c/||b_c||^2<tol^2`. New leaf kernel `axpby_col` (two per-column scalars) added to
  `multishift_block_kernels_claude.h`. Per-column dots = 2*nstack host-blocking `cublasZdotc`/iter
  (start-simple; fused block-dot = next opt). Driver `test_multishift_block_solve_claude.cu` (real
  overlap seed, free field, size=21): nstack=1 (reproduces solve_multishift) + nstack=4 + nstack=
  n_sites; PASS = max|block-single| < 1e-10 every column (expect ~machine, bit-identical).
  NOTE: per-pole freeze (`frozen[m]`) is the SAME freeze-converged already in `solve_multishift`
  (HMC-NaN fix); the block just adds the per-COLUMN extension (`colfrozen`). The old "freeze deferred
  to v2" line in the Decisions section is revised accordingly (it is now required, not deferred).
  **DONE / PASS (C6a+C6b, TITAN V, free field, size=11/npole=10):** C6a leaf kernels all
  `max|diff|=0`. C6b block-vs-single BIT-IDENTICAL (`max|block-single|=0`) at nstack=1/4/12.
  Timing (naive host-blocking dots): nstack=4 block 21.8ms vs serial 42.9ms = 1.97x; nstack=12 block
  51.8ms vs serial 130.9ms = **2.53x** (in the predicted 2-4x band; fused block-dot = next opt).
- **C6c** block `mult_ms`/`adj_ms`/`DDH_ms` in overlap (inherit block-ness from `solve_multishift_block`
  + block `A[m]` combine). Validate `D_m^{-1}`/`D_m^{-dag}` block vs single.
  **IMPLEMENTED (pending user compile/run).** `overlap_wmass_claude.h` (after `DDH_..._ms`, ~:482):
  `mult_/adj_/DDH_deviceAsyncLaunch_ms_block<int NSTACK>` (side-by-side; originals untouched; slot
  is=0). Each: inner sign solve -> `solve_multishift_block<N,NSTACK>` into `d_MemorySets[is]`
  `.d_Zblk_blk`/`.d_Yblk_blk`; residues `A[1..size-1]` -> `d_alm_blk` (free post-solve);
  `block_reduce_poles` fold; `(1/lambda_max) M_DW`/`M_DWH` via `mult_block`; `C`/`mass`/`(1/lmax)`
  folds via `axpy_uniform` (uniform complex scalar). adj uses `d_res` as the Ys0/RHS holder (`d_B`
  const in the solve). DDH = `mult_block` into `d_mp_blk`, `adj_block` into `d_mm_blk`, then the 3
  `axpy_uniform` folds (conj(m), 1+m, -(2Re m+|m|^2)). New kernels `axpy_uniform`/`block_reduce_poles`
  added to `multishift_block_kernels_claude.h`; `d_mp_blk`/`d_mm_blk` added to `DeviceMemorySetN`.
  Driver `test_overlap_block_claude.cu` (real op, free field, nstack=1/4/N_SITES): block vs
  per-column single `_ms`, expect bit-identical max|diff|<1e-10. In `tmp_claude.sh`.
  **DONE / PASS (TITAN V, free field, size=11):** mult/adj/DDH all `max|diff|=0` at nstack=1/4/12.
- **C6d** block outer CG: `MatPoly::solve` over `nstack` RHS (per-column scalars, iterate-to-slowest).
  Validate `op_Dmsq` block solve vs nstack single solves.
  **IMPLEMENTED (pending user compile/run).** `MatPoly::solve_block_cg<N,NSTACK,ApplyOp>(apply, d_X,
  d_B, tol)` (matpoly_claude.h, after the single `solve(CuC*,CuC*)`): plain block CG, operator via a
  callable `apply(d_res,d_v)` (the test passes `Dm.DDH_deviceAsyncLaunch_ms_block<NSTACK>` = op_Dmsq),
  per-column COMPLEX `al=mu/gam` (matches single `solve` -> bit-identical), per-column freeze, shared
  worst-column stop. New kernel `axpy_col_c` (per-column complex). Scratch from the PARENT set (see
  below): `d_ocg_p/q/r_blk` + `d_alc_blk` (added to `allocate_block`; NOT touched by the block
  operator's internals). Driver `test_block_outer_solve_claude.cu`: block vs nstack single
  `op.solve` (MatPoly{LinOpWrapper of DDH_ms}), nstack=1/4/N_SITES, PASS=max|diff|<1e-10 + reports
  block-vs-serial speedup. In `tmp_claude.sh`.
  **DONE / PASS (TITAN V, free field, size=11):** block op_Dmsq bit-identical (`max|diff|=0`) at
  nstack=1/4/12; speedup vs nstack-serial = 1.05x / 1.94x / **2.44x** (the sink-solve batch win).
  Parent-set swap left C6b/C6c still `max|diff|=0`.
- **C6e** batch the jj sink loop (`nstack = n_sites`); free-field re-run -> Vpp/Vmm unchanged + faster.
  **IN PROGRESS (`jj_corr_mrhs_claude.cu`, copy-edit of `jj_corr_claude.cu`, per user; target = jj_corr
  ONLY).** Scope = BLOCK ALL THREE Dsq operators (`op_Dsq`/D, `op_Dmsq`/Dm, `op_tilDmsq`/Dtil) wherever
  a per-site (tp, nstack=N_SITES) or per-link (sp, nstack=N_LINKS) SOURCE-solve loop occurs (the ~8
  loops at jj_corr_claude.cu :426/431/476/481/514/547/622-624/645-647; ylm N_ELL=3 + single phi'/tilphi
  solves left single). Pattern per loop: build host RHS block `hblk` (loop col n: K apply + `op_X.from_cpu`
  into `hblk+n*N`), H2D -> `d_Bblk`, `op_Xsq.solve_block_cg<N,NSTACK>(ab, d_Xblk, d_Bblk, TOL_OUTER)`
  with `ab` = the matching `{D,Dm,Dtil}.DDH_deviceAsyncLaunch_ms_block<NSTACK>`, D2H -> scatter to
  `psi_*`. Setup added: `NSTACK_TP=Comp::N_SITES`/`NSTACK_SP=Comp::N_LINKS`, `d_MemorySetParent.
  allocate_block(max,npole)`, `hblk`/`d_Bblk`/`d_Xblk`; freed at end. **DONE: tp(++) via op_Dmsq (the
  reviewable template).** REMAINING (identical pattern): tp(--) tilDmsq, sp(++) Dmsq, sp(--) tilDmsq,
  axial tp/sp (flavor->Dsq / parity->tilDmsq / else->Dmsq). Validate: free-field `corr.h5` tp Vpp/Vmm
  bit-identical to `jj_corr_claude.o`, + per-loop timing. (Bit-identical since block==per-column single.)
- **C6f** block conserved-current apply: `ConservedCurrent::apply_k`/`apply_k_dag` acting on an
  `nstack`-wide vector in one shot (the block analog of the C5 `apply_k_ms`/`apply_k_dag_ms`). Mirrors
  C6c: the same-RHS inner pole loops (Step 1 / Term A^dag) route through `solve_multishift_block`
  into `d_MemorySets[is].d_Zblk_blk`; the m-DEPENDENT Term B / Term B^dag and the COO matvecs become
  block kernels (`mult_block` + the `_block` axpys). Reuses the SAME `d_MemorySets` block scratch as
  the overlap block ops (the "retain extended vectors for K" decision). Side-by-side with `apply_k_ms`
  (originals untouched). Validate block `K`/`K^dag` vs per-column single `apply_k_ms`/`apply_k_dag_ms`
  (free field, bit-identical). NEEDED only if the estimator batches the `K` applies themselves;
  if `K` stays single-vector + batched sink solve, C6f is optional. [Added per user request 2026-06-06.]

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
  machine-precision floor issue). **UPDATE (C6b):** "all columns stay live; freeze deferred" was
  REVISED -- per-column freeze (`colfrozen`) is now implemented. It is the column-axis extension of
  the per-pole freeze-converged ALREADY in `solve_multishift` (the HMC-NaN fix), and is REQUIRED
  here: with the shared worst-column stop, an early-converged column kept live drives its search
  direction `p0_c -> 0` and hits `gam_c = <p0_c|A p0_c> -> 0/0` (same breakdown the per-pole freeze
  prevents -> would trip the `assert(false)`). Freezing it (no updates after its seed residual hits
  its floor) avoids that AND makes the block result bit-identical to the per-column single solve.

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
