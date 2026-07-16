# Blocking (mrhs) for the conn/disc/condensate measurements -- benchmark + plan

## Goal

Amortize the dominant cost (the outer $D_m^\dagger D_m$ solve, `op_Dmsq`) of the $Y_{\ell m}$ measurement drivers by
packing the many same-operator right-hand-sides per config into a block CG (`BlockedMat`, the Nf-flavor blocking
machinery). conn: block the $(l,m)\times$leg sources (option 2, pure speedup). disc / condensate: block the dilution
patterns + optionally `nhit` (option 1, cheap statistics -- deferred; `nhit=1` stays default).

## Machinery (existing, proven)

`includes/blocked_mat_claude.h`: `BlockedMat<N,NSTACK,Op>` + `BlockMemPool<N,NSTACK>` (RAII scratch). Deployed in
`jj_corr_block_t_claude.cu`, `hmc_Nfblocked_claude.cu`. Column-blocked layout `d_B[c*N+i]` (RHS $c$, site $i$); the
sparse $D_W$ matvec (`mult_block`) loads each CSR row ONCE and applies across all `NSTACK` columns -> that is the
bandwidth amortization. Method `blk.solve_sq_from_cpu(host, host, tol)` = outer per-column block CG over $D_m^\dagger D_m$,
inner Zolotarev multishift block; **bit-identical per column** to serial (validated below, `maxdiff=0`). `NSTACK` is
compile-time; total RHS processed in $\lceil \text{RHS}/\text{NSTACK}\rceil$ tiles, so `nhit` costs tiles, not memory.

## Benchmark (2026-07-14; `test_block_solve_bench_claude.cu`, GPU0 TITAN V, memBW 653 GB/s)

Per L, real massless config (L1 k2000 / L2 k1000 / L4 k100), pole count 11 ($npole=5$), tol $10^{-5}$. Times ONE
block `solve_sq` of width NB vs NB serial `op_Dmsq.solve`s on the same Z2 RHS. `lambda_min` = 0.218 (L1) /
0.094 (L2) / 0.062 (L4).

| L | N | NB=1 | NB=8 | NB=16 | NB=32 | NB=64 | best |
|---|---:|---:|---:|---:|---:|---:|---|
| **L1** | 3,072 | 1.05 | 1.98 | 2.09 | 2.16 | 2.14 | **~2.1x @ NB=16-32** |
| **L2** | 10,752 | 0.97 | **1.66** | 1.48 | 1.43 | 1.40 | **1.66x @ NB=8** |
| **L4** | 41,472 | **1.36** | 0.81 | 0.79 | 0.77 | 0.76 | blocking HURTS (NB>=8) |

Correctness: `maxdiff = 0.000e+00` at every (L, NB) -- block CG is bit-identical to serial. Pool memory tiny
(L1 NB=32 = 47 MB, L2 NB=8 = 41 MB, L4 NB=64 = 1274 MB) -- never the constraint.

## Insight

mrhs blocking amortizes the $D_W$ read ONLY when a single solve under-fills the GPU. Strongly true at **L1**
(tiny lattice -> 2.1x), marginal and decaying past NB=8 at **L2** (1.66x), and at **L4** the single solve already
saturates bandwidth, so blocking just enlarges the working set and **slows down** (0.76-0.81x). So blocking payoff
is INVERSE to lattice size.

vs MPS (~1.5x aggregate measured earlier): block wins at L1 (2.1x), edges L2 (1.66x), loses at L4. Blocking
saturates the GPU at L1/L2, so MPS on top adds little -> it is **block XOR MPS**, not both.

**Bonus (single measurement, verify):** at L4 the block-CG path at **NB=1** ran 1.36x faster than `op_Dmsq.solve`
(21.3 s vs 28.9 s) -- the block-CG implementation itself beats the MatPoly serial solve for the big lattice,
independent of blocking. If robust, L4 gets a free ~1.36x by routing through `solve_sq` at width 1.

## Recommendation -- per-L NSTACK

- **L1: NSTACK=16**, single blocked process (drop MPS) -> ~2.1x.
- **L2: NSTACK=8**, single blocked process -> ~1.66x.
- **L4: NSTACK=1** (no blocking; optionally adopt the block-CG path for the ~1.36x), single serial process.

Compile-time via `-DNSTACK_CLI` per L (like `-DN_REFINE_CLI`). Run mode: L1/L2 one blocked process/GPU; L4 one
serial process/GPU (MPS ~flat there since saturated).

## Implementation plan (deferred pending deflation/LMA -- see below)

1. conn `jj_local_ylm_scalar_conn_stoch_claude.cu` (option 2): build `BlockedMat<N,NSTACK,Fermion>` per operator
   (Dm forward/backward legs); collect the $(l,m)\times$leg sources $\{W_{lm}\eta\}$ (~128/config) into host blocks,
   `solve_sq_from_cpu` in tiles of NSTACK, unpack to per-$(l,m)$ sinks. Keep the serial path commented for A/B.
2. disc `jj_local_ylm_disc_stoch_claude.cu` + condensate `jj_scalar_condensate_eo_stoch_claude.cu` (option 1):
   block the dilution-pattern RHS (free speedup at L1/L2); `nhit>1` stacks as extra tiles (cheap statistics),
   left OFF by default.
3. Run scripts: swap the 3-wide MPS pool for one blocked process/GPU at L1/L2; keep serial at L4.

Expected conn speedup vs current: L1 ~2.1x, L2 ~1.66x, L4 ~1x (or ~1.36x via block-path NB=1).

## Open / next

- Re-check L2 + L4 with the added **NB=2, 4** (the original sweep jumped 1->8, missing the L4 crossover between
  NB=1=1.36x and NB=8=0.81x, and the L2 low end). `tmp_claude.sh` re-runs L2+L4 only (L1 settled). Confirms the L4
  optimum (likely NB=1 or 2) and whether the block-path NB=1 1.36x holds before routing L4 through `solve_sq`.
- **Deflation + LMA (exploring now):** L4's cost is CG-iteration-bound (small $\lambda_\min=0.062$), which blocking
  cannot fix but low-mode deflation can (fewer iterations/solve, independent of L). LMA (low-mode averaging) also
  boosts statistics cheaply on the low-mode-dominated massless correlators (esp. disc/condensate). This may SUPERSEDE
  the blocking calculus for L4 and change the run strategy -- so hold the L4 wiring until the deflation/LMA
  assessment (Grid `/mnt/baracuda_14/grid_claude/Grid/examples` model codes) is in.
