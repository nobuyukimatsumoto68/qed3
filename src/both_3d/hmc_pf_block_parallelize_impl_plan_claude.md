# HMC pseudo-fermion action-inversion parallelization (MRHS block) -- impl plan

**Reference (algorithm):** B. Jegerlehner, "Krylov space solvers for shifted linear systems",
hep-lat/9612014 (inner multishift) + standard block CG (outer). The block machinery
(`BlockedMat`/`BlockMemPool`) is already validated bit-identical per column (C6a-d, `blocked_mat_claude.h`).

## Context / motivation

Nf=6 HMC production is the bottleneck: each trajectory runs Nf/2 = 3 independent pseudo-fermion
inversions per fermion kick, currently done SERIALLY in the integrator
(`integrator.h::MinimumNorm2`, loop `for(pf : pfs)`). Each `pf->update_eta()` is an outer CG over
$D^\dagger D$ (multishift inner solve per iteration) -- the dominant fermion cost.

The three inversions share the SAME operator $D^\dagger D(U)$ and differ ONLY by the Gaussian RHS
$\phi_f$. That is a textbook multi-RHS (mrhs) problem. Instead of running them one-at-a-time (each
kernel narrow, launch-overhead / low-occupancy bound -- especially at small lattices), we FUSE the
Nf/2 columns into one $N\cdot(N_f/2)$ block and run a single block CG: every sparse mat-vec then
processes all flavors at once. This raises per-kernel occupancy and cuts launch count -- the same
lever that gave jj ~3x (C6f blocking).

**Scope (this pass): ACTION inversions only.** Heat-bath `gen` (block `adj`), `update_eta` (block
`solve_sq`), and `S()` go through the block. The FORCE (`precalc_grad` + per-link `grad`) stays
SERIAL per flavor -- it stores state in the shared `D` and would need new NSTACK x npole block-grad
code. The block manager is structured so the force can be added later as a drop-in (a later chunk).

**Interaction with MPS / run strategy (global design note).** The split is by **Nf**, not L:
- **Nf=2** (NSTACK=1, a single pseudo-fermion -> blocking is a no-op): stays on **MPS** process-packing
  to fill the GPU. NOT parallelized.
- **Nf=4 (NSTACK=2) and Nf=6 (NSTACK=3)**: use the **block**. Multiple RHS fused into one wider solve
  raises single-process occupancy so a single process can better fill the GPU (reducing / removing the
  need for MPS at these Nf).
All refinement levels are in scope, **including L=1**. The deciding measurement is a block-alone vs
current MPS-packing benchmark at representative L/Nf (see Verification).

## Approach chosen (decisions)

- Route **B, MRHS block** (not stream/thread-parallel): reuses validated `BlockedMat`, sidesteps the
  thread-safety surgery that privatizing `OverlapWMass`'s shared streams/handles/`d_MemorySets`/
  `is_precalc` would require.
- **Everything L/Nf-specific becomes a compile-time macro** (like the existing per-L split), and all
  seven target drivers collapse into ONE unified source. Nf joins L as a compile-time knob (chosen:
  "option 2" -- Nf baked in, separate binaries, differentiated by suffix). The unmodified serial
  originals continue to serve **Nf=2** (+ MPS).

## Operator equivalence (why this is correct)

- `PseudoFermion::update_eta()` = `Op_DHD.solve<N>(d_eta, d_phi, TOL_OUTER)`, where `Op_DHD` routes
  through `LinOpDHDWrapper::act` -> `D.DHD_deviceAsyncLaunch_ms` (`sparse_matrix_claude.h:80`).
  `BlockedMat::solve_sq` (`blocked_mat_claude.h:402`) is an outer CG over `DDH` = block
  `D^\dagger D` using the SAME `DHD_deviceAsyncLaunch_ms` seed matvec. At `NSTACK=1` it is documented
  bit-for-bit identical to `MatPoly::solve`; at `NSTACK>1` it is the same algorithm with a shared
  worst-column stop -> matches serial to tolerance.
- `PseudoFermion::gen` = `D.adj_deviceAsyncLaunch_ms(d_phi, d_xi)` then `update_eta`. Block: fill the
  Nf/2 Gaussian columns IN THE SAME ORDER as the serial `for(pf) pf->gen` (so the RNG stream matches),
  block `adj` -> phi, then block `solve_sq` -> eta.
- `PseudoFermion::S()` = `real(Op_DHD.dot(d_phi,d_eta))`, and `MatPoly::dot` is a plain
  `cublasZdotc` (`matpoly_claude.h:197`) = $\phi^\dagger\eta$. Block: sum per-column
  `cublasZdotc(N, phi_c, eta_c)` real parts.
- `M_mass`/`mass_coeff` are ALWAYS built in `OverlapWMass` (`overlap_wmass_claude.h:253`, massless =>
  `mass_coeff=0`), so the block `mult/adj/DDH` cover BOTH massless and Family-B massive drivers with
  no change.

## New headers (shared, in includes/)

### `includes/pseudofermion_Nfblocked_claude.h` -- `PseudoFermionBlock<Fermion, NSTACK>`
Owns the Nf/2 = NSTACK columns as contiguous device blocks + its own block pool/operator.
- members: `Fermion& D; BlockMemPool<N,NSTACK> pool(D.size-1,true); BlockedMat<N,NSTACK,Fermion> blk(D,pool);`
  `CuC* d_phi_blk; CuC* d_eta_blk;` (each `N*NSTACK`); a `cublasHandle_t` for `S()` dots (or reuse
  `blk.handle`).
- `gen(rng)`: for `c=0..NSTACK-1` draw `rng.fill_gaussian` into host xi, copy to column `c` of a
  device xi block (SAME order as serial); `blk.adj(d_phi_blk, d_xi_blk)`; then `update_eta()`.
- `update_eta()`: `blk.solve_sq(d_eta_blk, d_phi_blk, Comp::TOL_OUTER)`.
- `S()`: `sum_c real( cublasZdotc(N, phi_c, eta_c) )`.
- `eta_col(int f)` / `phi_col(int f)`: `return d_*_blk + (size_t)f*N;` (force accessor -- serial force
  reads column f). This accessor is the seam where the future block-force plugs in.
- ctor allocs `d_phi_blk`/`d_eta_blk`; also exposes `size_t pool_bytes()` for the startup memory check.

### `includes/integrator_Nfblocked_claude.h` -- `MinimumNorm2Block`
Structural copy of `MinimumNorm2` (`integrator.h`) but the fermion kicks take
`PseudoFermionBlock<Fermion,NSTACK>& bpf` instead of `std::vector<PseudoFermion>& pfs`. Each kick:
```
bpf.update_eta();                                  // ONE block solve -> all etas
for(int f=0; f<NSTACK; f++){
  fermion->precalc_grad_deviceAsyncLaunch(U, bpf.eta_col(f)); // serial force (unchanged)
  dSf.compute(U, bpf.eta_col(f), D);                          // Force::compute (gauge_ext.h:312)
  pi += <coeff>*tau * dSf;
}
```
Keep the same lambda / nsteps_inner gauge sweeps and the `-lambda / -(1-2lambda) / -2lambda` fermion
kick weights verbatim. (The 0th kick and the two per-outer-step kicks each replace the serial
`for(pf)` loop with `bpf.update_eta()` + the serial force loop above.)

### `includes/hmc_Nfblocked_claude.h` -- `HMC2Block`
Copy of `HMC2` (`hmc.h`) holding `PseudoFermionBlock<Fermion,NSTACK>& bpf` in place of `pfs`:
- `H()` = `0.5*pi.squared_norm() + Sg(U) + bpf.S()`.
- `run()`: `pi.gaussian(rng); Gauge U0(U); bpf.gen(rng); h0=H(); integrate(); h1=H();` accept/reject
  identical to `HMC2` (on reject `U=U0; fermion->update(U)`).

## Unified driver source -- `hmc_Nfblocked_claude.cu`

ONE source replaces all seven originals (they differed only in L, geometry path, mass/window/nsteps --
all of which are now macros or runtime). Drop L8. Compile-time macros:
- `-DLREFINE=<1|2|4>`  -> `Comp::N_REFINE` (sizes `Comp::N`).
- `-DNFPF=<4|6>`       -> Nf baked in; `NSTACK_PF = NFPF/2`. Still read `argv` Nf and
  `assert(Nf==NFPF)` so existing run scripts pass unchanged. (Nf=2 stays on the serial originals+MPS.)
- `-DGEOM_DIR="..."`   -> geometry directory string (local `../../geometry/` vs fermilab absolute);
  passed into the `geodesic.h` include path / `dir` var. Build-script set. NOT a toggle -- a parameter.
- `Nt` stays a hardcoded `constexpr Nt=128` (no macro).

Executables suffixed by the build script: `hmc_block_L2_Nf6.out`, `hmc_block_L1_Nf4.out`, ...

Toggle-free unification of the historical local-vs-fermilab differences:
- **argv**: always accept `[gsq Nf nu0 mass_re mass_im max_sec]` (extras default 0).
- **run loop**: `while(k<kmax && (max_sec==0 || elapsed<max_sec))` -- L1's kmax-scan and fermilab's
  wall-time budget are just the two limits, both active. `kmax` a large `constexpr` default.
- **output dir**: must byte-match each original for AUTO-RESUME. The L1 (local) and fermilab strings
  differ ONLY by the `"mRe"+mass.real()+"mIm"+mass.imag()` segment inserted between `nu0` and `nt`.
  One compile-time flag reproduces both exactly: `#ifndef DIR_NO_MASS` include the segment (fermilab,
  default); `-DDIR_NO_MASS` omits it (local L1) + `assert(mass==0)`. Local L1 build gets
  `-DGEOM_DIR="../../geometry/" -DDIR_NO_MASS`; fermilab gets just its `-DGEOM_DIR`.

Compile-time tuning (L-keyed; `nsteps` also on runtime Nf), per user 2026-07-09:
| L | tmax | nsteps (all Nf) | nsteps_inner | Zolotarev |
|---|------|-----------------|--------------|-----------|
| 1 | 1.9  | 10              | 100          | 21, 0.001 |
| 2 | 1.0  | 8               | 100          | 21, 0.001 |
| 4 | 1.0  | 8               | 100          | 21, 0.001 |
(This UPGRADES the stale massive L1/L2/heavy: 17->21 poles, default->0.001 window, L2 nsteps 6->8.)

Body edits vs an original driver (all mechanical):
1. `#include "blocked_mat_claude.h"` right after `#include "overlap_wmass_claude.h"`
   (block kernels already arrive via `matpoly_claude.h`).
2. swap `#include "integrator.h"` -> `#include "integrator_Nfblocked_claude.h"`;
   after `#include "hmc.h"` add `#include "pseudofermion_Nfblocked_claude.h"` +
   `#include "hmc_Nfblocked_claude.h"`.
3. `constexpr int NSTACK_PF = NFPF/2; assert(Nf==NFPF);`
4. replace the `pfs` vector with `PseudoFermionBlock<Fermion,NSTACK_PF> bpf(D);`
5. `MinimumNorm2 integrator(...)` -> `MinimumNorm2Block integrator(...)`;
   `HMC2 hmc(...,pfs,...)` -> `HMC2Block hmc(...,bpf,...)`.
6. startup print: block-pool bytes vs `cudaMemGetInfo` free/total, `assert(pool_bytes < free)`
   (answers the device-limit question). Pool bytes ~ `(13*N*NSTACK + 3*N*NSTACK*npole)*16`;
   worst target L4/Nf6/npole20 ~ 145 MB; L1 ~ a few MB. Comfortable on 12 GB.

The seven serial originals are LEFT UNTOUCHED (Nf=2 + rollback/A-B reference preserved).

## Ordered chunks

- **Chunk 1 -- 3 headers + unified source.** Write `includes/pseudofermion_Nfblocked_claude.h`,
  `includes/integrator_Nfblocked_claude.h`, `includes/hmc_Nfblocked_claude.h`, and
  `hmc_Nfblocked_claude.cu`. Compile-verify at e.g. `-DLREFINE=2 -DNFPF=6`. (Build is the
  USER's to run -- provide the exact `nvcc` line / a `tmp_claude.sh`; I read the log.)
- **Chunk 2 -- validation handoff (DONE, awaiting run).** Added a `-DSERIAL_REF` flag to the unified
  source that selects the ORIGINAL serial path (`pfs`+`HMC2`+`MinimumNorm2`) with byte-same params, so
  block-vs-serial compares cleanly FROM ONE SOURCE. `run_Nfblock_validate_claude.sh` (USER runs, GPU1)
  builds block + `-DSERIAL_REF` at L1/Nf6 (fast sanity) and **L4/Nf6 (the real target)**, massless,
  local geometry, runs each from a sibling scratch dir under `src/` (so `../../geometry/` resolves and
  production dirs in `src/both_3d` are untouched), and compares first-traj dH (expect ~1e-6) + per-traj
  wall-time (speedup). Tees to `run_Nfblock_validate_claude.log`; I read it back. SINGLE-USE (re-run =
  delete the two scratch dirs first).
- **Chunk 3 -- build targets / rollout** once Chunk 2 validates. Makefile or build-script targets that
  emit the suffixed binaries (`hmc_block_L{1,2,4}_Nf{4,6}.out`) with the right
  `-DLREFINE/-DNFPF/-DGEOM_DIR`, and update the run scripts to call them for Nf=4/6 (Nf=2 stays
  serial+MPS).

## Verification

- Bit-level: at `NSTACK=1` the block solve is documented bit-identical -> a Nf=2 sanity build
  (`-DNSTACK_PF=1`) must reproduce the serial trajectory exactly (RNG order preserved).
- Physics: Nf=6, identical start config+RNG -> `dH` matches serial to ~1e-6 for the first trajectory;
  acceptance + plaquette distribution match over a short run.
- Performance: report per-trajectory wall-time serial vs parallelized at L2 and L4, Nf=6, single GPU,
  MPS OFF (this is the block-alone-vs-MPS deciding number).
- Memory: startup print confirms pool bytes << free device memory at every target L.

## Result (2026-07-09)

- **Correctness PASS (L1/Nf6):** serial-vs-block first-traj `dH` identical to **16 digits**
  (`-1.269272270219517e-01`), then ~1e-11 chaotic divergence. Block reproduces serial physics.
- **Timing (L1/Nf6, GPU1):** serial ~193.8 s/traj -> block ~174.4 s/traj = **~1.1x (~10%)** only.
  The block packs the ~21 action solves/traj but NOT the ~63 force evals (Nf/2 x 21 kicks) or the 2000
  gauge sub-steps; the serial **force** (`precalc_grad` multishift + per-link `grad`) is the heavier
  fermion cost -> action-only blocking caps ~10%. L4 timing not yet captured (val L4 timed out; L1
  traj ~190s > the 900s per-run cap). **The real lever is the deferred FORCE blocking** (seam:
  `PseudoFermionBlock::eta_col`).

## Resolved decisions (2026-07-09)

- Route B (MRHS block), ACTION inversions only; force deferred (plug-in seam = `eta_col`/`phi_col`).
- Nf compile-time (`-DNFPF`), separate suffixed binaries; Nf=2 stays serial+MPS.
- ONE unified source; `-DLREFINE`, `-DNFPF`, `-DGEOM_DIR`; drop L8; `Nt=128` constexpr.
- Zolotarev 21/0.001 all L; tuning table above; toggle-free dir (always `mRe/mIm`, L1 loses auto-resume).
