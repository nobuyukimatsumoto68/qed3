# L=3 ensembles — impl plan (_claude, 2026-07-27)

## Goal
Add the **L=3** (spatial refinement `N_REFINE=3`) massless ensembles to fill the continuum/refinement
interpolation between the existing L1, L2, L4. Odd refinement — confirmed OK in the run code (see
Prereqs). Algorithm refs: Hasenbusch mass preconditioning (M. Hasenbusch, hep-lat/0107019); multishift
CG (B. Jegerlehner, hep-lat/9612014); Zolotarev sign approx (frozen window per L).

## Locked decisions (NM 2026-07-27)
- **at = 0.2 only** (no half-a_t for L3 for now).
- **gsq = {1.5, 3.0, 4.5}** (the campaign scales gsq ~ with L: L1{0.5,1,1.5}, L2{1,2,3}, L4{2,4,6}).
- **Nf = {2, 4, 6}** -> **9 ensembles**, massless (mRe=0).
- **KMAX = 800** to start (grow later; "start and see").
- Split affine + qed3 allocations; **SCC picks up overflow** (no hard allocation constraint).

## Prereqs
1. **n3 geometry** (`*_n3.dat` in `/project/affine/nmatsum/qed3/geometry/data/`): only n1/n2/n4/n8 exist.
   **NM is generating the N_REFINE=3 geometry** (via the geometry `get_geometry`/Makefile tooling). BLOCKS
   everything downstream. (If an "L even" assert lives in the geometry generator, NM catches it there --
   the RUN code has NO L/N_REFINE parity assert; the only even-ness assert is `Nf%2==0`, driver:224, fine.)
2. Driver builds at `-DLREF=3` (already parameterized). `at` is `-DAT_VAL` (default 0.2). dir3 auto =
   `...at0.200000...nt128L3_hb...`.

## Files to modify (all in production/, git-synced)
1. **`includes/hasenbusch_ladder_claude.h`** — add L==3 branches to ALL of: `hasenbusch_ladder` (aux
   masses), `hasenbusch_steps`, `hasenbusch_tau`, `hasenbusch_mg`, `hasenbusch_nforce`, `hasenbusch_naction`.
2. **`includes/frozen_window_claude.h`** — add L==3 (lmin,lmax). PROVISIONAL wide window first; FROZEN from
   the tuning-run measured extrema.
3. **`launch_redo_claude.sh`** (BOTH ~/ and production copies, keep identical) — 9 BUILDS entries
   `hmc_fermilab_redo_massless_L3_Nf{2,4,6}_claude.o 3 <Nf> 800 20` (KRNG 20; add -D flags if tuning wants
   finer steps). Build via a TARGETED script (like tmp_build_half_at) to avoid rebuilding RUNNING at=0.2
   bins -- the shared-header edits (#1/#2) touch SRC/includes -> mtime invalidates ALL bins -> ETXTBSY on
   running p4-p8/p17-18. So: targeted `tmp_build_L3_claude.sh` (builds the 9 L3 bins + touches prod bins).
4. **`run_wrapper_redo_claude.sh`** — add L3 slots (new indices, e.g. p19/p20/p21/p22/p23 -- 9 streams ->
   5 pair-slots, one solo). Type `ml` (massless at=0.2), like-cost MPS pairing.
5. `redo_ensembles_claude.txt` (blackboard) — add 9 L3 rows.
6. New `params_L3_claude.md` — record the tuned ladder/window/steps/nsteps/gsq/KMAX.

## L3 tuning constants (STARTING values -> confirm at the tuning smoke)
- **ladder**: start 3-stage `{0.0, 0.4, 1.0}` (mirror L4; L3 is better-conditioned so <=0.4 OK) with steps
  `{3,3,3}`; OR 2-stage `{0.0, 1.0}` with steps `{4,4}`. NM leans "{3,3,3} or {4,4}". Decide by acceptance.
- **tau** 1.0 ; **mg** ~100 (L4-like; gauge force grows with L) -- try 100, drop to 50/20 if cheap enough.
- **n_action** ~28 (between L2=25 and L4=31) ; **n_force** 11 (all L use 11).
- **frozen window**: interpolate ~ (0.02, 6-7) as the PROVISIONAL wide bracket; FREEZE from measured
  lambda_min/lambda_max after a short thermalization (watch "eval outside window").

## Cost / feasibility
L3 N_SITES=92 (L2=42, L4=162) -> per-traj ~2.5-4x L2 ~= **350-700 s/traj** (A100, MPS-packed, EST).
9 streams x KMAX 800 ~= **8-16 kcore-hr**. Normal avail: affine ~1.3 + qed3 ~2.0 = **~3.3 kcore-hr** ->
normal covers only ~20-40%. Plan: **bulk on opp (free) + normal slice + SCC overflow**. A timing smoke
(once geometry lands) pins the real per-traj cost.

## STATUS 2026-07-27
- **C0 DONE**: n3 geometry generated (NM) -> 10 `*_n3.dat` in geometry/data/ (matches n2). No L-parity assert.
- **C1 DONE**: L3 branches added to `hasenbusch_ladder_claude.h` (ladder {0,0.4,1.0}:36, steps {3,3,3}:53,
  tau 1.0:65, mg 100:75, nforce 11:87, naction 31:98) + `frozen_window_claude.h` (PROVISIONAL (0.015,8.0):26).
  L1/L2/L4 untouched. (These are shared includes -> a launcher build would rebuild RUNNING at=0.2 bins ->
  ETXTBSY; the test/L3 build scripts touch prod bins to prevent that.)
- **C3 few-traj test STAGED**: `/lustre2/affine/redo/tmp_L3_fewtraj_claude.sh` builds the L3 driver
  (-DLREF=3 Nf2 KMAX20) + submits a qos=test 30min COLD-START run (gsq3.0), tee'ing startup
  delta/admiss/alat/window + per-traj dH/acc, and creating configs for the force-norm probe. NM runs it.
- **PENDING**: force-norm probe = `both_3d/test_hasenbusch_tune_claude.cu` (reads a config -> per-frame force
  norms L2/Linf + Osborn Cost to size nsteps). It lives in both_3d (relative geom + both_3d includes) and
  reads an existing config -> run it AFTER the few-traj test makes L3 configs; may need L3 added to
  both_3d/includes (frozen_window/ladder) or adapting to production includes. Set up next.

## Ordered chunks
- **C0 (NM)**: generate n3 geometry. BLOCKS all.
- **C1**: L3 branches in hasenbusch_ladder + frozen_window (provisional). Files: the 2 includes.
- **C2**: `tmp_build_L3_claude.sh` -> build 9 L3 bins (-DLREF=3, KMAX 800, KRNG 20) + touch prod bins.
  Files: new script. (BUILDS entries added to launchers for the record.)
- **C3 tuning smoke** (qos=test 30min, 1 stream, wide window): confirm geometry loads, admiss/alat(L3),
  delta, MEASURE lambda_min/lambda_max, acceptance -> FREEZE window + pick ladder ({3,3,3} vs {4,4}) +
  aux mass. Files: frozen_window (final values), maybe ladder.
- **C4**: wrapper slots + blackboard rows; launch on opp+normal split. Files: wrapper, blackboard.
- **C5**: validate startup across all 9 (delta/admiss/acc/|dH|, no eval-outside-window); record params_L3.

## Open (settle at C3 tuning smoke)
1. ladder: {3,3,3} (3-stage, aux mass ~0.4) vs {4,4} (2-stage). 
2. frozen window (lmin,lmax) at L3 -- from measurement.
3. mg (100 vs cheaper) and n_action (~28) -- from delta/cost.
