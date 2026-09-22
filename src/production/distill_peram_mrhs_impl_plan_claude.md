# distill_peram MRHS-blocked + multi-source (impl plan)

Goal: accelerate the distillation perambulator generation and double the per-config statistics, WITHOUT
touching the currently-running `distill_peram_L1_claude.o` source. Two changes:
1. **mrhs block the $N_v$ mode solves** (the biggest win): the per-mode single-RHS `op_Dsq.solve` loop
   re-applies the (expensive Zolotarev sign-function) overlap operator $N_v$ times; batch the $N_v$ RHS into
   ONE `BlockedMat::solve_sq_from_cpu` block-CG. Expected ~2-3x (matches jj C6f 2.99x, HMC force ~2.5x).
2. **second source window at $t_0=N_t/2=64$**: only ~32 of 128 timeslices are used; a source's footprint is
   ~$\pm30$ (propagates BOTH directions), so exactly TWO non-interfering windows fit: $[0,32)$ and $[64,96)$.
   Doubles the independent statistics (~$\sqrt2$ error) at 2x solve cost, amortizing the expensive gauge gen.

## Files (VARIANT -- do NOT overwrite the running source)
- Copy `distill_peram_claude.cu` -> `distill_peram_mrhs_claude.cu`; edit the copy. Keep the single-RHS solve
  loop COMMENTED in place beneath the mrhs version (A/B + rollback). Reuse `../both_3d/includes/blocked_mat_claude.h`.
- Update the Makefile / provide a build command for the new target; handoff `tmp_claude.sh` + `*_claude.log`.
- Downstream: `distill_contract_claude.py` must read the multi-source h5 (chunk 3).

## Key API (from `blocked_mat_claude.h` + `jj_corr_mrhs_claude.cu`)
- `using Fermion = OverlapWMass<WilsonDirac>;` (as in the driver).
- `BlockedMat<N, NSTACK, Fermion> blk_Dsq(D);` created ONCE (ctor alloc heavy), reused every timeslice.
  `NSTACK` is COMPILE-TIME = $N_v = 2 N_\text{sites}$ (complete basis) -> a `constexpr` from `N_REFINE`
  (24 @L1, 84 @L2). Cap it like the jj code (`Comp::` constexpr).
- Block solve of $(D_{ov}^\dagger D_{ov})^{-1}$ over the NSTACK-wide host block via `solve_sq_from_cpu`
  (host-block staging `std::vector<Complex> hblk(N*NSTACK)`, device I/O hidden). Bit-identical per column to
  the single-RHS `op_Dsq.solve` (validated C6a-d in the header).

## Chunks
1. **mrhs block solve** (`distill_peram_claude.cu:461-475`, the `for(a) for(l<nv) op_DH; op_Dsq.solve; ...`):
   - Build the block RHS: for `l in 0..nv-1` embed $w_l(t)$ and apply `op_DH` (cheap, keep per-source) ->
     stage $D_{ov}^\dagger w_l$ into `hblk` column `l`.
   - ONE `blk_Dsq.solve_sq_from_cpu(hblk, tol)` -> $\psi_l = D_{ov}^{-1} w_l$ for all $l$.
   - Per `l`: read column `l` back to a Fermion, `op_oneMinusDdag` -> `dpsi`, extract sink slices $t'$,
     fill `Tau[ap][a].col(l)`, `Taup[ap][a].col(l)` (unchanged contraction). Comment the old loop below.
   - VALIDATION: the existing T2a/T2b/T2c reconstruction checks must still print ~solver tol; ALSO assert
     max|mrhs - single-RHS| < few*tol on the first config, and PRINT per-config wall-time + speedup vs a
     single-RHS timing run (feedback_benchmark_tests).
2. **multi-source**: `--tsrc-list 0,64` (or `--nsrc 2` -> {0, Nt/2}). Wrap the basis+peram build in a loop over
   `tsrc0` in the list (the $V(t)$ basis is per-timeslice, already all-$t$; only the solves/contraction repeat
   per window). Store per window.
3. **h5 format + analysis** (OPEN, see below): store `/peram/tau` as shape `(nsrc, twin, twin, Nv, Nv)` (+
   `tau_gw`), meta `tsrc_list`. `distill_contract_claude.load_peram` returns a list/extra axis; the correlator
   code treats the `nsrc` windows as independent source sets (translation-average within each, then combine
   -> effectively 2x configs). Keep backward-compat read for old single-window files (nsrc=1).

## Open questions for NM (resolve before chunk 3)
- h5 layout: extra leading `nsrc` axis on `tau`/`tau_gw` (compact) vs separate `/src0/`,`/src1/` groups
  (clearer, easier partial-resume). Recommend the `nsrc` axis.
- Combine the two windows as 2 independent "configs" (simplest, doubles ncfg) or average per config first?
- NSTACK = full $N_v$ (24@L1) in one block, or sub-block (e.g. 12) if register/scratch pressure hurts? Start full.

## Build / run
- `nvcc -arch=sm_70 -O2 -lcusolver -std=c++17 ... distill_peram_mrhs_claude.cu -o distill_peram_mrhs_L1_claude.o`
  (match the Makefile flags; include path for `blocked_mat_claude.h`). Handoff runs it stride-10 with
  `--tsrc-list 0,64`, tee to `distill_peram_mrhs_claude.log`; time vs the current .o on a few configs first.
- Do NOT touch the running `distill_peram_L1_claude.o` / its .cu. NO rm/kill in scripts.

Refs: multishift Jegerlehner hep-lat/9612014; distillation Peardon 0905.2160.
