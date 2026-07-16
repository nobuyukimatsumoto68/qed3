# Wire Nf-blocking (mrhs) into HasenbuschPF -- impl plan

## Goal

Pack the $N_f/2 = $ NSTACK Hasenbusch pseudofermion stacks into ONE $N\times$NSTACK block, driving the
per-frame ACTION solves (heatbath, eta, S) and the FORCE (Term A + Term B) through the validated
`BlockedMat` / `BlockedForce` (mrhs), instead of the current serial `for(pf)` loop in
`HasenbuschPF` + `FermionGroupLevel`. This is the deferred "C4" item, now on the two-operator split-pole path.

Only pays for $N_f\ge4$ (NSTACK $\ge2$). $N_f=2$ (NSTACK=1) stays serial -- nothing to pack.

## What already exists (reuse)

- `blocked_mat_claude.h` -- `BlockedMat<N,NSTACK,Op>`: block `adj`/`mult`/`solve_sq` (outer block CG over
  $D^\dagger D$) + `solve_multishift_block` (inner Zolotarev). **Mass-aware**: applies add $m_L=$
  `mass_coeff`$\cdot M_\text{mass}$ over NSTACK, so `set_frame_mass(i)` on the wrapped op composes (verified).
  Inner multishift + outer block CG (Jegerlehner hep-lat/9612014).
- `overlap_force_Nfblock_claude.h` -- `BlockedForce<N,NSTACK,Op>::compute()` = **Term A only** (block
  `precalc_grad`+`grad_l4`).
- `pseudofermion_Nfblocked_claude.h` -- `PseudoFermionBlock` = the single-frame (non-Hasenbusch) pattern to
  extend: `pool`+`blk`(+`bforce`), block buffers `d_phi_blk`/`d_eta_blk`, `gen`/`update_eta`/`S`/`eta_col`.
- `hmc_Nfblocked_claude.cu` -- the driver wiring for NSTACK (compile-time), block pool device-limit check.

## The gap to fill

`BlockedForce` has **no Term B (bilinear) block force**. The Hasenbusch ratio frames need
$-d[2\text{Re}\langle\phi_i|D_\text{ov}|\eta_i\rangle]/dU$ block-packed. Must add a block bilinear force
(mirror the serial `precalc_grad_bilinear_deviceAsyncLaunch_ms` + `compute`).

## Files to add / modify

- NEW `includes/pseudofermion_hasenbusch_block_claude.h` -- `HasenbuschPFBlock<Fermion,Force,NSTACK>`:
  per-frame block buffers `d_phi[i]`/`d_eta[i]`/`d_chi[i]` (each $N\cdot$NSTACK), TWO `BlockedMat` (action
  $D$ + force $D_f$) with their own pools (npole differ: $n_\text{act}$ vs $n_f$), block `gen` /
  `update_eta_frame` / `update_eta_force_frame` / `S` / `get_force_frames`. Mirrors `HasenbuschPF` frame
  logic (ratio + heavy, `set_frame_mass`, Term A + Term B) but every op is a block apply/solve.
- `includes/overlap_force_Nfblock_claude.h` -- add `BlockedForce::compute_bilinear(dSf, U, d_bra_blk, d_ket_blk)`
  (Term B, block).
- `includes/hmc_hasenbusch_ml_claude.h` -- a block `FermionGroupLevel` variant driving ONE `HasenbuschPFBlock`
  (no `for(pf)` loop); `HMCHasenbuschML::run` action-eta re-solve stays (block `update_eta`).
- `hmc_hasenbusch_claude.cu` -- for $N_f\ge4$: NSTACK = $N_f/2$ **compile-time** (see Q1), construct the
  blocked PF + block level; keep the serial path for $N_f=2$.

## Status (2026-07-14)

- **Chunk 1 DONE** -- block Term B in `includes/overlap_force_Nfblock_claude.h` (`precalc_bilinear` /
  `grad_bilinear_block` / `compute_bilinear`; bilinear grad reuses `grad_block` with bra in `eta_bra`).
- **Chunk 2+3 DONE** -- `includes/pseudofermion_hasenbusch_block_claude.h`: `HasenbuschPFBlock<Fermion,Force,NSTACK>`
  (two BlockedMats action D + force Df, per-frame N*NSTACK blocks, block gen/update_eta/update_eta_force/S +
  Term A/B force). Method-for-method mirror of the serial `HasenbuschPF`.
- **Chunk 4 DONE** -- NEW `hmc_hasenbusch_block_claude.cu` (`-DNF`, NSTACK=NF/2; one block PF in a size-1 pfs
  vector). **The ML integrator `hmc_hasenbusch_ml_claude.h` needed NO change** -- its `FermionGroupLevel` /
  `HMCHasenbuschML` templates only call methods the block PF also has. Serial `hmc_hasenbusch_claude.cu` kept.
- **VALIDATED (2026-07-14): DONE.** (a) COMPILE check `tmp_hb_block_build_claude.sh` PASS (NSTACK=1,2,3 x
  L1/L4 all OK). (b) PARITY `test_hasenbusch_block_validate_claude.cu` (same phi -> S + force vs the summed
  NSTACK serial `HasenbuschPF`, rel < 1e-6): **NSTACK=1 (Nf=2) PASS** and **NSTACK=3 (Nf=6, L2) PASS** -- the
  real block packing reproduces the serial. Handoffs: `tmp_hb_block_validate_claude.sh` (GPU-configurable).
- **Files (final)**: block Term B in `includes/overlap_force_Nfblock_claude.h`; `HasenbuschPFBlock` in
  `includes/pseudofermion_hasenbusch_block_claude.h`; unified `-DNF` driver `hmc_hasenbusch_block_claude.cu`
  (all Nf>=2, NSTACK=Nf/2; serial `hmc_hasenbusch_claude.cu` kept); ML integrator `hmc_hasenbusch_ml_claude.h`
  UNCHANGED (block PF has the same interface). Parity harness `test_hasenbusch_block_validate_claude.cu`.
- **Next (production)**: a short `dH~tau^2` / acceptance chain for a real Nf>=4 block run; then production.

## Ordered chunks

1. **Block Term B force** (Files: overlap_force_Nfblock_claude.h). Add `compute_bilinear`; unit-check vs the
   serial Term B on one config (block col 0 == serial).
2. **HasenbuschPFBlock** (Files: pseudofermion_hasenbusch_block_claude.h). Two BlockedMats + per-frame block
   buffers; block `gen`/`update_eta*`/`S`. Validate `S` + heatbath identity vs serial `HasenbuschPF`.
3. **Block force in the PF** (Files: pseudofermion_hasenbusch_block_claude.h). `get_force_frames` = block
   Term A (`bforce.compute`) + block Term B (`compute_bilinear`) per frame. Validate force vs serial.
4. **Block level + driver** (Files: hmc_hasenbusch_ml_claude.h, hmc_hasenbusch_claude.cu). Block
   `FermionGroupLevel`; driver NSTACK wiring + device-limit check. dH-validate a short chain vs serial.

## Decisions (NM 2026-07-14, "proceed with blocking")

- Multiple pseudofermions ARE needed ($N_f/2$ independent $\phi$ raise the det power; a common factor only
  rescales $\phi$ -> still 1 power). RHMC (one PF, rational $(D^\dagger D)^{-N_f/2}$) is the alternative but a
  nested-rational-on-overlap -> heavier; blocking chosen.
- **Q1/Q3 (revised, NM "unified script") = ONE blocked driver for ALL $N_f$.** Copy `hmc_hasenbusch_claude.cu`
  -> NEW `hmc_hasenbusch_block_claude.cu`, built with `-DNF=` ($N_f$ constexpr, NSTACK$=N_f/2$). It handles
  $N_f=2$ too via NSTACK=1 (no serial fallback inside it -- uniform block path). The original serial
  `hmc_hasenbusch_claude.cu` is KEPT UNCHANGED (reference / current running ensembles); the blocked driver
  becomes the go-to going forward. (NSTACK is a template param, so $N_f$ is fixed at build via `-DNF`.)

## Open questions

1. **Compile-time NSTACK.** `BlockedMat<N,NSTACK,Op>` needs NSTACK = $N_f/2$ at COMPILE time; the driver
   currently reads $N_f$ from `argv[2]`. Options: (a) a dedicated blocked driver with $N_f$ `constexpr` /
   `-DNF=` (like `hmc_Nfblocked`); (b) runtime `if(Nf==4)...else if(Nf==6)...` template dispatch. Propose (a).
2. **Two block pools (action $n_\text{act}$ + force $n_f$)** roughly double the block scratch
   ($16 N\cdot$NSTACK $+ 3 N\cdot$NSTACK$\cdot$npole per op). Confirm it fits at the target $N_f$ (4/6) and
   the biggest L (L4, $N=41472$). Device-limit check at startup (reuse `block_bytes`).
3. **$N_f=2$ fallback**: keep the serial `HasenbuschPF` path (NSTACK=1 blocking is a no-op). Blocked path is
   $N_f\ge4$ only. OK?
4. **Term B RHS packing**: the ratio-frame bilinear bra is $\phi_i$ (block), ket is $\eta_i$ (block) -- both
   already per-frame blocks, so `compute_bilinear` takes two $N\cdot$NSTACK blocks. Confirm no extra seam.
