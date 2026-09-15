# Adapt the mixing + flavor-GEVP analysis to the nsrc=2 (`distill_Nv24_v2`) data -- impl plan

## Goal
Run the F$^2$-$\sigma^2$ glueball-mixing analysis and the flavor-cross GEVP on the higher-statistics
**nsrc=2** perambulators (`distill_Nv24_v2`, source windows at `tsrc_list=[0,64]`, maximally separated on
the $N_t=128$ torus, so nearly independent). Use BOTH source windows per config (average), roughly doubling
the source statistics per gauge config. Validate each `_nsrc2` script against its v1 parent by checking that
**window 0 on configs common to v1 reproduces v1 to machine precision** (same `tsrc0=0`, only a different
output dir) -- v1 and v2 are otherwise different config subsets, so full correlators differ.

## Naming
Keep the existing analysis tag `_v2` and append the data tag `_nsrc2`: `<existing>_nsrc2_claude.py`.
The `_v2` = corrected-FS/triple-subtraction ANALYSIS version; `_nsrc2` = the nsrc=2 DATA.

## Data facts (verified)
- v2 peram `peram/tau` shape `(2, 32, 32, 24, 24)`; `meta/tsrc_list = [0, 64]`; `meta/nsrc = 2`; `twin=32`.
- v2 configs so far: 233 (stride-4, run still going -> ~1000). v1: 400 configs, single window `tsrc0=0`.
- `dc.load_peram_windows(k)` -> `(V, [(tsrc0_s, tau_s, taugw_s), ...])`; `dc.load_peram(k)` = window 0
  (backward-compat, unchanged). NVDIR selects the dir: `distill_Nv24` (v1) vs `distill_Nv24_v2` (v2).

## Shared helper (one edit to a shared file, ADD only)
Add `make_config_win(k, w)` to `fs_gevp_point_claude.py`, a sibling of `make_config` that builds `AblkS`/
`AblkSt` from source window `w` via `dc.load_peram_windows` (window's own `tsrc0`). `make_config` is left
untouched. `make_config_win(k, 0)` on a v1 file is IDENTICAL to `make_config(k)` (window 0) -- the validation
hook. Returns `(AblkS, AblkSt, twin, nsite, U, tau, tsrc0_w)`.

## Correlator combination over windows (POST-CONTRACTION ONLY)
CRITICAL: $\tau$ is NEVER averaged across windows. Each window's propagator blocks are built from that
window's own $\tau_w$ (`make_config_win(k,w)`), and the ENTIRE diagram contraction is done separately within
window `w` -- this already translation-averages over the source times $s$ inside the window -- yielding a
fully-contracted correlator $C_w(dt)$ (any global object, F$^2$ over $N_t$, indexed at absolute times relative
to that window's `tsrc0_w`). ONLY THEN are the finished correlators averaged at fixed separation:
$$
C(dt) = \frac{1}{n_\text{win}} \sum_w C_w(dt).
$$
No cross-window object is ever formed (never `tau[0]` and `tau[1]` in one contraction). This is exactly the
within-window source translation-average, extended to the second (well-separated) window. One sample per gauge
config, so the jackknife/binning is unchanged.

## Chunks

### Chunk 1 -- helper + the validated O_2m cross + VALIDATION
Files: EDIT `fs_gevp_point_claude.py` (add `make_config_win`); NEW `f2_sigma2_cross_o2m_v2_nsrc2_claude.py`
(copy of `f2_sigma2_cross_o2m_v2_claude.py`; loop windows via `make_config_win`, index F$^2$ per `tsrc0_w`).
VALIDATE: (a) on v1 data, `make_config_win(k,0)` path == `make_config(k)` path (machine precision);
(b) run on v2 (NVDIR=distill_Nv24_v2), report max|S/N| and effmass vs the v1 numbers (max|S/N|=3.73,
$a_t m\sim0.30$-$0.40$). Deliverable: does v2 sharpen the cross?

### Chunk 2 -- remaining mixing scripts
Files: NEW `f2_sigma2_cross_v2_nsrc2_claude.py` (local $\sigma^2_{00}$), `f2_sigma2_cross_o1m_v2_nsrc2_claude.py`
(coincident), `f2_o2m_gevp_v2_nsrc2_claude.py` (coupled 2x2 {F$^2$,O_2m} GEVP). Same window loop. Compare to v1.

### Chunk 3 -- flavor cache builder + P+/P- GEVP on v2
Files: NEW `sigma2_flavorgeom_full_v2_nsrc2_claude.py` (build the 9x9 flavor$\times$geom cache with both
windows -> `sigma2_flavor_cache_claude/sigma2_flavorgeom_FULL_..._nsrc2_..._claude.npy`); NEW
`sigma2_Peven_6x6_gevp_nsrc2_claude.py`, `sigma2_Podd_3x3_gevp_nsrc2_claude.py` (read the v2 cache).
VALIDATE: v2 window-0 cache on common configs == v1 cache (per-config, machine precision), then compare the
0.46/0.62/~0.83/~0.92 tower and errors v2 vs v1.

## Open questions
- None outstanding (scope = mixing + flavor GEVP; suffix = `_nsrc2`; use both windows averaged; validate
  window-0 vs v1 on common configs).

Refs: distillation Peardon 0905.2160; corrected FS/triple-subtraction `fs_diag_corr_v2_claude.py`;
mixing CP arXiv:1603.05582; parent plan `f2_sigma2_mixing_redo_impl_plan_claude.md`,
`ps_fs_flavor_cross_impl_plan_claude.md`.
