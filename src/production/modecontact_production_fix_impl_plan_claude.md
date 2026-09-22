# Production fix: mode-space contact in the P+ sigma^2 GEVP chain, and the L2 (2,2) state

## Goal
Apply the MODE_CONTACT fix (mode-space GW contact `-1/2 I_Nv` on the peram, not position-space `-1/2 I_2Ns`) to the
PRODUCTION P+ sigma^2 flavor x geometry GEVP chain, rebuild the L2 truncated (Nv=24) cache with it, and show the
L2 spectrum with the single-meson leak removed -- in particular the (2,2) single-meson excitation (Delta=4) and the
two-meson, which the leak previously buried under a spurious m_PS ground.

## Mechanism (settled)
The "leak" is the position-space contact breaking GW anti-hermiticity of M = D_ov^{-1}-1/2 under truncation
(revives the tadpole D_S ~ #removed modes). Mode-space contact keeps P M P anti-hermitian. Full derivation:
`gw_antiherm_exclusion_mechanism_claude.md`. Verified free + interacting L2 (C_22: m_PS -> two-meson).

## Design
- The contact flag already exists: `fs_gevp_point_claude.py` `MODE_CONTACT` (env; both `AblkS` bodies). The
  builders route through `G.make_config` / `G.make_config_win` -> AblkS, so `MODE_CONTACT=1` in the env flows in.
- The BLOCKER is the on-disk cache: builder skips rebuild if the file exists, and the filename has no contact tag,
  so a rebuild would either be skipped or silently overwrite the (buggy) old cache.
- FIX: tag the cache filename with the contact mode -- `_mc1` when MODE_CONTACT=1, "" otherwise (backward compatible;
  the existing position-contact caches are untouched and still found by the default path). Readers pick the tag from
  the SAME env, so `MODE_CONTACT=1` selects the fixed cache end to end.

## Chunks

### Chunk 1 -- cache tag in the builder + P+ readers
Files: `sigma2_flavorgeom_full_v2_nsrc2_claude.py` (builder), `sigma2_Peven_6x6_gevp_nsrc2_claude.py`,
`sigma2_Peven_partialhankel_claude.py`, `sigma2_Peven_partialhankel_fit_claude.py`,
`sigma2_Peven_partialhankel_scan_claude.py` (readers).
- Add `MC = int(os.environ.get("MODE_CONTACT","0"))` and `MCTAG = "_mc1" if MC else ""`.
- Builder cache name: `..._nsrc2_d%d%s_claude.npy % (..., SPLIT, MCTAG)`.
- Reader globs: append `%s` MCTAG before `_claude.npy` so MODE_CONTACT=1 finds only the `_mc1` cache (and a fallback
  message if none yet). Keep the default (MC=0) glob exactly as-is.
- Add an optional `NCFG` cap to the builder (env) if not present, so a reduced-cfg fixed cache can be built fast for
  a first look (the `_mc1` tag makes the reader pick it regardless of cfg count).

### Chunk 2 -- rebuild the L2 fixed cache (ANALYSIS, I run)
`MODE_CONTACT=1 ENS=<L2 g1.0> LREF=2 NVDIR=distill_Nv24_v2 NCFG=<~150 first> python3 sigma2_flavorgeom_full_v2_nsrc2_claude.py`
-> `sigma2_flavorgeom_FULL_..._L2_<ncfg>cfg_nsrc2_d1_mc1_claude.npy`. Background (heavy 9x9 at L2). Full 400 after.

### Chunk 3 -- L2 P+ GEVP with the fix -> the (2,2) state
`MODE_CONTACT=1 ENS=<L2> LREF=2 ... python3 sigma2_Peven_6x6_gevp_nsrc2_claude.py` (and partial-hankel).
Expect: ground no longer m_PS; the physical 0++ tower -- (2,2) single-meson excitation (Delta=4) and the two-meson
(2m_PS) -- resolved. Plot with m_PS / 2m_PS reference lines; SendUserFile.

## Open question
- (2,2): reached via the sigma3-odd O_A kernel; the s2/O2m/O1m geometry basis overlaps it. Confirm the fixed 6x6
  shows it as a distinct level (vs needing O_A explicitly). If the geometry basis is insufficient, add O_A.
