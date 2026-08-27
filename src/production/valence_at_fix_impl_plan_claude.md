# Valence $a_t$ fix -- derive $a_t$ from the ensemble, not hardcoded 0.2

## Problem

The valence measurement drivers hardcode `const double at = 0.2;` and pass it into the
temporal Wilson-Dirac operator (`DiracExt`), where it sets the temporal hopping
$\kappa_t = (1/\nu_0)\, A_n / \bar\ell / a_t$ (`includes/dirac_ext.h::set_kappa_t`).

For the half-$a_t$ ensembles (`at0.100000` in the dir name) the sea gauge configs were
generated with $a_t = 0.1$ (HMC uses `-DAT_VAL=0.1`), but the valence operator was built
with $a_t = 0.2$ -- a sea/valence temporal-spacing mismatch. The at0.2 ensembles are correct.

The affected at0.1 valence results (conn + disc, 6 ensembles: L1 Nf2/4/6 g1.0, L2 Nf2/4/6 g2.0)
have been archived into `data_.../archive_at0.2/`. This plan fixes the drivers and the at0.1
results get recomputed with $a_t = 0.1$.

## Fix

Derive $a_t$ from the ensemble dir name (token `atD.DDDDDD`, uniquely "at"+digit in the name),
so it can never silently mismatch the ensemble again. Add a `--at <x>` flag that OVERRIDES the
derived value (needed for free-field `U=1` runs, which have no ens-dir; default there is 0.2).

Resolution order (in `main`, after `ParseArgs`):
- if `--at` given (>= 0): use it.
- else if `ens_dir` non-empty: `at = at_from_ensdir(ens_dir)` (assert it parsed).
- else (free field): `at = 0.2`.
Print the resolved `at` prominently in the run header.

Helper (static free function, added near `seed_from_string`):
```
// extract a_t from an ensemble dir name: the unique 'at' immediately followed by a digit,
// e.g. "...gsq2.000000at0.100000nu0..." -> 0.1 .  Returns -1.0 if not found.
static double at_from_ensdir(const std::string& s);
```

## Files

- `jj_local_ylm_scalar_conn_stoch_claude.cu` (conn: VECTOR+AXIAL+scalar) -- at at line 262.
  Files: add helper, add `--at` long_opt + default sentinel, resolve at in main, use in ctor+assert.
- `jj_local_ylm_disc_stoch_claude.cu` (disc) -- at at line 233. Files: same edits.
- `jj_local_ylm_scalar_conn_stoch_fnal_claude.cu` (fnal conn) -- at at line 265. Files: same edits
  (for consistency; at0.1 recompute runs the LOCAL conn driver, but keep fnal in sync).

## Out of scope (flagged, not touched)

- `jj_scalar_condensate_eo_stoch_claude.cu:238` -- same latent bug, but condensate at0.1 was
  not archived / not part of this recompute. Leave unless NM asks.
- glue drivers (`glue_f2*`, `glue2_msm*`, `glue_therm*`) -- $a_t$ there is a pure-gauge
  analysis-time normalization; glue at0.1 GEVP already handled via `run_glue_gevp_at01` (AT=0.1).

## Rerun (MEASUREMENT -- user runs via handoff scripts)

After the fix builds: rerun conn + stride-2 disc for the 6 at0.1 ensembles only. The existing
run scripts pass `--gsq/--Nf/--nu0`; they will get `--at` derived automatically now (no flag
needed) OR we pass `--at` explicitly per ensemble. Handoff scripts to be prepared after code fix
is validated.

## Open questions

- none blocking; proceeding with auto-derive + `--at` override.
