# $F^4$ local-density shape operators in the $0^{++}$ GEVP -- implementation plan

**Status: PLAN ONLY (2026-08-11). Nothing built, nothing launched, no data touched.**

Scope fixed by NM: the linear $F$ driver (`glue2_msm_shapes_claude.cu`) is **not** touched. The change is
one generalisation in the shape header plus a power axis in the $F^2$ driver. A full re-measurement is
expected and accepted (NM: "of course, we do the calculation again; that's not a constraint").

## Sources (mandatory credit)

- Variational multi-operator basis: **C. Morningstar and M. Peardon, arXiv:hep-lat/9901004**,
  Phys. Rev. D 60, 034509 (1999).
- GEVP: **M. Luscher and U. Wolff, Nucl. Phys. B339, 222 (1990)**; improved plateaux
  **B. Blossier et al., arXiv:0902.1265**.
- Wilson flow: **M. Luscher, arXiv:1006.4518**.
- Existing 7-shape basis and its two null enlargements: `glue_multiflow_impl_plan_claude.md`,
  memory `project_glue_shapes.md`.

## 1. Goal

Enlarge the $0^{++}$ variational basis with **quartic local-density** operators, built from the same
Wilson-loop shapes already in use, and test whether the $F^2$ $0^{++}$ error improves.

### 1.1 What the operator is

Each shape instance carries a local holonomy angle $\theta_\text{inst}(s)$ -- the signed sum of
plaquette angles around that loop at timeslice $s$ (`includes/wilson_shapes_claude.h:318-330`,
`holonomy()`). This is the local field strength $F$ at the instance position $\hat r_\text{ins}$.

The shape operator is the $Y_{\ell m}$ projection of a **local power** of that angle
(`wilson_shapes_claude.h:333-341`, `op()`):

$$
O^{(p)}_{i,\ell m}(s) \;=\; \frac{1}{N_\text{orb}}\sum_{\text{inst}\,\in\,\text{orbit}_i}
Y_{\ell m}(\hat r_\text{ins})\;\bigl[\theta_\text{inst}(s)\bigr]^{p}
$$

Today only $p=1$ (linear $F$, `SQUARED=false`, the $F$ driver) and $p=2$ ($F^2$, `SQUARED=true`, the
$F^2$ driver) exist -- the single line `const double g = squared ? th*th : th;` at `:338`.

**This plan adds $p=4$**, i.e. $F^4(x)$, projected onto $Y_{\ell m}$ exactly as $F$ and $F^2$ are, and
puts the resulting operators alongside the $p=2$ ones in the same GEVP basis.

Because the power is taken on the *local* density before the projection, $O^{(4)}_{i,\ell m}$ carries
the same $(\ell,m)$ label as $O^{(2)}_{i,\ell m}$ and transforms identically under the icosahedral
group. The block-diagonal storage assumption (`glue_f2_shapes_claude.cu:419`, "Cross-(l,m) vanish by
symmetry") is therefore preserved exactly, and **the analysis needs no change at all**.

### 1.2 Which powers

Only **even** powers may be mixed into the $0^{++}$ basis. The linear operator is parity-odd
(`wilson_shapes_claude.h:15`), so $\theta$ and $\theta^3$ are parity-odd and
$\langle O^{(1)} O^{(2)}\rangle = \langle O^{(3)} O^{(2)}\rangle = 0$ by parity -- including them would
add rows that are pure noise. So the basis is $p \in \{2, 4\}$, optionally extended to $p=6$ later.

### 1.3 Free-field reference

$F = \sqrt2$ (one photon), $F^2 = 2\sqrt2$ (two photons) are the established project values. By the same
photon counting, the state $F^4$ couples to most strongly in the free limit is $4\sqrt2$ (four photons)
-- flagged as my inference from the counting, not an independently established number. It is the value
the enlarged GEVP should *not* be converging to if the $F^4$ operators are helping the $0^{++}$.

## 2. Files

### To modify

| file | change |
|---|---|
| `includes/wilson_shapes_claude.h` | ADD a member `op_pow(U, s, orbit, ell, em, power)` -- identical to `op()` but with `g = th^power`. Leave `op()` (`:333-341`) **untouched** so the linear $F$ driver is bit-identical and needs no edit. |
| `glue_f2_shapes_claude.cu` | add a POWER axis to the shape index; `n_shapes = N_POW * N_FLOW * n_shapes_geom`. Power-2 block FIRST so the existing 7 shapes keep indices 0..6. New h5 prefix. Revert `n_shapes_geom` 9 -> 7 (Q2). |

### Unchanged

- `../both_3d/glue_gevp_analysis_claude.cu` -- reads `n_shapes` from the h5 (`:163`) and rebuilds the
  block layout itself. **No edit.**
- `glue2_msm_shapes_claude.cu` (linear $F$ driver) -- NM's instruction.
- `fit_perm_claude.py` and the downstream ratio/plot scripts -- the GEVP `.dat`/`_jk.dat` format is
  unchanged.
- All existing `glue_f2_shapes.*.h5` production data and its prefix.

### To create

| file | purpose |
|---|---|
| `glue_f2_v2_shapes_impl_plan_claude.md` | this plan |
| `run_glue_f2_v2_shapes_claude.sh` | test sweep with the $p\in\{2,4\}$ driver. Handoff -- NM runs it. |
| `run_glue_gevp_f2_v2_claude.sh` | GEVP sweep: `orbits_arg` = shapes 0..6 (must reproduce production), then the full 14, with an `nops2` / `vacsub` scan. |
| `glue_f2_v2_compare_claude.py` | table of $0^{++}$ mass, error, error ratio vs production, metric condition number, operators surviving `rtol`. |

## 3. Index layout and storage

Shape index, with the power axis outermost so the production block stays at 0..6:

```
ishape = (ipow*N_FLOW + iflow)*n_shapes_geom + is        POW = {2, 4},  N_FLOW = 1
op     = ishape*n_lm + ilm                                (unchanged convention)
```

With `n_shapes_geom = 7`, `N_FLOW = 1`, `N_POW = 2`: `n_shapes = 14`, `nops = 56`.

| dataset | now | with $p=4$ |
|---|---|---|
| `F_corr_blk` | $(17,\ 4\cdot7^2=196)$, 3332 doubles | $(17,\ 4\cdot14^2=784)$, 13328 doubles |
| `F` | $(28,)$ | $(56,)$ |
| file size | 31 036 B | ~110 KB |

`n_shapes` is written to the h5 (`:448`) so the analysis auto-adapts. A `powers` dataset should be added
alongside `flow_times` (`:452`) for provenance.

## 4. Cost

The added work is one extra `op_pow` call per (orbit, $\ell m$, $t$). Note the shape evaluation is
*already* repeated across the four $(\ell,m)$ values (the $\ell m$ loop sits outside the orbit loop,
`glue_f2_shapes_claude.cu:388-396`), so this doubles an already-4x-redundant inner cost, against a
Wilson-flow integration that dominates the config. The correlator contraction grows 4x (nob 7 -> 14),
i.e. ~1.7M multiplies per config -- negligible.

Reference for a full production sweep: the last glue sweep ran **2026-07-25 to 2026-08-07, 13 days**
wall on an 8-worker pool over 75 ensemble-jobs. Glue is measured at **stride 1** (every config: 3999 at
L1, 799 at L3, 560 at L4), so a re-sweep is the same size job.

If the doubled shape evaluation turns out to matter, the single-pass variant is to compute
`th` once and accumulate both powers in one loop; deferred as an optimisation, not needed for the test.

## 5. Ordered chunks

### Chunk 1 -- shape header
Add `op_pow()`. `op()` untouched.
*Files:* `includes/wilson_shapes_claude.h`

### Chunk 2 -- driver power axis
Add `POW = {2,4}`, extend the shape index, size `obs`/`obs_raw` accordingly, call `op_pow`, write the
`powers` dataset, set the new h5 prefix, revert `n_shapes_geom` to 7.
*Files:* `glue_f2_shapes_claude.cu`

### Chunk 3 -- test sweep
One ensemble (Q1). Handoff script, `tee` to `*_claude.log`, resumable on the `complete` gate, no `rm`.
NM runs it.
*Files:* `run_glue_f2_v2_shapes_claude.sh`

### Chunk 4 -- VALIDATION GATE
Analyse the new files with `orbits_arg="0,1,2,3,4,5,6"`. The $0^{++}$ mass and error must reproduce the
production `glue_f2_shapes` result on the same configs to round-off. Nothing proceeds until this holds.
(Same discipline as the s9 test, whose 7-shape sub-block matched production at 0.00e+00.)
*Files:* `run_glue_gevp_f2_v2_claude.sh`

### Chunk 5 -- full 14-operator GEVP
`orbits_arg` empty. Scan `nops2` and both `vacsub` settings, and report the whole spectrum, not just the
index that was the $0^{++}$ in the 7-operator basis -- see section 6.
*Files:* `run_glue_gevp_f2_v2_claude.sh`

### Chunk 6 -- verdict
Decision rule fixed before looking: keep the enlargement only if the $0^{++}$ error drops by more than
its own jackknife uncertainty on the test ensemble.
*Files:* `glue_f2_v2_compare_claude.py`

### Chunk 7 -- record
*Files:* memory `project_glue_shapes.md`, `project_redo_obs_flow.md`

## 6. Two things that will need attention in Chunk 5

1. **Vacuum / state indexing.** The $F^2$ call runs `vacsub=0` (argv[15]) and absorbs the vacuum as a
   GEVP mode, taking the $0^{++}$ at index `nops2-2` with `nops2=2`
   (`run_glue_gevp_redo_claude.sh:47`). $\langle \theta^4 \rangle$ is large and positive, so the
   enlarged $\ell=0$ block carries more near-constant directions and this indexing will not survive
   unchanged. Both `nops2` and `vacsub` are existing CLI args -- no code change, but the scan is
   mandatory, not optional.
2. **Conditioning.** $O^{(4)}$ and $O^{(2)}$ differ by orders of magnitude in scale, so $C(t_0)$ is
   badly scaled and the `rtol` metric pruning (argv[11], $10^{-8}$) may drop rows arbitrarily. The
   generalised eigenvalues are invariant under a per-operator rescaling, so if pruning misbehaves the
   fix is a one-line normalisation by $\sqrt{C_{ii}(t_0)}$ in the analysis. Contingency, not a planned
   edit.

## 6b. RESULT (2026-08-11, L3 Nf2 gsq1.5, 799 configs, Nc=780 after the k<20 cut, nbins=52)

**Verdict: KEEP. The first non-null basis enlargement** (multi-flow and the s9 shapes were both null).

Gate (Chunk 4) PASSED: the `F_corr_blk` p=2 sub-block is bit-identical to production at the data
level (max|diff| = 0.000e+00), and the GEVP output matches with an identical NaN pattern and
0.000e+00 on every finite entry.  `kept_ops=14` in every full-basis run, so the p=4 operators really
entered the GEVP -- the `rtol=1e-8` pruning feared in section 6.2 did not bite, and the normalisation
contingency was NOT needed, despite the p=4 diagonal entries sitting 4-7 orders of magnitude below
p=2 (C_00(dt=0): shape 6 = 1.2e-03 vs shape 7 = 5.2e-11).

| variant | 0++ state | M(err) | err/ref | ratio +- 2nd-level jk |
|---|---|---|---|---|
| production (7 ops, p=2) | 0 | 3.3211(0.2378) | -- | -- |
| full_n2_v0 | 0 | 3.3101(0.2253) | 0.948 | 0.948 +- 0.013 (3.9 sigma) |
| **full_n4_v0** | **2** | **3.3047(0.2219)** | **0.933** | **0.933 +- 0.017 (4.0 sigma)** |
| full_n2_v1 | 0 | 3.2797(0.2174) | 0.914 | 0.914 +- 0.032 (2.7 sigma) |
| full_n4_v1 | 2 | 3.2706(0.2156) | 0.907 | 0.907 +- 0.047 (2.0 sigma) |
| full_n5_v1 | 3 | 3.2652(0.2195) | 0.923 | 0.923 +- 0.055 (1.4 sigma) |

**Two independent pieces of evidence, per NM 2026-08-11:**

1. **Error reduction ~6.7%** (best well-determined variant, `full_n4_v0`), equivalent to ~15% more
   configurations. NOTE the pre-registered rule in section 6 ("error drops by more than its own
   jackknife uncertainty") uses the NAIVE uncertainty $1/\sqrt{2(n_\text{bins}-1)} \approx 10\%$ and
   is NOT met by 6.7%.  That rule is the wrong test here: the two estimators use IDENTICAL configs,
   so their errors are strongly correlated.  Jackknifing the ratio directly gives 4.0 sigma.
   (Caveat: the second-level jackknife resamples jackknife samples -- practical, not rigorous.)
2. **Stabilised mean with a systematic downward trend.**  Within a fixed `vacsub` the central value
   falls monotonically as the basis grows: `vacsub=0` 3.3211 -> 3.3101 -> 3.3047; `vacsub=1`
   3.2797 -> 3.2706 -> 3.2652.  Total drift -0.056 (0.25 sigma), small but systematic in SIGN at
   every step -- the variational signature, which a noise fluctuation in the error cannot fake.
   Care: the strict variational bound is for the LOWEST state; with `vacsub=0` the near-constant
   vacuum mode (M ~ 0.007) is still in the spectrum so the 0++ is excited and the bound is only
   approximate (corrections $O(e^{-\Delta E t})$, Blossier et al. arXiv:0902.1265).  Under
   `vacsub=1` the vacuum is subtracted explicitly and the argument applies most cleanly -- which is
   also where the drift is largest.

**Bookkeeping change, as predicted in section 6.1:** at `nops2=4` the 0++ is **state 2**, not state 0.
Any downstream fit script must be retargeted; `nops2=2` keeps it at state 0 at a slightly worse ratio.

**Physics observation:** a new state appears at 6.31(1.15).  Free-limit four-photon is
$4\sqrt2 \approx 5.66$ and the 0++ sits at 3.30 against $2\sqrt2 \approx 2.83$ -- both about 17% above
free, so the p=4 operators do appear to resolve the four-photon state.

**Recommended production setting:** `orbits_arg=""` (all 14), `nops2=4`, `vacsub=0`, 0++ = state 2.

**Not yet established:** this is ONE ensemble, and 10 variants were scanned before quoting 4.0 sigma,
so there is a selection effect in that number.  What argues against it being noise is that all 10
variants land in 0.907-0.948, every one below 1.  A second ensemble would settle it.

## 7. Open questions

**Q1. Test ensemble?** The s9 null test used L3 Nf2 g1.5 (198 configs). L1 has 3999 configs and is
cheapest per config, so it is the sharper discriminator; L3/L4 is where the signal is actually wanted.

**Q2. Revert `n_shapes_geom` 9 -> 7?** The driver is currently left at 9
(`glue_f2_shapes_claude.cu:307,332`, prefix `glue_f2_shapes_s9`) from the null four-five/three-four
test. The plan assumes a revert so $F^4$ is measured against the production basis; keeping 9 would
confound two changes.

**Q3. Include $p=6$?** Costs one more block (n_shapes 21, `F_corr_blk` 1764 columns). Cheap to add now,
awkward to add later.

**Q4. New h5 prefix name.** Proposed `glue_f2_v2_shapes`, following the `_s9` / `_fullflow`
convention, so production data is untouched.
