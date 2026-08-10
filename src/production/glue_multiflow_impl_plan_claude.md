# Multi-flow (multi-smearing) glueball operator basis -- impl plan

**Method source (mandatory citation):** the multi-smearing-level variational basis for glueballs is
C. Morningstar & M. Peardon, *"The glueball spectrum from an anisotropic lattice study"*,
Phys. Rev. D **60**, 034509 (1999), **arXiv:hep-lat/9901004** (same shapes measured at several
smearing levels, then a GEVP over the enlarged basis). Smearing here = Wilson/gradient flow,
M. Lüscher, **arXiv:1006.4518**. GEVP: Lüscher & Wolff, Nucl. Phys. B **339**, 222 (1990).
Put both citations in the driver comment blocks.

## Goal

Enlarge the glueball operator basis from **7 -> 28 operators per $(l,m)$ block** by measuring the
SAME 7 shapes at **4 flow times** $t_\text{flow} = 0.5, 1.0, 2.0, 3.0$ (fixed $dt = 0.01$, the
current value). Each flow time is a different smearing radius $\approx \sqrt{4 t}$, so the operators
have genuinely different overlaps with the state -- unlike the current 7 shapes, which all sit at one
smoothing scale and are strongly correlated.

Target channels: **F $\ell=2$ at L4** and **$F^2$ $0^{++}$** (currently nan / 7-11% errors).
F $\ell=1$ is already 2.5-4% and has little to gain.

Cost is near-free in flow integration: the flow is applied **in place**, so one trajectory
$0 \to 0.5 \to 1.0 \to 2.0 \to 3.0$ (300 steps total vs 200 now) yields all 4 levels; only the (cheap)
shape-operator evaluations multiply by 4.

## Key design decision: reuse the `n_shapes` axis

The analysis binary auto-detects the basis from the h5 (`n_lm`, `n_shapes`) and indexes shape
operators as `op = ishape*n_lm + ilm`. So we do NOT need a new axis or any analysis-side change:
store the flow level in the shape index,

```
op = (iflow*n_shapes_geom + ishape)*n_lm + ilm ,   n_shapes_written = N_FLOW * 7 = 28
```

and write `n_shapes = 28`. `glue_gevp_analysis_claude.o` then builds a 28x28 per-m GEVP unchanged.
Its existing `orbits_arg` (iorbit subset, e.g. `"0,3"`) lets us select flow subsets at ANALYSIS time
-> we can scan basis size (7 / 14 / 21 / 28) without re-measuring.

## Files to modify

- `glue2_msm_shapes_claude.cu` -- F (linear) driver: segmented flow + iflow loop + n_shapes/h5 fields.
- `glue_f2_shapes_claude.cu` -- F^2 (squared) driver: same change (identical flow block at :282-284,345-346).
- `run_glue_shapes_multiflow_claude.sh` -- NEW sweep script (copy of `run_glue_shapes_redo_claude.sh`,
  new binaries `*_mf_L<L>_claude.o`, so the existing single-flow sweep/scripts stay intact).

**Output prefix MUST change** (`glue_msm_shapes_mf`, `glue_f2_shapes_mf`): the existing single-flow
data is complete and good, and the per-config "complete" gate would otherwise skip everything.
No `rm` anywhere; the old data is untouched.

## Chunks

### Chunk 1 -- segmented flow + iflow loop in the F driver
Files: `glue2_msm_shapes_claude.cu`
- Replace the single `Flow flow(&SW, 2.0, 200)` with fixed `FLOW_DT = 0.01` and
  `FLOW_T[4] = {0.5, 1.0, 2.0, 3.0}`; build 4 segment flows with lengths `{0.5, 0.5, 1.0, 1.0}` and
  steps `{50, 50, 100, 100}` (dt preserved exactly).
- Wrap the existing per-orbit measurement in `for(iflow=0..3){ seg[iflow](Uflow); <measure> }`,
  applying each segment IN PLACE to the same `Uflow` (cumulative flow times).
- Consolidate orbits->shapes as now, then place into `op = (iflow*7 + ishape)*n_lm + ilm`.
- h5: write `n_shapes = 28`, plus provenance datasets `n_flow = 4` and `flow_times = {0.5,1,2,3}`.
- `H5PREFIX = "glue_msm_shapes_mf"`.

### Chunk 2 -- same for the F^2 driver
Files: `glue_f2_shapes_claude.cu`  (mirror of Chunk 1; `H5PREFIX = "glue_f2_shapes_mf"`)

### Chunk 3 -- sweep script + build
Files: `run_glue_shapes_multiflow_claude.sh`
- Copy of the redo sweep; `for L in 1 2 3 4`; binaries `glue2_msm_shapes_mf_L<L>_claude.o`,
  `glue_f2_shapes_mf_L<L>_claude.o`; stride 1; `NWORK=16`; CPU-only so it runs concurrent with the
  GPU disc job; complete-gated/resumable.

### Chunk 4 -- validation (BEFORE the full sweep)
Files: (none -- run + compare)
- Run the mf driver on ONE ensemble, few configs.
- **Correctness check:** the `iflow=2` sub-block ($t=2.0$) must reproduce the existing single-flow
  `glue_msm_shapes` values to round-off -- same dt, same total steps, in-place segmentation.
- Sanity: correlators at larger $t_\text{flow}$ are smoother / smaller-error at large $dt$.

### Chunk 5 -- basis scan + decide
Files: `glue_gevp_*` analysis calls (no code change; use `orbits_arg`)
- On a high-stat ensemble (L3 Nf2 g1.5, 780 cfg) and a hard one (L4 Nf4 g4, 387 cfg), run the per-m
  GEVP with flow subsets: `{2}` (= current baseline), `{1,2}`, `{0,1,2}`, `{0,1,2,3}`.
- Compare fitted mass + error for F $\ell=1$, F $\ell=2$, $F^2$ $0^{++}$.
- Accept the basis size that minimizes the error WITHOUT destabilizing the fit (28 ops from ~400
  configs may be over-parameterized; `rtol` pruning is the guard).

## Cost / disk

- Flow integration 200 -> 300 steps (1.5x) + 4x operator evaluations -> roughly **2x** the current
  sweep, i.e. a few hours at `NWORK=16` (glue is ~6-9 s/config now).
- h5 per config: dominated by `F_corr_blk` = nsep x n_lm*n_shapes^2, which scales as **n_shapes^2**,
  so 7 -> 28 is **16x**, NOT 4x:  17 x 8 x 28^2 x 8 B = **~0.85 MB/cfg** (vs 57 KB now).
  Over ~93k configs ~ **79 GB** (4.3 TB free) -- affordable, but it is the main argument for
  restricting the sweep to L2/L3/L4 (open question 2) and/or dropping a flow level.

## Open questions

1. ~~Flow times~~ **RESOLVED**: $\{0.5, 1.0, 2.0, 3.0\}$ at fixed $dt=0.01$ (NM). At L1 (20 faces)
   $t=3.0$ is likely over-smeared -- harmless, droppable at analysis time via `orbits_arg`.
2. ~~Which ensembles to re-sweep?~~ **RESOLVED: ALL 42** (NM). Disk (79 GB of 4.3 TB) is not the
   constraint. Including L1 keeps the basis UNIFORM across L -- otherwise a 7-op L1 vs 28-op L2/3/4
   would inject a methodology change into the ratio-vs-1/L^2 plots, i.e. a systematic that looks like
   L-dependence. L1 also has the most configs (3999) = the cleanest testbed for how much the enlarged
   basis actually buys. The t=3.0 over-smearing worry at L1 is handled at ANALYSIS time (`orbits_arg`
   flow-level subset), so nothing is lost by measuring it.
3. ~~Keep `SQUARED`/face_sign conventions identical~~ **RESOLVED**: yes, unchanged.
4. **SCOPE (NM 2026-08-07): flow levels ONLY.** No new shapes, no cross-product operators -- see the
   appendix for why the cross-products were dropped.

## Appendix -- why instance-wise cross-products were DROPPED (analysis, 2026-08-07)

Considered adding local bilinears $\theta_1\theta_2$ (same anchor) instead of only $\theta^2$. It
turns out to add ~nothing, because the face-anchored shapes are LINEARLY DEPENDENT by construction.
All four are $w=+1$ face sums over $f_0$ and its 3 edge-neighbours $n_1,n_2,n_3$:

- $T=\{f_0\}$, $R=\{n_1,n_2,n_3\}$, $S=\{f_0,n_1,n_2,n_3\}$, $H_k=\{f_0,n_i,n_j\}$ (3 placements)

`holonomy()` is a linear weighted face sum with common `face_signs`, so $\theta_S=\theta_T+\theta_R$
EXACTLY, hence

$$\theta_T\theta_R = \tfrac12\left(\theta_S^2-\theta_T^2-\theta_R^2\right)$$

i.e. the $T\cdot R$ cross-product is already a linear combination of operators we measure today; the
GEVP depends only on the SPAN, so it adds nothing (and makes the metric exactly singular).

More generally the face-anchored shapes involve only **4 elementary fluxes** $\{t,a,b,c\}$, so the
local bilinear space is $\dim = 10$, and the part surviving the equal-weight orbit consolidation
(which symmetrizes over the 3 neighbours) is only 4-dimensional, spanned by
$t^2,\; t(a{+}b{+}c),\; a^2{+}b^2{+}c^2,\; ab{+}bc{+}ca$ -- and the four CURRENT operators
$T^2, R^2, S^2, \sum_k H_k^2$ already span it (the $4\times4$ change-of-basis has $\det=4$).

What WOULD enlarge the span (kept for later, NOT in scope now):
(a) shapes reaching past the first ring (2-ring contour, annulus, scaled triangle) = more elementary
fluxes; (b) DE-CONSOLIDATING the orbits -- at $L\ge2$ the 3 $H$-placements live in different Ih
orbits (`wilson_shapes_claude.h:157`) and are summed away; they are already computed and discarded.
