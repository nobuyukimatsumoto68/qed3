# Larger 0++ glueball operator basis via icosahedral-orbit Wilson loops

## Goal / physics

Enlarge the variational basis for the scalar ($0^{++}$) glueball GEVP by adding **distinct spatial
Wilson-loop shapes** (connected AND disconnected), each mapped over the icosahedral group into an
**orbit**, instead of leaning on near-redundant **flow times** (whose collinearity gives an
ill-conditioned $C(t_0)$ metric -- the exact failure that forced the `rtol` pruning / NaN in the current
F$^2$ GEVP). Distinct shapes add genuine rank; flow copies do not.

Works already at **L=1** (icosahedron): beyond the single-triangle plaquette we can enlarge the spatial
loop to cover two triangles (the "rectangle" counterpart), and use disconnected shapes. For **L$\ge$2**
each icosahedral face is subdivided, so the same shapes admit many more orbit patterns depending on the
base point (inside vs. straddling the icosahedral triangle).

## Locked-in decisions (from discussion 2026-07-04)

1. **Drop $E^2$ (temporal / electric).** Too noisy. Keep only magnetic $B^2 = F_{12}^2$ = spatial
   Wilson-loop flux. => the F$^2$ shape driver uses spatial loops only (`NCH=1` equivalent).
2. **Target quantum number = $I_h$ trivial irrep $A_g$** (= continuum $0^{++}$, P-even, T-even).
   The magnetic flux $\theta = F_{12}$ is a pseudoscalar (P-odd) and **T-odd**; $\theta^2$ is **T-even**;
   the scalar glueball sits in $A_g$.
3. **Continuum $\ell$ is not a good quantum number** on the discretized sphere -- only the $I_h$ irrep is.
   $Y_{lm}$ weighting is optional; any $\ell$ that subduces to $A_g$ (l=0, l=6, ...) overlaps the target,
   and higher-$\ell$ operators overlap the lower-$\ell$ (ground) states. The GEVP disentangles the mixing.
   (Same $I_h$ subduction the ylm-tower work saw: $l{=}3\to T_2+G$, etc.)
4. **Each distinct ORBIT = one operator.** Do NOT merge two different orbits even if they share the same
   shape type (same #faces, etc.) -- merging would require an inter-orbit relative normalization we want
   to avoid.
5. **Redundancy is fine on the first pass.** Keep possibly-redundant operators; if the GEVP degenerates,
   THEN sum or drop a subset. Do not over-engineer independence up front.
6. **Reuse existing orbit code**: `s2ising/groups.h` -- `orbits[ix][ig]` = site index that vertex `ix`
   maps to under $I_h$ element `ig` (NIh=120), built from `Ih.rotation(ig, ms, mt)`. A loop labelled by
   its site sequence maps under $g$ by `site -> orbits[site][ig]`; the orbit is the dedup of the NIh images.
7. **No index bridge needed (Q5 RESOLVED).** We build the orbit table NATIVELY in the qed3 lattice: take
   the $I_h$ rotation matrices and apply them to our own `base.sites` (VE unit vectors), matching each
   image to a site -> `orbits[ix][ig]` in qed3 labelling. Reuse only the MATH from `s2ising/groups.h`:
   - `Rotation`: $R(\theta) = \exp(-\sum_i \theta_i J_i)$ (Eigen `.exp()`), $\theta$ = axis $\times$ angle.
   - generators $m_s$ = $\pi$-rotation about an edge-midpoint axis (order 2, $=R(\hat n_{\rm edge}\pi)$),
     $m_t$ = $2\pi/3$ about a face-center axis (order 3), built from the base icosahedron's edge / face.
   - full group $I_h = \langle m_s, m_t, -\mathbb{1}\rangle$ -> generate ALL 120 matrices by CLOSURE (BFS
     under matrix mult, dedup at tol); the multiplication-table FILE is NOT needed.
   - RUNTIME CHECK: every element must permute `base.sites` (each $g\,x$ matches a site to tol); catches
     any icosahedron-orientation mismatch between the generators and the qed3 lattice.
   Axes for $m_s,m_t$: use the base icosahedron (the 12 corner vertices; same C2/C3 axes at every L).

## Operator construction pipeline

1. **Shape** = a set of one or more spatial faces/loops, specified as ordered oriented link sequences
   anchored at a base site (or base sites, for disconnected shapes).
2. **Loop holonomy** (gauge-invariant magnetic flux):
   $$\theta_{\rm loop} = \sum_{\ell \in \partial(\rm loop)} \pm\, a_\ell \;=\; \sum_{\triangle \in \rm loop} \theta_{\rm sp}(\triangle)$$
   (interior shared edges cancel; $a_\ell$ = link angle, sign from orientation).
3. **Shape observable** $O_{\rm shape}$ (T-even, $A_g$-projectable) -- **OPEN (Q1)** for disconnected shapes:
   connected loop $\to \theta_{\rm loop}^2$; disconnected $\{A,B\}$ $\to (\theta_A+\theta_B)^2$ vs
   $\theta_A\,\theta_B$ (different operators; decide).
4. **Orbit + operator**: apply `orbits[.][ig]`, ig$=0..$NIh$-1$, to the shape's site(s) $\to$ NIh images;
   dedup $\to$ orbit. Operator = orbit-sum (optionally $Y_{lm}$-weighted at the base-point centroid) of
   $O_{\rm shape}$. Each orbit is a separate operator (decision 4).
5. **Driver**: fill `obs[op][t]` from the shape-operator list, then the existing correlator / h5 / GEVP
   machinery is unchanged (nops is generic).

## RESOLVED 2026-07-04b -- figure-8 + twist (five shapes)

- Q1 RESOLVED: disconnected shape = **figure-8 SUM** ($\theta_1+\theta_2$), NOT the product. Reason: the
  linear operator must be parity-ODD (no square), which the sum gives and the product (parity-even,
  squared-only) does not. Renamed disc -> figure-8. Two vertex-sharing triangles; insertion = shared vertex.
- **Twist**: each two-triangle shape has a relative-orientation weight $\sigma=\pm1$ between the two
  triangles -> holonomy $\Phi_0+\sigma\Phi_1$, $\Phi(f)=\mathrm{face\_sign}[f]\theta(f)$. Untwisted
  $\sigma{=}{+}$, twisted $\sigma{=}{-}$. Independent info: linear $\Phi_0\pm\Phi_1$; squared spans
  $\{\Phi_0^2{+}\Phi_1^2,\ \Phi_0\Phi_1\}$.
- **FIVE shapes**: triangle, rect, twisted-rect, figure-8, twisted-figure-8. Instance carries per-face
  weights `w[]`; canonicalized (sort faces, global sign w[0]=+1) for dedup. Orbit counts:
  L1 = 1 each (5 total); L2 = 2/2/2/5/5 (16); L4 = 5/6/6/16/16 (49). Both drivers wired + BUILD
  (f2 nops=n_orb*16, msm nops=n_orb*15). CAVEAT: if op-count changes, clear old `_shapes_*` h5 dir
  (resume checks only the complete flag, not nops).

## Two-trace / disc (product) -- DROPPED 2026-07-04 (redundant with twist)

Considered explicit disconnected two-trace operators $\Phi_0\Phi_1$ (T-even, F^2-only). But since we keep
BOTH twisted and untwisted squared operators at the same insertion point,
$O_{\rm disc} = \tfrac14(O_{\rm untw^2}-O_{\rm tw^2})$ EXACTLY -> zero added rank, singular metric.
The twist already spans the two-trace. So NO separate disc operators; the five connected shapes stand.

## Open questions (to resolve together)

- **Q1** Disconnected-shape operator: $(\theta_A+\theta_B)^2$ (sum-then-square) vs $\theta_A\theta_B$
  (product)? Both are T-even; they are different operators / different overlaps. Start with which?
- **Q2** Initial shape list. Proposal: (a) single triangle [= current $B^2$]; (b) two adjacent triangles
  (rhombus loop); (c) disconnected triangle pairs at separation $r$ (a few shells). Confirm / extend.
- **Q3** $Y_{lm}$ weighting: start with $l=0$ only ($A_g$), or also include higher $\ell$ (e.g. $l=6$)
  that subduce to $A_g$? (Higher $\ell$ = extra operators overlapping the target.)
- **Q4** Base-point definition for a multi-site / disconnected shape (which site anchors the $Y_{lm}$ /
  orbit label -- centroid, first site, or per-piece?).
- **Q5** Location of the s2ising$\leftrightarrow$s2n_simp index-mapping code (USER to provide).
- **Q6** Per-operator normalization -- any fixed choice is fine (GEVP is invariant under per-operator
  rescaling); confirm we do not need cross-operator norms (consistent with decision 4).
- **Q7** Keep flow times in the basis or replace? Proposal: replace with shapes; optionally keep only
  $t_{\rm flow}=0$ (unsmeared) + one modest flow as 1-2 extra operators, not a stack.

## Implementation chunks (draft -- to firm up after Q1-Q7)

1. **DONE (2026-07-04)** -- `includes/icos_orbits_claude.h` + `test_icos_orbits_claude.cu`. Group closes
   to 120; every element permutes sites to ~1e-15 at L=1,2,4 (orbit histograms = 12 verts / 30 edges /
   face orbits, textbook). No alignment needed. Original chunk-1 spec:
   **Native orbit generator** (UNBLOCKED, Q5 resolved) -- `includes/icos_orbits_claude.h`:
   port `Rotation`; build $m_s,m_t$ from the qed3 base icosahedron edge/face; close $\langle m_s,m_t,-\mathbb1\rangle$
   to 120 matrices; apply to `base.sites` -> `orbits[ix][ig]`; runtime permutation check. Self-contained
   (no s2ising link, no mult-table file). Reference: `s2ising/groups.h` (math only).
2. **DONE (2026-07-04)** -- `includes/wilson_shapes_claude.h` + `test_wilson_shapes_claude.cu`.
   GENERIC: Instance = set of faces (a loop); orbit generated by applying the 120 elements to a base
   instance (`orbit_of`); `orbits_from(candidates)` partitions; `holonomy = sum face_sign*theta(face)`
   (interior edges cancel), `centroid`, `op(orbit,ell,em,squared)`. Shapes = enumerators of base
   instances: `triangles()` (1 face), `rectangles()` (link-adjacent face pair). Orbit counts:
   L1 tri1/rect1, L2 tri2/rect2, L4 tri5/rect6. face_signs are mixed +-1 -> face_sign applied for a
   clean linear Y_lm moment (harmless when squared). Instance carries an INSERTION POINT r_ins (Y_lm
   eval point) that rotates with the group: triangle -> face centroid, rectangle -> shared-link
   midpoint (per NM). Adding a shape = one more enumerator.
   Original chunk-2 spec below.
3. **Operator generator** -- enumerate orbits per shape, emit the operator list (each = set of site
   sequences + optional $Y_{lm}$ weight). Files: same header / driver.
4. **DONE (2026-07-04)** -- `glue_f2_shapes_claude.cu` (SQUARED, ell 0-3, `_shapes_f2/`) +
   `glue2_msm_shapes_claude.cu` (LINEAR, ell 1-3, `_shapes_msm/`). Single flow tmax=1.0/100 (Q7).
   obs[op][t] = shp.op(Uflow, t, orbits[iorbit], ell, em, squared), op = iorbit*n_lm + ilm;
   orbits = orbits_from(triangles()) + orbits_from(rectangles()). nops = n_orbits*n_lm (driver prints
   n_orbits). h5 output. Both BUILD (nvcc+h5).
5. **DONE** -- analysis = existing `glue_gevp_analysis_claude.cu` UNCHANGED: run with
   N_FLOW=n_orbits, NCH=1, drop_l0 (0 squared / 1 linear); opset=tile(lm_set,n_orbits), zero_noise via
   op%n_lm all still valid.
6. **L=1 smoke DONE** (11 cfg): f2_shapes nops=32 / msm_shapes nops=30; GEVP finite (linear +ve plateau
   ~0.4 even at 11 cfg). TODO: run scripts (full ensembles), then L>=2 (richer basis), then disconnected.

## What each of us brings

- USER: the s2ising$\leftrightarrow$s2n_simp index-mapping code (Q5); answers to Q1-Q4, Q7.
- CLAUDE: refine chunks, prototype the shape/holonomy + orbit generation once the bridge + Q's land.

## References
- Orbit code: `s2ising/groups.h` (`Orbits`, `orbits[ix][ig]`, NIh=120, `Ih.rotation`).
- Current F$^2$: `glue_f2_claude.cu`, `includes/action_ext_claude.h` (density_Ylm_spatial = $B^2$).
- QN / $I_h$ subduction context: [[project-ylm-tower]] (ylm_tower_degeneracy_claude.md).
- $F^2$ definition: qed3_v2-6.pdf Eqs. IV.24-35.
