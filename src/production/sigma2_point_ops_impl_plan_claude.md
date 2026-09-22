# Point (unaveraged) sigma^2 operators for the two-meson GEVP -- impl plan

## Goal / physics
The two-meson mode is the NOISY one in the L2 tower (ground ~0.35 rock-solid, m1 ~0.7 marginal). Our whole
basis -- $\sigma^2_{00}$ (double site-sum), $O_{2m}$ (antipodal sum), $O_{1m}$ (coincident sum) -- is l=0
SPATIALLY-SUMMED, so it can only reach the two mesons at zero relative momentum / a few relative separations.
POINT operators (two $\sigma$'s pinned to FIXED sites, no $\sum_x$) add the RELATIVE-POSITION variational
directions the summed basis can't reach, which is exactly the two-meson's internal structure. Idea + method
from NM & the stress-tensor agent (12 degree-5-vertex point ops in her O_H GEVP).

## Mechanism (verified against fs_gevp_point_claude.py)
A point-$\sigma^2$ operator = `op_vspec` returning the two bilinear vertices weighted by ONE-HOT (delta)
vectors at the fixed sites $x_1,x_2$ instead of `dual*Y00`. `perm_contrib_folded` appends the weight vector on
each vertex's site-letter and the einsum SUMS the letter; a one-hot weight collapses that sum to the fixed
site. So `make_config`/`AblkS`/`perm_contrib_folded` are UNCHANGED -- only the weight vectors:
```
def op_vspec_point(x1, x2, letters, nsite):     # point sigma^2 at fixed sites (x1,x2), unit weight
    l0, l1 = letters
    e1 = np.zeros(nsite); e1[x1] = 1.0
    e2 = np.zeros(nsite); e2[x2] = 1.0
    return [(l0, e1, False), (l1, e2, False)]
```
Antipode handling stays via the (False) flag; for an antipodal point op just set x2 = P(x1).

## Sites (from the stress-tensor agent)
The 12 degree-5 (icosahedral) vertices `five=[i for i in range(nsite) if len(nns[i])==5]` (L1: all 12; L2:
12 of 42). Symmetry-equivalent -> no privileged site. Two ways to use them:
- (a) TRUE point ops: individual fixed pairs (noisy; the GEVP recovers signal -- her warning #1).
- (b) ORBIT-AVERAGED "shell" ops: sum a fixed relative geometry over the 12-vertex orbit -> orientation-
  averaged fixed-separation operator with good S/N (generalizes $O_{2m}/O_{1m}$ to INTERMEDIATE separations).
Relative geometries to span the wavefunction: coincident (r=0), nearest-neighbor, a mid shell, antipodal
(r=max). $O_{2m}$/$O_{1m}$ are the endpoints already; the NEW value is the intermediate shells.

## Warnings (from the stress-tensor agent -- MUST heed)
- Individual point ops fan into noise at large t -> rely on the GEVP and/or the orbit average.
- NORMALIZATION: none needed -- the GEVP whitens per-operator scale; just symmetrize $0.5(C+C^T)$.
- HANKEL: do NOT add a Hankel offset ladder to a multi-point basis (over-determines, rank<=Nv, rebase errors
  blow up 5-10x). Use PLAIN off=[0] and a SMALL NKEEP. (Our two-meson production uses Dt=[0,2,4] -- must drop
  it for the point-enriched basis.)
- CONTACT: coincident sites (x1=x2) carry the -1/2 I self-contraction -- same care as $O_{1m}$.

## Chunks
### Chunk 1 -- fixed-site vspec + validation (free/L1)
Add `op_vspec_point` to `fs_gevp_point_claude.py` (new function; op_vspec untouched). Validate on the free
peram: (a) a point op's correlator is finite/real; (b) a point op summed over ALL site pairs with dual*Y00
weights REPRODUCES the corresponding summed op (sanity that the one-hot path == the weighted path); (c)
`five_vertices(nns)` picks 12 sites. Files: EDIT `fs_gevp_point_claude.py`; NEW `sigma2_point_test_claude.py`.

### Chunk 2 -- point-enriched GEVP on L1 (clean, decisive test)
Build a cache with the summed P+ ops PLUS a handful of point/shell ops (coincident, nn, mid, antipodal at the
12 vertices; orbit-averaged first). GEVP with PLAIN off=[0], NKEEP=2-3, symmetrized. Question: does adding
the relative-position ops sharpen the two-meson (0.62 at L1) or add a resolved excited two-meson that the
summed basis misses? L1 is the clean control. Files: NEW `sigma2_point_gevp_claude.py`.

### Chunk 3 -- L2 (the target)
Same on the L2 Nf2 g1 (400 cfg, the noisy tower). Does the point/shell enrichment recover the two-meson m1?
Compare to the summed-only reb3 result (m1 ~0.7 marginal). Files: same, ENS/LREF switch.

## Open questions for NM (resolve before/at chunk 2)
1. TRUE point ops (no average -- your explicit ask, matches the stress-tensor agent) vs ORBIT-AVERAGED shells
   (better S/N)? Recommend: test orbit-averaged shells FIRST (cleaner signal to see if relative-position helps
   at all), then true point ops. Which do you want prioritized?
2. Which relative separations/shells (how many intermediate)?
3. Confirm dropping Hankel (plain off=[0]) for the point basis per the warning -- OK to deviate from the
   Dt=[0,2,4] production for this sub-study?

Refs: stress-tensor agent's `t00_ham_pole_gevp_claude.py` (mask_star, five_vertices, corr_matrix, gevp);
distillation Peardon 0905.2160; two-meson machinery `fs_gevp_point_claude.py`,
`sigma2_flavorgeom_full_v2_claude.py`; ps_fs_flavor_cross_impl_plan_claude.md.
