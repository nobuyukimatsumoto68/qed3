# Glueball GEVP upgrades: multi-flow basis, vacuum subtraction, irrep projection

## Goal

Improve the glueball spectroscopy pipeline
(`glue2_claude.cu` -> `F_corr.k`/`F.k` -> GEVP in `glue_ylms3.ipynb`).

This plan covers four of the techniques discussed:

- **Item 1 (vacuum subtraction).** The $l=0$ ($0^{++}$) operator has
  $\langle O\rangle \neq 0$ and mixes with the vacuum. The connected correlator
  must be $C_{ij}(\Delta t)=\langle O_i(t)O_j(t+\Delta t)\rangle-\langle O_i\rangle\langle O_j\rangle$.
  The disconnected subtraction is currently commented out in the notebook
  (`conn = mean_corr # -prod1pt`). Re-enable it, done per jackknife sample.
- **Item 2 (normalization bug).** `glue2_claude.cu:476` reads
  `cdt_avg(i,j) / Comp::Nt;` -- a no-op (result discarded). Must be `/=`.
  Harmless for GEVP eigenvalue ratios, but the time-origin average is never
  normalized, so `F_corr` and `F` (which IS divided by Nt) live on
  inconsistent scales -- which breaks the vacuum subtraction of item 1.
  Fixing the `/=` puts both on the $\langle\cdots\rangle$ scale.
- **Item 7 (multi-smearing / multi-flow variational basis).** Currently one
  Wilson-flow time only. Measure the same $Y_{lm}$ operators at several
  cumulative flow times -> richer variational basis -> earlier/cleaner GEVP
  plateaus. Gradient-flow analog of multi-level APE/HYP smearing.
  Ref: Morningstar & Peardon, PRD 60 (1999) 034509.
- **Item 6 (lattice-symmetry irrep projection).** NEXT STEP, separate copy.
  See assessment at the bottom.

GEVP references (for the notebook side, not implemented here yet):
Luscher & Wolff, Nucl. Phys. B339 (1990) 222;
Blossier, Della Morte, von Hippel, Mendes, Sommer, JHEP 0904 (2009) 094,
arXiv:0902.1265.

## Files

- `glue2_msm_claude.cu` -- copy of `glue2_claude.cu`; implements items 2 + 7.
- `glue_ylms3_claude.ipynb` -- copy of `glue_ylms3.ipynb`; implements item 1
  and consumes the larger (multi-flow) operator set from item 7.
- `glue2_irrep_claude.cu` -- copy of `glue2_msm_claude.cu`; scaffold for item 6.

## Operator layout (items 7 + analysis contract)

Fixed $(\ell,m)$ channel list (matches the notebook `opset`), 16 channels:

```
(0,0)
(1,-1),(1,0),(1,1)
(2,-2),(2,-1),(2,0),(2,1),(2,2)
(3,-3),(3,-2),(3,-1),(3,0),(3,1),(3,2),(3,3)
```

Operators are indexed `op = iflow*n_lm + ilm`, so `nops = N_FLOW * n_lm`.
Default `N_FLOW = 3` cumulative flow checkpoints at flow times
`{1.0, 2.0, 3.0}` (each `FLOW_INCR = 1.0` over `FLOW_NSTEP = 100` integrator
steps; cumulative -- the running `Uflow` is advanced, not restarted).
These three constants are editable knobs at the top of the k-loop.

Output unchanged in format: `F_corr.k` is a `Nt x (nops*nops)` text matrix
(one row per `dt`), `F.k` is the length-`nops` 1-pt vector. With the item-2 fix
both are normalized by `Nt`.

## Implementation chunks

### Chunk A -- `glue2_msm_claude.cu` (items 2 + 7)   [DONE-pending-user-build]
Files: `glue2_msm_claude.cu`
- Replace the ~16 named `flow_plaq_avg_*` vectors + `obs_ptrs` with a
  `lm_set` table and a `std::vector<std::vector<double>> obs(nops, ...)`
  filled by a flow-checkpoint loop (cumulative `flow(Uflow)`).
- Fix `cdt_avg(i,j) /= Comp::Nt;`.
- `F`/`F_corr` writers iterate `obs[op]`.

### Chunk B -- `glue_ylms3_claude.ipynb` (item 1)   [DONE-pending-user-run]
Files: `glue_ylms3_claude.ipynb`
- `nops`/`opset` rebuilt as `N_FLOW` tiles of the 16 `lm` channels;
  `get_lmlm` returns the flow index too.
- Connected correlator: subtract `kron(mean1pt, mean1pt)` built from the
  jackknife-resampled `F1pts` (same deleted slice as `corrs`), inside the
  resampling functions so errors are correct. Only the $l=0\times l=0$ block
  is materially affected; $l\geq1$ 1-pt averages to zero.
- Point the data `label` at the new `data_*` dir produced by `glue2_msm`.

### Chunk C -- `glue2_irrep_claude.cu` (item 6 scaffold)   [assessment only]
Files: `glue2_irrep_claude.cu`
See assessment below; not implemented in this step.

## Item 6 assessment: is it a massive modification?

Short answer: **no, it is moderate, not massive** -- and almost none of the
work is in `glue2`. Reasoning:

The operator is already a spatial weighted sum of plaquette angles,
`O = (1/N_face) sum_face theta(face) w(face)`, with the weight currently
`w(face) = Ylm_real(l, m, r_face)`. An icosahedral-irrep operator is just a
**different weight vector** `w^{Gamma,k}(face)` transforming in irrep `Gamma`.
So the kernel change in `glue2` is tiny: add a
`plaquette_angle_avg_weight(s, w)` that dots the plaquette angles with an
arbitrary precomputed per-face weight table, and feed it irrep weights instead
of (or in addition to) `Ylm_real`.

The real work is producing the irrep weight tables, which is offline group
theory:
1. Realize the icosahedral group on the faces (60 rotations as a permutation
   rep of `lattice.faces`, or as rotations acting on `dual_sites`).
2. Build the projectors $P^\Gamma=(d_\Gamma/|G|)\sum_g \chi^\Gamma(g)^*\,g$
   and extract an orthonormal weight basis per irrep
   ($A,T_1,T_2,G,H$ of $I$; add parity for $I_h$).

Useful simplification for the current $l\leq3$ data: under the icosahedral
group, $l=0\to A$, $l=1\to T_1$, $l=2\to H$ are each a SINGLE irrep, so the
existing $Y_{lm}$ operators are already irrep-pure for $l\leq2$ (no $m$-mixing,
no cross-$l$ mixing within measured $l$). Only $l=3$ reduces reducibly,
$l=3 \to T_2 \oplus G$ (3+4). So with present data, item 6's concrete content
is: split the 7 $l=3$ operators into a 3-dim $T_2$ block and a 4-dim $G$ block
via the Wigner-D action of the 60 rotations on $Y_{3m}$ -- a moderate Python
(or one-off C++) computation, and it can even be done as a **linear
recombination in the notebook** without rerunning `glue2`. Cross-$l$ aliasing
(e.g. $H$ in $l=2,4,5,6$) only matters once higher $l$ are measured.

Recommended scope for item 6 when we get to it:
- (small) notebook recombination of the measured $l=3$ block into $T_2 + G$,
  to remove the only genuine within-$l$ mixing in the current basis; OR
- (moderate) generate full icosahedral-irrep weight tables and add
  `plaquette_angle_avg_weight` to `glue2_irrep_claude.cu`, enabling higher-$l$
  and aliasing-correct channels.

## Open questions

- Flow-radii set for item 7: default cumulative `{1,2,3}` (`FLOW_INCR=1.0`,
  `N_FLOW=3`). Confirm or adjust the radii / count.
- Item 6 route: notebook-side $l=3$ split (quick) vs. full weight-table
  machinery in `glue2_irrep` (general). Decide before starting Chunk C.
