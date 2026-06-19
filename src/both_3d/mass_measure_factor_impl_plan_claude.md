# Mass setting for $L>1$ via the site-area measure factor — write-up + implementation plan

**Status (2026-06-19):** supersedes `mass_scaling_ssb_claude.md` (in the
`project_qed3/` root, now marked OBSOLETE). Pre-implementation write-up for review.
No code written yet. (Lives in `qed3/src/both_3d/` so it is git-tracked with the code.)

## Source (mandatory citation)

P. A. Boyle, R. C. Brower, G. T. Fleming, E. Katz, N. Matsumoto, R. Misra,
*"Studying QED$_3$ with radial quantization on the lattice: Free limit,"*
arXiv:2510.03085 (FERMILAB-PUB-25-0712-T). The free-limit operator there fixes the
mass normalization used here. Key equations:

- **(IV.7)** $C^{(\triangle,t)} = \dfrac{1}{\bar a_s a_t}\,A_\triangle\,a_t\,
  \sigma^a e_a^\mu\overset{\leftrightarrow}{\nabla}^S_\mu\,(1+O(a^2))
  = \dfrac{A_\triangle}{\bar a_s}\,D_\text{cont}$,
  with $D_\text{cont}\equiv\sigma^a e_a^\mu\overset{\leftrightarrow}{\nabla}^S_\mu$
  the massless curved-space Dirac operator (the $a_t$ cancels).
- **(IV.11)** the lattice operator carries the local volume factor $\sqrt{g(x)}$.
- **(IV.12)** generalized eigenproblem $D_\text{lat}\psi=\bar\lambda\,
  (\overline{\delta V})^{-1}\delta V\,\psi$, with $\delta V\equiv\text{diag}(A_y a_t)$
  (per-site dual-cell area $A_y$); divide by the mean to recover the continuum spectrum.

## Why the old scheme was wrong (two reasons)

1. **Varying site measure.** On the refined simplicial $S^2$ the dual-cell areas
   vary site to site. The mass enters weighted by that local measure, so a uniform
   bare $\hat m\,\mathbb{1}$ (current code: `D_ov + m`) does not correspond to a
   uniform physical mass. The old `mean_ell` factor is only a global average.
2. **Order of limits.** SSB is $\Sigma_0=\lim_{m\to0}\lim_{\bar a_s\to0}
   \langle\bar\psi\psi\rangle$ — continuum first at fixed physical $m$, then $m\to0$.
   The schedule must keep a single physical $m$ as $\bar a_s\to0$.

This is a curved-space study (not flat-space QCD); the flat-space uniform-measure /
$Z_m$ intuition does not carry over. The free limit (arXiv:2510.03085) is where $m$
is defined.

## Conventions

- $m$ is the **physical** mass in $R=1$ units (a single number used at **every** $L$);
  it is *not* a conventional $am_q$ bare lattice mass.
- $\bar a_s$, $a_t$ are also in $R$ units. Limits: $\bar a_s\to0$ first, then $m\to0$.

## The prescription

The massive operator should mirror (IV.7): the mass rides the **same** measure factor
$A_\text{tri}/\bar a_s$ that multiplies $D_\text{cont}$,
$$
\frac{A_\text{tri}}{\bar a_s}\big(D_\text{cont}+m\big)
\qquad\Longrightarrow\qquad
\text{add }\ m\,\frac{A_\text{tri}}{\bar a_s}\ \ (\times\ \text{spinor }\mathbb{1}),
$$
i.e. a piece $\propto\mathbb{1}$ (a literal $\propto D_\text{cont}$ piece would merely
rescale the kinetic term, not add a mass). Implemented as a **diagonal matrix** in
site space (uniform at $L=1$ by icosahedral symmetry, site-varying for $L>1$):
$$
m_L \equiv \text{diag}\!\left(m\,\frac{A_\text{tri}}{\bar a_s}\right),
\qquad
D_m \equiv D_\text{ov} + m_L .
$$

### Reverse-engineering the physical $m$ from the $L=1$ runs

At $L=1$ the measure is uniform, and the previous runs implemented $D_\text{ov}+m_1\mathbb{1}$
with bare values $m_1\in\{0.01,0.05,0.10,0.20\}$. Identifying $m_1$ with the uniform
$L=1$ diagonal entry,
$$
m_1 = m\,\frac{A_\text{tri}^{(1)}}{\bar a_s^{(1)}}
\qquad\Longrightarrow\qquad
m = m_1\,\frac{\bar a_s^{(1)}}{A_\text{tri}^{(1)}},
$$
with $A_\text{tri}^{(1)},\bar a_s^{(1)}$ the $L=1$ geometric values. These physical $m$
are then reused unchanged at $L=2,4$, with the $L$- and site-dependence carried entirely
by $A_\text{tri}/\bar a_s$. Consistency: $A_\text{tri}\sim\bar a_s^2\Rightarrow
A_\text{tri}/\bar a_s\sim\bar a_s\to0$, so at fixed physical $m$ the mass term vanishes
in the continuum (correct for the $\bar a_s\to0$-first protocol), and the global average
of $m_L$ reduces to roughly the old `mean_ell` scaling — but now per-site-exact.

## Code grounding (what exists today)

- Mass is a scalar `Complex mass` in `includes/overlap_wmass_claude.h:148`, applied as
  `cplx(mass)` axpy in the operator (`:373` $D_\text{ov}v+mv$, `:405` $D^\dagger_\text{ov}v+m^*v$,
  and the `_ms` variants `:433,:456`).
- **Scalar-mass squared-operator identity** at `:462-503`:
  $(D+m)^\dagger(D+m)v=(1+m^*)(D+m)v+(1+m)(D+m)^\dagger v-(2\,\mathrm{Re}\,m+|m|^2)v$,
  valid because $m$ is scalar and $D^\dagger D=D+D^\dagger$ (GW). It lets the inner
  Zolotarev multishift solve be **shared** between $(D+m)v$ and $(D+m)^\dagger v$ (same
  input $v$) — see the shared-seed solve `:568-586`. Inner-product/force routines
  `:633-692` also assume scalar mass.
- Geometry already exposes the per-site measure: `includes/dirac_ext.h` has
  `lattice.dual_areas[ix]` ($A_y$), `lattice.mean_dual_area`, `lattice.mean_ell`, and a
  `volume_matrix(...)` helper (`:350`) that builds `diag(dual_areas[ix]/mean_dual_area)`
  — i.e. the $\delta V/\overline{\delta V}$ of (IV.12).

## Files to modify

- `includes/overlap_wmass_claude.h` — replace scalar-mass axpy with a per-site diagonal
  multiply by a device array `d_mL[ix]` in `mult/adj` (`:373,:405`) and the `_ms` variants
  (`:433,:456`); rework the squared-operator path (`:462-503`) and the mass-dependent
  inner products (`:633-692`); add storage + upload for `m_L`.
- `includes/dirac_ext.h` (or the lattice/base geometry) — provide the per-site measure
  array used to build `m_L` (`dual_areas`, `mean_ell`, `mean_dual_area` already here).
- `hmc_fermilab_wmass_L{2,4}_claude.cu` — `Fermion D(DW, mass, npole)` now takes the
  **physical** $m$ and builds the diagonal `m_L` from lattice geometry; the CLI
  `mass_re/mass_im` becomes physical $m$. (Mirror in any measurement/`glue*` code that
  builds $D_\text{ov}+m$.)

## Ordered implementation chunks

1. **Pin numbers / decisions** (resolve open questions below): exact per-site measure,
   exact normalization, $L=1$ constants → physical $m$ values. *Files: this doc.*
2. **Geometry plumbing**: build `m_L` (host) from `dual_areas`/`mean_ell`, upload device
   array. *Files: `dirac_ext.h`, `overlap_wmass_claude.h`.*
3. **Operator apply**: diagonal multiply in `mult/adj` + `_ms`. *Files: `overlap_wmass_claude.h`.*
4. **Normal operator**: rework `:462-503` for diagonal $M$ (see Q4) + mass inner products
   `:633-692`. *Files: `overlap_wmass_claude.h`.*
5. **Driver/CLI**: pass physical $m$, construct `m_L`. *Files: `hmc_fermilab_wmass_L{2,4}_claude.cu`.*
6. **Validation**: $L=1$ must reproduce the previous run bit-for-physics (uniform
   $m_L=m_1$); spot-check free-limit spectrum / a known observable at $L=2$.

## Open questions (resolve before coding)

1. **Per-site measure $A_\text{tri}$ = dual-cell area $A_y$?** The natural per-site
   diagonal measure is `lattice.dual_areas[ix]` ($A_y$), matching $\delta V$ in (IV.12)
   and the existing `volume_matrix`. Confirm this is what "$A_\text{tri}$" means (vs. a
   triangle *face* area $A_\triangle=4\pi/20$ at $L=1$, which differs from the dual area
   $A_y=4\pi/12$).
2. **Exact normalization of the diagonal entry.** Is it literally
   $m\cdot A_y/\bar a_s$ (with $\bar a_s=$ `mean_ell`), or the dimensionless
   $m\cdot A_y/\overline{A}$ (`dual_areas/mean_dual_area`, as `volume_matrix` builds),
   or $m\,(A_y/\overline A)(\overline A/\bar a_s)$? Pin the precise constant — it sets
   the overall mass scale.
3. **$L=1$ constants → physical $m$.** Provide/confirm $A_\text{tri}^{(1)}$ and
   $\bar a_s^{(1)}$ (and whether $A_\text{tri}^{(1)}$ is the dual area or face area), so
   the four physical $m=m_1\,\bar a_s^{(1)}/A_\text{tri}^{(1)}$ are fixed.
4. **Overlap normalization.** $D_\text{ov}$ is GW-normalized (eigenvalues on the unit
   circle), not in the $\frac{A_\triangle}{\bar a_s}D_\text{cont}$ normalization of (IV.7).
   Confirm that adding $m_L$ directly to $D_\text{ov}$ is correct (vs. needing a
   compensating factor so $m_L$ matches the low-mode normalization of $D_\text{ov}$).
5. **Solver cost (important).** Diagonal $M$ breaks the scalar-$m$ identity at `:462`:
   $(D+M)^\dagger(D+M)=(D+D^\dagger)+MD+D^\dagger M+M^2$, and the $D^\dagger M$ term needs
   $D^\dagger$ applied to $Mv$ (a different vector), so the inner Zolotarev solve can no
   longer be shared between the two applies — roughly $2\times$ the fermion-operator cost
   in the CG, hitting L4 hardest. Decide: (a) accept the cost via composition, or
   (b) reformulate the pseudofermion action to avoid the $D^\dagger M$ cross term.
6. **Force unchanged?** $m_L$ is gauge-independent, so $\partial S/\partial U$ from the
   mass term is zero; the HMC force only changes through the modified solves. Confirm no
   explicit mass-derivative term is needed.
