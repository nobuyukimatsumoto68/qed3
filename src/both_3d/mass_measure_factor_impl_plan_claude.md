# Mass setting for $L>1$ via the site-area measure factor — write-up + implementation plan

**Status (2026-06-19):** supersedes `mass_scaling_ssb_claude.md` (in the
`project_qed3/` root, now marked OBSOLETE). (Lives in `qed3/src/both_3d/` so it is
git-tracked with the code.)

**IMPLEMENTATION STATE (2026-06-19, local Claude) -- see the chunk list below for detail:**
- DONE: operator path (`includes/overlap_wmass_claude.h`) -- `mult`/`adj`/`DHD`(+`_ms`)
  and ALL HMC force variants (`grad` + `grad_l1`/`_l2`/`_l4`) fold the diagonal `m_L`.
- DONE: `BlockedMat` (`includes/blocked_mat_claude.h`) `mult`/`adj`/`DDH` -- the parallel
  mrhs operator used by JJ/condensate measurement; was the one stale scalar-mass path
  (audit `mass_measure_audit_handoff_claude.md`), now diagonal.
- DONE + VERIFIED (2026-06-19 rerun): `test_diag_mass_l1_claude.cu` (operator obsolete-vs-production +
  force-vs-FD + `BlockedMat`-vs-MatPoly), nontrivial gaussian gauge; runner `tmp_claude.sh` (10 phases:
  all 4 grad variants default/l1/l2/l4 at L=1 & L=2, default+l4 at L=4). ALL PASS, 0 errors.
  BlockedMat(NSTACK=1) vs MatPoly `_ms` = 0.0 (bit-identical). L>1 massive force-vs-FD ~1e-4..1e-5
  (solver-tol/eps limited; the L=1 obsolete-vs-production force = 1e-16 is the machine-precision proof).
- TODO: driver/CLI (`hmc_fermilab_wmass_L{2,4}_claude.cu`) physical-`m` plumbing; the audit's
  extra numerical checks (adjointness L=2, Ward identity); a separate L>2 massive HMC test `.cu`
  (real + imag m). Action checklist: `mass_diag_l1_task_claude.md`.

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

The massive operator mirrors (IV.7): the mass rides the **same** measure factor that
multiplies $D_\text{cont}$. Note IV.7's prefactor symbol $A_\triangle$ is the *triangle
(dual) area* and is **not** the per-site dual-cell area $A_y$; the quantity that belongs
in the per-site **action** mass term is $A_y$ (the dual-cell area). So the added piece is
$$
m\,\frac{A_y}{\bar a_s}\ \ (\times\ \text{spinor }\mathbb{1})
$$
— $\propto\mathbb{1}$ (a literal $\propto D_\text{cont}$ piece would merely rescale the
kinetic term, not add a mass). Implemented as a **diagonal matrix** in site space
(uniform at $L=1$ by icosahedral symmetry, site-varying for $L>1$):
$$
m_L \equiv \text{diag}\!\left(m\,\frac{A_y}{\bar a_s}\right),
\qquad
D_m \equiv D_\text{ov} + m_L ,
$$
with $A_y=$ `lattice.dual_areas[ix]` and $\bar a_s=$ `lattice.mean_ell` (literally
$m\,A_y/\bar a_s$, **not** normalized by $\overline A$). This rides the same $A_y/\bar a_s$
measure as the kinetic hoppings (IV.2): spatial $\kappa_{y_1y_2}=2A_{y_1y_2}/(\bar a_s\ell)$
(`dirac_simp.h:360`) and temporal $\kappa'_y=A_y/(\bar a_s a_t)$ (`dirac_ext.h:452`).

### Reverse-engineering the physical $m$ from the $L=1$ runs

At $L=1$ the measure is uniform, and the previous runs implemented $D_\text{ov}+m_1\mathbb{1}$
with bare values $m_1\in\{0.01,0.05,0.10,0.20\}$. Identifying $m_1$ with the uniform
$L=1$ diagonal entry,
$$
m_1 = m\,\frac{A_y^{(1)}}{\bar a_s^{(1)}}
\qquad\Longrightarrow\qquad
m = m_1\,\frac{\bar a_s^{(1)}}{A_y^{(1)}},
$$
with $A_y^{(1)},\bar a_s^{(1)}$ the $L=1$ geometric values (uniform over sites by
icosahedral symmetry). These physical $m$ are then reused unchanged at $L=2,4$, with the
$L$- and site-dependence carried entirely by $A_y/\bar a_s$. Consistency:
$A_y\sim\bar a_s^2\Rightarrow A_y/\bar a_s\sim\bar a_s\to0$, so at fixed physical $m$ the mass term vanishes
in the continuum (correct for the $\bar a_s\to0$-first protocol), and the global average
of $m_L$ reduces to roughly the old `mean_ell` scaling — but now per-site-exact.

### Geometry constants per $L$ (from `test_diag_mass_l1_claude.cu` startup print; $R=1$ units)

Saved 2026-06-19 (local Claude) from the test's `# L=.. check: ... mean_dual_area=.. mean_ell=..`
line. $\bar A=$ `mean_dual_area` $=4\pi/n_\text{sites}$ exactly (dual cells tile the sphere);
$\bar a_s=$ `mean_ell` is triangulation-dependent (NOT analytic — read from the run).

| $L$ | $n_\text{sites}$ | $\bar A$ = mean\_dual\_area = $4\pi/n_\text{sites}$ | $\bar a_s$ = mean\_ell | $\bar A/\bar a_s$ | $\bar a_s/\bar A$ |
|---|---|---|---|---|---|
| 1 | 12  | 1.047197551196597    | 1.107148717794090   | 0.945850845840348 | 1.0572491470487 |
| 2 | 42  | 0.2991993003418845   | 0.5909464448075018  | 0.506305271773566 | 1.97509300366762 |
| 4 | 162 | 0.0775701889775253   | 0.2994744726728923  | 0.259021038705536 | 3.860690255114117 |
| 8 | 642 | 0.0195737860036747   | needs a build (not in the test) | — | — |

At $L=1$ the dual cell is uniform by icosahedral symmetry, so the per-site
$A_y^{(1)}=\bar A^{(1)}=1.047197551196597$.

### CANDIDATE physical $m$ (reverse-engineered from the old $L=1$ bare $m_1$) — REMOTE AGENT TO VERIFY

To reproduce the old $L=1$ physics, choose physical $m$ so the (uniform) $L=1$
$m_L=\texttt{mass\_coeff}=m\,(\bar A/\bar a_s)^{(1)}$ equals the old bare $m_1$:
$$
m = m_1\,\frac{\bar a_s^{(1)}}{\bar A^{(1)}} = m_1 \times 1.0572491470487 .
$$

| old bare $m_1$ | candidate physical $m$ |
|---|---|
| 0.01 | 0.010572491470487 |
| 0.05 | 0.052862457352435 |
| 0.10 | 0.105724914704870 |
| 0.20 | 0.211449829409740 |

**These are CANDIDATES computed locally from the $L=1$ geometry above. The REMOTE AGENT must
verify**: (i) recompute $\bar a_s^{(1)}/\bar A^{(1)}$, (ii) confirm the direction
$m=m_1\,\bar a_s/\bar A$ (NOT $\bar A/\bar a_s$), and (iii) decide whether to reproduce the old
$L=1$ ensembles or define $m$ fresh. The CLI `mass_re/mass_im` is now this physical $m$; `dir3`
encodes it, so these land in NEW checkpoint dirs (see the driver/CLI TODO).

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

1. **[DONE] Pin numbers / decisions**: see "Resolved decisions" below (per-site measure $A_y$,
   normalization $m A_y/\bar a_s$, $L=1$ reverse-engineering, no extra solver cost, force unchanged).
2. **[DONE] Geometry plumbing**: `m_L = mass_coeff * M_mass`, `M_mass = volume_matrix(1)`,
   `mass_coeff = m * mean_dual_area/mean_ell`, built in the `OverlapWMass` ctor. *Files: `overlap_wmass_claude.h`.*
3. **[DONE] Operator apply**: `apply_mL`/`apply_mLdag`; diagonal add in `mult`/`adj` + `_ms`. *Files: `overlap_wmass_claude.h`.*
4. **[DONE] Normal operator + force**: `DHD`(+`_ms`) diagonal identity $(1{+}M^*)D{+}D^\dagger(1{+}M){+}|M|^2$;
   HMC force `grad`+`grad_l1/l2/l4` fold $(1{+}M^*)$ via `M_mass`. ALSO the mrhs `BlockedMat::mult/adj/DDH`
   (`blocked_mat_claude.h`) -- the JJ/condensate measurement operator (audit-found stale scalar). *Files: `overlap_wmass_claude.h`, `blocked_mat_claude.h`.*
5. **[TODO] Driver/CLI**: pass physical $m$ via `mass_re/mass_im`. *Files: `hmc_fermilab_wmass_L{2,4}_claude.cu`.*
6. **[PARTIAL] Validation**: $L=1$ operator obsolete-vs-production + force-vs-FD + BlockedMat-vs-MatPoly
   all PASS in `test_diag_mass_l1_claude.cu` (default grad + `-DGRAD_L4`). PENDING: a separate L>2 massive
   HMC test (real + imag m), the audit's adjointness/Ward-identity checks, and a free-limit spectrum spot-check.

## Resolved decisions (2026-06-19, NM)

1. **Per-site measure = $A_y$ (dual-cell area) = `lattice.dual_areas[ix]`.** IV.7's
   $A_\triangle$ (triangle/dual area) is a *different* geometric object from $A_y$, but
   the quantity that belongs in the per-site **action** mass term is $A_y$ — so the code
   strategy (use `dual_areas`) is correct.
2. **Normalization: literally $m\,A_y/\bar a_s$** with $\bar a_s=$ `mean_ell`. NOT the
   dimensionless $A_y/\overline A$ form.
3. **$L=1$ reverse-engineering confirmed:** $m=m_1\,\bar a_s^{(1)}/A_y^{(1)}$, with
   $A_y^{(1)},\bar a_s^{(1)}$ the (uniform) $L=1$ dual area and `mean_ell`.
4. **Overlap normalization OK** — adding $m_L$ directly to $D_\text{ov}$ is on the same
   $A_y/\bar a_s$ footing as the kinetic hoppings (IV.2), verified against
   `dirac_simp.h:360` ($\kappa=2A_{y_1y_2}/(\bar a_s\ell)$) and `dirac_ext.h:452`
   ($\kappa'_y=A_y/(\bar a_s a_t)$). No compensating factor needed.
5. **Solver: NO extra cost** (corrects an earlier mistaken "2x"). The squared operator
   keeps a linear-combination form (analog of the scalar `:462` identity), using GW
   $D^\dagger D=D+D^\dagger$ and $M=M^\dagger$ (real diagonal):
   $$D_m^\dagger D_m=(D+M)^\dagger(D+M)=(1+M)\,D + D^\dagger(1+M) + M^2.$$
   Computed as $(1+M)(D_\text{ov}v) + D_\text{ov}^\dagger\big[(1+M)v\big] + M^2 v$ — ONE
   forward apply (on $v$) + ONE dagger apply (on $(1+M)v$) + pointwise $M$ ops. Same
   overlap-apply count as the scalar case ⇒ **no doubling of CG cost**. (The $D^\dagger M$
   cross term is absorbed by feeding $(1+M)v$ to the single $D^\dagger$ apply; an overlap
   apply costs the same for any input vector.) Implementation: generalize `:462-503` from
   scalar `mass` to this diagonal form (`mult`/`adj` give $D_\text{ov}v$, $D_\text{ov}^\dagger u$;
   the $(1+M)\cdot$ and $M^2\cdot$ are pointwise).
6. **Force unchanged** — $m_L$ is gauge-independent, $\partial S/\partial U$ from the mass
   term is zero; the HMC force changes only through the modified solves. No explicit
   mass-derivative term.

**All open questions resolved → ready to implement (chunks above).**

## Implementation knowledge & validation strategy (2026-06-19)

**Geometry access.** `OverlapWMass` holds `const WilsonDirac& DW`; `DW` is a `DiracExt`
with `Base& lattice` (`dirac_ext.h:8`). Measure data: `lattice.dual_areas[ix]` ($A_y$),
`lattice.mean_dual_area` ($\bar A$), `lattice.mean_ell` ($\bar a_s$) — declared in
`s2n_simp.h:18,38,39` (also `s2n_dual.h`).

**Build $m_L$ via `volume_matrix` (mimic the GEVP code).** `DW.volume_matrix(elem, pow)`
(`dirac_ext.h:350`) fills a COO with $\text{diag}\big((A_y/\bar A)^{\text{pow}}\big)$
over all sites (NS=2 spin comps) $\times N_t$ slices. GEVP usage
(`saved_scripts/eig.cu:263-265`):
```
COO<N> gmfourth;  DW.volume_matrix(gmfourth.en, -0.5);  gmfourth.do_it();
// used as a matrix in MatPoly: Op.push_back(coeff, {&gmfourth, ...}); Op.on_gpu<N>(out,in);
```
For the mass use **pow = 1** ($\to A_y/\bar A$) with scalar coefficient
$$
\texttt{mass\_coeff}=m_\text{phys}\,\frac{\bar A}{\bar a_s}
= m_\text{phys}\,\frac{\texttt{mean\_dual\_area}}{\texttt{mean\_ell}},
$$
so $\texttt{mass\_coeff}\times\texttt{volume\_matrix}(1)=\text{diag}(m\,A_y/\bar a_s)=m_L$.
At $L=1$, `volume_matrix(1)` is the identity, so $m_L=\texttt{mass\_coeff}\cdot\mathbb{1}$.

**COO/CSR/MatPoly** live in `includes/matpoly.h` (`push_back` :81, `on_gpu` :106,
`from_cpu` :90) and `includes/sparse_matrix.h` (`COO`/`CSR`, `do_it`).

**Edit sites in `overlap_wmass_claude.h`** (scalar `mass` $\to$ diagonal $m_L$):
- ctor `:173-192` — reinterpret `mass` as physical $m$ (real); build member
  `COO<N> M_mass` from `volume_matrix(.,1)`+`do_it()`, store `mass_coeff`, alloc scratch.
- `mult_deviceAsyncLaunch` `:373` — replace `Taxpy_gen(d_res, cplx(mass), d_xi)` with
  $+\,m_L v$ (apply `M_mass` with coeff `mass_coeff` to `d_xi` via MatPoly $\to$ scratch, add).
- `adj_deviceAsyncLaunch` `:405` — same ($M$ real $\Rightarrow m_L^\dagger=m_L$).
- squared op `:462-503` — $D_m^\dagger D_m=(1+M)D+D^\dagger(1+M)+M^2$ (decision 5):
  $(1+M)(D_\text{ov}v)+D_\text{ov}^\dagger[(1+M)v]+M^2v$.
- mass inner products `:633-692` and `_ms` variants `:494-503` — same substitution.

**Validation (L=1 reduction to scalar).**
1. Freeze reference: `cp overlap_wmass_claude.h overlap_wmass_obsolete_claude.h`. To let
   one test TU include BOTH, wrap the obsolete copy's file-scope structs (`Zolotarev` :7,
   the rational struct ~:117, `OverlapWMass`) in `namespace obsolete { ... }` placed
   AFTER the `#include`s (after line 4) to EOF (`#pragma once` is per-path, fine) $\to$
   `obsolete::OverlapWMass` (scalar) vs production `OverlapWMass` (diagonal).
2. Edit production header per chunks.
3. `test_diag_mass_l1_claude.cu` at $L=1$: build `Base`/`DW` as in the driver; set
   $c=m\cdot\texttt{mean\_dual\_area}/\texttt{mean\_ell}$; construct `obsolete::OverlapWMass`
   with scalar `mass=c` and production `OverlapWMass` with physical $m$; apply `mult`
   (then `adj`, then squared op) to the **same** random vector; compare
   $\|\Delta\|/\|\cdot\|\sim10^{-14}$. At $L=1$ both equal $D_\text{ov}+c$.

Action checklist for the local run: see `mass_diag_l1_task_claude.md`.
