# Impl plan — displaced (link) sp current "disp" in jj_local_deter

Add a third current type **disp** (displaced / link current), parallel to **loc** (site $\sigma^a$,
already done) and **exact** (full overlap $K$). **sp-only.** Deterministic trace, continuum (and
lattice) propagator, $L=1,2$. New comparison notebook `comp_loc_disp_jj_claude.ipynb` (loc vs disp).

Source: **qed3int_v2-12.pdf** (LaTeX `qed3int_v2(5)/main.tex`), Eqs. (3.33)-(3.35), (3.45)-(3.46),
and the spatially-projected estimator `eq:spatially_projected_JJ`.

## Equations

The Wilson point-split link current kernel (Eq. 3.33-3.35), two-component:
$$
{\cal W}^{wz}_{xy}=\frac{1}{\lambda_M}\big(\delta_{xw}\delta_{yz}{\cal C}_{wz}-\delta_{xz}\delta_{yw}{\cal C}_{zw}\big),
\qquad {\cal C}_{xy}=-\kappa_{xy}\,{\cal P}_{xy}\,U_{xy},
\qquad {\cal P}_{xy}=\tfrac12(1-e^a_{xy}(x)\sigma_a)\,\Omega_{xy}.
$$

**Displaced current (this task).** The genuine link current with the per-link coupling $\kappa$ divided out:
$$
j^{\rm disp}(x,y)=\bar\psi(x)\,\frac{{\cal W}^{wz}_{xy}}{\kappa_{xy}}\,\psi(y)
\quad\Longleftrightarrow\quad
\frac{{\cal C}_{xy}}{\kappa_{xy}}=-{\cal P}_{xy}U_{xy}=-\tfrac12(1-e^a_{xy}\sigma_a)\Omega_{xy}U_{xy}.
$$
This is the **older / un-reduced form of the loc current**: it keeps the full Wilson Dirac structure
$(-r\sigma_0+e^a\sigma_a)$, the spin connection $\Omega$, and the gauge link $U$ on the **link** (legs
displaced to two distinct sites $x\neq y$) -- it is NOT collapsed to the on-site $\sigma^a$ of loc.
Per Eq. (3.46), $J^{wz}_V\simeq\sqrt{C_{J_V}}\kappa^{(0)}_{wz}e^a_{wz}j^a_V$; dividing the lattice
operator by $\kappa$ removes the spatially-dependent $\kappa^{(0)}_{wz}$ coupling.

**Implementation shortcut (exact, faithful).** `DiracExt::d_coo_format(coo, U, {t,lk})` already builds
the spatial-link entries $=i\,{\cal C}_{xy}=i\lambda_M\,{\cal W}^{wz}_{xy}$ (the same `0.5*kappa*(-r*sigma0
+gamma)*i*exp(i u_sp)*Omega` used by the Dirac operator and `build_W`). So
$$
\texttt{disp entries} = \texttt{d\_coo\_format}/\kappa_{il}
= i\lambda_M\,({\cal W}^{wz}/\kappa)
$$
i.e. the requested $W/\kappa$ up to the overall constant $i\lambda_M$ (absorbed into $C_J$; irrelevant
for shape/ratio). Reusing `d_coo_format` guarantees the exact $W$ (sign/$\Omega$/$U$/$r$) without
re-deriving by hand.

## Estimator + weight (Eq. `spatially_projected_JJ`)

$$
\sum_{(n,n')\in\triangle}\frac{\ell_{nn'}\ell^*_{nn'}}{\kappa^{(0)2}_{nn'}}\,J_V^{nn'}(t_0)J_V^{nn'}(t)
\simeq C_J(\delta^{ab}-e^a_3 e^b_3)\,j_V^a j_V^b\,A_\triangle .
$$
With $j^{\rm disp}=J_V/\kappa$ the $1/\kappa^2$ is already in the operator, so the per-link weight is
$\ell_{nn'}\ell^*_{nn'}=A_{nn'}=$ `link_volume[il]` (sums to $4\pi$; the dual link area, both triangles).
Insertion-diagonal over spatial links:
$$
C^{\rm disp}(t_0\!\to\!t)=\sum_{il} \texttt{link\_volume}[il]\,\mathrm{tr}\big[W_d(il,t_0)\,P\,W_d(il,t)\,P\big],
\quad W_d=\texttt{d\_coo\_format}/\kappa .
$$
`write_corr` folds the global $1/(4\pi)$ (same as the s1/s2/s3 channels). (Consistency check: the OLD
loc-sp used `w_sp=link_volume/kappa^2` with a $\kappa\gamma$ current; disp moves the $1/\kappa^2$ into
the operator and uses `link_volume`, same net normalization.)

## Files / chunks

**Chunk 1 -- SEPARATE program `jj_disp_deter_claude.cu` (DONE).** (User: separate .cu, not folded into
jj_local_deter.) Mirrors jj_local_deter scaffolding (includes, CLI, P-load `Dm_inv`, atomic h5), but
computes ONLY the displaced sp link current:
- `build_W_disp(en, DW, U, t, lk)`: `DW.d_coo_format(coo,U,{t,lk})` then `/DW.bd.kappa[il]`,
  `COOEntry{v(CuC),i,j}` -> `Ent{i,j,Complex(c.v.x,c.v.y)}`.
- main: link loop, `w=link_volume[il]`, insertion-diagonal `tr_WPWP`, output
  `data_<ESNID>/corr_deter_disp[_<tag>]_L<L>/corr.<k>.h5`, key `h0/t0_b/sp/{Vpp,Vmm}` + `h0/disc/sp/J`.
  Distinct dir from loc (`corr_deter_local*`) -- caches separately.

**Chunk 2 -- run `tmp_disp_claude.sh`** (lattice L=1,2: compile `jj_disp_deter_L<L>.o`, run `--n-t0 2`
reading `prop_deter_L<L>`). No stale `rm` needed (new program, fresh dirs).

**Chunk 3 -- notebook `comp_loc_disp_jj_claude.ipynb`** (lattice, `tag=''`).

**Chunk 2 -- notebook `comp_loc_disp_jj_claude.ipynb`.**
- Load loc $G_s=$ `s1`+`s2` and `disp` from `corr_deter_local_<tag>_L{1,2}`.
- Plot both vs the Eq.(4.28) $G_s$ shape (each amplitude-matched), and the ratio disp/loc and each/CFT,
  for $L=1,2$ -- they share the continuum limit, differing by lattice (operator) artifacts.

## Why insertion-diagonal over links reproduces Eq. (4.28)  [flagged]

The estimator `eq:spatially_projected_JJ` is itself a sum of **same-link** products
$J_V^{nn'}(t_0)\,J_V^{nn'}(t)$ (identical superscript $nn'$ on both times) over the edges of a triangle:
$$
\sum_{(n,n')\in\triangle}\frac{\ell_{nn'}\ell^*_{nn'}}{\kappa^{(0)2}_{nn'}}J_V^{nn'}(t_0)J_V^{nn'}(t)
\simeq C_J(\delta^{ab}-e^a_3 e^b_3)\,j_V^a(x_\triangle,t_0)j_V^b(x_\triangle,t)\,A_\triangle .
$$
No cross-link ($nn'\neq mm'$) terms appear. The reason the **diagonal** edge sum still builds the full
isotropic projector is the triangle identity `eq:triangle_formula`
$\sum_{(n,n')\in\triangle}\ell_{nn'}\ell^*_{nn'}e^a_{nn'}e^b_{nn'}=A_\triangle(\delta^{ab}-e^a_3e^b_3)$:
with $J_V^{nn'}\simeq\sqrt{C_{J_V}}\kappa^{(0)}_{nn'}e^a_{nn'}j^a$ (Eq. 3.46), the per-edge diagonal terms
$\propto e^a_{nn'}e^b_{nn'}$ sum to $(\delta^{ab}-e_3e_3)$. Summing over all triangles (equivalently over
all links once, with $A_{nn'}=$`link_volume` $=\ell\ell^*$ covering both triangles) gives
$\int(\delta-e_3e_3)\langle j^a j^b\rangle = G_s$ = **Eq. (4.28)**.

So my contraction is exactly that diagonal-over-links sum: $C^{\rm disp}(t)=\sum_{il}$`link_volume`$[il]\,
\mathrm{tr}[W_d(il,t_0)P\,W_d(il,t)P]$, with $W_d=W/\kappa$ (the $1/\kappa^2$ folded in). The displaced
$W_d$ has entries on **both** endpoints of the link, so each `tr_WPWP` already pulls in the cross-site
propagators $G(n,n')$ — the diagonal is in the **link index**, not in the sites. I believe this
reproduces Eq. (4.28) up to lattice artifacts; the comparison notebook checks it numerically.

## Addendum (2026-06-11) -- TEMPORAL displaced current (tp channel)

Goal: give `disp` its OWN temporal current so the $G_s/G_t$ ratio can be formed within disp (it should
$\to-(D-1)=-2$ as $L\to\infty$, with the lattice $Z$ canceling between disp's own sp and tp). Until now
disp had sp only, so the ratio had to borrow loc's tp (mixed normalization, ratio not $\to-2$).

**Kernel.** The Wilson TEMPORAL point-split current, exact analog of the spatial one. The temporal
hopping `DiracExt::d_coo_format(coo, U, {t, ix})` (the `std::pair<int,Idx>` overload, dirac_ext.h:406)
builds $i\,{\cal C}^t$ with every entry scaled by $\kappa_t[ix]$ and the Dirac $i$:
$$
{\cal C}^t_{(t,ix),(t+1,ix)} = -\tfrac12\,\mathrm{sgn}\,\kappa_t[ix]\,(-r\sigma_0 - \sigma_3)\,e^{-i\,u_{\rm tp}(t,ix)} ,
$$
($\mathrm{sgn}=-1$ at $t=N_t-1$ for antiperiodicity; the backward block is its partner). The displaced
temporal current divides out $\kappa_t$ (Eq. 3.46, temporal direction), exactly like sp divides $\kappa$:
$$
W_d^t = {\cal C}^t/\kappa_t \;=\; \texttt{d\_coo\_format}\,/\,(i\,\kappa_t[ix]) .
$$
So `build_W_disp_t(en, DW, U, t, ix)` = `DW.d_coo_format(coo, U, {t,ix})` then $\times\,1/(i\,\kappa_t[ix])$.
$W_d^t$ touches the two timeslices $t,t+1$ at the SAME spatial site $ix$ (cf. sp: two sites, one slice).

**Weight (full-sum mode).** Temporal current is a sum over SITES, so the per-insertion weight is the dual
SITE area `dual_areas[ix]` (sums to $4\pi$) -- same convention as loc's `s3` tp weight `w_site`. (sp uses
`link_volume`, the dual LINK area.) Insertion-diagonal:
$$
C^{\rm disp,tp}(t_0\!\to\!t)=\sum_{ix}\texttt{dual\_areas}[ix]\,\mathrm{tr}\big[W_d^t(ix,t_0)\,P\,W_d^t(ix,t)\,P\big].
$$

**Single-insertion mode (`--ins i`).** site $i$, weight $1$ (raw), matching exact1 tp (site $i$). Output
key `h0/t0_b/tp/Vpp` (+ `h0/disc/tp/J`) alongside the existing `sp`. disp1 then supplies BOTH channels,
so `comp_threesome_jj_1` can plot disp1's own $G_s/G_t\to-2$ and disp1$_{tp}$ vs exact1$_{tp}$ at site $i$.

**Note on the ratio.** sp divides $\kappa$ (spatial link coupling), tp divides $\kappa_t$ (temporal); these
differ, but that is the principled bare-current normalization (Eq. 3.46 removes each direction's coupling).
Any residual of $G_s/G_t$ from $-2$ at finite $L$ is the cutoff effect we are measuring.

## Addendum (2026-06-11) -- unify the bare kernel as a ConservedCurrent member `build_W_bare`

Design decision (user): the current kernel $W$ categorically belongs to $K$; the disp check is just
exercising one of $K$'s routine kernels.  So the BARE kernel lives as a member of `ConservedCurrent`,
not as a free function and not lifted to a separate header.  $\kappa$ is treated as a relative weight
(lattice $\to$ continuum), factored OUT of the bare kernel.

- **Conserved kernel `build_W` is UNCHANGED** ($W=-C/\lambda_M$, $\kappa$-in; it is $\partial D_W/\partial\theta$,
  the Noether kernel; feeds `apply_k`).  Shared with meson_pq/disc -- do not touch.
- **New member `ConservedCurrent::build_W_bare`** (two overloads, spatial `{s,BaseLink}` and temporal
  `{t,Idx}`): same `d_coo_format` $=iC$, but divided by $(i\kappa)$ instead of $\times i/\lambda_M$, giving
  $\tilde W = C/\kappa = -P U$ (spatial) / $C^t/\kappa_t$ (temporal).  HOST `std::vector<COOEntry>` output
  (for dense-propagator traces), NOT the GPU `COO<N>` of `build_W`.  $\kappa$ from `D.DW.bd.kappa[il]`
  (spatial) / `D.DW.kappa_t[ix]` (temporal).
- **disp** constructs a `ConservedCurrent` (its ctor `cudaMalloc`s the multishift scratch
  $\sim(2\,\text{size}{+}2)N\cdot16$ B $\approx$ 29 MB at $L{=}4$, and needs an `OverlapWMass` -- no
  `update()`/solve, since `build_W_bare` only uses `D.DW`).  Accepted minor cost; $L{=}4$ is rarely run
  post-debug.  `build_W_disp`/`build_W_disp_t` become thin wrappers: call `kop.build_W_bare(coo,U,el)`
  then repack `COOEntry`$\to$`Ent` (the host type adaptation stays in disp; the $W$ DEFINITION is in $K$).
- **exact** is unchanged: it builds the full conserved $K$ via `apply_k` (the $\kappa$-in `build_W`) and
  divides the assembled $K$ by $\kappa$ (the same weight, pulled out) -- identical by linearity.

Pure refactor: disp/exact numerics are unchanged; only the home of the $\kappa$-division moves into $K$.

### Naming refinement (2026-06-11, later): `_ov_kappa` normalized kernel + single $\kappa$ lookup

$\kappa$ is "just a tedious factor"; name the normalized objects explicitly and let the linearity of $K$
in $W$ carry the normalization through:
- `ConservedCurrent::insertion_kappa(el)` -- the ONE place that owns the $\kappa$ lookup (spatial:
  `kappa[link]`; temporal: `kappa_t[site]`).
- `ConservedCurrent::build_W_ov_kappa(en,U,el)` (renamed from `build_W_bare`) $=W/\kappa=C/\kappa$, using
  `insertion_kappa`.  disp consumes this.
- **K_ov_kappa** $=K/\kappa$: because $K$ is linear in the single insertion $W$, dividing the assembled
  dense $K$ by `kop.insertion_kappa(el)` IS $K$ built with $W_{ov\kappa}$ -- exact, no parallel `apply_k`.
  jj_exact_diag_deter_free does exactly this (one scalar per channel, $\kappa$ from `kop`).
(If a non-dense K_ov_kappa operator is ever needed, a `cfg_ov_kappa` flag scaling `operator()` by
$1/$`insertion_kappa` reuses the same lookup -- not needed now.)

## Resolved

- **Q1 (propagator for the comparison) = lattice** (`tag=''`, loc `corr_deter_local_L{1,2}` vs disp
  `corr_deter_disp_L{1,2}`). Notebook keeps a `tag` switch (continuum one edit away).
- **Q3 (contraction) = insertion-diagonal over links is correct** (user-confirmed); reproduces
  Eq.(4.28) via the triangle identity, as derived above. No rework.
- **Separate program** (not folded into jj_local_deter); output key `sp` under `corr_deter_disp*`.
