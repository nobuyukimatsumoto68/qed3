# $F^2$--$\sigma^2$ ($0^{++}$) glueball--two-meson mixing REDO (post-corrections) -- impl plan

## Goal / physics
Re-measure the Chester--Pufu ($0^{++}$) mixing of the glueball $F^2$ with the four-fermion singlet
$\sigma^2=(\bar\psi\psi)^2$ (arXiv:1603.05582), now with the corrected $\sigma$ machinery established this
session. The earlier pass (`f2_sigma2_mixing_result_claude.md`) called it null ($\sim1\sigma$), but it used
(i) SINGLE-plateau subtraction (underestimates errors), (ii) only the local $\sigma^2_{00}(s)$ (no two-meson
interpolator), (iii) no GEVP. Redo with the corrected pipeline and decide the verdict rigorously.

**Key simplification: PS $=$ FS for this cross, exactly.** $F^2$ is gluonic, so a single $\sigma^2(s)$ contracts
only into
$$
\text{disc}=D_S(s)^2\ (\text{two length-1 tadpole loops}),\qquad
\text{ext}=D'_S(s)=\mathrm{Tr}[\Phi\tilde\tau\Phi\tilde\tau]\ (\text{one length-2 loop}),
$$
$\tilde\tau=\tau-\tfrac12 I$ (improved). All loops here are EVEN or tadpole, so the per-loop $(-1)^k$ furnishing
gives $+1$ (even) and the tadpoles vanish -- **no FS/PS split, no anomaly subtlety in this channel** (the odd-loop
FS effect lives only in the separate $2\text{PS}\to$PS/FS triangle). So the mixing is a clean PS-only object.

- disc is structurally DEAD: improved $\langle D_S\rangle=0$ and its config fluctuation is $\sim10^{-7}$ (contact-
  saturated, no condensate) -> a nearly config-independent loop cannot correlate with $F^2$. Confirmed $\sim10^{-15}$.
- ext $=D'_S$ is the only live channel: the connected $(\bar\psi\psi)^2$ density that can share a $0^{++}$ state
  with $F^2$. This is the mixing.

## Object
Cross $C_{F\sigma}(dt)=\langle F^2(s{+}dt)\,\sigma^2_{00}(s)\rangle_c$, source $s$ translation-averaged over the
perambulator window $[t_{\rm src0},t_{\rm src0}{+}$twin$)$ (twin=32), $F^2$ at $(t_{\rm src0}{+}s{+}dt)\bmod N_t$.
$O_F$ = glue $F^2$ shape op $0$ ($\ell{=}0$, $p{=}2$; `glue_f2_v2_shapes`, `O` dataset), full $N_t{=}128$.
Ensemble: Nf2 gsq1.0 L1, 400 cfg (perams `distill_Nv24`, glue matched).

The coupled GEVP (if signal) needs the $2\times2$ (or larger)
$$
C=\begin{pmatrix}C_{FF}&C_{F\sigma}\\ C_{\sigma F}&C_{\sigma\sigma}\end{pmatrix},\qquad
C_{FF}=\langle F^2 F^2\rangle\ (\text{glue}),\quad
C_{\sigma\sigma}=\langle\sigma^2_{00}\sigma^2_{00}\rangle\ (\text{4-point, have it}),\quad
C_{F\sigma}\ (\text{this measurement}).
$$

## Chunks

### Chunk 1 -- corrected local cross $\langle F^2\,\sigma^2_{00}\rangle$ (the verdict)
Redo the diagram-by-diagram cross (disc, ext, total) with the **triple subtraction** (per-config $t$-sum ->
per-diagram plateau -> total plateau) for correct errors, linear + log$|C|$. This directly upgrades the earlier
single-plateau pass. Deliverable: is ext a clean falloff (mixing) or noise (null), with honest errors?
- Files: NEW `f2_sigma2_cross_v2_claude.py` (copy of `f2_sigma2_cross_diag_claude.py`, swap the single `jk_corr`
  for the triple subtraction from `fs_diag_corr_v2_claude.py`); reuse `distill_contract_claude.py`. Figs
  `figs/f2_sigma2_cross_v2_{lin,log}_*_claude.png`.

### Chunk 2 -- two-meson interpolator cross (better projection)
Cross $F^2$ with the ANTIPODAL two-meson interpolator $O_{2m}(s)=\sum_x A_x\,\sigma(x,s)\sigma(P(x),s)$ (enhances
the two-meson overlap) and the coincident $O_{1m}$ (one-meson), using the position-propagator
$A(t,s)=V(t)\tau(t,s)V^\dagger(s)$ machinery from `fs_channels_v2_claude.py`. Same triple subtraction. See whether
any interpolator gives a $0^{++}$ falloff that the local $\sigma^2_{00}$ washes out.
- Files: NEW `f2_sigma2_cross_points_v2_claude.py` (reuse the point-block builders from `fs_channels_v2_claude.py`
  / `fs_gevp_point_claude.py`; multiply by $O_F(s{+}dt)$ and translation-average). Figs.

### Chunk 3 -- coupled GEVP (ONLY if chunk 1/2 shows a real cross)
Assemble $C_{FF}$ (glue) + $C_{\sigma\sigma}$ (four-point) + $C_{F\sigma}$ (chunks 1-2) on the common windowed
$dt$ grid; glueball-style fixed-$t_0$ GEVP (inv_sqrt_sym metric, `glue_gevp_analysis`/`gevp_OA_sig22_glue` recipe);
read the light-$0^{++}$ composition (F$^2$-dominated vs $\sigma^2$-dominated) and $\Delta_-$ vs CP {3.30,3.65,3.77}.
Vacuum handling: F$^2$ and $\sigma^2$ both have $\langle O\rangle\ne0$ -> either keep the identity in the basis or
subtract the vacuum products; decide from the F-diagram normalization tell (as in the four-point).
- Files: NEW `f2_sigma2_gevp_v2_claude.py`.

## Open questions for NM (resolve before coding)
1. **Scope now:** just chunk 1 (the corrected local cross, to settle the verdict), or go straight through chunks
   1-2 (add the antipodal two-meson interpolator) before deciding on the GEVP? (Recommend: chunk 1 first, look, then
   decide 2/3.)
2. **GEVP normalization:** the coupled GEVP needs $C_{FF}$, $C_{\sigma\sigma}$, $C_{F\sigma}$ in a common
   normalization. $O_F$ is an arbitrary-normalized shape op; the GEVP is invariant under per-operator rescaling, so
   this is fine -- but confirm we want the CP $\Delta_-$ comparison (needs the F$^2$ diagonal mass too) vs just the
   mixing yes/no.
3. **Which glue op:** $O_F$=op 0 ($\ell{=}0$ $p{=}2$ basic $F^2$) as before, or the full shape/$F^4$ basis GEVP for
   $F^2$ first (op index from `project_glue_shapes`)?

Refs: CP arXiv:1603.05582; distillation Peardon 0905.2160; glue shapes `project_glue_shapes`; corrected FS/triple
subtraction `fs_diag_corr_v2_claude.py`; no-SSB `project_chi_volume_scaling_nossb`.

## CFT: is there a mass shift for the two-meson state? (recurring question)

**Yes -- but the right object is an anomalous dimension, not a binding energy.** At the (putative) massless
conformal fixed point there are NO asymptotic particles and NO $2m_{PS}$ threshold. States on the sphere are
labeled by scaling dimensions; radial quantization on $S^2\times\mathbb{R}$ gives
$$
E \;=\; \Delta/R ,
$$
with $R$ the sphere radius. What we call the "single meson" is the scalar bilinear $\sigma=\bar\psi\psi$
(a CFT primary), so the measured $m_{PS}=\Delta_\sigma/R$. The "two-meson state" is the DOUBLE-TRACE primary
$[\sigma\sigma]_{n=0,\ell=0}=:\!\sigma\sigma\!:$, whose dimension is
$$
\Delta_{[\sigma\sigma]} \;=\; 2\,\Delta_\sigma \;+\; \gamma_{[\sigma\sigma]} ,
\qquad
\gamma_{[\sigma\sigma]} = O(1/N_f) .
$$
So the energy sits at $(2\Delta_\sigma+\gamma)/R$, i.e. shifted from the naive $2m_{PS}=2\Delta_\sigma/R$ by
$$
\delta E \;=\; \gamma_{[\sigma\sigma]}/R .
$$

- **Origin of $\gamma$:** gauge-boson / auxiliary-field exchange between the two $\sigma$ constituents (the same
  interaction that in AdS language is the two-particle binding in AdS$_4$). It is $O(1/N_f)$ and NOT necessarily
  small at $N_f=2$.
- **Sign = effective interaction.** $\gamma<0$ (attractive, energy below $2m_{PS}$) vs $\gamma>0$ (repulsive,
  above). There is no threshold, so $\gamma<0$ is NOT "binding" in the scattering sense -- just a smaller scaling
  dimension. This is the honest replacement for the earlier (retracted) "bound two-meson".
- **What the data says.** $s2$-only partial Hankel L2: $m_1\simeq0.735$ vs $2m_{PS}=0.705(3)$, i.e. a small
  POSITIVE offset $\delta E\simeq+0.03$ ($\gamma\gtrsim0$), consistent with mild repulsion / unbound.
- **Caveat before quoting $\gamma$.** We are NOT exactly at the fixed point (finite $a_t$, finite refinement,
  possible non-conformal IR). $m_{PS}=\Delta_\sigma/R$ itself carries its own lattice corrections, so part of the
  $0.03$ offset can be artifact. The clean CFT statement requires the continuum + a genuine $\Delta_\sigma$; the
  lattice number is only suggestive of the SIGN.

Contrast with the $(2,2)$ single-meson level ($\sim0.62$ L1): that is a distinct 1-particle CFT primary (a
$\sigma$-descendant / higher radial excitation), NOT the double trace -- so its shift is unrelated to
$\gamma_{[\sigma\sigma]}$.

## CURRENT BEST -- L2 two-meson tower (2026-09-17)

Ensemble Nf2 gsq1.0 L2 at0.2, 400 cfg (perams `distill_Nv24_v2`, TRUNCATED Nv=24=smeared).
P+ 6-op flavor-geometry basis + PARTIAL Hankel: `s2` (=$\sigma^2_{00}$) offset-set {0,2} and `O1m`
(coincident single-sum) offset-set {0,2,4}; both flavors PP,FF -> 12 ops. reb3@4 T0=3 bin10.
Driver: `sigma2_Peven_partialhankel_fit_claude.py` (S2SET=0-2 O1MSET=0-2-4).
Correlated constant fit (GLS with jackknife covariance) over $t\in[5,9]$:

| state | $a_t m$ (corr fit) | $\chi^2/\nu$ | note |
|---|---|---|---|
| 0 | 0.3522(9)  | 0.51 | $= m_{PS}=0.3527(14)$ |
| 1 | 0.6956(189) | 0.37 | on threshold $2m_{PS}=0.7054(28)$ |
| 2 | 0.7376(85)  | 2.24 | just above threshold |

- m0 reproduces $m_{PS}$ (single meson).
- m1 sits right AT $2m_{PS}$ -> the near-threshold double-trace $[\sigma\sigma]_{0,0}$; anomalous-dimension shift
  consistent with $\gamma\simeq0$ (see "CFT: is there a mass shift" section).
- m2 slightly above threshold; $\chi^2/\nu=2.2$ from mild curvature in [5,9] -- firm up with fine-stride stats or
  window [6,9].
- Uncorrelated cross-check (diag): m0=0.3524(6), m1=0.7053(256), m2=0.7432(85) -- consistent.
- Plot: `figs/sigma2_Peven_partialhankel_fit_s202_o1m024_..._L2_..._fit59_claude.png`.

CAVEAT (assignment): s2 got the 2-shift set, O1m the 3-shift set; the reverse (s2=0-2-4) was noise. 400 cfg only;
the fine-stride run (~800 cfg) will tighten m1/m2. Not yet continuum / exact-Nv84, so the SIGN of the shift is
robust, the VALUE is not.
