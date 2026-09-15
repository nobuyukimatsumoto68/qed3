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
