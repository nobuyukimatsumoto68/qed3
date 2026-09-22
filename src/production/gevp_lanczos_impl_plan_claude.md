# Coupled asymmetric GEVP + Lanczos (polyiterated) for the sigma\sigma - F^2 (0++) mixing

Plan for the coupled GEVP that extracts the $0^{++}$ spectrum and the $F^2$-$\sigma\sigma$ mixing ($\Delta_-$ vs
Chester-Pufu 1603.05582 {3.30,3.65,3.77}) from the EXACT distillation four-point (`distill_contract_claude.py`)
+ the glue $F^2$. Implements NM's `asymm_gevp.pdf` idea: an ASYMMETRIC GEVP (distinct source/sink bases, because
$\sigma_{FS}$ is non-Hermitian) and a Lanczos-type "polyiterated" enlargement (time-shifted operator copies ->
block-Hankel correlator; GEVP<->Lanczos, Wagman 2406.20009).

## Operators (source/sink basis)
$O = \{\,F^2,\ \sigma_{PS}\sigma_{PS},\ \sigma_{FS}\sigma_{FS},\ \sigma_{FS}\sigma_{PS}\,\}$ at a timeslice, $\ell{=}0$.
$\sigma_{PS}=\eta^\dagger\xi+\xi^\dagger\eta$ (parity-even), $\sigma_{FS}=\eta^\dagger\xi-\xi^\dagger(1-D_{ov}^\dagger)\eta$
(GW, NON-Hermitian). Non-Hermiticity -> the correlator matrix is generically ASYMMETRIC $C_{ab}\ne C_{ba}$.

## Correlator matrix $C_{ab}(t)=\langle O_a(t)\,O_b(0)\rangle_c$ (vacuum-subtracted)
- **$\sigma\sigma$-$\sigma\sigma$ block** (distillation four-point): each $\sigma_a\sigma_a$ operator's four-point is
  $G_4=\langle S_4\rangle+\langle\tilde S_4\rangle$ (5.5). $S_4$ = all legs $\tau$ (same for every channel).
  $\tilde S_4$ = legs per vertex: a vertex at timeslice $t$ (sink operator $O_a$) uses $a$'s $\tilde S$
  ($\tau$ if PS, $-\tau'$ if FS); a vertex at $0$ (source $O_b$) uses $b$'s. So `compute_diags` is generalized to
  TWO leg-providers keyed by the SINK timeslice of each leg (furnishing sits at the sink): `legs_t` (sink at $t$)
  and `legs_0` (sink at $0$). Then $C_{ab}=G_{10}[\tau,\tau] + G_{10}[\text{leg}_a@t,\ \text{leg}_b@0]$.
  - PS·PS = $2G_{10}[\tau]$; FS·FS = $G_{10}[\tau]+G_{10}[-\tau',-\tau']$; PS·FS = $G_{10}[\tau]+G_{10}[\tau@t,-\tau'@0]$;
    FS·PS = $G_{10}[\tau]+G_{10}[-\tau'@t,\tau@0]$ -> PS·FS $\ne$ FS·PS (asymmetry).
  - The $\sigma_{FS}\sigma_{PS}$ operator (mixed within one timeslice) = DEFERRED to 1b (its two same-timeslice
    vertices are different types -> per-vertex, not per-timeslice, leg keying).
- **$F^2$-$F^2$**: from the existing glue ($O_F(t)$, `glue_f2_v2_shapes` op 0), same configs.
- **$F^2$-$\sigma\sigma$ cross**: $\langle F^2(t)\sigma\sigma(0)\rangle$ AND $\langle\sigma\sigma(t)F^2(0)\rangle$
  (both, since asymmetric) -- gluonic $O_F$ times the equal-time scalar $[D_S^2+D'_S]$ (exact, de-noised).
- **VACUUM SUBTRACTION**: $C_{ab}^c(dt)=\langle O_a(t)O_b(0)\rangle-\langle O_a\rangle\langle O_b\rangle$
  ($0^{++}$ = vacuum quantum numbers, so $\langle O_a\rangle\ne0$: the $J,H,I,F$ pieces). Config-average each.

## Asymmetric GEVP
$C(t)v_R=\lambda(t,t_0)\,C(t_0)v_R$ (right), $w_L^\dagger C(t)=\lambda\,w_L^\dagger C(t_0)$ (left), distinct because
$C\ne C^\dagger$. Solve `scipy.linalg.eig(C(t), C(t0))`. $\lambda_n=e^{-E_n(t-t_0)}$ ->
$E_n(t)=\log[\lambda_n(t,t_0)/\lambda_n(t{+}1,t_0)]$ (or $-\log\lambda/(t{-}t_0)$). Light branch = lowest $E$; its
eigenvector's load on $F^2$ vs $\sigma\sigma$ = the mixing; $\Delta_-$ from its plateau. Watch (diagram breakdown)
whether the light eigenvector loads onto the genuinely-connected $A/B$ direction.

## Lanczos / polyiterated enlargement (NM asymm_gevp.pdf right panel)
Dress operators by the transfer matrix $O_i^{(m)}=O_iT^m$ -> $C_{i\kappa^{(m)}}(t)=C_{i\kappa}(t+m)$. Enlarged
index $I=(i,\kappa^{(1)},\dots)$ gives a BLOCK-HANKEL correlator
$C_{I\tilde I}(t)=\big[\,C(t+p+q)\,\big]_{p,q=0..M}$ built ONLY from measured $C$ at consecutive times (no new
measurement). Do the asymmetric GEVP in the polyiterated basis -> more states / earlier plateaus. Exact
distillation (tiny errors to $dt\sim29$) makes the Hankel structure clean (Lanczos/Prony's noise weakness is
mitigated). Ref Wagman 2406.20009.

## Chunks
1. **$\sigma\sigma$ 2x2 asymm GEVP** (this chunk): generalize `compute_diags` (two sink-keyed leg-providers);
   build vacuum-subtracted $2\times2$ $\{\sigma_{PS}^2,\sigma_{FS}^2\}$ $C(dt)$ on the 400-cfg Nf2 sets; asymm GEVP;
   effective energies + eigenvector content. Files: `distill_gevp_claude.py`.
2. **Full 4x4 coupled GEVP**: add $F^2$ (glue) + $F^2$-$\sigma\sigma$ cross + the $\sigma_{FS}\sigma_{PS}$ operator;
   vacuum-subtract; $\Delta_-$ vs CP.
3. **Lanczos/polyiterated**: block-Hankel enlargement; asymm GEVP in the polyiterated basis; compare spectra/
   plateaus; check the $A/B$-carried fermionic state resolves earlier.

## Validation
- Hermitian limit check: the PS·PS-only (Hermitian) GEVP eigenvalues match a plain effmass.
- $C(t_0)$ well-conditioned (drop near-null directions); eigenvalues real+positive within errors for the
  Hermitian sub-block; light $E$ stable in $t_0$.
- Lanczos vs standard: same energies where both plateau; Lanczos should extend the plateau / add excited states.
- Jackknife over configs for all energies + $\Delta_-$.
