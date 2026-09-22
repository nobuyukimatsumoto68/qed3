# Perambulator regeneration -- consumer changelog (2026-09-17)

Author: "Fin: Two-meson" session. One-page reference for EVERY consumer of the distillation
perambulators (T_00 / axial-VSH / sigma^2 / single-meson / any `load_peram`). Two changes are landing.
Companion derivation: `fs_gw_collapse_v_agent_two_meson_claude.md`.

## Change 1 -- SYMMETRIZED distillation basis (`-DBASIS_SYM=1`)

The basis $V(t)$ is now the low modes of the SYMMETRIZED Wilson normal operator
$$
M = D_W^\dagger D_W + D_W D_W^\dagger
$$
NOT the one-sided $D_W^\dagger D_W$. Rationale: the one-sided operator's eigenbasis does not carry the
theory symmetry (parity / sigma3 maps $D_W^\dagger D_W$ eigenmodes to $D_W D_W^\dagger$ eigenmodes -- a
DIFFERENT operator), so a TRUNCATED basis ($N_v < 2 N_s$) spans a non-symmetry-invariant subspace and breaks
symmetry-protected zeros. Concretely it leaked the single meson $m_\text{PS}$ into $\sigma^2$ at L2
$N_v = 24$-of-$84$. The symmetrized operator is invariant under the swap of the two normal ops, so its
low-mode truncation is symmetry-invariant.

Impact on consumers:
- COMPLETE basis (L1 $N_v=24=2\times12$; L2 $N_v=84=2\times42$; free): symmetrized $\equiv$ one-sided
  (IDENTICAL, verified $M$ and $M^\dagger$ commute at complete rank). So all COMPLETE-basis results are
  UNCHANGED. (T2a on the new L1 `_sym`: $|V\tau V^\dagger - D_\text{ov}^{-1}| = 1.7\times10^{-5}$, exact.)
- TRUNCATED basis (L2 $N_v=24$-of-$84$): the new basis DIFFERS and is the whole point -- any truncated-L2
  number may shift. Re-run and compare. Known truncated-L2 consumers: T_00 pole-GEVP $\Delta=3.16(23)$
  (Fin: Stress); higher-ell tangential VSH $\Phi_3,\Psi_2,\Psi_3$ + $\Psi_l\sim\Phi_{l+1}$ (Fin: axial sp);
  sigma^2 P+ GEVP (this session).
- Flag `BASIS_SYM` lives in BOTH `distill_peram_claude.cu` (~line 390) AND `distill_peram_mrhs_claude.cu`
  (main eig ~line 494 + gauge-check eig ~line 527). The regen builds the mrhs file.

## Change 2 -- the DAGGER: `/peram/tau_gw` is BUGGY. DO NOT USE IT.

`peram/tau_gw` $= V^\dagger (1 - D_\text{ov}^\dagger) D_\text{ov}^{-1} V$ is built by applying
$(1 - D_\text{ov}^\dagger)$ to the FORWARD inverse. Via GW ($D_\text{ov}^{-1} + D_\text{ov}^{-\dagger} = 1$)
that equals $1 + D_\text{ov}^{-1} - D_\text{ov}^\dagger$ -- it carries the BARE $D_\text{ov}^\dagger$
(Dirac adjoint, $O(1)$), which is NOT the adjoint inverse. That bare term is a spurious $O(1)$ artifact.

RULE for all consumers:
- Use `tau` for the forward leg $\langle\xi\,\eta^\dagger\rangle = \tau = D_\text{ov}^{-1}$.
- For the adjoint / FS-furnished / backward leg use $\delta - \tau = 1 - \tau$ (= T2b-verified
  $\bar\tau = D_\text{ov}^{-\dagger}$), NEVER `tau_gw`, NEVER a literal $(1 - D_\text{ov}^\dagger)$ apply on
  a forward inverse.
- Consequence for scalars: the FS leg collapses (GW) to plain $\tau$, so $\sigma_\text{FS} = \sigma_\text{PS}$
  (smeared: $\Box\tau\Box$, box LEFT of the furnishing, Eq 1.51). Use $\tau$ only.

Backward compatibility: the `_sym` perams STILL contain a `/peram/tau_gw` dataset (same buggy definition) so
existing `load_peram` readers do NOT break -- it simply must remain UNUSED. (It was left in rather than
stripped because unpacking it costs nothing and removing it would `KeyError` older readers; it is NOT a
second solve, just one mat-vec, so it is not on the critical path.)

Consumer status (grepped / confirmed 2026-09-17):
- Fin: Stress tensor -- unpacks `tau_gw` but never uses it; contractions on `tau` / `dtau` only. CLEAN.
- Fin: axial sp -- forward `AblkS` + conj-transpose only; no `tau_gw`, no FS-furnished axial. CLEAN (no-op).
- `distill_contract_claude.py` `four_point` FS path STILL references `-taugw` (line 121/193/288) --
  PRODUCTION FLAG: its FS.FS $\ne$ PS.PS (~0.57x on free L1) because of this; fix = adjoint leg $\delta-\tau$,
  or just use PS (FS==PS). Not used by the P+ flavorgeom GEVP (that uses `AblkS=tau` + flavfac = correct).

## The swap

Nothing is repointed yet. New perams go to a SEPARATE dir `data_<ens>/distill_Nv24_sym/` (suffix `_sym`);
the existing `_v2` perams are untouched. Running now: Nf2 at0.2, L1 gsq0.5 first, then L2 gsq1.0, $N_v=24$,
stride 4, nsrc2 fused. When the `_sym` perams land and validate, consumers repoint from `_v2` to `_sym` and
re-run the truncated-L2 comparisons. A COMPLETE $N_v=84$ `_sym` L2 (for real magnetic GEVP / ground truth)
is a possible later heavier pass (at complete rank `BASIS_SYM` is moot).
