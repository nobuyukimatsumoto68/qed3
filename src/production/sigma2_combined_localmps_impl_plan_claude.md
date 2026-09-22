# Add the local $m_{PS}$ single-meson operator to the combined GEVP basis

## Goal

The combined basis has 3 shell single-meson ops + 3 $\sigma^2$ two-meson ops, and $m_{PS}$ only
appears as the emergent GEVP ground state (no dedicated interpolator). Add the LOCAL $m_{PS}$
operator -- the single $\sigma_{PS}$ (scalar density) meson, vertex $\Gamma=\text{diag}(w Y_{00})$ --
so the GEVP has a clean, high-statistics $m_{PS}$ interpolator and can orthogonalize it out of the
$(2,2)$ / two-meson states sharply.

Mode-space vertex (distillation elemental):
$$
\Phi_\text{loc}(t) = U(t)^\dagger \,\text{diag}(w Y_{00})\, U(t),\qquad U=V\ (2N_s\times N_v).
$$
Its kernel for the two-meson triangle is $K_\text{loc}=U\Phi_\text{loc}U^\dagger$ (= distilled local
vertex). This is exactly the single-$\sigma_{00}$ operator used in the earlier
$\langle\sigma_{00}|\sigma^2\rangle$ triangle check (`sigma2_mPS_gevp` $C_{12}$), so the cross
machinery is validated.

## Basis order (7 ops)

`[local_mPS, shell l1/2, shell l3/2 (2,2), shell l5/2, s2, O2m, O1m]`

Treat `local_mPS` as one more "meson vertex" alongside the shell projectors: the meson x meson
block, meson x $\sigma^2$ cross, and $\sigma^2\times\sigma^2$ block all generalize with the vertex
list `Phi_list = [Phi_loc, Ps_l1/2, Ps_l3/2, Ps_l5/2]`, kernels `K = U Phi U^H`.

## Files

1. `sigma2_combined_localmps_claude.py` (NEW; variant of `sigma2_combined_permode_hankel_claude.py`,
   original untouched). New `one_config_lm` computing the 7x7 tensor; caches to a `_lm`-tagged npy;
   reuses `hankel_permode` / `staged_project` / `rebased_effmass_fixed`.

## Chunks

- Chunk 1 (Files: `sigma2_combined_localmps_claude.py`): `one_config_lm` -- add `Phi_loc`, generalize
  the meson-vertex loop to include it (diagonal + cross-shell in mode space; cross-$\sigma^2$ via the
  blob with `K_loc`). Cache-tagged `_lm`.
- Chunk 2: run 200 cfg (~90 min cross-blob recompute; the machine was loaded, may be faster now),
  then GEVP + plot. No-Hankel first (matches the current baseline), NKEEP=4 or 5.

## Cost

~90 min for the 200-cfg recompute (the 7th op adds a triangle to the blob). After it caches, all
Hankel/window/NKEEP variations are instant.

## Open question

- Keep the $\ell 1/2$ shell too, or replace it with local $m_{PS}$? Both mainly overlap $m_{PS}$, so
  keeping both risks near-collinearity in the GEVP metric. Default plan: KEEP both (7 ops), NKEEP=5,
  and check conditioning; drop $\ell 1/2$ if collinear. Confirm with NM.

## Refs

- Chester, Pufu arXiv:1603.05582 (physics context). Block-Hankel/GPOF Aubin-Orginos arXiv:1010.0202.
- Distillation Peardon et al. arXiv:0905.2160.
