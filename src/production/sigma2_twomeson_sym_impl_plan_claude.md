# Two-meson-only GEVP on the symmetrized-basis perambulators

## Goal (NM 2026-09-22)
Redo the two-meson-only GEVP (`{s2, O2m, O1m}`, PP flavor, MODE_CONTACT=1) on the new
`distill_Nv24_sym` perams (L2 gsq1.0, truncated 24-of-84), and compare with the `distill_Nv24_v2`
result (no-Hankel NKEEP=3: ground state drifting 0.78->0.74 by t=8-12; partial Hankel O2m,O1m [0,2]
reb3@4: ~0.72-0.73 at t=5-7). Reference $2m_{PS}=0.7054$.

## Files
1. `sigma2_twomeson_gevp_claude.py` (NEW) -- lean builder of the 3x3 $\sigma^2\times\sigma^2$ PP block
   per config (same `perm_contrib_folded` + `FFPP` path as `sigma2_combined_gevp_claude.py`, NO shell
   blob), multiprocessing over configs (4 procs x OMP=1), cache keyed by NVDIR. Then per-mode Hankel
   (`hs.hankel_permode`) + reb + fixed-$t_0$ GEVP; runs the two settings above and plots both.

## Chunks
- Chunk 1 (Files: `sigma2_twomeson_gevp_claude.py`): builder + GEVP + plots; run L2 sym 200 cfg.
- Chunk 2: same driver on `distill_Nv24_v2` (same k list) for a same-config A/B if wanted.

## Refs
Block-Hankel/GPOF Aubin-Orginos arXiv:1010.0202; distillation Peardon et al. arXiv:0905.2160.
