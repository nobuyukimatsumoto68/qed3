# Per-mode block-Hankel for the combined (2,2)-shell + $\sigma^2$ GEVP

## Goal

The combined GEVP (`sigma2_combined_gevp_claude.py`) currently runs with a single common
offset list applied to every operator (default `OFFSETS=[0]` = no Hankel). Give each operator
its OWN block-Hankel offset list ("per-mode Hankel") and see whether the four-state resolution
(m_PS 0.35 / (2,2) 0.66 / two-meson 0.76 / excited 0.83) tightens.

Requested per-operator offsets (basis order 0..5):

| idx | operator          | offsets   |
|-----|-------------------|-----------|
| 0   | shell $\ell 1/2$  | `[0,4]`   |
| 1   | shell $\ell 3/2$ (2,2) | `[0,3]` |
| 2   | shell $\ell 5/2$  | `[0,2,4]` |
| 3   | s2                | OPEN (see question) |
| 4   | O2m               | `[0,2]`   |
| 5   | O1m               | `[0,2]`   |

Augmented matrix dimension $D=\sum_i |\text{offs}_i|$; the block at $((i,p),(j,q))$ is $C(t+p+q)_{ij}$.

Method/refs: block-Hankel / GPOF, Aubin-Orginos arXiv:1010.0202 (already cited in `hankel_rebase_scan_claude.py`).

## Files

1. `hankel_rebase_scan_claude.py` -- ADD `hankel_permode(Cts, offs_list)` (new function, leaves `hankel_off` untouched).
2. `sigma2_combined_permode_hankel_claude.py` (new) -- computes/caches the raw per-config
   correlator tensor `allC` (6x6xDTMAX) via `CB.one_config`, then applies per-mode Hankel,
   `staged_project` (reb `REBT@NKEEP`), fixed-$t_0$ GEVP, plot. Caches `allC` to `.npy` so
   Hankel-knob iteration is instant (the expensive cross-blob recompute runs once).

## Chunks

- Chunk 1 (Files: `hankel_rebase_scan_claude.py`): `hankel_permode`. Row index enumerates
  `(i, p) for i in ops for p in offs_list[i]`; `Big[t, r, c] = Cts[t + p + q][i, j]`;
  `tmax = twin - 2*max(all offsets)`.
- Chunk 2 (Files: `sigma2_combined_permode_hankel_claude.py`): caching driver + per-mode GEVP + plot.

## Open question

- s2 (idx 3) offset: not specified. Candidates `[0]` (no Hankel, leave the cleanest op alone)
  or `[0,2]` (match the other two-meson ops). Ask NM.
