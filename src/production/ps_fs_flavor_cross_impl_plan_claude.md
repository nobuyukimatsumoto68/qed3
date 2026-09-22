# Flavor-cross matrix $\{\sigma_{PS}^2,\sigma_{FS}^2,\sigma_{FS}\sigma_{PS}\}$ (FS $\ne$ PS) -- impl plan

## Goal
The $0^{++}$ two-meson flavor basis (Chester-Pufu, minus $F^2$): $\{$PP, FF, FP$\}$ = $\{\sigma_{PS}^2,
\sigma_{FS}^2,\sigma_{FS}\sigma_{PS}\}$, each on the point geometry $\{\sigma^2_{00},O_{2m},O_{1m}\}$.
The diagonal PP$=$FF (proven PS$=$FS for the 4-point), but the CROSS $\langle$FF$\,$PP$\rangle$ and the mixed
FP genuinely differ -- a loop threading both FS and PS vertices breaks the GW loop-by-loop no-op. _v1 (Nf2
gsq1.0 L1, 400 cfg).

## Mixed-furnishing contraction rule (derived)
$\sigma=\eta^\dagger S\xi+\xi^\dagger\tilde S\eta$, $S=1$, $\tilde S=1$ (PS) / $-(1-D_{ov}^\dagger)$ (FS). A
closed fermion loop has an $S$-version (all forward improved legs $\tilde\tau$) and an $\tilde S$-version. By
GW the FS furnishing collapses the backward-improved leg to $-\tilde\tau$, so each FS vertex in a loop flips
one leg's sign:
$$
(\tilde S\text{-version of a loop}) = (-1)^{n_{FS}}\,(\text{forward improved loop}),\qquad
n_{FS}=\#\{\text{FS vertices in the loop}\}.
$$
So each closed loop contributes a **flavor factor**
$$
\boxed{\;(1+(-1)^{n_{FS}})\;}\times(\text{forward improved loop}),
$$
and loops multiply. Checks:
- PP (all PS, $n_{FS}{=}0$): factor $2$ per loop $\Rightarrow 2^{\#\text{loops}}$ = the existing `fs_channels_v2`
  NCYC. (reproduces PS.)
- FF (all FS, $n_{FS}{=}k$=loop length): $(1+(-1)^k)$ -> even loops $2$, odd loops $0$; the 4-point's odd
  (3-)cycle always pairs a 1-cycle tadpole (also $0$) -> FF $=$ PP. (reproduces PS$=$FS.)
- **CROSS $\langle$FF$\,$PP$\rangle$** (sink 2 FS, source 2 PS): diagram $E=C_S^2$ = two loops each (1 sink FS
  + 1 source PS) -> $n_{FS}{=}1$ per loop -> factor $0$ -> **$E$ VANISHES**; diagram $B$ (nested, 1 loop, 2 FS)
  -> $n_{FS}{=}2$ -> survives. So the cross loses the $E$ two-meson piece -> genuinely $\ne$ PP$\,$PP.

**Implementation:** identical to `fs_channels_v2` (forward improved `AblkS`, the 10 diagrams via `PERMS` cycle
decomposition), but replace the uniform `NCYC=2^{#cyc}` with a per-perm flavor factor
`FLAVFAC[ip] = prod over cycles c of (1 + (-1)^{ sum_{v in c} fsmask[v] })`, where `fsmask` over the 4
vertices (0,1 sink; 2,3 source) is set by the flavor operator pair. Vertex-flavor patterns:
PP=[PS,PS], FF=[FS,FS], FP=[FS,PS]. All legs forward improved $\tilde\tau$ (the $(-1)$ signs live in FLAVFAC).

## Chunks
### Chunk 1 -- flavor machinery + validation (sigma^2_00 geometry)
Build the 3x3 flavor matrix $\{$PP,FF,FP$\}$ at $\sigma^2_{00}$ with FLAVFAC. VALIDATE: (a) PP reproduces
`fs_channels_v2` C00; (b) **FF $=$ PP to machine precision** (PS=FS); (c) the cross FF-PP $\ne$ PP and $E$
vanishes. Files: NEW `sigma2_flavor_cross_v2_claude.py` (reuses `fs_gevp_point` PERMS/perm_contrib_folded/
op_vspec + the smear-free forward AblkS; config-parallel cache).

### Chunk 2 -- full point set x flavor (9 ops) + GEVP
Extend to $\{\sigma^2_{00},O_{2m},O_{1m}\}\times\{$PP,FF,FP$\}$; flavor GEVP + Hankel; does FF/FP resolve a
state PP misses? Parity note: $\sigma_{FS}\sigma_{PS}$ is parity-odd in continuum -> FP overlaps parity-odd
states, only $0^{++}$ at finite $a$.

## CHUNK 1 RESULT (2026-09-10) -- DONE, flavor cross is real AND useful
- `sigma2_flavor_cross_v2_claude.py` (FLAVFAC rule) + `sigma2_flavor_hankel_claude.py`. Validated: PP reproduces
  `fs_channels_v2` (up to a larger s0s range, omax=0 vs 1), **FF $=$ PP to $\sim10^{-4}$** (PS=FS), cross
  $\langle$FF PP$\rangle$ distinct (negative, $\sim3\%$ of PP, $E=C_S^2$ vanishes as predicted).
- KEY: PP and FF overlap the SAME $0^{++}$ tower with nearly ORTHOGONAL amplitudes ($\langle$PP FF$\rangle$ tiny
  vs $\langle$PP PP$\rangle$) -> $\{$PP,FF,FP$\}$ is a genuine variational basis (NOT redundant like $N_v$-smear).
- Flavor Hankel+rebase (Dt=[0,2,4] reb3@4 T0=3) at $\sigma^2_{00}$: state0=state1$\approx0.46$ (ground, doubly),
  **state2$\approx0.62$ = two-meson ($2m_{PS}$)**, clean small errors. So the flavor cross extracts the two-meson
  comparably to the geometry basis, via the FS furnishing. Fig `figs/sigma2_flavor_hankel_*`.

## CHUNK 2 RESULT (2026-09-10) -- DONE. FS-PS difference = single-meson FILTER (key finding)
Drivers: `sigma2_flavor_geom_v2_claude.py` (one flavor x 3-geom), `sigma2_flavorgeom_full_v2_claude.py` (full 9-op).
- **FP (PS-FS) on {s2_00,O_2m,O_1m}** production Hankel+rebase: state0~0.46 (ground), state1~0.62 (two-meson) --
  SAME spectrum as PP (overlaps the same tower). Fig `figs/sigma2_flavorgeom_FP_*`.
- **Full 9-op flavor x geometry** GEVP: state0=state1~0.46 (ground, doubly) + state2~0.62 (two-meson). NO new
  3rd state -- combining both variational axes reconfirms the 2-state spectrum (ground + two-meson). Fig
  `figs/sigma2_flavorgeom_full_*`.
- **FF-PP cross** <sigma_FS^2 sigma_PS^2>: small, NEGATIVE, ~-3 to -5% of PP-PP (E=C_S^2 vanishes as derived).
- **KEY: (PP-FF) ANTISYMMETRIC channel** (FS-PS difference), Hankel+rebase on the geom basis: the one-meson-rich
  0.46 ground CANCELS (PP and FF overlap it equally), so **state0~0.62 = two-meson becomes the GROUND**, plus a
  NEW **state1~0.92** (higher, clean at t=3-7). So the FS furnishing is a FILTER for single-meson contamination
  -- it projects out the 0.46 one-meson tower and isolates the two-meson + higher. Directly addresses the
  single-meson contamination of sigma^2. Fig `figs/sigma2_antisym_PPmFF_*`. Full cache
  `sigma2_flavor_cache_claude/sigma2_flavorgeom_FULL_*_d1_claude.npy` (ncfg,9,9,DT; idx=flavor*3+geom).

## Open questions / NEXT
- The ~0.92 antisym state1: identify (two-meson excited? (2,2)?). Run at other gsq/Nf/L. The antisym-as-filter
  is the main deliverable -- could sharpen the bound-state m1-2m_PS by removing the 0.46 contamination.
- FP parity-odd in continuum -> overlaps parity-odd states, $0^{++}$ only at finite $a$ (interpretation caveat).
Refs: qed3int_v3-4 Ch.5 (Eq 5.1-5.5); GW IV.17; distillation Peardon 0905.2160; CP 1603.05582.
