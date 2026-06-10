# Continuum (Eq C.28) vs lattice $D_{\rm ov}^{-1}$ — comparison diagnosis (free, L=1)

> Summary: $a_t$ identical (not a mismatch); BOTH sides homogeneous (icos symmetry intact on the
> correct lattice file `prop_deter_L1`; continuum exact). Differences = UV $S^2$ mode tower (small
> $dt$) + antiperiodic temporal BC (large $dt$); effective masses agree to ~few %. Compare ONLY
> SU(2)-frame-invariants ($\lVert P\rVert_F$/det). RETRACTED: the earlier "7% mass deficit" AND the
> "~3% real vertex scatter" (the latter was from analyzing the STALE `prop_exact_L1` file).

Comparison of the continuum propagator `cont_prop_L1/Dinv.0.h5` against the other agent's exact
lattice `data_free_vmRe0.000000vmIm0.000000/prop_exact_L1/Dinv.0.h5` (generator
`jj_propagator_deter_claude.cu`), free field, L=1, massless, north-pole same-site $(1,1)$ block
$G(dt)$ vs temporal separation $dt$ (slices), $\tau=a_t\,dt$.

## $a_t$ is the SAME in both — not a naive mismatch

- Lattice generator: `jj_propagator_deter_claude.cu:252` `const double at = 0.2;`.
- My continuum: default `at = 0.2` (`cont_prop_eigbasis_claude.cpp:536`).

So both use $a_t=0.2$, $N_t=128$, $r=1$, $M_5=-1$, unit-sphere geometry. Lattice anisotropy
$\nu_0=\nu_1=1.0$ (default, `jj_propagator_deter_claude.cu:219,225`); my continuum assumes the
isotropic $c=1$ radial metric $ds^2=d\tau^2+d\Omega^2$ that $\nu_0=1$ is meant to realize.
**A code-level $a_t$ (or $\nu$) mismatch is ruled out.**

## RETRACTION of the earlier "7% finite-spacing mass deficit"

An earlier version of this doc claimed the lattice ground mass is ~5-7% below the continuum
$a_t\lambda_0=0.2$. **That was an overclaim and is retracted.** It came from a LOCAL effective-mass
slope $m_{\rm eff}(dt)=\ln[G(dt)/G(dt{+}1)]$ which here is contaminated by two effects that have
nothing to do with a clean ground-state mass (below). There is NO demonstrated clean mass deficit;
the masses are consistent to within a few % but cannot be cleanly extracted at this coarse corner.

## What actually differs (all gauge-invariant-confirmed)

North-pole same-site temporal $G(dt)$, then site/gauge checks. $m_{\rm eff}$ values bracket
$0.19$-$0.23$ (cosh-fit 0.23 UV-biased; local plateau 0.19; circle-distance 0.196) — that SPREAD is
the contamination, not a measurement.

1. **Antiperiodic temporal BC (dominates the tail).** The lattice lives on a TIME CIRCLE
   ($N_t=128$): $G(dt)$ falls to a minimum near $dt=N_t/2=64$ and climbs back to ~0.12 by $dt=127$
   (decays with CIRCLE distance + a backward image). My continuum is on the infinite line (no image,
   the "ABC not necessary" choice) and decays monotonically. So beyond $dt\sim30$-40 they are
   structurally different objects; a naive $+$image overshoots the far side ~15x, a $-$image zeros
   $dt{=}64$, so the fermion antiperiodic sign structure must be DERIVED not guessed (not done).

2. **UV mode tower (dominates small $dt$).** At $dt{=}1$ continuum$=1.98\gg$lattice$=0.30$. The
   continuum carries the full $S^2$ angular tower $\sum_n\frac{n+1}{4\pi}e^{-(n+1)a_t dt}$ (all $n$ to
   40, full weight); the lattice UV is REGULATED (12-site spectrum cutoff + dispersion). In the clean
   tail window $dt=7$-16 the continuum $m_{\rm eff}$ sits ABOVE the lattice (e.g. $dt{=}10$: 0.256 vs
   0.224), converging as the tower dies — the lattice is the better-behaved one there.

3. **Difference is the SPATIAL ($S^2$) spectrum (user's diagnosis, confirmed).** At
   $\theta{=}\theta'{=}0$ the temporal correlator is purely a sum over $S^2$ angular modes; $a_t$ is a
   shared scale, so the effective-mass-curve difference is set by the $S^2$ eigenvalues+residues. The
   ground angular mode (mass) agrees to a few %; higher modes differ (12-site discretization + cutoff)
   — consistent with the paper's statement that L=1 reproduces only $n\lesssim6$.

## CORRECTION: I used the WRONG lattice file; the "vertex scatter" was retracted

An earlier version of this doc claimed the lattice has a ~3% site-to-site scatter (a "real
asymmetry, not gauge"). **WRONG, retracted.** Cause: I analyzed the STALE
`data_free_.../prop_exact_L1/Dinv.0.h5` (leftover from the old, since-deleted `jj_propagator_exact`
generator), NOT the live `data_free_.../prop_deter_L1/Dinv.0.h5` (written by
`jj_propagator_deter_claude.cu`, the file the other agent uses).

On the CORRECT file `prop_deter_L1`, using the (right) layout $r=N_x t+N_S\,\text{site}+\text{spin}$
(spinor adjacent, $2i{+}a$): the same-site block is **purely $c_3\sigma_3$**, with $c_3=|P_{00}|$
**constant across all 12 sites to 3.7e-15**, $\lVert P\rVert_F$ to 2e-15, off-diagonal machine-zero.
**Icosahedral symmetry is intact.** This matches the other agent's authoritative note
`propagator_frame_dependence_claude.md` exactly.

## Frame dependence is SU(2), not U(1) (other agent's `propagator_frame_dependence_claude.md`)

Each site carries its own local spinor frame (from link angles `alpha` + spin connection `Omega`);
the symmetry maps $n\to n'$ together with an $SU(2)$ rotation, so the block transforms COVARIANTLY
$P(n')=U_n P(n) U_n^\dagger$. Therefore individual components $|P_{ab}|$ are NOT invariant; compare
only $\lVert P\rVert_F$, $\det P$, $\mathrm{tr}(P^\dagger P)$, or eigenvalues. (The same-site,
time-separated block is special: purely $\sigma_3$ = the global time axis, so it has NO frame
freedom and is constant. Off-SITE blocks are the genuine SU(2) case.) Continuum is likewise
homogeneous (same-site $\lVert P\rVert_F$ const across sites to ~1e-13).

## Clean gauge-invariant comparison (prop_deter_L1 vs continuum, $\lVert P\rVert_F$)

| $dt$ | $\lVert P\rVert_F^{\rm lat}$ | $\lVert P\rVert_F^{\rm cont}$ | $m_{\rm eff}^{\rm lat}$ | $m_{\rm eff}^{\rm cont}$ |
|---|---|---|---|---|
| 8  | 8.16e-3 | 3.57e-2 | 0.265 | 0.290 |
| 16 | 1.25e-3 | 4.99e-3 | 0.211 | 0.215 |
| 30 | — | — | 0.194 | 0.201 |
| 40 | — | — | 0.191 | 0.200 |

Same physics as before survives: effective masses agree to ~few % (lat->0.19, cont->0.20); the
small-$dt$ difference is the UV $S^2$ mode tower (continuum heavier), the large-$dt$ difference is the
antiperiodic temporal BC (lattice on a circle). Overall normalization differs (lat $\lVert P\rVert_F$
~4-17x smaller = the $1/\bar a_s a_t$ scale, irrelevant). NO spurious asymmetry, NO 7% mass deficit.

## Status / next
- $a_t$ identical (0.2), $\nu_0=1$ both — no parameter mismatch.
- Both sides homogeneous (icos symmetry intact on the correct lattice file; continuum exact).
- Differences are physical: UV $S^2$ tower (small $dt$) + antiperiodic BC (large $dt$); masses agree
  to few %. Compare INVARIANTS ($\lVert P\rVert_F$/det), per the SU(2) frame point.
- **NEXT:** L=2 lattice $D_{\rm ov}^{-1}$ (`prop_deter_L2`) — UV tower fills in, gauge-invariant
  magnitudes track closer. Use `prop_deter_L<L>`, NOT `prop_exact_L<L>` (stale).
