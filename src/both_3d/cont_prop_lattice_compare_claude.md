# Continuum (Eq C.28) vs lattice $D_{\rm ov}^{-1}$ — temporal-decay discrepancy diagnosis

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

## Effective-mass evidence

$m_{\rm eff}(dt)=\ln[G(dt)/G(dt{+}1)]$ (lowest-state mass per slice):

| $dt$ | lattice $G$ | continuum $G$ | $m_{\rm eff}^{\rm lat}$ | $m_{\rm eff}^{\rm cont}$ |
|---|---|---|---|---|
| 1  | 3.01e-1 | 1.98e0  | 0.601 | 1.394 |
| 5  | 5.55e-2 | 7.33e-2 | 0.281 | 0.401 |
| 10 | 1.58e-2 | 1.44e-2 | 0.224 | 0.256 |
| 20 | 1.95e-3 | 1.51e-3 | 0.198 | 0.207 |
| 39 | 4.94e-5 | 3.26e-5 | 0.1905 | 0.2001 |
| ~60 | (BC-far) | | ~0.186 | ~0.200 |

Two distinct effects:

1. **Small $dt$ (UV):** continuum $\gg$ lattice ($G(1)$: 1.98 vs 0.30, factor ~6.6). The continuum
   carries the FULL operator tower $\sum_{n}(n{+}1)e^{-(n+1)a_t dt}/4\pi$ — at $a_t dt=0.2$ this sums
   to $\approx x/(1-x)^2/4\pi$, $x=e^{-0.2}$, a near-divergent UV peak. The lattice has a UV
   REGULATOR (finite mode count + lattice dispersion caps high frequencies), so its short-distance
   value is finite and much smaller. **Expected, not a bug.**

2. **Large $dt$ (IR mass):** continuum plateaus at EXACTLY $m=0.2000 = a_t\lambda_0$
   ($\lambda_0=n{+}|m|{+}\tfrac12=1$, the $\Delta=1$ primary), by construction. The lattice plateaus
   ~5-7% LOWER, $m^{\rm lat}\approx0.186$. Ratio $m^{\rm lat}/m^{\rm cont}\approx0.93$.

## Interpretation — finite-lattice artifacts, not a formula/$a_t$-mapping error

The continuum result is the $a_t\to0,\ L\to\infty$ limit. The $\approx7\%$ IR-mass deficit at
$(L=1,\ a_t=0.2)$ is consistent with finite-lattice artifacts on the COARSEST setup:

- **Spatial ($L=1$):** the 12-site icosahedral $S^2$ Dirac gap $\lambda_0^{\rm lat}$ is below the
  continuum $\lambda_0=1$. $m^{\rm lat}=0.186\Rightarrow\lambda_0^{\rm lat}\approx0.93$ (a few-% gap
  deficit on 12 points is very plausible).
- **Temporal ($a_t=0.2$):** the lattice dispersion gives $E_{\rm lat}=a_t\lambda_0(1-O((a_t\lambda_0)^2))$,
  a $\sim$0.7% shift for a simple $\cosh$; the overlap construction may enlarge it.
- Note the paper's `prop_ov` validation (Fig., $T=12,L_t=168$) used $a_t=12/168=0.0714$ — nearly
  $3\times$ FINER in time than $a_t=0.2$ here, so its temporal artifacts were much smaller.

## Decisive tests to settle continuum-vs-lattice (proposed)

1. **Spatial refinement:** have the other agent generate the lattice $D_{\rm ov}^{-1}$ at L=2 (and
   L=4). If $m^{\rm lat}$ rises toward 0.2 as $L$ grows, the deficit is spatial discretization and
   the continuum is vindicated. (My continuum mass is $L$-independent: 0.2 at all $L$.)
2. **Temporal refinement:** regenerate BOTH at $a_t=0.0714$ ($T=12,L_t=168$, the paper's setup). If
   the lattice plateau moves to $0.0714$ and the curves overlap, the residual at $a_t=0.2$ is
   temporal dispersion.
3. The continuum is the LIMIT; do NOT expect exact agreement at $(L=1,a_t=0.2)$, which is the
   coarsest corner.

## Caveat (genuinely open)
If, after spatial+temporal refinement, the lattice mass does NOT converge to $a_t\lambda_0$, then the
$\tau=a_t\,dt$ mapping or an overall scale in the continuum is wrong and must be revisited — the
user's $a_t$ instinct would then point at a scale convention (e.g. the anisotropy $\nu$ / speed-of-
light tuning relating the temporal and spatial scales), not a literal $a_t$ value.
