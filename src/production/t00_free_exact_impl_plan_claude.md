# Free-limit exact overlap propagator at L3/L4 (third lattice point for $\Delta_T$) -- impl plan

## Goal (NM 2026-09-22)
Confirm the free $\Delta_T=3$: we have only L1 (12 sites, $m_{\ell=2}\approx0.470$) and L2 (42 sites, 0.5410) vs $3/R=0.567$.
Need L3/L4 free energies without a GPU peram: build the exact free $D_{\rm ov}$ in python (temporal momentum space,
$N_t=128$), invert it, and evaluate the $\ell=2$ $O_H$ and $\ell=0$ $\sigma$ correlators directly (no distillation).

## Method
- Free field: $D_W$ is block-diagonal in temporal momentum $p_t$. Per $p_t$: $D_W(p_t)=$ spatial part ($W$ hop with $r$ =
  the production Wilson $r$, plus the Wilson diagonal / $M_5$) + temporal hop $\propto\sigma_3$ (Wilson temporal term,
  $a_t$ anisotropy). Conventions read off from `includes/dirac_simp.h` / `dirac_ext.h` / `overlap.h` (NOT guessed).
- $D_{\rm ov}(p_t)=1+ H/\sqrt{H^\dagger H}$-type form exactly as in overlap.h, evaluated by eigen/SVD of the $2N_s\times2N_s$ block.
- $G(t)=\frac1{N_t}\sum_{p_t}e^{ip_tt}D_{\rm ov}^{-1}(p_t)$ -> position-space blocks $G(t,0)$.
- Correlators (same as distillation contractions but with $\Phi\to$ full vertex): $C_{\ell=2}(t)=-2\,{\rm Re\,Tr}[W^{(m)}G(t)W^{(m)}G(-t)]$,
  $C_\sigma(t)=-{\rm Tr}[\Sigma G(t)\Sigma G(-t)]$, $\Sigma={\rm diag}(A_xY_{00})$.
- VALIDATION (mandatory before L3/L4): reproduce peram numbers L1 $m_{\ell2}\to0.470$ trend / L2 plateau 0.5410, and
  $m_\sigma$ L1 0.379 / L2 0.393 (log-effmass, same dt).

## Files
- NEW `t00_free_exact_claude.py` (all of the above; caches $G(t)$ per L in t00_ham_cache_claude/).
- reuse read-only: `geom_hopping_claude.py`, `t00_stress_ham_interacting_claude.build_W_gauge`, `t00_ell2_claude.real_ylm`,
  `t00_dump_links_claude.o` (link tables n3/n4, compile with -DN_REFINE_CLI=3/4), `distill_contract_claude.dual_areas_from_mesh`
  (n4 has no dual files -> fall back to link-volume based site areas from the link table kappa if needed).

## Chunks
- A. Read conventions; build $D_W(p_t)$, $D_{\rm ov}$, $G(t)$; validate on L1 vs peram (`data_free/distill_Nv24`).  Files: driver.
- B. Validate on L2 (Nv84).  Files: driver.
- C. Run L3, L4; tabulate $m_{\ell2}$, $m_\sigma$, $\Delta_T=2m_{\ell2}/m_\sigma$ vs $1/N_{\rm sites}$; extrapolate.  Files: driver + md/memory update.

## Open questions
- exact $M_5$, $r$, $a_t$ used for the free peram generation (read from the distill driver source / peram meta).
- site-area weights for $\sigma$ at n4 (no dual mesh file).

## RESULTS (2026-09-22)
- Chunk A/B VALIDATION: exact $G(t,t')$ vs peram blocks $V\tau V^\dagger$: worst $|{\rm diff}|=1.4\times10^{-5}$ (L1), $1.2\times10^{-5}$ (L2)
  = the peram solve tolerance ($10^{-5}$). Conventions (dirac_ext.h temporal hop, kappa_t=dual_area/mean_ell/at, M5=-1, r=1,
  antiperiodic $p=(2n+1)\pi/N_t$, $D_{\rm ov}=1+D_W(D_W^\dagger D_W)^{-1/2}$) are exactly right. GW check $|G_0+G_0^\dagger-1|\sim10^{-16}$.
- Exact correlators (a_t=0.2, N_t=128, plateau dt~30-45; caches t00_ham_cache_claude/free_exact_{G,corr}_n*_nt128_at0.2):

  | L | N | $m_\sigma$ ($\ell=0$, $2E_1$) | $m_{\ell=2}$ ($E_1+E_2$) | $\Delta_T=2m_{\ell2}/m_\sigma$ |
  |---|---|---|---|---|
  | 1 | 12 | 0.3784 | 0.4672 (still falling @dt30) | 2.469 |
  | 2 | 42 | 0.39171 | 0.53332 | 2.7230 |
  | 3 | 92 | 0.39479 | 0.56310 | 2.8526 |
  | 4 | 162 | 0.39590 | 0.57409 | 2.9002 |

  Extrapolation in $1/N$ (or mean_ell$^2$): linear L3-L4 2.9626, linear L2-L4 2.9619, quadratic 2.9630 -> **2.963(1)**,
  NOT 3: a 1.2% deficit that is $N$-independent => not a spatial artifact. Plot figs/t00_free_exact_DeltaT_vs_N_claude.png.
- NOTE the peram-based numbers were precision-floor-limited: true L2 $\ell=2$ asymptote 0.5333 (not 0.5410); L1 $\ell=0$
  keeps falling below 0.563. At L3/L4 the $\ell=0$ $O_H$ correlator -> CONSTANT ($m_{\rm eff}\to0$): $\int T_{00}=H$ is
  conserved, its connected 2pt is $t$-independent -- direct confirmation that $\ell=0$ is the charge, not a state.
- Lattice single-particle energies at L4: $E_1=0.19795$, $E_2=0.37614$, ratio 1.900 (continuum 2). Both $m_\sigma$ (0.396 vs
  2/R=0.378) and $m_{\ell2}$ (0.574 vs 0.567) OVERSHOOT the nominal $R=1/0.189$.
- HYPOTHESIS: temporal discretization at FIXED $a_t=0.2$ (nonlinear dispersion $m_{\rm lat}=h(a_tE)$, $E_2/E_1\ne2$).
  Test running: L2-L4 at $a_t=0.1$, $N_t=256$ (log t00_free_exact_at0p1_claude.log).
- **$a_t$ TEST DONE ($N_t=256$, plateau dt~95):** $\Delta_T$ at $a_t=0.1$: L2 2.7372, L3 2.8743, L4 2.9242 -> $1/N\to0$: **2.990**
  (vs 2.963 at $a_t=0.2$). Deficit from 3: 0.0374 -> 0.0102, ratio 3.67 = $O(a_t^2)$ (4 expected). Quadratic-in-$a_t$
  extrapolation of the two spatial-continuum intercepts: **2.999**. => free $\Delta_T=3$ CONFIRMED (spatial $O(a^2)$ +
  temporal $O(a_t^2)$ artifacts, both extrapolate away). $E_2/E_1$ at L4: 1.900 ($a_t$=0.2) -> 1.924 ($a_t$=0.1).
  Plot figs/t00_free_exact_DeltaT_vs_N_claude.png (both $a_t$). Caches free_exact_*_n{2,3,4}_nt256_at0.1.
- Implication for the interacting analysis: ratios $2m_{\ell2}/m_\sigma$ (or $/m_{\rm axial}$) at fixed $a_t$ carry a ~1.2%
  ($a_t$=0.2) / 0.3% ($a_t$=0.1) temporal bias in the free limit; the L1/L2 spatial artifact is far larger (-0.28/-0.10).
