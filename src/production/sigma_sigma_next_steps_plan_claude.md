# $\sigma\sigma$ - 0++ mixing: next-steps plan (resume note)

Session 2026-09-05 wrap-up. Benchmark established: connected $\sigma\sigma$ (PS) $\to a_t m\approx1.0$ (lattice),
a two-$\sigma$ scattering state (see `sigma_connected_gevp_benchmark_claude.md`). Ordered plan below.
Ensemble for development: Nf2 $g^2{=}0.5$ L1, $a_t{=}0.2$, massless, 400 cfg, distillation $N_v{=}24$ (exact).

## Step 1 -- Re-derive the four-point coefficients for $\sigma = \mathrm{PS}-\tfrac12$

The A--J weights $\{4,2,4,4,2,1,4,1,1,1\}$ were derived for the RAW PS operator. With
$$
O=\sigma^2=(\mathrm{PS}-\tfrac12)^2 = \mathrm{PS}^2 - \mathrm{PS} + \tfrac14 ,
$$
the full (non-connected) $\langle\sigma\sigma\sigma\sigma\rangle$ mixes $\langle\mathrm{PS}^2\mathrm{PS}^2\rangle$,
$\langle\mathrm{PS}^2\mathrm{PS}\rangle$, $\langle\mathrm{PS}\,\mathrm{PS}\rangle$ and a constant, so the diagram content /
coefficients change. Subtracting $\tfrac12$ from the equal-time $\tau(s,s)$ legs generates SOME cross-terms; verify
whether that already equals the correct $\sigma=\mathrm{PS}-\tfrac12$ four-point or is only an approximation. The
**F-diagram normalization** is where a discrepancy would show. Write the derivation with equations here.
Expect the same qualitative result (scattering $\sim1.0$) as a consistency check.

## Step 2 -- Main study: does $\sigma(x)\sigma(y)$ mix with 0++ via the disconnected same-timeslice piece?

Hypothesis (NM): the 0++ in $\sigma\sigma$ arises from the DISCONNECTED part of $\sigma\sigma$ on the SAME timeslice;
if so the $\sigma\sigma$ channel contains BOTH the 0++ AND the scattering state.

The disconnected same-timeslice operator is exactly the **F diagram** $D'_S(s)\,D'_S(t)$, i.e. the local scalar-
density two-point $\langle\,\delta\!:\!\sigma^2\!:(t)\ \delta\!:\!\sigma^2\!:(0)\,\rangle$ -- the natural FERMIONIC 0++
candidate. TEST: build F's CONNECTED correlator (subtract $\langle D'_S\rangle=-0.00985$ from each leg) and compare
its mass to the 0++ $a_t m_{F^2}=0.616$.
- F carries $0.616$ -> $\sigma\sigma$ holds both 0++ and scattering (hypothesis confirmed).
- F stays noise (as in the raw data: $|\text{sig}|\sim10^{-8}$, effmass $\approx0$) -> the 0++ is essentially
  gluonic, coupling to $\sigma\sigma$ only weakly off-diagonal -> the gluonic-mixing GEVP is essential.

## Step 3 -- FS channel study (parallel to PS)

Repeat Steps 1-2 for the FS scalar
$$
\sigma_\text{FS} = \eta^\dagger\xi - \xi^\dagger(1-D_\text{ov}^\dagger)\eta \quad\text{(indefinite, GW-furnished)} ,
$$
using the FURNISHED perambulator $\tau' = V^\dagger(1-D_\text{ov}^\dagger)D_\text{ov}^{-1}V$ (`tau_gw`), with
$\mathrm{FS}\cdot\mathrm{FS} = G_{10}[\tau] + G_{10}[-\tau']$.
- **OPEN PROBLEM (carried over):** the FS contact is UNSOLVED -- the forward-$\tau'$ equal-time contact
  $\ne 1/2$ (the analytic GW contact applies to the backward/collapsed leg, not the forward $\tau'$). So the FS
  one-point must be determined FROM DATA (measure $\langle D_S^{FS}\rangle$, $\langle D'^{FS}_S\rangle$ with the
  furnished leg, as in `one_point_claude.py` for PS), then subtracted -- the data-driven one-point approach that
  worked for PS ($\langle D_S\rangle$ pure contact, $\langle D'_S\rangle=-0.00985$).
- Then the same connected cumulant $\langle\sigma_\text{FS}^4\rangle_c$ and the same disconnected-same-timeslice
  (F-type) test for a fermionic 0++.
- The coupled program ultimately wants both channels: $\{F^2, \mathrm{FS}{\cdot}\mathrm{FS}, \mathrm{PS}{\cdot}\mathrm{PS},
  \mathrm{FS}{\cdot}\mathrm{PS}\}$.

## Step 4 -- GEVP

Only after Steps 1-3 confirm the operator content:
1. **Lanczos idea** first (asymmetric GEVP + polyiterated / block-Hankel, `asymm_gevp.pdf` right panel;
   Wagman 2406.20009). Plan stub: `gevp_lanczos_impl_plan_claude.md`.
2. **Mixing with gluonic operators** last: coupled $\{F^2,\sigma\sigma\}$ GEVP, the $\sigma\sigma$ level should
   land near the benchmark $a_t m\approx1.0$; resolve whether the 0++ ($0.616$) is a distinct mixed level.
   Compare $\Delta_-$ vs CP {3.30, 3.65, 3.77}.

## Key facts to carry in (units + scales, this ensemble, $a_t m$ lattice)

- Units: redo records are PHYSICAL $m=\mathrm{acosh}/a_t$; distillation is LATTICE $a_t m$. Divide fermion effmass
  by $a_t$ ($\times5$) to compare, or quote physical. axial $\ell{=}1$: $a_t m_A=0.357$ ($=$ redo's $1.787\times a_t$).
- single $\sigma$ $a_t m\sim0.35$ ($\Delta/\Delta_A\approx0.98$); free two-$\sigma$ $2m_\sigma\sim0.67$;
  0++ $a_t m_{F^2}=0.616$ ($\Delta/\Delta_A\approx1.72\sim2.0$); connected $\sigma\sigma\sim1.0$.
- Vacuum coeff in the cumulant $=\tfrac12\langle O\rangle^2$ (verified); free-two-$\sigma$ coeff $2\langle C_S\rangle^2$
  (flavor factor to be derived rigorously).

## Files (src/production)

Analysis: `two_meson_connected_claude.py` (cumulant), `one_point_claude.py`, `two_meson_op1sub_claude.py`,
`diag_fit_claude.py`, `two_meson_summed_claude.py`, `diag_effmass_claude.py`, `distill_contract_claude.py`.
Notes: `sigma_connected_gevp_benchmark_claude.md`, `sigma_onepoint_subtraction_note_claude.md`, this file.
Production sweep running (Nf2 g0.5/1.0/1.5 done @400cfg, on Nf4, ~7d).
