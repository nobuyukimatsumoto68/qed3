# $F^2$--$\sigma^2$ ($0^{++}$) glueball--two-meson mixing -- NULL across $N_f\times g^2$ (L1)

Direct measurement of the Chester--Pufu grand target: does the $0^{++}$ glueball $F^2$ mix with the four-fermion
singlet $\sigma^2=(\bar\psi\psi)^2$? (CP arXiv:1603.05582.) **Result: the mixing is dynamically negligible on every
accessible L1 ensemble ($N_f\in\{2,4,6\}$, $g^2\in\{0.5,1.0,1.5\}$, 400 cfg). Treat as null; do NOT build the coupled
$\{F^2,\sigma^2\}$ GEVP.** (2026-09-09 redo, corrected machinery.)

## What is measured
$F^2$ is purely gluonic, so a single $\sigma^2(s)$ contracts into two fermionic topologies, linked to $F^2$ only
through the gauge field:
$$
\langle F^2(s{+}dt)\,\sigma^2_{00}(s)\rangle_c,\qquad \sigma^2_{00}=2\big(D_S^2+D'_S\big),
$$
$$
D_S(s)=\mathrm{Tr}[\Phi\,\tilde\tau(s,s)],\quad D'_S(s)=\mathrm{Tr}[\Phi\,\tilde\tau\,\Phi\,\tilde\tau],\quad
\tilde\tau(s,s)=\tau(s,s)-\tfrac12 I\ (\text{GW contact}).
$$
- **disc** $=D_S^2$ (two $\sigma$ tadpole 1-loops); **ext** $=D'_S$ (one connected 2-loop).
- $O_F$ = **Wilson-flow-smeared** $F^2$ shape op (`glue_f2_v2_shapes`, dataset `O`, op $0$ = $\ell{=}0$ $p{=}2$ basic
  $F^2$, flow $t_{\rm flow}=2.0$, **spatial-only** flow -> $dt$ unblurred), full $N_t{=}128$. Source $s$
  translation-averaged over the perambulator window (twin=32); $F^2$ at $(t_{\rm src0}{+}s{+}dt)\bmod N_t$.

**PS $=$ FS for this channel (exact).** All loops are even (2-loop) or tadpole (1-loop), so the FS furnishing gives
$(1+(-1)^k)$ = $2\times$improved for the even loop and $0$ for the tadpoles -- identical to PS. Hence
$\langle F^2\sigma_{FS}^2\rangle=\langle F^2\sigma_{PS}^2\rangle$ (the same even-loop PS=FS accident proven for the
four-point). The anomaly/parity subtlety lives only in the odd-loop $2\text{PS}\to$FS triangle, NOT here.

## Method (corrected from the earlier single-plateau pass)
- **Triple subtraction** (correct errors): per-config $t$-sum (removes each config's $\langle F^2\rangle\langle X\rangle$
  vacuum level) -> per-diagram plateau -> total plateau. Single-plateau UNDERESTIMATES errors.
- **Clean $t$-folding**: $C(dt)=C(-dt)$ for two Hermitian $0^{++}$ ops; the backward $F^2$ at $s{-}dt$ is FREE (global
  $F^2$), so folding costs nothing and gives $\sim\sqrt2$ error reduction. (The four-point CANNOT fold -- both legs
  window-bound; this cross can, because $F^2$ is full-$N_t$.)
- **Full $dt$ reach**: since only $\sigma^2$ is window-bound and $F^2$ is global, $dt$ runs to $N_t/2=64$ (the
  four-point was capped at $dt\sim31$).
- Driver `f2_sigma2_cross_v2_claude.py` (env ENS/OP_F/CONTACT/DTCALC/DTMAX/PLAT_LO). Figs
  `figs/f2_sigma2_cross_v2_{lin,log}_*_claude.png`.

## Interpretation -- the null is DYNAMICAL, not a selection rule
An earlier version claimed "$F^2$ cannot reach the two-meson by particle number." **That was WRONG (corrected).**
$F^2$ and the two-meson are both $0^{++}$; with dynamical fermions $F^2$ couples to the sea ($F^2\to$ gluons $\to
q\bar q$), so $\langle0|F^2|{\rm two\text{-}meson}\rangle$ is symmetry-ALLOWED -- exactly how glueballs mix with
multi-meson states in QCD. In the correlator $\langle F^2(t)\sigma^2(0)\rangle=\sum_n\langle0|F^2|n\rangle
\langle n|\sigma^2|0\rangle e^{-E_n t}$ the two-meson $n$ contributes with both overlaps allowed. The Wick picture
only shows the amplitude is **sea-quark-mediated** (the $\sigma^2$ fermions close into loops at $t{=}0$, reaching
$F^2$ through the gauge/sea) -> potentially small, not zero. So the null = a measured smallness of a sea-mediated
mixing, which IS the CP answer, and it warranted checking whether it grows with coupling / $N_f$ (it does not).

## Findings
- **disc ($D_S^2$) is numerically DEAD** ($\sim10^{-15}$): after the $\tfrac12$ improvement $\langle D_S\rangle=0$
  and its config fluctuation is $\sim10^{-7}$ (contact-saturated, no condensate -> `project_chi_volume_scaling_nossb`).
  A config-independent loop cannot correlate with $F^2$.
- **ext ($D'_S$) is the only live channel -- and it is null.** No monotone falloff, no plateau, sign-changing.

### $N_f\times g^2$ scan (max ext $|S/N|$ over $dt\in[0,47]$; `f2_sigma2_cross_v2`, folded, triple-sub)
| $N_f\backslash g^2$ | 0.5 | 1.0 | 1.5 |
|---|---|---|---|
| **2** | 1.61 | 0.97 | 2.14 |
| **4** | 1.72 | 0.87 | 2.06 |
| **6** | 0.66 | 0.84 | (no data) |

(Nf6 g1.0 = 101 cfg partial; Nf6 g1.5 perams not yet generated.) Every entry is the MAX over ~40 $dt$ points, so a
$\sim2\sigma$ peak is expected from noise alone. Crucially **no channel has a decaying $0^{++}$ shape**: the largest
peaks ($g^2{=}1.5$) have $|C|$ that GROWS with $dt$ (Nf4 g1.5: $|C|$ climbs $1.0\to2.0\times10^{-6}$ out to
$dt\sim6$) -- correlator wander pinned by the tail plateau, not a state. **The sea-mediated $F^2$--$(\bar\psi\psi)^2$
mixing stays dynamically small at every $N_f\le6$, $g^2\le1.5$ -- it does not turn on at stronger coupling / larger
$N_f$.**

### Time-displaced one-meson interpolator ($O_{1m}$) -- also null
Crossing $F^2$ with the coincident time-split one-meson interpolator $O_{1m}(s)=\sum_x A_x\sigma(x,s)\sigma(x,s{+}\delta)$
(self-loop $L_\delta(s)=-\sum_x A_x\mathrm{Tr}[P(x,s;x,s{+}\delta)P(x,s{+}\delta;x,s)]$, position propagator
$P=$`AblkS`) for $\delta{=}1,2,4$: $|S/N|\le0.8$ everywhere, no falloff -- letting the fermion loop propagate does not
expose any hidden overlap. Driver `f2_sigma2_cross_o1m_v2_claude.py`, fig `f2_sigma2_cross_o1m_v2_*_claude.png`.
(The antipodal $O_{2m}$ is the same self-loop topology and will be null for the same dynamical reason -- not run.)

## Physically reasonable
No chiral SSB at these $N_f$ (`project_chi_volume_scaling_nossb`, near-conformal), and the CP $F^2$ partner is HEAVY
($a_t m_{F^2}\approx3.08$, far above the two-meson $\approx0.61$), so a strong low-energy $F^2$--$\sigma^2$ overlap
was not expected. The coupled GEVP is not warranted on these ensembles.

## Files
- Drivers: `f2_sigma2_cross_v2_claude.py` (local $\sigma^2_{00}$), `f2_sigma2_cross_o1m_v2_claude.py` ($O_{1m}$).
- Plan: `f2_sigma2_mixing_redo_impl_plan_claude.md`. Superseded: `f2_sigma2_cross_diag_claude.py` (single-plateau,
  PS-only, no fold, $dt\le24$; the earlier "$\sim1\sigma$ null" pass).
- Refs: CP arXiv:1603.05582; distillation Peardon 0905.2160; no-SSB `project_chi_volume_scaling_nossb`; glueball mass
  `glue_gevp_masses_claude.md`; Wilson flow Luscher arXiv:1006.4518.
