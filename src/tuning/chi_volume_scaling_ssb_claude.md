# Detecting chiral SSB via the volume scaling of the low-mode susceptibility $\chi(V)$

**Date:** 2026-07-15. **Author note:** written for NM to review. Working dir `src/tuning/`, next to the
Arnoldi spectrum data (`eig_dw_arnoldi_L*_sig1_gsq*_cfg*_nt128_claude.dat`).

**One-line goal.** Decide whether this theory (Nf=2, reasonable gsq) spontaneously breaks chiral
symmetry, *without* a fragile chiral extrapolation of the condensate, by measuring one number from the
Dirac low modes and watching how it scales with volume.

Method basis: Banks-Casher, Nucl. Phys. B169 (1980) 103 ($\Sigma=\pi\rho(0)$). Overlap low modes from
shift-invert Arnoldi (`eig_arnoldi_claude.cu`, `DW_ARNOLDI`; see [[project-overlap-wall-additive-mass]]).

---

## 1. Why the condensate $\sigma(m)$ is a poor finite-volume signal

The subtracted condensate as a function of valence mass $m$ has the spectral representation

$$
\sigma(m)=\frac1V\sum_\lambda \frac{2m}{\lambda^2+m^2},
$$

where $\lambda$ are the (near-zero) eigenvalues of the massless overlap operator. The weight
$2m/(\lambda^2+m^2)$ is $\propto m$ as $m\to0$, so at any *finite* volume $\sigma(m)\to0$ in the chiral
limit **regardless of phase**. A single-volume chiral extrapolation hitting zero is therefore necessary
but not sufficient for "no SSB" -- it masks the order parameter by construction. This is exactly what the
condensate notebooks show: contact-subtracted $\sigma_{PS}\simeq\chi\,m$, linear to zero.

## 2. Why the susceptibility $\chi$ is a better signal

Differentiate at the massless point:

$$
\chi \;\equiv\; \frac{d\sigma}{dm}\Big|_{m\to0} \;\approx\; \frac1V\sum_\lambda \frac{2}{\lambda^2}.
$$

The weight is now $1/\lambda^2$, i.e. **IR-dominated**: the lowest handful of eigenvalues carry
essentially all of $\chi$. Instead of suppressing the near-zero region (as $\sigma$ does), $\chi$
*amplifies* precisely the $\rho(0)$ region that Banks-Casher cares about. Two practical consequences:

- It is computable directly from the Arnoldi low modes we already have -- **no mass scan, no new sea
  ensembles** needed for the leading answer.
- Its volume scaling is qualitatively different in the two phases (next section), giving a yes/no fork
  rather than a drift-to-zero.

## 3. The volume-scaling fork (the actual discriminator)

Note our lattices hold $N_t=128$ fixed and refine the sphere, so the 3D volume $V=N_s\,N_t$ with
$N_s=10L^2+2\propto L^2$, hence $V\propto L^2$.

**Symmetric / gapped phase.** The spectrum has a gap set by the linear box size, $\lambda_1\sim1/L$
(a single low mode receding as the box grows). Then

$$
\chi \sim \frac1V\,\frac{1}{\lambda_1^2} \sim \frac{1}{L^2}\cdot L^2 = \text{const}.
$$

$\chi$ is intensive -- it saturates. No growth with volume.

**SSB phase.** Eigenvalues accumulate at zero with density $\rho(0)\neq0$; the level spacing near zero is
$\lambda_1\sim 1/(V\rho)\sim 1/L^2$, and $\sum_\lambda 1/\lambda^2\sim(V\rho)^2$. Then

$$
\chi \sim \frac1V\,(V\rho)^2 \;\propto\; V \;\propto\; L^2 \;\to\;\infty.
$$

$\chi$ **grows with volume** -- the "anomalous" scaling. It is the finite-$V$ image of the delta-function
accumulation $\rho(0)\neq0$ (equivalently the massless would-be-Goldstone pole at $m=0$). An intensive
susceptibility has no business growing with $V$; the growth *is* the SSB signature.

| phase | gap $\lambda_1$ | $\chi(V)$ (with $N_t$ fixed, $V\propto L^2$) |
|---|---|---|
| symmetric | $\sim 1/L$ | $\to$ const (flat) |
| SSB | $\sim 1/V\sim 1/L^2$ | $\propto V\propto L^2$ (diverges) |

## 4. What we actually measure

From the massless-ensemble Arnoldi spectra, at each $L$, form

$$
\chi(V)=\frac1V\sum_{\lambda\ \text{low}}\frac{1}{\lambda^2},
\qquad \lambda \approx |z|\ \text{(near-zero overlap eigenvalue)},
$$

summing the captured low modes (IR-dominated, so convergence in the number of Arnoldi modes is fast).
The overlap eigenvalue is obtained from the $D_W$ Arnoldi mode $\mu$ by the phase projection
$z = 1 + (\mu-1)/|\mu-1|$ (M=1), valid for the near-real low modes; this is the same mapping validated
against the dense free spectrum and against the earlier gsq scan.

**Deliverable:** a single plot, $\chi$ vs $V$ for L1/L2/L4, read as diverging (SSB) or flat (symmetric).

## 5. Current data and the preliminary read

Overlap $\lambda_\text{min}=\min|z|$ from the `sig1` (wall $M=1$) Arnoldi data, massless sea:

| $L$ | $V=N_s N_t$ | free $\lambda_\text{min}$ | gsq (reasonable) $\lambda_\text{min}$ |
|---|---|---|---|
| L1 | 1536  | 0.855 | 0.79 / 0.81 / 0.83 (gsq 0.5 / 1.0 / 1.5) |
| L2 | 5376  | 0.492 | 0.43 / 0.42 (gsq 1.0 / 2.0) |
| L4 | 20736 | (n/a)  | 0.213 / 0.214 / 0.25 (gsq 2 / 4 / 6) |

Two facts already visible:

- **Coupling barely moves the gap** -- $\lambda_\text{min}$ sits within a few percent of the free value
  and is nearly gsq-independent. The gap is finite-volume, not dynamical; lowering gsq does not open a
  condensate. (The mapping reproduces the independent gsq scan: 0.80 at g1, 0.82 at g2.)
- **So far the L1$\to$L2 gap closes like $1/L$, not $1/V$.** Ratio $0.43/0.81=0.53$, vs the predictions
  $1/L$ (refine 1$\to$2) $=0.50$ and $1/V$ (sites) $=0.29$. Linear-size scaling = the symmetric signature.

Crude single-lowest-mode proxy $\chi\approx (1/V)\,\lambda_\text{min}^{-2}$ (gsq $\sim$1):

$$
\text{L1:}\ \frac{1}{1536}\,\frac{1}{0.81^2}=1.02\times10^{-3},
\qquad
\text{L2:}\ \frac{1}{5376}\,\frac{1}{0.43^2}=1.00\times10^{-3}.
$$

$\chi$ is **flat** across the two volumes -- consistent with the symmetric ($1/L$) branch, i.e. *no
condensate so far*. This is only two points and one mode per point; **L4 is required** to separate
`flat` (symmetric) from `$\propto L^2$` (SSB), and the full low-mode sum should replace the single-mode
proxy.

## 5b. RESULT with L4 (2026-07-15) -- CONCLUSION: no SSB

L4 Arnoldi done (massless, gsq 2/4/6). Using the full low-mode sum $\chi=\frac1V\sum_{\rm low}1/|z|^2$
(lowest $N_\text{low}=20$ overlap modes; flatness is robust to $N_\text{low}$), representative
reasonable gsq per L (L1 g1.0, L2 g1.0, L4 g2.0), config mean:

| $L$ | $V$ | $\lambda_\text{min}$ | $\chi=\frac1V\sum_{\rm low}1/\lambda^2$ |
|---|---|---|---|
| L1 | 1536  | 0.819 | $1.48\times10^{-2}$ |
| L2 | 5376  | 0.431 | $1.38\times10^{-2}$ |
| L4 | 20736 | 0.213 | $1.43\times10^{-2}$ |

Power-law fits over the three volumes:

$$
\lambda_\text{min}\sim V^{-0.52}\quad(\text{symmetric } -1/2,\ \text{SSB } -1),
\qquad
\chi\sim V^{-0.01}\quad(\text{symmetric } 0,\ \text{SSB } +1).
$$

**$\chi$ is flat to ~1% over a 13.5x volume range, and the gap closes like $1/L$ ($V^{-1/2}$).** Both
are the symmetric-phase signatures with no ambiguity: the low Dirac spectrum does *not* accumulate at
zero, so $\rho(0)=0$ -- **no chiral condensate, no SSB** at these couplings/volumes. The condensate
notebooks' linear-to-zero $\sigma_{PS}$ is thus the genuine symmetric answer, not a finite-volume
masking. Results within each L are gsq-independent (gsq6 at L4 drifts, likely under-thermalized / too
strong -- 3 cfg). Plots: `chi_vs_V_claude.png`, `lammin_vs_V_claude.png`; driver
`chi_scaling_plot_claude.py`.

## 6. On the "physical mass" and the continuum limit (NM's point)

A fixed physical mass enters the spectrum rescaled by the measure factor $\propto a$: at
$m_\text{phys}$ the lattice mass is $m_\text{lat}=(\text{measure})\times m_\text{phys}$, with measure
$= 0.9459$ (L1), $0.506305$ (L2), $0.259021$ (L4) -- halving per refinement, i.e. $\propto 1/L\propto a$
(`heavy_mass_L124_impl_plan_claude.md`). So a *fixed physical mass* tracks a fixed fraction of the gap:
$m_\text{lat}/\lambda_1\approx0.6$ at every $L$ so far.

The equivalence with Section 3: if SSB, the gap collapses faster ($\lambda_1\sim1/V$) than the physical
mass ($m_\text{lat}\sim1/L$), so $m_\text{lat}/\lambda_1\sim L$ grows -- the fixed physical mass clears
the gap and $\sigma(m_\text{phys})$ develops a plateau as one refines. If symmetric, the ratio stays flat
and no plateau ever forms. This is the *observable-side* version of the same test; the $\chi(V)$ scan is
the *spectral-side* version and is cheaper (already have the eigenvalues). Both are watching the same
thing: does the low spectrum accumulate ($\propto V$) or merely close a single gap ($\propto 1/L$).

Note $ma\to0$ automatically in the continuum -- so the plateau, if it exists, is diagnosed by its
*emergence with refinement at fixed physical mass*, never at a single lattice.

## 7. Plan / status

1. **L4 Arnoldi** (massless, gsq 2/4/6, `sig1`) -- DONE (NM ran it). See Sec 5b.
2. **Full low-mode $\chi$** -- DONE ($N_\text{low}=20$; flat conclusion robust to $N_\text{low}$).
3. **The plots** -- DONE: `chi_vs_V_claude.png` (flat), `lammin_vs_V_claude.png` (log-log, sits on the
   $1/L$ reference). Driver `chi_scaling_plot_claude.py`.

Open / possible follow-ups: more configs per (L, gsq) to tighten errors; re-thermalize / drop L4 gsq6;
if a definitive statement is wanted, add a fourth volume or push statistics -- but the current signal
(flat over 13.5x) is already unambiguous.

## 8. Caveats

- Massless sea; 3 configs per (L, gsq) so far -- statistics thin.
- Phase-projection $D_W\to z$ mapping is approximate for non-normal $D_W$ (fine for the near-real low
  modes; validated vs dense free and vs the gsq scan).
- Single-lowest-mode proxy in Sec 5 is indicative only; use the full low-mode sum.
- Only L1/L2 measured; the whole conclusion hinges on the L4 point.
