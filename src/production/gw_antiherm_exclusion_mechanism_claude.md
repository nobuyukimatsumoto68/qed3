# Why $\sigma^2$ does not couple to a single meson: GW anti-hermiticity of the normal-ordered propagator

Author: "Fin: Two-meson" session, with the "{1,1,1,1}" agent. 2026-09-17/18.
Supersedes the earlier "$\sigma_3$-hermiticity selection rule" attribution (`sigma2_single_meson_exclusion_claude.md`):
the interacting mechanism is **Ginsparg-Wilson (GW) anti-hermiticity of the normal-ordered propagator**, not
$\sigma_3$-hermiticity. This note gives the identity, why it holds **configuration-by-configuration** in the
interacting theory, the role of time-reflection ($T$) symmetry, the free-vs-interacting difference, and how it
unifies the truncation-robustness of the cross with the $C_{22}$ "leak" (a contact-subtraction bug, not physics).

## 1. Setup

Zero-momentum $\ell=0$ scalar density $\sigma_{00}(t) = \sum_x A_x Y_{00}\,\bar\psi_x(t)\psi_x(t)$, with the
Hermitian vertex $\Phi = \mathrm{diag}(A_x Y_{00})$ (real, diagonal). The four-fermion $0^{++}$ operator is
$\sigma^2_{00}$. The object of interest is the cross correlator
$$
C_{12}(t) = \langle \sigma_{00}(t)\,\sigma^2_{00}(0)\rangle ,
$$
whose **connected** part is a single closed fermion loop through three vertices (the triangle: one sink vertex at
time $t$, two source vertices at $0$). If $C_{12}=0$ then $\sigma^2$ has no single-meson ($\bar\psi\psi$-tower)
overlap and the $\sigma^2$ correlator is a clean two-meson interpolator. This is the exclusion.

The massless overlap propagator is $G = D_\text{ov}^{-1}$. The GW relation (qed3_v2-6.pdf, Eq. IV.17) is
$$
D_\text{ov}^{-1} + D_\text{ov}^{-\dagger} = 1 ,
$$
equivalently $\mathrm{Re}\,G = \tfrac12\,\mathbb 1$ as an operator: the equal-point diagonal of $G$ is exactly
$1/2$ (the ultralocal GW contact). Subtracting it is the normal ordering $\sigma \to \sigma - \tfrac12$.

## 2. The normal-ordered propagator is anti-hermitian (GW)

Define
$$
M := D_\text{ov}^{-1} - \tfrac12\,\mathbb 1 = \tfrac12\left(D_\text{ov}^{-1} - D_\text{ov}^{-\dagger}\right).
$$
Then, using GW,
$$
M^\dagger = D_\text{ov}^{-\dagger} - \tfrac12 = (1 - D_\text{ov}^{-1}) - \tfrac12 = -\left(D_\text{ov}^{-1} - \tfrac12\right) = -M .
$$
So $M$ is **anti-hermitian**, $M^\dagger = -M$. This is a per-configuration operator identity (GW holds on every
gauge background), and it is the whole engine. In the distillation basis it is the $T2b$ identity
$\bar\tau = \delta - \tau$ (backward $=$ $\delta$ minus forward perambulator), verified at build time.

Note this is NOT $\sigma_3$-hermiticity. $\sigma_3$-hermiticity would be $\Gamma D_\text{ov}\Gamma = D_\text{ov}^\dagger$
for some $\Gamma$ (a similarity to the adjoint); GW anti-hermiticity is the additive statement $G + G^\dagger = 1$.
The interacting theory has the latter but not the former (Section 5).

## 3. Odd loops are purely imaginary; $T$-symmetry lets us take the real part

The connected triangle is a closed loop of three $M$'s with Hermitian vertices $\Phi$. Absorbing
$\Phi = \Phi^{1/2}\Phi^{1/2}$ into similarity factors, the loop value is
$$
T = \mathrm{Tr}\big[(\Phi M)^3\big]\ \text{-type} \;=\; \mathrm{Tr}\big[\widetilde M^3\big],
\qquad \widetilde M := \Phi^{1/2} M\,\Phi^{1/2},\quad \widetilde M^\dagger = -\widetilde M .
$$
For an anti-hermitian $\widetilde M$ and a loop of $n$ factors,
$$
T^* = \mathrm{Tr}\big[(\widetilde M^\dagger)^n\big] = (-1)^n\,\mathrm{Tr}\big[\widetilde M^n\big] = (-1)^n\,T .
$$
- **Odd loop** ($n=3$, the cross triangle): $T^* = -T$, so $T$ is **purely imaginary**, $\mathrm{Re}\,T = 0$.
- **Even loop** ($n=4$, the $\sigma^2$-$\sigma^2$ two-meson): $T^* = +T$, so $T$ is **real** and generically nonzero
  — the two-meson survives, as it must.

Why the real part is the physical answer (and why this is config-by-config):

- **$T$-symmetry.** The gauge action and measure are real, so the ensemble is invariant under $U \to U^*$
  (a time-reflection / CP). Under it the correlator's imaginary part is odd, hence
  $\langle \sigma_{00}\sigma^2\rangle$ is **real**: $\mathrm{Im}$ cancels in the average.
- Because the true expectation is real and the per-config $\mathrm{Im}$ is $T$-odd, replacing the estimator by its
  **real part is exact** (and variance-reducing), not an approximation. This is the `.real` taken in the contraction
  code.
- Once only $\mathrm{Re}$ matters, GW ($M^\dagger = -M$) forces $\mathrm{Re}\,T = 0$ on **every configuration**. No
  gauge average is needed for the (real part of the) result; the exclusion is exact configuration-by-configuration.

Two symmetries, distinct jobs: **$T$ drops the imaginary part; GW kills the remaining real part.** GW is a per-config
operator identity, so the zero is exact per config.

## 4. The two orientations are equal, not conjugate

The triangle has two loop orientations, $T_1$ (sink $\to$ src$_1$ $\to$ src$_2$) and $T_2$ (sink $\to$ src$_2$ $\to$
src$_1$). They are related by relabeling the two symmetric source vertices, so they are **equal as complex numbers**,
not complex conjugates:
$$
T_1 = T_2 \quad(\text{verified } |T_1-T_2|\sim 10^{-21},\ |T_1-\overline{T_2}|\sim 10^{-7}).
$$
Hence $C_{12} = T_1 + T_2 = 2\,T_1$, still purely-imaginary-dominated. Reality is **not** produced by summing the two
orientations; it is the explicit `.real`. (A common "term $+$ conjugate $= 2\,\mathrm{Re}$" pattern does NOT apply
here.)

## 5. Free vs interacting: the imaginary part

Writing $M = iH$ with $H = -iM$ Hermitian, $\mathrm{eig}(M) = i\mu_k$ ($\mu_k$ real), and
$$
\mathrm{Im}\,\mathrm{Tr}(M^3) = -\,\mathrm{Tr}(H^3) = -\sum_k \mu_k^3 .
$$
$\mathrm{Re} = 0$ always (Section 2). $\mathrm{Im}$ vanishes **iff** the spectrum $\{\mu_k\}$ is symmetric under
$\mu \to -\mu$, i.e. there is an extra reflection $\Gamma M \Gamma^{-1} = -M$ — equivalently a
$\sigma_3$/$\gamma$-hermiticity $\Gamma D_\text{ov}\Gamma = D_\text{ov}^\dagger$.

- **Free theory:** the reflection symmetry holds, the spectrum is $\pm$-paired, $\sum \mu_k^3 = 0$, so
  $\mathrm{Im} = 0$ as well. The free triangle vanishes entirely (indeed site-by-site).
- **Interacting theory:** the gauge field **breaks** the reflection ($\sum \mu_k \ne 0$, spectrum not $\pm$-paired),
  so $\mathrm{Im} \ne 0$. But this never matters: $T$-symmetry already licenses dropping $\mathrm{Im}$, and GW kills
  $\mathrm{Re}$. This is exactly the sense in which "we do not have $\sigma_3$-hermiticity" interacting — we lose the
  identity that would have killed $\mathrm{Im}$, but not the one (GW) that kills $\mathrm{Re}$.

## 6. Truncation robustness, and the $C_{22}$ "leak" that was a contact bug

**Cross is truncation-robust.** Under distillation truncation the propagator is $P G P$ with
$P = V V^\dagger$ a Hermitian projector onto the kept modes. With the contact subtracted **consistently in mode space**
the truncated normal-ordered propagator is $P M P$, and
$$
(P M P)^\dagger = P M^\dagger P = -\,P M P ,
$$
still anti-hermitian. So $\mathrm{Re}\,T = 0$ survives truncation: the cross exclusion holds at any $N_v$ (verified to
$N_v = 6$).

**The $C_{22}$ collapse was a bug, not physics.** The production contraction subtracted the contact as $-\tfrac12\,\mathbb 1$
on the **full** $2N_s$ space, i.e. it used $P G P - \tfrac12\,\mathbb 1$. That is **not** anti-hermitian unless $P=\mathbb 1$:
$$
(P G P - \tfrac12)^\dagger = P G^\dagger P - \tfrac12 = P(1-G)P - \tfrac12 = P - PGP - \tfrac12 \ne -(PGP - \tfrac12)\ \ (P\ne\mathbb 1).
$$
The mismatch $-\tfrac12(\,\mathbb 1 - P)$ reintroduces the single-$\sigma$ tadpole $D_S = \mathrm{Tr}[\Phi M]$ (which grows
linearly in the number of removed modes), feeding single-meson disconnected diagrams into $\langle\sigma^2\sigma^2\rangle$
and collapsing the $C_{22}$ ground to $m_\text{PS}$. The fix is to subtract the contact **in mode space**,
$\tau(t,t) \to \tau(t,t) - \tfrac12\,\mathbb 1_{N_v}$ before $U(\cdots)U^\dagger$ (equivalently $-\tfrac12 P$), which keeps
$P M P$ anti-hermitian; the GW contact is diagonal $=1/2$ in the distillation basis (since $\mathrm{Re}\,G = \tfrac12$ in
any orthonormal basis), so mode-space subtraction is exact under slicing. Flag `MODE_CONTACT=1` in
`fs_gevp_point_claude.py`. This needs **no** larger basis ($N_v=84$) and no basis change.

## 7. Numerical verification (free and interacting $L1$, single config)

Anti-hermiticity and the purely-imaginary triangle (`antiherm_check_claude.py`):

| quantity | free | interacting (Nf2 g0.5) |
|---|---|---|
| $\lVert M + M^\dagger\rVert/\lVert M\rVert$ (eq-time block) | $5\times10^{-10}$ | $1.3\times10^{-4}$ (= solve tol) |
| off-diag GW $\lVert \tau(a,b)^\dagger + \tau(b,a)\rVert/\lVert\tau\rVert$ | $10^{-7}$ | $\sim10^{-4}$ |
| triangle $T_1$: $\mathrm{Re}$ | $\sim10^{-12}$ | $\sim10^{-10}$ (zero) |
| triangle $T_1$: $\mathrm{Im}$ | $\sim10^{-18}$ (zero) | $\sim10^{-7}$ (nonzero) |

Spectrum reflection (`M_spectrum_check_claude.py`, $M=\tau(0,0)-\tfrac12$):

| quantity | free | interacting |
|---|---|---|
| $\sum_k \mu_k$ | $8\times10^{-15}$ | $-3.3\times10^{-3}$ |
| $\sum_k \mu_k^3 = -\mathrm{Im}\,\mathrm{Tr}(M^3)$ | $4\times10^{-16}$ (zero) | $-8.0\times10^{-5}$ (nonzero) |
| reflection asym $\max_k\lvert\mu_k+\mu_{N-1-k}\rvert$ | $8\times10^{-15}$ | $4.7\times10^{-3}$ |
| $\mathrm{Tr}(M^3)$ | $9\!\times\!10^{-12} - 4\!\times\!10^{-16} i$ | $-1\!\times\!10^{-7} + 8.0\!\times\!10^{-5} i$ |

Contact / tadpole (`tadpole_trunc_check_claude.py`), $D_S=\mathrm{Tr}[\Phi(\tau_{tt}-\tfrac12)]$ vs contact space:

| $N_v$ kept | $D_S$ (position-space $-\tfrac12\mathbb 1$, the bug) | $D_S$ (mode-space $-\tfrac12\mathbb 1_{N_v}$, the fix) |
|---|---|---|
| all (24) | $\sim10^{-10}$ | $\sim10^{-10}$ |
| 18 | $-0.886$ | $\sim10^{-10}$ |
| 12 | $-1.772$ | $\sim10^{-10}$ |
| 6 | $-2.659$ | $\sim10^{-11}$ |

Cross robustness (`sigma2_cross_trunc_check_claude.py`): $|C_{12}|/|C_{11}| \sim 10^{-7}$ at $N_v = 24,18,12,6$.

$C_{22}$ leak vs fix (`sigma2_leak_check_claude.py`, free $L1$, ground effmass at $t=15$; $m_\sigma=0.378$, $2m_\sigma=0.756$):
complete $N_v24$ $\to 0.659$; truncated $N_v12$ position-contact $\to 0.385$ ($=m_\sigma$, the leak); truncated $N_v12$
mode-contact $\to 0.659$ (recovers the complete two-meson curve).

## 8. Consequences

- $\sigma^2_{PS}$ is a clean $(\bar\psi\psi)^2$ / two-meson interpolator: it has zero single-meson overlap, exactly and
  config-by-config, at any lattice and any (mode-space-consistent) truncation. The mechanism is GW $+$ $T$, not
  $\sigma_3$-hermiticity, not completeness, not a basis symmetry.
- The $L2$ $N_v=24$ (truncated) "leak" is fixed by `MODE_CONTACT=1` on the existing perambulators — no $N_v=84$, and
  the symmetrized-basis regen is not required for it (it remains desirable only for cleanliness).
- Even/odd rule: an $n$-vertex zero-momentum $\ell=0$ loop of the normal-ordered propagator has
  $T^* = (-1)^n T$. Odd $\Rightarrow$ $\mathrm{Re}=0$ (single-meson-excluding). Even $\Rightarrow$ real, survives. This
  is why the 3-point cross vanishes while the 4-point two-meson does not.

## 9. Files

- `antiherm_check_claude.py` -- $M^\dagger=-M$, off-diagonal GW, triangle $\mathrm{Re}$ vs $\mathrm{Im}$.
- `M_spectrum_check_claude.py` -- $\mu\to-\mu$ reflection, $\sum\mu^3$, free vs interacting.
- `triangle_split_claude.py` -- the two orientations $T_1,T_2$ separately (equal, not conjugate) + site-by-site tensor.
- `tadpole_trunc_check_claude.py` -- $D_S$ revival, position vs mode contact.
- `sigma2_cross_trunc_check_claude.py` -- $|C_{12}|/|C_{11}|$ vs $N_v$ (truncation-robust).
- `sigma2_leak_check_claude.py` -- $C_{11}$ / $C_{22}$-only / GEVP; `NVKEEP` truncation, `MODE_CONTACT` fix.
- `fs_gevp_point_claude.py` -- `MODE_CONTACT` flag on `AblkS`; `distill_contract_claude.py` -- `NVKEEP` in the loader.
- Companions: `sigma2_single_meson_exclusion_claude.md` ({1,1,1,1}, to be updated to this mechanism),
  `fs_gw_collapse_v_agent_two_meson_claude.md`. Source: GW Eq. IV.17, `qed3_v2-6.pdf`.
