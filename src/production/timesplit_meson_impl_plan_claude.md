# Time-split mesonic operators for the $0^{++}$ GEVP — implementation plan

## Goal

Enlarge the free $0^{++}$ GEVP with genuinely new interpolators (not GPOF/Hankel recombinations):
place $\bar\psi$ and $\psi$ of the mesonic operators on **neighboring timeslices**, inserting the
one-timestep perambulator $\tau(t,t{+}1)$ as the non-local kernel. Motivation: the equal-time kernel
$\tilde\tau=\tau(t,t)-\tfrac12$ needs the ad-hoc GW **contact** subtraction; the off-diagonal leg
$\tau(t,t{+}1)$ is **contact-free** and still $\sigma_3$-mixed, so it plays $O_A$/$O_{22}$'s role cleanly.

Target: resolve vacuum, $(1,1,1,1)=2E_1\approx0.52$, and the two-meson $2m_\sigma\approx0.756$ (L1)/$0.786$ (L2)
as clean plateaus. Prior finding: $\{1,\sigma^2,O_{22}\}$ gives vac+$2E_1$ but the two-meson does not
plateau (it is heavier than $2E_1$; $\sigma^2$ is dominated by $2E_1$). Time-split gives the extra views.

## Definition of $Q_2$ (the $\lambda=2$ spectral projector) — as requested

$Q_2=V^\dagger\,\Pi_{\lambda=2}\,V$ is the physical shell-$2$ projector, extracted from the perambulator's
**time dependence**. NOTE ON ATTRIBUTION: the distillation vectors $V$, the smearing $\square=VV^\dagger$,
and the perambulator $\tau=V^\dagger D^{-1}V$ are from Peardon et al. 0905.2160. The step below —
diagonalizing the time-averaged perambulator $K(\Delta t)$ and clustering its eigenvalues into shells to
build $Q_2$ — is NOT in 0905.2160; it is a project-specific construction (a one-quark variational read of
the perambulator, `sigma22_lattice_claude.py`), standard GEVP spirit applied to the 1-quark $\tau$.
Concretely (`shell_projector` in the scripts, $\Delta t_\text{dec}=3$):

1. Time-averaged one-step-block perambulator (Nv$\times$Nv):
$$
K(\Delta t)=\frac{1}{N_s}\sum_s \tau(s{+}\Delta t,\,s)\;=\;\sum_{\text{shell }\ell}Q_\ell\,e^{-E_\ell\,\Delta t},
$$
the spectral sum over shells with degenerate energies $E_\ell$ (free single-fermion levels).
2. Diagonalize $K=R\,\mathrm{diag}(\mu)\,R^{-1}$; set $E=-\ln|\mu|/\Delta t$; sort ascending in $E$.
3. **Cluster** by $E$ (gap $<0.02$): the clusters are the shells, with degeneracies
$$
\dim = 4,\,8,\,12,\dots \;=\; \lambda=1\,(j{=}\tfrac12),\ \lambda=2\,(j{=}\tfrac32),\ \lambda=3\,(j{=}\tfrac52),\dots
$$
(Verified: `[4, 8, 12]` at L1.) The **second** cluster (8-fold, $j=\tfrac32$) is $\lambda=2$.
4. Non-orthogonal spectral projector from the right/left eigenvectors of $K$:
$$
Q_2 \;=\; R[:,\,\mathrm{sel}]\;R^{-1}[\mathrm{sel},\,:],\qquad \mathrm{sel}=\text{the }\lambda{=}2\text{ cluster indices}.
$$
$Q_2^2=Q_2$ (idempotent), $\mathrm{tr}\,Q_2=8$. It sharply projects the intermediate quark line onto the
$\lambda=2$ shell. ($Q_2$ is $\sigma_3$-EVEN by itself and decouples from $\sigma^2$; composing with the
non-local $\tilde\tau$ — or here $\tau(t,t{+}1)$ — makes it $\sigma_3$-mixed so it couples.)

## Operators

Spatial $\ell=0$ (Y00) smear: $\Phi(a)=V^\dagger(a)\,\mathrm{diag}(w_{00})\,V(a)$, $w_{00}=\text{dual}\cdot Y_{00}$.

Existing (equal-time), vertex matrices at time $a$:
- $O_A(a)=\bar\psi\,[\Phi\,\tilde\tau]\,\psi$,  vertex $P_A(a)=\Phi(a)\,\tilde\tau(a)$, $\tilde\tau(a)=\tau(a,a)-\tfrac12 I$.
- $O_{22}(a)=\bar\psi\,[Q_2\,\tilde\tau]\,\psi$, vertex $P_Q(a)=Q_2\,\tilde\tau(a)$.
- $\sigma^2=\sigma_{00}^2$ (two-meson), $\mathbb 1$ (vacuum).

New (time-split), $\bar\psi$ at $a$, $\psi$ at $a{+}1$, kernel $=\tau(a,a{+}1)$ (contact-free):
$$
\boxed{\;O_A^{s}(a)=\bar\psi_i(a)\,[\Phi(a)\,\tau(a,a{+}1)]_{ij}\,\psi_j(a{+}1),\qquad
O_{22}^{s}(a)=\bar\psi_i(a)\,[Q_2\,\tau(a,a{+}1)]_{ij}\,\psi_j(a{+}1)\;}
$$
i.e. vertex $M_A(a)=\Phi(a)\,\tau(a,a{+}1)$, $M_Q(a)=Q_2\,\tau(a,a{+}1)$.

**OPEN CHOICE (flagged for NM):** whether $\Phi$ (the $\ell=0$ smear) should be evaluated at $a$ (as
above) or symmetrized between $a$ and $a{+}1$. Default = at $a$; trivial to change.

## Contraction (split$\times$split)

Source meson spans $(s,s{+}1)$, sink $(t,t{+}1)$, lag $dt=t-s$:
$$
\langle O^{s}(t)\,O^{s}(0)\rangle
= \frac1{N_s}\sum_s -\,\mathrm{Tr}\!\big[\,M(t)\,\tau(t{+}1,s{+}1)\,M(s)^\dagger\,\tau(s,t)\,\big],
$$
derived from Wick: $\psi(t{+}1)\!-\!\bar\psi(s{+}1)\to\tau(t{+}1,s{+}1)$, $\psi(s)\!-\!\bar\psi(t)\to\tau(s,t)$.
With $M(a)=K\,\tau(a,a{+}1)$: $-\mathrm{Tr}[K\,\tau(t,t{+}1)\,\tau(t{+}1,s{+}1)\,\tau(s,s{+}1)^\dagger K^\dagger\,\tau(s,t)]$
($K=\Phi$ or $Q_2$, both Hermitian). All $\tau$ available; perambulators do NOT compose (keep separate).

Cross-correlators needed for the full basis $\{1,\sigma^2,O_A,O_{22},O_A^{s},O_{22}^{s}\}$ (6$\times$6):
split$\times$equal-time (one vertex $M$, one vertex $P$; legs $\tau(t{+}1,s)$/$\tau(s,t)$ mixed) and
split$\times\sigma^2$ (one $M$ vertex into the two-loop $\sigma^2$ triangle). One-points $\langle O^{s}\rangle$
for the identity row. These follow the same Wick rules; will enumerate in code with explicit leg times.

## Plan / chunks

1. `Q2` + vertices + split$\times$split correlator; sanity: $O_{22}^{s}$ alone effmass (expect a clean pair).
2. All cross-correlators; assemble 6$\times$6; GEVP (T0 scan); look for vac / $2E_1$ / two-meson plateaus.
3. Compare additional-vs-replacement (empirical, per NM): if split ops dominate, drop the equal-time ones.

Files: `gevp_timesplit_meson_free_claude.py` (new). Data: `data_free/distill_Nv24` (L1), `distill_Nv84` (L2).
Ref: distillation Peardon 0905.2160; GPOF (rejected here) Aubin-Orginos 1010.0202.
