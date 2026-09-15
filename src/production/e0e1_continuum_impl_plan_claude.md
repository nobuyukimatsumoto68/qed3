# Continuum enumeration: which bilinear combinations create the $E_0+E_1$ meson

## Goal

Using the analytic free $S^2$ Dirac eigenspinors (`free_wavefunctions_claude.py`, qed3_v2-6.pdf
App C.1), find which mode combinations $(m,n,\iota_3)$ a local meson bilinear must couple in order
to create the single-meson state at energy $E_0+E_1$ (one member in the $\lambda=1$ shell, the other
in $\lambda=2$), and NOT the lower $2E_0$ ground.

This is the pure-continuum precursor to the lattice $(2,2)$/excited-meson operator. Frame transport
$R(x)$ (geometry/connection) is deferred: here everything is in the continuum $(\hat\theta,\hat\phi)$
frame of the analytic spinors.

## Spectrum / notation

Spatial Dirac eigenvalue (Eq C.17): $\lambda = n+|m|+\tfrac12$, so total angular momentum
$j=\lambda-\tfrac12=n+|m|$.

- $\lambda=1$ shell ($j=\tfrac12$): 4 modes $(m=\pm\tfrac12,\,n=0,\,\iota_3=\pm1)$.
- $\lambda=2$ shell ($j=\tfrac32$): 8 modes $(m=\pm\tfrac12,\pm\tfrac32;\ n=1$ for $|m|=\tfrac12$, $n=0$ for $|m|=\tfrac32;\ \iota_3=\pm1)$.

Single-fermion energy in continuum units $E\propto\lambda$ (I set $E_k=\lambda_k=n+|m|+\tfrac12$, so
$\lambda=1\to E=1$, $\lambda=2\to E=2$). This is NOT the lattice dispersion (which compresses
$E_1/E_0=0.26/0.189=1.38$ instead of $2$); only the mode-SELECTION logic is being checked here, so
$\lambda$-units are enough. Meson energies in these units:

$$
2E_0 \leftrightarrow 2,\qquad E_0+E_1 \leftrightarrow 3,\qquad 2E_1 \leftrightarrow 4 .
$$

## Free meson correlator in the mode basis

For a local bilinear $O=\sum_x A(x)\,\bar\psi(x)\,\Gamma(x)\,\psi(x)$ with
$\Gamma(x)=Y_{\ell M}(\hat x)\,\Sigma$ ($\Sigma$ a $2\times2$ spin matrix), the vertex matrix element
between modes $a,b$ is

$$
\Gamma^{\ell M}_{ab} \;=\; \int d\Omega\; \psi_a^\dagger(\hat x)\,\Sigma\,\psi_b(\hat x)\;Y_{\ell M}(\hat x).
$$

The free two-point (single fermion loop) is

$$
C(t) \;=\; \sum_{a,b} P_\ell(a,b)\; e^{-(E_a+E_b)\,t},
\qquad P_\ell(a,b)=\sum_{M=-\ell}^{\ell}\big|\Gamma^{\ell M}_{ab}\big|^2 ,
$$

so the effective mass plateaus at the SMALLEST $E_a+E_b$ over pairs with $P_\ell(a,b)\neq0$.
$P_\ell$ (summed over $M$) is the rotationally-invariant $\ell$-channel strength of the pair.

## Selection rules (what to expect)

Azimuthal: $\Gamma^{\ell M}_{ab}\neq0$ needs $m_a=m_b+M$ (from $e^{i(m_b+M-m_a)\phi}$), and spin matrix
$\Sigma$ can shift the effective $m$ by the ladder structure of the two spinor components.

Angular-momentum triangle (dominant rule): coupling $j_a\otimes j_b$ contains $\ell$ iff
$|j_a-j_b|\le\ell\le j_a+j_b$.

- $\ell=0,1$ vertices connect $\lambda=1\!\leftrightarrow\!\lambda=1$ ($j=\tfrac12\otimes\tfrac12=0\oplus1$)
  $\Rightarrow$ lowest surviving pair is $2E_0$. So a naive $Y_{1M}$ scalar sees $2E_0$, not $E_0+E_1$.
- $\ell=2$ vertex: $\tfrac12\otimes\tfrac12$ has max $\ell=1<2$, so the $\lambda=1\!\times\!\lambda=1$
  ($2E_0$) channel is FORBIDDEN. The lowest allowed is $\tfrac12\otimes\tfrac32$ ($1\le2\le2$)
  $=\lambda=1\times\lambda=2=E_0+E_1$. **So an $\ell=2$ density has its ground at $E_0+E_1$.**

The role of $\Sigma$ ($\sigma_0,\sigma_3,\sigma_\pm$): the two-component $\iota_3$ (chirality-like)
structure can additionally kill the shell-diagonal ($2E_0$) piece of an $\ell=1$ vertex, which would
promote an $\ell=1$ operator's ground to $E_0+E_1$ as well. This is exactly the "combination in the
bilinear" being checked. The script enumerates all $\Sigma\in\{\sigma_0,\sigma_1,\sigma_2,\sigma_3\}$
and $\ell\in\{0,1,2\}$ and reports, per $(\Sigma,\ell)$: the minimal surviving $E_a+E_b$, whether the
$2E_0$ channel survives, and the explicit contributing $(m,n,\iota_3)$ pairs at $E_0+E_1$.

## Files

- `e0e1_modes_continuum_claude.py` — enumeration + continuum effmass, uses
  `free_wavefunctions_claude.py`.

## Refs

qed3_v2-6.pdf App C.1 (Eq C.10, C.17, C.18); Peardon distillation 0905.2160 (basis only, not used here).
