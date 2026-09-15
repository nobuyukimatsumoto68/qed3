# Quantum numbers of the fermion bilinear: $(\Delta_\text{tot},\ell_\text{tot},m_\text{tot})$ by mode pair

Continuum classification of the single-flavor fermion bilinear on $S^2$, looping over every ordered pair of
free Dirac eigenmodes, **validated against the full lattice overlap propagator** $\tau$. Wavefunctions:
`free_wavefunctions_claude.py` (qed3_v2-6.pdf App C.1). Scripts: `sigma_pair_quantum_numbers_claude.py`,
`sigma_ell_validate_propagator_claude.py`.

## NOTE on terminology (two DIFFERENT axes -- do not conflate)

- **This document = the spinor-vertex axis:** $\Gamma=\mathbb 1$ ($\rho=u^*u+d^*d$) vs $\Gamma=\sigma_3$
  ($\rho=u^*u-d^*d$), within ONE flavor. $\Gamma=\mathbb 1$ is the code's $\sigma_{00}$ vertex (identity in
  spin) = the scalar $\sigma$.
- **The paper's PS/FS (qed3int_v3-4.pdf Eq 5.1-5.3) = a FLAVOR axis:** $\sigma_{\rm PS/FS}=\eta^\dagger\xi+
  \xi^\dagger\tilde S\eta$ with $\xi,\eta$ the TWO FLAVORS of two-component overlap fermions, $\tilde S=+1$
  (PS) or $-(1-D_{\rm ov}^\dagger)$ (FS). PS vs FS = the SIGN of the flavor cross-term; $(1-D_{\rm ov}^\dagger)$
  is the O($a$) furnishing. **"PS and FS have the same dimension" is the flavor degeneracy at $m=0$**
  ($\xi,\eta$ share the massless propagator) -- it is NOT the spinor $\Gamma=\mathbb 1/\sigma_3$ split below.
  To connect to PS/FS one redoes the bilinear as the two-flavor $\eta^\dagger\xi$ object; the $\ell$/energy
  content here is flavor-blind and transfers unchanged.

## Setup and the $\iota_3$ treatment

Mode $\psi_{m,n,\iota_3}$: $\lambda=n+|m|+\tfrac12$, $j=n+|m|$. $\iota_3=\pm1$ is the **sign of the Dirac
eigenvalue** ($D\psi=i\iota_3\lambda\psi$); the two $\iota_3$ at fixed $(m,n)$ share the same UPPER component
and have OPPOSITE lower: $\psi_{m,n,\pm}=(u;\pm L)$, $L=\iota_m i(-1)^n\xi(-\iota_m z)$. Meson = antifermion
$\bar a=(m_1,n_1,\iota_1)$ $\times$ fermion $a=(m_2,n_2,\iota_2)$; density $\rho=\psi_{\bar a}^\dagger\Gamma\psi_a$.

- $\Delta_\text{tot}=\lambda_1+\lambda_2$ (physical $\varepsilon_{\lambda_1}+\varepsilon_{\lambda_2}$;
  $\varepsilon_1=E_0=0.189$, $\varepsilon_2=E_1\approx0.26$).
- $m_\text{tot}=m_2-m_1$ (exact).
- $\ell_\text{tot}\in\{|j_1-j_2|,\dots,j_1+j_2\}$, with the **parity selected by $\Gamma$ through $\iota_3$**.

Because $d_{\iota_3}=\iota_3 L$, pairing OPPOSITE $\iota_3$ flips the lower-component sign, turning
$u^*u+d^*d$ into $u^*u-d^*d$. So the two must NOT be summed together (that incoherent-$\iota_3$ sum was the
bug in the first draft, which manufactured a spurious $\Gamma=\mathbb 1$ signal at odd $\ell$).

## Parity split (validated)

| $(\lambda_1,\lambda_2)$ | $\Delta$ | $\Gamma=\mathbb 1$ ($u^*u{+}d^*d$, same-$\iota_3$) | $\Gamma=\sigma_3$ ($u^*u{-}d^*d$, opp-$\iota_3$) |
|---|---|---|---|
| $(1,1)$ | 2 | $\ell=0$ | $\ell=1$ |
| $(1,2)$ | 3 | $\ell=1$ | $\ell=2$ |
| $(2,2)$ | 4 | $\ell=0,2$ | $\ell=1,3$ |

Rule: $\Gamma=\mathbb 1\Rightarrow\ell=|j_1-j_2|,+2,\dots$ ; $\Gamma=\sigma_3\Rightarrow$ the opposite parity.

## Spectra of the two spinor channels

**$\Gamma=\mathbb 1$ (the scalar $\sigma$):**
- $\ell=0$ tower (diagonal shell traces $\sum_{a\in\lambda}\bar\psi_a\psi_a$): $2E_0$, $2E_1$ (the $(2,2)$), $2E_2$.
- $\ell=1$: $E_0+E_1$ (adjacent shells) -- **ground $E_0+E_1$, NOT $2E_0$**.
- $\ell=2$: $D{=}4$ ($2E_1$ / $E_0+E_2$).

**$\Gamma=\sigma_3$:**
- no $\ell=0$.
- $\ell=1$: $2E_0$; $\ell=2$: $E_0+E_1$; $\ell=3$: $2E_1$.

A genuine spinor-parity degeneracy: $\Gamma=\mathbb 1$ at $\ell=0$ and $\Gamma=\sigma_3$ at $\ell=1$ sit at the
same dimension $2E_0$ (validated $0.378$ vs $0.379$). This is a fact about these two spinor channels; it is
NOT the flavor PS$=$FS statement above.

## Full-propagator validation (FREE L1, lattice $\tau$; `sigma_ell_validate_propagator_claude.py`)

$C=-\mathrm{Tr}[\Phi(t)\tau(t,s)\Phi(s)\tau(s,t)]$, $\Phi=V^\dagger[Y_{\ell m}\otimes\Gamma_\text{spin}]V$,
$m$-summed. Effmass at $dt=30$:

| channel | prediction | $\tau$ effmass | |
|---|---|---|---|
| $\Gamma{=}\mathbb 1$, $\ell=0$ | $2E_0=0.378$ | **0.378** | ✓ |
| $\Gamma{=}\mathbb 1$, $\ell=1$ | $E_0{+}E_1\approx0.449$ | **0.443** | ✓ |
| $\Gamma{=}\mathbb 1$, $\ell=2$ | $\approx0.51$ | 0.488 (falling) | ✓ |
| $\Gamma{=}\sigma_3$, $\ell=1$ | $2E_0=0.378$ | **0.379** | ✓ |
| $\Gamma{=}\sigma_3$, $\ell=2$ | $E_0{+}E_1\approx0.449$ | 0.456 | ✓ |
| $\Gamma{=}\sigma_3$, $\ell=0$ | none | noise (no state) | ✓ |

Every channel matches.

**FREE L2** (lattice $\tau$, $2E_0=0.393$; correlators hit solve-noise past $dt\sim15$, so quote $dt=15$):

| channel | prediction | $\tau$ effmass ($dt{=}15$ / $30$) | |
|---|---|---|---|
| $\Gamma{=}\mathbb 1$, $\ell=0$ | $2E_0=0.393$ | **0.395** / 0.394 | ✓ |
| $\Gamma{=}\mathbb 1$, $\ell=1$ | $E_0{+}E_1$ | **0.538** / 0.548 | ✓ |
| $\Gamma{=}\mathbb 1$, $\ell=2$ | $D{=}4$ | 0.624 / 0.577 (falling) | ✓ |
| $\Gamma{=}\sigma_3$, $\ell=1$ | $2E_0=0.393$ | **0.393** / 0.394 | ✓ |
| $\Gamma{=}\sigma_3$, $\ell=2$ | $E_0{+}E_1$ | **0.537** / 0.546 | ✓ |
| $\Gamma{=}\sigma_3$, $\ell=0$ | none | noise (no state) | ✓ |

L2 confirms L1: $2E_0=0.393$ ($\Gamma{=}\mathbb 1$ $\ell0$ $=$ $\Gamma{=}\sigma_3$ $\ell1$); $E_0{+}E_1\approx0.54$
($\Gamma{=}\mathbb 1$ $\ell1$ $=$ $\Gamma{=}\sigma_3$ $\ell2$); the spinor $2E_0$ degeneracy holds at both refinements.

### The $\ell=0$ scalar tower: ground AND excited (the $(2,2)$)

The validation tables above show each channel's GROUND only, so the $\Gamma{=}\mathbb 1$ $\ell=0$ row is
$2E_0$; the **excited $\ell=0$ scalars are buried in that same correlator** and must be isolated by the
diagonal shell projectors $Q_\lambda=V^\dagger\Pi_\lambda V$ (`sigma22_lattice_claude.py`). The $\ell=0$
scalar tower (each state = the shell trace $\sum_{a\in\lambda}\bar\psi_a\psi_a$), effmass ($dt{=}15$ / $30$):

| $\ell=0$ scalar state | operator | L1 | L2 |
|---|---|---|---|
| $2E_0$ (ground $\sigma$) | $Q_{\lambda1}$ ($\sim\sigma_{00}$) | 0.382 / 0.378 | 0.392 / 0.393 |
| $2E_1$ = **$(2,2)$, 1st excited** | $Q_{\lambda2}$ | 0.566 / 0.527 | 0.676 / 0.698 |
| $2E_2$ (2nd excited) | $Q_{\lambda3}$ | 0.687 / 0.622 | 0.823 / 0.798 |

This is the $0^{++}$ single-meson ladder. The $(2,2)=2E_1$ excited scalar -- the state $\sigma_{00}^2$ mixes
with (diagram A) -- is the $Q_{\lambda2}$ row: L1 $\approx0.53$--$0.57$, L2 $\approx0.68$--$0.70$ (the free
lattice $2E_1$ shifts with refinement; still descending on a single config). The shell degeneracies are
$4,8,12$ ($\lambda=1,2,3$) at both L.

## The $O_A$ operator (definition) and the open mechanism question

$$O_A(t)=\sum_{x,y}A(x)\,Y_{00}\;\bar\psi(x,t)\,\tilde\tau(x,y;t)\,\psi(y,t),\qquad
  \tilde\tau \equiv D_\text{ov}^{-1}-\tfrac12\ \ (\text{equal-time, GW contact removed}).$$

- **Bilocal $\ell=0$ scalar**: antifermion at $x$ (weighted by the $\ell=0$ area measure $A(x)Y_{00}$),
  fermion at $y$, connected by the equal-time (same timeslice $t$) contact-subtracted propagator
  $\tilde\tau(x,y)$. Contrast the LOCAL $\sigma_{00}=\sum_x A(x)Y_{00}\,\bar\psi\psi$ (fermion=antifermion point).
- **Mode vertex** (complete distillation basis): $\Phi_A=\Phi_{00}\,\tilde\tau$, with
  $\Phi_{00}=V^\dagger\mathrm{diag}(A\,Y_{00})V$ and $\tilde\tau=\tau(t,t)-\tfrac12 I$ (equal-time perambulator,
  contact-subtracted). In code: `PA[a] = Phi[a] @ (tau[a,a] - 0.5*I)`.
- **$\sigma_3$-parity**: $\Phi_{00}$ is $\sigma_3$-even; $\tilde\tau$ is $\sigma_3$-MIXED,
  $\sigma_3\tilde\tau\sigma_3=-\tilde\tau-1$, i.e. even part $-\tfrac12$ (a c-number) + odd part $D_\text{ov}^{-1}$.
  So $O_A=\Phi_{00}\tilde\tau$ is $\sigma_3$-mixed, hence parity-orthogonal to the whole $\sigma_3$-even
  $\sigma_{00}/Q_\lambda$ tower.

**Open puzzle (to resolve):** the lattice $\langle O_A O_A\rangle$ effmass sits near $2E_1$ ($k{=}1$) before
noise, and $\langle O_A\,\sigma^2\rangle\neq0$, but I have NOT derived this from the continuum. The even part
of $\tilde\tau$ is $-\tfrac12$, so $O_A$ appears to contain $-\tfrac12\,\sigma_{00}$, which should carry a
$k{=}0$ ($2E_0$) component — yet the two-point shows no clean $2E_0$. Resolving this (decompose $O_A$ in the
exact $\sigma_3$-odd $\ell=0$ eigenbasis; track where the even-part $2E_0$ overlap goes) is the next task.
$O_A$'s coupling to $k{=}1$ is at present a LATTICE observation, NOT a continuum result.

## Consequence for the GEVP

The $0^{++}$ (scalar $\ell=0$) sector has distinct states at $2E_0$ ($\sigma_{00}$), $2E_1$ (the $(2,2)$,
$=Q_{\lambda2}$), $2E_2$ ($Q_{\lambda3}$), and the two-meson $2m_\sigma$ ($\sigma_{00}^2$). $E_0+E_1$ is
$\ell\ge1$ and does NOT enter $0^{++}$. Pick one operator per $(\Delta,\ell)$ target.

Refs: qed3_v2-6.pdf App C.1; qed3int_v3-4.pdf Eq 5.1-5.3 (PS/FS flavor operators);
project_cont_prop (reconstructed free propagator).
