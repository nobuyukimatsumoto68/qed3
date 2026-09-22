# Flavor SSB signatures in radial quantization ($S^2\times\mathbb{R}$)

Sources: P. Hasenfratz, hep-lat/9802007 (exact GW Ward identities, no tuning/mixing/current
renormalization); S. Chandrasekharan, hep-lat/9805015 (GW chiral symmetry, lattice condensate and
susceptibility identities); H. Leutwyler, PLB 189 (1987) 197 and P. Hasenfratz, F. Niedermayer,
Z. Phys. B 92 (1993) 91, hep-lat/9212022 (rotor / $\delta$-regime spectrum of Goldstone bosons in a
finite spatial volume); S. Chester, S. Pufu, arXiv:1603.03771 (large-$N$ QED3 bilinear dimensions);
R. Pisarski, PRD 29 (1984) 2423 and T. Appelquist, D. Nash, L. Wijewardhana, PRL 60 (1988) 2575 (QED3
flavor SSB pattern).

## 1. Setup

$N$ two-component flavors ($N = 2N_f$ if $N_f$ counts the $(\xi,\eta)$ pairs), global flavor group
$G = SU(N)$. The parity-even mass $m\,\bar\psi\,\mathrm{diag}(1_{N/2},-1_{N/2})\psi$ is in the
**adjoint** of $G$; the $G$-singlet $\bar\psi\psi$ is parity-odd. The candidate SSB pattern is

$$
SU(N) \;\to\; H = S\big(U(N/2)\times U(N/2)\big), \qquad \dim(G/H) = N^2/2 \ \text{Goldstones},
$$

order parameter $\Sigma = \langle\bar\psi\,\mathrm{diag}(1,-1)\psi\rangle$, Goldstone decay constant
$F^2$ (mass dimension 1). QED3 has one scale, so $F^2 = c_F\, g^2$, $\Sigma = c_\Sigma\, g^4$ with
$c_F, c_\Sigma$ possibly exponentially small (Miransky-type onset).

Radial dictionary (only exact in a CFT): the Hamiltonian on $S^2$ of radius $R$ has

$$
E_n(R)\,R = \Delta_n \quad (\text{CFT}),
$$

$R$-independent up to $R^{-\omega}$ from the leading irrelevant operator; with the coupling as the
only scale, "$R$" enters as $g^2 R$, and the QED3 IR fixed point is $g^2R\to\infty$. In the SSB case
the theory is not conformal and $E_n(R)R$ is only an **effective** dimension; its $R$-dependence is
the signal.

## 2. Unbroken (CFT) spectrum

- Unique $G$-singlet ground state.
- Adjoint scalar $\bar\psi T^a\psi$ ($\ell=0$ ground state of the connected scalar-density channel):

$$
\Delta_{\rm adj} = 2 - \frac{64}{3\pi^2 N} + O(1/N^2)
\quad(N=4:\ 1.46,\ \ N=8:\ 1.73,\ \ N=12:\ 1.82),
$$

  singlet (parity-odd) $\Delta_{\rm sing} = 2 + \frac{128}{3\pi^2 N}$.
- All $SU(N)$ currents (vector and "axial" alike) conserved: $\Delta_J = 2$ exactly, tower on $S^2$
  $E_\ell R = \ell + 1$ for $\ell = 1,2,3,\dots$ ($\ell=0$ is the charge, annihilates the vacuum).
  Ratios $\ell 2/\ell 1 = 3/2$, $\ell 3/\ell 1 = 2$.
- Descendants of a scalar primary: $E R = \Delta + \ell$.

## 3. Broken (SSB) spectrum on the finite sphere

Finite volume restores the symmetry through the zero mode; leading order in the IR-free
$G/H$ sigma model (Leutwyler; Hasenfratz-Niedermayer), moment of inertia
$\Theta = F^2\cdot 4\pi R^2$ (up to the $O(1)$ normalization of $F$), valid for $F^2 R \gg 1$:

**(a) Rotor tower** (zero mode of the Goldstone field). States are $G$-irreps $\lambda$ that contain an
$H$-singlet; the lowest is the **adjoint**, i.e. the quantum numbers of the order parameter
$\bar\psi T^a\psi$:

$$
E_\lambda = \frac{C_2(\lambda)}{2\Theta},\qquad
E_{\rm adj}\,R = \frac{N}{8\pi F^2 R}\ \propto\ \frac{1}{g^2 R}\ \to 0 .
$$

This is the finite-$R$ form of "degenerate vacua": the effective dimension of the adjoint scalar
channel flows to **zero** as $1/R$ instead of saturating at $\Delta_{\rm adj}\in(1.4,2)$.
NLO corrections are $O\!\big(1/(F^2 R)\big)$ relative.

**(b) One-Goldstone states with $j\ge1$ on $S^2$.** The shift symmetry forbids the curvature coupling
$\xi\mathcal{R}\pi^2$, so the Goldstone is minimally coupled:

$$
E_j\,R = \sqrt{j(j+1)} = 1.414,\ 2.449,\ 3.464,\dots \qquad (\text{plus rotor offsets } \propto 1/R),
$$

vs. $\ell+1 = 2,3,4$ for the conserved-current tower. The coset ("axial") currents create these
states, $\langle 0|A^a_\mu|\pi\rangle = F p_\mu$, so the axial $\ell$-tower moves from
$E_\ell R = \ell+1$ to $\sqrt{\ell(\ell+1)}$: ratios $\ell2/\ell1 = \sqrt3 = 1.73$ (vs $1.5$),
$\ell3/\ell1 = \sqrt6 = 2.45$ (vs $2$). The unbroken ($H$) currents stay at $\ell+1$.
The $\ell=0$ current channel is null in **both** cases (the finite-volume ground state is a
$G$-singlet, $Q^a|0\rangle = 0$), so dropping $\ell=0$ from the current towers loses nothing.

**(c) Non-Goldstone states** (parity-odd singlet "sigma", non-Goldstone mesons, glue $0^{++}$): gapped
at a physical mass $M \propto g^2$, so

$$
E R = M R \ \propto\ g^2 R \ \to \infty \quad\text{linearly}.
$$

Summary of the $R$ (equivalently $g^2R$) dependence of $ER$:

| channel                          | CFT (no SSB)              | SSB                                  |
|----------------------------------|---------------------------|--------------------------------------|
| adjoint scalar, $\ell=0$         | $\to \Delta_{\rm adj}$    | $\propto 1/(g^2R) \to 0$ (rotor)     |
| coset (axial) current, $\ell\ge1$| $\ell+1$                  | $\sqrt{\ell(\ell+1)}$ (Goldstone)    |
| unbroken current, $\ell\ge1$     | $\ell+1$                  | $\ell+1$                             |
| massive states                   | const                     | $\propto g^2R \to\infty$             |

With $m\neq0$: $M_\pi^2 = 2m\Sigma/F^2$; the rotor picture holds for $M_\pi R\ll1$, and for
$M_\pi R\gg1$ the Goldstone levels become $E R = R\sqrt{M_\pi^2 + j(j+1)/R^2}$.

## 4. Link to the GW Ward identity (Hasenfratz, Chandrasekharan)

With overlap fermions the flavor Ward identity is exact at finite $a$ with the modified densities
$P^a = \bar\psi\gamma_5 T^a(1-\tfrac{a}{2}D)\psi$, $S = \bar\psi\,\mathrm{diag}(1,-1)(1-\tfrac{a}{2}D)\psi$,
no $Z_A$, no additive mass shift, and the contact term fixed by the GW relation (your "contact = 1/2").
Summed over the whole $S^2\times[0,T]$ (the divergence of $A^a_\mu$ integrates to zero):

$$
2m\,\chi_P(m,R) = -\Sigma(m,R), \qquad
\chi_P = \frac{1}{V}\int d^3x\,d^3y\,\langle P^a(x)P^a(y)\rangle .
$$

Insert the $\ell=0$ spectral decomposition, $P^a_0(\tau) = R^2\!\int_{S^2}\! d\Omega\, P^a(\tau,\Omega)$,
$C_0(\tau) = \langle P^a_0(\tau)P^a_0(0)\rangle = \sum_n |\langle n|P^a_0|0\rangle|^2 e^{-E_n|\tau|}$,
$T\to\infty$:

$$
-\Sigma(m,R) = \frac{2m}{4\pi R^2}\sum_n \frac{2\,|\langle n|P^a_0|0\rangle|^2}{E_n(m,R)} .
$$

This is the exact-lattice version of "SSB $\Leftrightarrow$ Goldstone pole $\Leftrightarrow$ $\chi_P\sim\Sigma/m$":
the $1/m$ behavior can only come from terms whose $E_n\to0$. On the finite sphere:

- **SSB**: the sum is saturated by the rotor state, $E_{\rm adj}\propto 1/(F^2R^2)$ and
  $|\langle {\rm adj}|P^a_0|0\rangle|^2 \propto \Sigma^2 (4\pi R^2)^2$, giving at $m\to0$
  $\chi_P \propto \Sigma^2 F^2 R^4$ (growing) and $\Sigma(m,R)$ a scaling function of
  $m\Sigma\,\Theta\cdot 4\pi R^2$.
- **CFT**: $E_n \sim \Delta/R$, $|\langle n|P^a_0|0\rangle|^2\sim R^{2(2-\Delta_{\rm adj})}$, so
  $\chi_P\sim R^{\,3-2\Delta_{\rm adj}}$ and $\Sigma(m)\sim m^{\Delta_{\rm adj}/(3-\Delta_{\rm adj})}$ at
  $R\to\infty$.

Caveat on the susceptibility test alone: a CFT with $\Delta_{\rm adj} < 3/2$ (marginally the
$N=4$ large-$N$ value) also has a $\chi_P$ that grows with $R$, only more slowly than the SSB
$R^4$. The spectrum (Sec. 3) discriminates unambiguously: $E_{\rm adj}R\to0$ vs
$E_{\rm adj}R\to\Delta_{\rm adj}$.

## 5. Reading off from the existing channels

- Scalar-density connected $\ell=0$ ground state: recorded $E R\sim 1.6$, flat between $a_t=0.2$ and
  $0.1$, i.e. an $O(1)$ constant in the range of $\Delta_{\rm adj}$, not a rotor level.
- Axial tower ratios recorded as $\ell2/\ell1\to3/2$, $\ell3/\ell1\to2$: the conserved-current
  column ($\ell+1$), not the Goldstone column ($\sqrt3$, $\sqrt6$).
- $V/A\sim1$: all $SU(N)$ currents degenerate at $\Delta=2$, as for unbroken $G$.

All three sit in the "CFT" column at the $g^2R$ reached so far. Because $F^2 R = c_F\, g^2 R$ with
$c_F$ possibly small, the SSB column is only reachable once $c_F\,g^2R \gg 1$; a monotone decrease of
the $\ell=0$ adjoint-scalar $ER$ with $g^2R$ would be the first sign.

## Open questions

1. $N$ convention: does the runtime `Nf` count two-component flavors or $(\xi,\eta)$ pairs? (Sets
   which $\Delta_{\rm adj}(N)$ to compare the $\ell=0$ scalar against.)
2. Is the runtime `gsq` $g^2 a$ or $g^2 R$? Needed to place the L1-L4 x gsq grid on the $g^2R$ axis
   used above.
