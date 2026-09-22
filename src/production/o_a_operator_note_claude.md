# What $O_A$ is, and why it resolves the $(2,2)$ excited meson

## Definition

$$
O_A(t) \;=\; \sum_{x,y} w(x)\;\bar\psi(x,t)\,\tilde\tau(x,y;t)\,\psi(y,t),
\qquad w(x)=A(x)\,Y_{00},\quad \tilde\tau = D_\text{ov}^{-1}-\tfrac12 .
$$

It is a **bilocal scalar density**: the antifermion sits at $x$, the fermion at $y$, and the two are
connected by the improved (contact-subtracted) equal-time propagator $\tilde\tau(x,y)$. Contrast with
the ordinary local scalar
$$
\sigma_{00}(t)=\sum_x w(x)\,\bar\psi(x,t)\,\psi(x,t)\quad(\text{fermion and antifermion at the same point}).
$$

In the complete distillation basis the vertex is remarkably simple (since $V V^\dagger=1$):
$$
\Phi_A = V^\dagger\,\big[\operatorname{diag}(w)\,\tilde\tau\big]\,V \;=\; \Phi_{00}\,\tilde\tau ,
$$
i.e. the ordinary $\ell=0$ wall vertex $\Phi_{00}$ times the equal-time perambulator $\tilde\tau$.

## What it creates: orthogonal to $(1,1)$, overlaps $(2,2)$

The local $\sigma_{00}=\bar\psi\psi$ makes the fermion pair at a single point -> its lowest overlap is
the **ground** meson $(1,1)$ = both fermions in the $\lambda=1$ Dirac shell (Eq C.17,
$\lambda=n+|m|+\tfrac12$), $a_t m = 2E_0$.

$O_A$ spreads the pair over the sphere through $\tilde\tau$, giving it a spatial wavefunction with a
node. Measured:
- $\langle O_A(t)O_A(0)\rangle$ effmass $\to 2E_1$ = the **$(2,2)$ excited meson** (L1: 0.56 clean).
- $\langle O_A\,\sigma_{00}\rangle \approx 0$ -> **$O_A$ is orthogonal to the $(1,1)$ ground.**

So $O_A$ is effectively "the $\ell=0$ scalar with the ground meson projected out": its lowest state is
the first radial excitation $(2,2)$.

Why the orthogonality is exact (not just small): the free triangle/loop obeys $\sigma_3$-hermiticity
($D\sigma_3=-\sigma_3 D$, Eq C.16 $\Rightarrow \sigma_3 G\sigma_3=-G$). Any $\sigma_3$-even vertex
$\Gamma$ (identity, $\Omega$-transport, currents) gives $\mathrm{Tr}[\Gamma GGG]=0$, so $\sigma_{00}$
cannot mix with $\sigma_{00}^2$. The propagator kernel is $\sigma_3$-**mixed**,
$\sigma_3\tilde\tau\sigma_3=-\tilde\tau-1$, which is exactly what lets $O_A$ (i) overlap $(2,2)$ and
(ii) mix with $\sigma_{00}^2$ ($C_{2A}\neq0$) -- the two things every local/geometric operator lacked.

## Role in the two-meson GEVP

$\sigma_{00}^2$ (the two-meson operator) is contaminated by the $(2,2)$ single meson through diagram A
(its non-local equal-time loop is literally $\tilde\tau$). Because $\sigma_{00}$ has $C_{12}=0$, it can
never subtract $(2,2)$ out. $O_A$ does: it overlaps $(2,2)$ and mixes with $\sigma_{00}^2$, so the
connected $2\times2$ GEVP $\{\sigma_{00}^2,O_A\}$ splits them -- level 0 $\to (2,2)$, level 1 $\to$
two-meson $2m_\sigma$. Confirmed L1 ($0.756$) and L2 ($0.786$).

## Seeing $(2,2)$ cleanly

The cleanest $(2,2)$ is just the $O_A$ two-point $\langle O_A O_A\rangle$ (a single-meson correlator):
L1 plateaus at $0.56$. At L2 it decays faster and the single free configuration runs into
solve-tolerance noise past $t\sim15$ (no plateau) -- proliferation does not help (shown $K=1,2,3$
identical). To sharpen $(2,2)$: (a) statistics -- the interacting ensembles average over gauge configs;
(b) a dedicated $\lambda=2$-shell projector operator as a second $(2,2)$ interpolator.

## Refs
`qed3_v2-6.pdf` App C (Eq C.16-C.17, free $S^2$ Dirac spectrum); `qed3int_v3-4.pdf` Eq (5.5) diagram A.
