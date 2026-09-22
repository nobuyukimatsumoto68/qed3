# Collision contact and the single-$\sigma$ admixture in $O=\sigma^2$ (derivation)

Working note for the two-scalar $\sigma\sigma$ vs 0++ study. Goal: understand why the summed two-meson correlator
plateaus at the single-$\sigma$ mass ($a_t m\approx0.35$), and what operator gives the clean two-$\sigma$ / 0++.
Companion: `sigma_ps_minus_half_coeffs_claude.md` (diagram coefficients), `sigma_connected_gevp_benchmark_claude.md`.

## 1. Operators

Local scalar density $\rho(x)=\bar\psi S\psi(x)$ ($S=1$ for PS). Its contact VEV is a c-number,
$\langle\rho\rangle=\tfrac12$ (from $D_\text{ov}^{-1}=\tfrac12\delta + \text{prop}$; sign details in the coeffs note).
Single-density normal-ordering (level 1):
$$
\sigma(x) = \rho(x) - \langle\rho\rangle = \rho(x) - \tfrac12 .
$$
$\ell{=}0$-projected density and the two-meson (composite) operator:
$$
\Sigma(s) = \sum_x w_x\,\sigma(x,s),\qquad w_x = A_x Y_{00};\qquad O(s) = \Sigma(s)^2 .
$$

## 2. Wick building blocks (distillation)

$$
D_S(s)=\mathrm{Tr}[\Phi\,\tau(s,s)]\ \ (\text{self-contraction / tadpole}),\qquad
D'_S(s)=\mathrm{Tr}[\Phi\,\tau(s,s)\,\Phi\,\tau(s,s)]\ \ (\text{equal-time } x\!-\!y \text{ loop}),
$$
$$
C_S(s,t)=\mathrm{Tr}[M],\ \ M=\Phi\tau(s,t)\Phi\tau(t,s)\ \ (\text{single-}\sigma\text{ propagator}),\quad
T_S=\mathrm{Tr}[M^2],\ V_S,\ S_S\ (\text{connected}).
$$
$\Phi = V^\dagger W_0 V$. The 10 diagrams A--J are traces of these; weights $\{4,2,4,4,2,1,4,1,1,1\}$.

## 3. Two distinct subtractions

There are **two independent** contact removals, not one:

**(L1) Tadpole / self-contraction** -- $\sigma = \rho - \tfrac12$. Removes the VEV of a *single* bilinear:
$D_S \to D_S - c_0$, $c_0=\tfrac12\mathrm{Tr}\,\Phi$. Diagrammatically kills C,D,G,H,I,J (they carry $\delta D_S\approx0$).
Connected legs ($D'_S,V_S,S_S,C_S,T_S$) stay **raw**. This is a genuine c-number, verified $\langle D_S-c_0\rangle=0$.

**(L2) Collision contact** -- the equal-time $\tau(s,s)$ that links the **two different** $\sigma$'s of $O(s)=\Sigma^2$
at the same timeslice. This is NOT a self-contraction; L1 leaves it in place. It is the source of the single-$\sigma$
admixture (Sec 4).

## 4. Why $\sigma^2$ carries single-$\sigma$: the collision collapse

Diagram A $= -S_S = -\mathrm{Tr}[\Phi(s)\,\tau(s,s)\,\Phi(s)\,\tau(s,t)\,\Phi(t)\,\tau(t,t)\,\Phi(t)\,\tau(t,s)]$.
At the contact value $\tau(s,s)\to\tfrac12 I$, the colliding leg collapses:
$$
\Phi(s)\,\tfrac12 I\,\Phi(s) = \tfrac12\,\Phi(s)^2\ \Longrightarrow\
-\tfrac12\,\mathrm{Tr}[\Phi(s)^2\,\tau(s,t)\,\Phi(t)\tau(t,t)\Phi(t)\,\tau(t,s)] ,
$$
a **single $s\!\to\!t$ propagation** dressed by a local vertex -- i.e. one $\sigma$ is eaten into a constant, the
other propagates as a single meson. Hence A $\to$ single-$\sigma$ ($a_t m\approx0.35$) at large $dt$. Formally:
$\langle 0|O|\sigma\rangle \ne 0$ through this equal-time collision. $O=\Sigma^2$ is **not** a pure two-particle
operator; its lightest overlap is the single $\sigma$, which dominates the tail.

## 5. Removing the single-$\sigma$: exact vs ultralocal

The clean two-particle interpolator is the point-split composite
$$
:\!\Sigma^2\!:(s) \;=\; \lim_{x\to y}\Big[\Sigma_x(s)\,\Sigma_y(s) - \langle \Sigma_x(s)\,\Sigma_y(s)\rangle\Big],
$$
i.e. subtract the **full** equal-time contraction $\langle\Sigma(x)\Sigma(y)\rangle$ (the collision) between the two
$\sigma$'s. That contraction has an ultralocal part ($\tfrac12$) and a finite coincidence part:
$$
\langle\Sigma(x)\Sigma(y)\rangle = \underbrace{\tfrac12\,\delta_{xy}}_{\text{ultralocal}} + \underbrace{(\text{coincident propagation})}_{\text{finite}} .
$$
- **Ultralocal subtraction** $\tau(s,s)\to\tau(s,s)-\tfrac12 I$ in the colliding legs (A, $V_S$, $S_S$): removes only
  the $\tfrac12\delta$. This is the "$\tau-\tfrac12 I$ everywhere" prescription.
- **Full subtraction**: also removes the finite coincident-propagation overlap; equivalent to projecting out the
  single-$\sigma$ state $\Sigma$, i.e. $O \to O - c\,\Sigma$ with $c=\langle 0|O|\sigma\rangle/\langle 0|\Sigma|\sigma\rangle$.

## 6. Numerical status (Nf2 $g^2{=}0.5$ L1, 400 cfg)

A+B+E summed effmass (weights $2\{4,2,2\}$):

| operator | large-$dt$ plateau | reading |
|---|---|---|
| raw $\sigma^2$ (L1 only) | $\mathbf{0.35}$ | single-$\sigma$ (collision kept) |
| ultralocal $:\!\sigma^2\!:$ (L1+L2 $\tfrac12 I$) | $\mathbf{0.50}$ | single-$\sigma$ mostly gone, not complete |

Per diagram after the ultralocal L2 (raw effmass): **B $=0.69$, E $=0.72$ -- clean two-$\sigma$, above 0.6**;
**A $=0.53$ -- still below 0.6** (residual single-$\sigma$ from the finite coincident part). So the summed 0.50 is
dragged down by A alone; the $\tfrac12 I$ removed the ultralocal collision but not the finite coincidence overlap.

## 7. Consequences / open points

- B ($T_S$) and E ($C_S^2$) are already clean two-$\sigma$ (no equal-time legs to collide) -- the plateau *does*
  come above 0.6 for them.
- A ($S_S$) has the residual single-$\sigma$. The **full** collision subtraction (Sec 5) or, equivalently, the
  **GEVP** with the single-$\sigma$ operator $\Sigma$ in the basis, is what removes it exactly.
- Reflection on the "awkwardness": the ultralocal $\tfrac12 I$ is a *partial* single-$\sigma$ projection; there is no
  reason it should fully clean A. The principled objects are (i) $\Sigma$ (single meson), (ii) $:\!\Sigma^2\!:$
  (two meson), (iii) $F^2$ (0++), fed to a GEVP that orthogonalizes them.

**Open question for NM:** is there a closed-form "full collision" subtraction ($O - c\,\Sigma$ with a computable
$c$ from the contact structure), or do we let the GEVP do the single-$\sigma$ projection? The former would give a
clean single-operator two-$\sigma$ effmass; the latter is the robust route.

## 8. Files

`two_meson_cnumber_claude.py` (L1 raw), `diag_effmass_claude.py` (L2 ultralocal $\tau-\tfrac12 I$),
`two_meson_normord_ABE_sum_claude.py` (raw vs ultralocal A+B+E), `two_meson_normord_diageff_claude.py`
(per-diagram, collision-removed).
