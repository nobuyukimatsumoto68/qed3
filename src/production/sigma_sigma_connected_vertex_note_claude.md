# The two-vertex connected $\sigma\sigma$ operator (distillation vs the stochastic estimator)

Short theory note for the $\sigma^2$-$F^2$ (0++) mixing program (Chester-Pufu, arXiv:1603.05582). It records a
concrete operator-definition point found while validating chunk 3 of exact distillation
(`distill_contract_claude.py`, `distillation_impl_plan_claude.md`): the **connected** part of the
two-scalar operator needs **two** area-weighted vertices, whereas the stochastic loop estimator in
`jj_sigma_loops_stoch_claude.cu` uses only one. It is a constant rescaling at L1 but a genuinely different
operator at L2+.

## 1. The scalar operator

The zero-momentum ($\ell=0$) scalar density at timeslice $t$ is a site sum of the local bilinear with the
mesh area weight:

$$
\sigma(t) = \sum_x A_x\,\bar\psi(t,x)\,\psi(t,x),
\qquad
A_x = \text{dual\_area}_x \cdot Y_{00},
\qquad
Y_{00} = \frac{1}{\sqrt{4\pi}} .
$$

This matches `valence_claude.h::mult_Ylm_real(0,0)` (multiply each site by $\text{dual\_area}_x\,Y_{00}$, both
spins) followed by `accumulate_loop_raw` (sum $\bar\eta\,W\,\phi$, no $1/N_\text{sites}$). Define the
**vertex matrix** $W_0 = \mathrm{diag}_x(A_x)\otimes\mathbb{1}_\text{spin}$ on the $2N_s$ spinor space at a
timeslice.

## 2. Wick contraction of $\sigma\sigma$

The 0++ four-fermion operator is $\sigma\sigma(t) = \sigma(t)^2$:

$$
\sigma(t)^2 = \sum_{x,y} A_x A_y\,\bar\psi\psi(t,x)\,\bar\psi\psi(t,y).
$$

Contracting the four fermion fields gives two pieces:

- **Disconnected** (each density self-closes):
$$
\langle\bar\psi\psi(x)\rangle\langle\bar\psi\psi(y)\rangle
= \mathrm{tr}\,G(x,x)\,\mathrm{tr}\,G(y,y)
\;\longrightarrow\; D_S^2,
$$
two separate one-point loops, $D_S(t) = \sum_x A_x\,\mathrm{tr}\,G(t,x;t,x)$.

- **Connected** (cross-contraction of the two densities):
$$
\langle\bar\psi(x)\psi(y)\rangle\langle\bar\psi(y)\psi(x)\rangle
= -\,\mathrm{tr}\,[\,G(x,y)\,G(y,x)\,]
\;\longrightarrow\;
D'_S(t) = -\sum_{x,y} A_x A_y\,\mathrm{tr}\,[\,G(t,x;t,y)\,G(t,y;t,x)\,].
$$

The connected term carries **two** area weights, $A_x$ **and** $A_y$: one at each point where the single
fermion loop meets a $\sigma$. This is the whole point below.

## 3. Distillation (mode) language

With the perambulator $\tau_{kl}(t',t) = w_k^\dagger(t')\,D_{ov}^{-1}\,w_l(t)$ and
$\Phi(t) = V(t)^\dagger W_0 V(t)$ (an $N_v\times N_v$ vertex matrix), the propagator restricted to a
timeslice is $G|_{t,t} = V(t)\,\tau(t,t)\,V(t)^\dagger$, and:

$$
D_S(t) = \mathrm{Tr}\,[\,\Phi(t)\,\tau(t,t)\,],
\qquad
D'_{S,\text{phys}}(t) = \mathrm{Tr}\,[\,\Phi(t)\,\tau(t,t)\,\Phi(t)\,\tau(t,t)\,].
$$

The connected $D'_{S,\text{phys}}$ has **two $\Phi$'s** — one per physical vertex — exactly mirroring the two
area weights $A_x, A_y$ in Section 2.

## 4. What the stochastic estimator computes

The `jj_sigma_loops` extended-loop estimator builds

$$
D'_{S,\text{est}}(t) \;=\; \bar\eta\;W_0\;D_{ov}^{-1}\;\Pi_t\;D_{ov}^{-1}\;\eta ,
$$

where $\Pi_t$ is the class/timeslice projector (`project_to_class`, a bare zeroing of off-class timeslices,
**no** area weight). Two facts remove the second vertex:

1. the noise closure $E[\eta\,\eta^\dagger] = \mathbb{1}$ puts **no** area weight at the source end;
2. the middle $\Pi_t$ is a **bare** projector, not $W_0$.

So only one $W_0$ survives. In mode language $V^\dagger V = \mathbb{1}$ collapses the middle
$D_{ov}^{-1}\Pi_t D_{ov}^{-1}$ to $\tau\,\tau$, giving

$$
D'_{S,\text{est}}(t) = \mathrm{Tr}\,[\,\Phi(t)\,\tau(t,t)\,\tau(t,t)\,] \qquad\text{— one }\Phi.
$$

This is a legitimate object ($\mathrm{tr}[W_0 G G]$, a self-energy-like insertion) but it is **not** the
connected part of $(\sum_x A_x\bar\psi\psi)^2$.

## 5. Consequence — constant at L1, genuine at L2+

$$
\frac{D'_{S,\text{phys}}}{D'_{S,\text{est}}}
= \frac{\mathrm{Tr}[\Phi\,\tau\,\Phi\,\tau]}{\mathrm{Tr}[\Phi\,\tau\,\tau]} .
$$

- **L1**: the refinement-1 mesh is the regular icosahedron, which is vertex-transitive, so **all $A_x$ are
  equal** (verified: `dual_areas` uniform to $10^{-15}$, $= 4\pi/12$). Then $W_0 = c\,\mathbb{1}$ and
  $\Phi = c\,\mathbb{1}$, so the missing $A_y$ is a pure constant:
$$
D'_{S,\text{phys}} = w_{00}\;D'_{S,\text{est}},
\qquad
w_{00} = \frac{4\pi/12}{\sqrt{4\pi}} = 0.2954089752 .
$$
  `distill_contract_claude.py` confirms the ratio $= 0.295409$ to all printed digits. An overall constant
  does not change "is the mixing nonzero", so the **L1 mixing check stayed valid**.

- **L2+**: the mesh is **not** vertex-transitive, $A_x$ varies site to site, so $\Phi \neq c\,\mathbb{1}$ and
  $\tau\,\Phi\,\tau \neq c\,\tau\,\tau$ — the middle $\Phi$ reweights the distillation modes unequally. The
  estimator then computes a genuinely different operator, not a rescaled one. **The exact coupled GEVP must
  use $\mathrm{Tr}[\Phi\,\tau\,\Phi\,\tau]$**, which distillation provides directly (no extra solves).

## 6. It is a vertex COUNT, not a weighting error (`mult_Ylm_real` already carries $A$)

`mult_Ylm_real` folds the area $A_x=\text{dual\_area}_x$ **together** with $Y_{\ell m}$ (Section, the
primitive). So wherever the $Y_{\ell m}$ vertex is applied, the area comes with it automatically -- the
weighting practice is complete **at each vertex where it is applied**. The $D'_S$ issue is therefore not a
wrong weight but a **miscount of vertices**: the connected $\sigma\sigma$ has TWO vertices, and the estimator
applies `mult_Ylm_real` only ONCE.

In `jj_sigma_loops_stoch_claude.cu` the extended-loop step is (schematically)
```
project_to_class(phiP, phi, ...);        // bare timeslice projector -- NO mult_Ylm_real here
op_Dsq.solve(chi1, D^dag phiP);          // second solve
Gamma = chi1;  Gamma.mult_Ylm_real(...); // SINK vertex only  (one A Y_lm)
```
The middle point `phiP` gets a bare `project_to_class` with **no** `mult_Ylm_real`, so that second vertex has
neither its $Y_{00}$ nor its area $A$. The stochastic fix is exactly the missing insertion:
```
project_to_class(phiP, phi, ...);
phiP.mult_Ylm_real(0, 0, base);          // <-- MIDDLE vertex: the l=0 sigma vertex, carries A Y_00
op_Dsq.solve(chi1, D^dag phiP);
```
which turns $\mathrm{Tr}[\Phi\tau^2]$ into $\mathrm{Tr}[\Phi\tau\Phi\tau]$. (Both $\sigma$'s are $\ell=0$, so
the middle vertex is $Y_{00}$, independent of the sink $\ell m$ used for the $F^2_{\ell m}$ cross.)

This is why the disc and conn drivers are fine -- the number of `mult_Ylm_real` calls matches the number of
vertices:
- **disc one-point** $D_S=\mathrm{tr}[W_0 G]$: ONE vertex, ONE `mult_Ylm_real`. Correct.
- **single-meson conn** $C_S=\mathrm{tr}[W_0 G W_0 G]$: TWO vertices, `mult_Ylm_real` at source
  (`jj_local_ylm_scalar_conn_stoch_claude.cu:471`) AND sink (`:481/:483`). Correct (two $A$ factors).
- **$D'_S$ extended loop**: TWO vertices but only ONE `mult_Ylm_real` (sink) -> missing the middle $A Y_{00}$.

Distillation needs no such fix: $\Phi=V^\dagger W_0 V$ carries the area at EVERY vertex, so
$\mathrm{Tr}[\Phi\tau\Phi\tau]$ has $A$ at both automatically.

**Governing rule.** The number of area factors $A$ equals the number of independent site sums $\sum_x$, which
equals the number of vertices: exactly **one** $A_x$ (one `mult_Ylm_real`, or one $\Phi$) per $\sum_x$. A
connected loop through $n$ vertices carries $n$ sums and $n$ area factors; a disconnected product carries one
$A$ per factor on its OWN independent index. So:
- $D_S=\sum_x A_x\mathrm{tr}\,G(x,x)$: one $\sum_x$, one $A$.
- $D_S^2=(\sum_x A_x\ldots)(\sum_y A_y\ldots)$: two INDEPENDENT sums, one $A$ each ($A_xA_y$, never $A_y^2$).
- $D'_S=\sum_{x,y}A_xA_y\mathrm{tr}[G(x,y)G(y,x)]$: two vertices, two $A$'s.
- $C_S=\sum_{x,y}A_xA_y\ldots$ (source+sink), $E=C_S^2$: four independent sums, four $A$'s.

**Corollary for the factorizable four-point diagrams (chunk 4).** Build $E=C_S^2$, $F=D'_S D'_S$, $J=D_S^4$,
etc. as products of EXACT, independently per-vertex-$A$-weighted distilled factors ($\Phi$ at each vertex on
its own site) -- no shared-index $A^2$, and no square bias (distillation is exact, so no distinct-hit trick is
needed). The existing stochastic $C_S$ is already area-summed (its per-site structure is spent) and noisy, so
$C_S$ and $E$ are REBUILT from the perambulators, where each vertex is correctly weighted $A_x$ at its site.

## 7. FS channel and the remaining pin

The parity-indefinite $\sigma_\text{FS}$ density has a furnished vertex $(1-D_{ov}^\dagger)$. Its connected
$\sigma_\text{FS}\sigma_\text{FS}$ has the **same** one-vs-two-vertex structure (two $\Phi$'s in the physical
operator), **plus** a placement question: where the furnishing $(1-D_{ov}^\dagger)$ sits across the two
vertices, i.e. which legs are forward $\tau$, which carry $\tau_\text{gw}=V^\dagger(1-D_{ov}^\dagger)D_{ov}^{-1}V$,
and which collapse via the GW identity (IV.17) of `qed3_v2-6.pdf`,

$$
D_{ov}^{-\dagger} = \mathbb{1} - D_{ov}^{-1},
\qquad
(1-D_{ov}^\dagger)\,D_{ov}^{-\dagger} = -\,D_{ov}^{-1}.
$$

That FS four-point bookkeeping is pinned against v3-4 Ch.5 in **chunk 4** (four-point + coupled GEVP).

## 8. Bottom line

Chunk 3 both **validated** the perambulator machinery (T3a: distilled $D_S,D_S^{1mD}$ equal the stochastic
loops to 5 digits) **and surfaced** that the correct connected two-scalar operator is the two-vertex trace
$\mathrm{Tr}[\Phi\,\tau\,\Phi\,\tau]$. Distillation gives it for free; the stochastic driver's single-vertex
form is a constant off at L1 and wrong at L2+.
