# Why the lattice propagator $P(n{@}0\to n{@}t)$ "varies with $n$" at L=1 — it's the spinor frame, not a symmetry violation

## TL;DR

The $L=1$ lattice is the icosahedron, which is **vertex-transitive** — all 12 sites are equivalent, so
icosahedral symmetry *does* hold. But the propagator $P$ is a **2-component spinor matrix**, and each site
carries its **own local spinor frame** (fixed by the link angles `alpha` and the spin connection `Omega`).
The symmetry maps site $n\to n'$ **together with an $SU(2)$ rotation of that local frame**, so $P$ transforms
**covariantly, not invariantly**:

$$ P(n'{@}0\to n'{@}t) = U_n\, P(n{@}0\to n{@}t)\, U_n^\dagger ,\qquad U_n\in SU(2). $$

Therefore the **individual spinor components** (and their absolute values $|P_{ab}|$) vary with $n$, while the
**frame-invariant** quantities — block Frobenius norm $\|P\|_F$, $\det P$, $\mathrm{tr}(P^\dagger P)$, the
eigenvalues — are site-independent. Any per-$n$ wiggle in the *components* is expected; it is not a broken
symmetry and (for the same-site block) it is roundoff. Compare invariants, not raw entries.

## What was observed and the resolution

The check was on the **same-site, time-separated** block $P(n{@}0\to n{@}t)$ (source at site $n$, time 0;
sink at the *same* site $n$, time $t$), and "it varies with $n$" was seen by taking absolute values of the
entries. Measured on `data_free_vmRe0.000000vmIm0.000000/prop_deter_L1/Dinv.0.h5` (`Dm_inv`, $N_t=128$,
12 sites), $t=1$, across all 12 sites:

| quantity | mean over $n$ | spread $(\max-\min)/\text{mean}$ |
|---|---|---|
| $\|P[0,0]\|=\|P[1,1]\|$  (the $\sigma_3$ coefficient) | $1.167\times10^{-1}$ | $3.7\times10^{-15}$ |
| $\|P[0,1]\|,\ \|P[1,0]\|$  ($\sigma_1,\sigma_2$ parts) | $2.3\times10^{-16}$ | $\sim 2.6$ |
| $\|c_0\|$ (identity part) | $3.4\times10^{-10}$ | — |
| $\|\text{block}\|_F$ | (const) | $2.2\times10^{-15}$ |
| $\|\det\|$ | (const) | $3.6\times10^{-15}$ |

Reading this:

- The **same-site** block is **purely $c_3\,\sigma_3$** — the *temporal* spinor structure — and its coefficient
  $c_3=|P[0,0]|=|P[1,1]|$ is **constant across all 12 sites to $2\times10^{-15}$** (roundoff). The $\sigma_3$
  direction is the global time axis, common to every site, so this block has **no frame freedom at all**.
- The $\sigma_1,\sigma_2$ and identity parts are **$\sim10^{-16}$ — machine zero**. Their relative spread
  "$\approx 2.6$" looks large only because it is $(\max-\min)/\text{mean}$ on numbers that are all $\sim10^{-16}$
  (noise / noise). **That is the "variation": roundoff on null components.** Taking $|\cdot|$ does not help —
  the entries are numerically zero.
- $\|\text{block}\|_F$ and $\det$ are constant to roundoff $\Rightarrow$ **icosahedral symmetry intact.**

So at $t\ge1$, $P(n{@}0\to n{@}t)\propto\sigma_3$ for every $n$ with the *same* coefficient; nothing physical
varies with $n$.

## Same-site vs off-site (the genuine frame case)

The $SU(2)$-frame caveat is real, but it shows up in the **off-site** blocks $P(n{@}0\to m{@}t)$ with $m\neq n$
(and in the spatial Dirac hops generally). There the point-split direction $\hat e^a_{nm}\sigma^a$ between two
sites genuinely rotates as you move around the icosahedron, so:

- the individual components $P_{ab}$, and even $|P_{ab}|$, **do** vary from one (oriented) link to another;
- but $\|P(n{@}0\to m{@}t)\|_F$ over the set of links related by the symmetry **is** constant.

This is exactly the spin-1/2 covariance $P\to U P U^\dagger$: norms/traces/determinants invariant, components not.

## The right diagnostic

To test the lattice symmetry, compare **frame-invariant scalars** built from the $2\times2$ spinor block:

- $\|P\|_F=\big(\sum_{ab}|P_{ab}|^2\big)^{1/2}$, or $\mathrm{tr}(P^\dagger P)$,
- $\det P$, or the eigenvalues of $P$ (or of $P^\dagger P$).

Do **not** compare a single entry $P_{ab}$ or its absolute value $|P_{ab}|$ across sites/links — those carry
the frame and are not invariant. If an invariant varies by more than $\sim10^{-13}$ (well above the LU
roundoff; the propagator's own self-check reports $\|D_mD_m^{-1}-I\|_F\sim10^{-13}$), *then* it would signal a
real effect; here they are constant to $\sim10^{-15}$.

## Code references (where the local frame comes from)

- `includes/dirac_simp.h:227` — `gamma(ix,iy) = cos(alpha)*sigma[1] + sin(alpha)*sigma[2]`: the per-link
  $\hat e\cdot\sigma$ direction set by the link angle `alpha` (read from `geometry/data/alpha_n*.dat`).
- `includes/dirac_simp.h` — `Omega(ix,iy) = cos(om/2)*sigma[0] - i*sin(om/2)*sigma[3]`: the spin connection
  (parallel transport between the two sites' frames), `om` from `omega_n*.dat`.
- The dense `Dov`/`Dm_inv` are written by `jj_propagator_deter_claude.cu`; layout = flat index
  $= N_x\,t + N_S\,\text{site} + \text{spin}$, $N_x=N_S\,N_\text{sites}$ (time is the OUTER block).
  Same-site block at $(n,t)$: rows $N_x t + N_S n + \{0,1\}$, cols $N_S n + \{0,1\}$.

(Reproduce the table: load `Dm_inv`, slice the $2\times2$ same-site block for each $n$, print $|P[0,0]|$ and
$\|\text{block}\|_F$ across $n$.)
