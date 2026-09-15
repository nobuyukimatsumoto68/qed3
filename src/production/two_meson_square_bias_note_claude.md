# Squaring loop factors: is $\langle\langle\cdot\rangle\rangle^2$ a valid finite-hit estimator?

Focused note on the disagreement about diagram $\mathfrak{A}_J$, Eq. (5.14) of `qed3int_v3-2.pdf`.

**The claim under dispute (NM):** given the two-hit disconnected loop data for $D_S(t)$, Eq. (5.15),
we can estimate $\mathfrak{A}_J$ by
$$
\mathfrak{A}_J \;\overset{?}{\approx}\; \bar D_S^{(2)}(t)^2\,\bar D_S^{(2)}(0)^2,
$$
where $\bar D_S^{(2)}$ is the two-hit average. This note argues that expression is **biased** at
finite hits, quantifies the bias, and gives the zero-extra-cost fix that reuses the very same two hits.

---

## 1. Setup and notation

Fix a gauge configuration $U$. A single stochastic "hit" $h$ draws a $Z_2\times Z_2$ noise source
$\phi^{(h)}$ (time- and spin-diluted, interval $N_t/4$, per Sec. 5) and produces a single-hit loop
estimate. Write the elementary loop factor at timeslice $\tau$ as
$$
X^{(h)}(\tau)\;\equiv\;\phi_\tau^{(h)\dagger}\,S\,D_{ov}^{-1}\,\phi_\tau^{(h)} .
$$
The stochastic (noise) expectation of one hit reproduces the exact loop **on this configuration**:
$$
E_\eta\!\big[X^{(h)}(\tau)\big]\;=\;\mathcal{D}(\tau)\;\equiv\;\mathrm{tr}\big[\Pi_\tau\,S\,D_{ov}^{-1}\big],
\qquad
\mathrm{Var}_\eta\!\big[X^{(h)}(\tau)\big]\;\equiv\;v(\tau).
$$
Here $\mathcal{D}(\tau)$ is a deterministic number (no noise left); $v(\tau)$ is the per-hit stochastic
variance. Different hits are i.i.d.: $E_\eta[X^{(h)}X^{(h')}] = \mathcal{D}^2$ for $h\ne h'$.

Two-hit average and its square:
$$
\bar X(\tau)\;=\;\tfrac12\big(X^{(1)}(\tau)+X^{(2)}(\tau)\big),
\qquad
\bar D_S^{(2)}(\tau)\;\equiv\;\bar X(\tau).
$$

**Two distinct expectations must not be conflated.**
- $E_\eta[\,\cdot\,]$ = average over stochastic noise at fixed $U$ (removes noise, leaves the exact
  loop $\mathcal{D}$).
- $\langle\,\cdot\,\rangle_g$ = average over gauge configurations $U$ (the physics).

The double bracket $\langle\langle\cdot\rangle\rangle$ in the writeup is the **noise** average $E_\eta$
(finite number of hits in practice). The physical correlator is $\langle \mathfrak{A}\rangle_g$. The bias
below lives in $E_\eta$ and therefore is **not** removed by taking more configurations.

---

## 2. The elementary identity: square of a mean is biased

For any single stochastic quantity $X$ with $E_\eta[\bar X]=\mathcal{D}$,
$$
E_\eta\!\big[\bar X^2\big]
\;=\;\big(E_\eta[\bar X]\big)^2+\mathrm{Var}_\eta[\bar X]
\;=\;\mathcal{D}^2+\frac{v}{H},
$$
using $\mathrm{Var}_\eta[\bar X]=v/H$ for $H$ i.i.d. hits. So the **square of the hit-average
overestimates $\mathcal{D}^2$ by $v/H$** — a strictly positive bias.

Because this term is an $E_\eta$-level object, averaging over gauge configs gives
$$
\big\langle E_\eta[\bar X^2]\big\rangle_g
=\big\langle \mathcal{D}^2\big\rangle_g+\frac1H\big\langle v\big\rangle_g,
$$
i.e. the bias $\tfrac1H\langle v\rangle_g$ **persists at any number of configurations**. Only
increasing the number of hits $H$ shrinks it (as $1/H$); it never vanishes at finite $H$.

### The unbiased replacement, using the SAME hits

The distinct-hit (off-diagonal) product is unbiased for $\mathcal{D}^2$:
$$
\widehat{\mathcal{D}^2}
\;\equiv\;\frac{1}{H(H-1)}\sum_{h\ne h'}X^{(h)}X^{(h')}
\;=\;\frac{H\,\bar X^2-\overline{X^2}}{H-1},
\qquad
\overline{X^2}\equiv\frac1H\sum_h \big(X^{(h)}\big)^2 ,
$$
with $E_\eta[\widehat{\mathcal{D}^2}]=\mathcal{D}^2$ exactly. The two written forms are algebraically
identical (expand $\sum_{h\ne h'}=(\sum_h X^{(h)})^2-\sum_h (X^{(h)})^2$). For $H=2$ it collapses to
$$
\widehat{\mathcal{D}^2}\Big|_{H=2}\;=\;X^{(1)}X^{(2)} .
$$
No new solves: the two hits are the ones already computed for $D_S$. The **only** change is
multiply-hit-1-by-hit-2 instead of average-then-square. The difference between the two prescriptions is
exactly the removable self-variance:
$$
\bar X^2-X^{(1)}X^{(2)}=\tfrac14\big(X^{(1)}-X^{(2)}\big)^2\;\ge 0 ,
$$
whose noise expectation is $v/2$ — precisely the bias of $\bar X^2$ at $H=2$.

---

## 3. Diagram $\mathfrak{A}_J$, Eq. (5.14), explicitly

Eq. (5.14) reads
$$
\mathfrak{A}_J=\big\langle\!\big\langle \phi_0^\dagger S D_{ov}^{-1}\phi_0\big\rangle\!\big\rangle^2
\cdot
\big\langle\!\big\langle \tilde\phi_t^\dagger S D_{ov}^{-1}\tilde\phi_t\big\rangle\!\big\rangle^2 ,
$$
i.e. $\bar X_0^2\,\bar Y_t^2$ with $X_0^{(h)}=\phi_0^{(h)\dagger}SD_{ov}^{-1}\phi_0^{(h)}$ (mean
$\mathcal{D}_0$, variance $v_0$) and $Y_t^{(k)}=\tilde\phi_t^{(k)\dagger}SD_{ov}^{-1}\tilde\phi_t^{(k)}$
(mean $\mathcal{D}_t$, variance $v_t$). The two factors use **independent** sources
$\phi_0\perp\tilde\phi_t$ (that is what the tilde denotes), so the $t{=}0$ block and the $t$ block do
**not** cross-contaminate. Hence
$$
E_\eta\!\big[\mathfrak{A}_J^{\text{literal}}\big]
=E_\eta[\bar X_0^2]\;E_\eta[\bar Y_t^2]
=\Big(\mathcal{D}_0^2+\tfrac{v_0}{H}\Big)\Big(\mathcal{D}_t^2+\tfrac{v_t}{H}\Big).
$$
Expanding, the physical target is $\mathcal{D}_0^2\mathcal{D}_t^2$, and the excess is
$$
\text{bias}\big(\mathfrak{A}_J^{\text{literal}}\big)
=\frac{\mathcal{D}_0^2\,v_t+\mathcal{D}_t^2\,v_0}{H}
+\frac{v_0 v_t}{H^2}\;>\;0 .
$$

**Unbiased version (same two hits, no extra solves):** debias each square separately,
$$
\mathfrak{A}_J^{\text{unb}}
=\Big[\tfrac{1}{H(H-1)}\!\sum_{h\ne h'}X_0^{(h)}X_0^{(h')}\Big]
\Big[\tfrac{1}{H(H-1)}\!\sum_{k\ne k'}Y_t^{(k)}Y_t^{(k')}\Big]
\;\xrightarrow{H=2}\;
X_0^{(1)}X_0^{(2)}\;Y_t^{(1)}Y_t^{(2)} ,
$$
with $E_\eta[\mathfrak{A}_J^{\text{unb}}]=\mathcal{D}_0^2\mathcal{D}_t^2$ exactly, for any $H\ge2$. The
cross-block factorization is clean because $\phi_0\perp\tilde\phi_t$; the only thing that had to be fixed
is each **within-block square**.

---

## 4. Does the bias matter numerically?

The bias is $\sim v/H$ per squared factor, relative to a signal $\sim\mathcal{D}^2$. So the relative
bias per squared loop is
$$
\frac{\text{bias}}{\text{signal}}\;\sim\;\frac{1}{H}\,\frac{v}{\mathcal{D}^2}
=\frac{1}{H}\,\frac{\mathrm{Var}_\eta[X]}{\big(E_\eta[X]\big)^2}
=\frac{1}{H\,\mathrm{SNR}^2},
$$
where $\mathrm{SNR}=\mathcal{D}/\sqrt v$ is the per-hit stochastic signal-to-noise of the loop.

For a **disconnected** scalar loop this is exactly the regime where the naive square fails: single-hit
disc loops are notoriously noisy, $v\gtrsim\mathcal{D}^2$ (SNR $\lesssim 1$ per hit), so at $H=2$ the
bias on each squared factor is $\gtrsim \tfrac12\mathcal{D}^2$ — **order the signal itself**, and
$\mathfrak{A}_J$ has two such squared factors. This is not a small $O(a^2)$-style correction we can drop;
it is a leading contamination. (For a factor with SNR $\gg1$ per hit the bias would be negligible and the
distinction academic — but disc loops are the opposite limit.)

---

## 5. Where this applies among the ten diagrams

The bias is present **iff a factor is squared**, i.e. wherever a $\langle\langle\cdot\rangle\rangle^2$
appears:
- $\mathfrak{A}_E$ (5.10): $\langle\langle\phi_0^\dagger SD_{ov}^{-1}\Pi_t SD_{ov}^{-1}\phi_0\rangle\rangle^2$
  — square of the meson correlator $C_S$. Debias: $C_S^{(1)}C_S^{(2)}$.
- $\mathfrak{A}_H$ (5.13): $\langle\langle\phi_0^\dagger SD_{ov}^{-1}\phi_0\rangle\rangle^2\cdot(\dots)$
  — the $D_S(0)$ factor is squared; debias it. The lone $D'_S(t)$ factor is linear (unbiased as is).
- $\mathfrak{A}_I$: time-reverse of $H$ — same treatment.
- $\mathfrak{A}_J$ (5.14): both $D_S(0)$ and $D_S(t)$ squared — this note.

The remaining diagrams have **no** squared factor and are unbiased exactly as written:
- $\mathfrak{A}_F$ (5.11) $=\langle\langle\cdots\phi_0\rangle\rangle\cdot\langle\langle\cdots\tilde\phi_t\rangle\rangle$
  — product of two **distinct linear** factors with independent sources; $E_\eta=\mathcal{D}'_0\mathcal{D}'_t$.
- $\mathfrak{A}_G$ (5.12) — three distinct linear factors with independent sources.
- $\mathfrak{A}_C,\mathfrak{A}_D$ (5.9) — a linear loop times an independent nested chain.
- $\mathfrak{A}_A,\mathfrak{A}_B$ (5.7-5.8) — a single connected chain, one source; no product of stochastic
  estimates at all.

So exactly the four squared diagrams $\{E,H,I,J\}$ need the distinct-hit product; this is what the memory
line "distinct-hit for squares E/H/I/J" means. It does **not** mean "acquire a second independent source";
it means "pair the two hits you already have."

---

## 6. Summary of the disagreement

| quantity | estimator | $E_\eta$ | verdict |
|---|---|---|---|
| $\mathcal{D}^2$ | $\bar X^2$ (square of average) | $\mathcal{D}^2+v/H$ | **biased** ($+v/H$) |
| $\mathcal{D}^2$ | $X^{(1)}X^{(2)}$ (distinct-hit, $H=2$) | $\mathcal{D}^2$ | unbiased |
| $\mathfrak{A}_J$ | $\bar D_S^{(2)}(0)^2\,\bar D_S^{(2)}(t)^2$ (NM's proposal) | $(\mathcal{D}_0^2+\tfrac{v_0}{H})(\mathcal{D}_t^2+\tfrac{v_t}{H})$ | **biased** |
| $\mathfrak{A}_J$ | $X_0^{(1)}X_0^{(2)}\,Y_t^{(1)}Y_t^{(2)}$ | $\mathcal{D}_0^2\mathcal{D}_t^2$ | unbiased |

- **Agreement:** no extra sources or solves are needed; the two hits already computed for $D_S$ suffice.
- **Disagreement:** the plain $\bar D_S^{(2)}(0)^2\bar D_S^{(2)}(t)^2$ (Eq. 5.14 read at finite hits) is
  biased by the per-loop self-variance $v/H$; for noisy disc loops this is order the signal. The
  distinct-hit product $X^{(1)}X^{(2)}$ removes it exactly, at zero extra cost, from the same data.

**If you still disagree**, the crux is Sec. 2: whether $E_\eta[\bar X^2]=\mathcal{D}^2$ or
$\mathcal{D}^2+v/H$. The $+v/H$ is just $\mathrm{Var}_\eta[\bar X]\ge0$; it is zero only if the loop
carries no stochastic variance (exact) or $H\to\infty$. That single line decides it.

---

*Refs: one-end/stochastic bias in loop products — Foster & Michael hep-lat/9810021; McNeile & Michael
(one-end trick) hep-lat/0603007; all-to-all/dilution Foley et al. hep-lat/0505023.*
