# Symmetries of the exact lattice current correlators -- all mass cases (Sec. 3.3.5)

Working note to fill in Section 3.3.5 of `qed3int_v2-14.pdf` (source `qed3int_v2(7)/main.tex`).
It studies the time-reflection ($t$-) symmetry and the reality of the gauge-averaged
current correlators for the three mass setups (massless, flavor-breaking $m_F$, parity-breaking
$m_P$), for both the vector and the axial channels, and summarizes the independent time regimes.

Supersedes / extends the massless-vector-only note `jj_reality_tsymmetry_claude.md` (which was
written against the v2-11 numbering, where the $t$-reversal relation was Eq. (3.58)); here it is
Eq. (3.68).

Equation numbers below are those of `qed3int_v2-14.pdf`.

---

## 0. Setup and conventions

Per gauge configuration $U$, the connected vector and the axial estimators are (massless,
Eqs. (3.52), (3.55); the disconnected vector piece is (3.53)):
$$
C^{nn'}_{V_{++},c}(t_0,t;U)=\operatorname{tr}\!\big[K^{nn'}(t_0)\,P\,K^{nn'}(t)\,P\big],\qquad P\equiv D_m^{-1},
$$
$$
C^{nn'}_{A_{+-}}(t_0,t;U)=\operatorname{tr}\!\big[K^{nn'}(t_0)\,(1-D_{\rm ov}^\dagger)\,\bar P\,K^{nn'\dagger}(t)\,(1-D_{\rm ov})\,P\big],
$$
with $P,\bar P$ the relevant propagators (Sec. 3.3.2-3.3.4). The kernel $K^{nn'}(t)=\partial D_{\rm ov}/\partial\theta_{nn'}(t)$
sits on a link of timeslice $t$; the temporal kernel $K^{t,t+1}(n)$ is analogous. We suppress the
$nn'$ / link labels and the projection (tp, sp, $Y_{lm}$) below: every relation here holds for any
*diagonal* projection (same kernel at both ends).

After the gauge average $\langle\cdot\rangle_U$, translation invariance in $t$ is restored, so
$$
\langle C(t_0,t;U)\rangle_U \;=\; g(\Delta t),\qquad \Delta t\equiv t-t_0 \pmod{N_t}.
$$
Time reflection acts as $\Delta t\to-\Delta t\equiv N_t-\Delta t$. An *even* function therefore has
fundamental domain $\Delta t\in[0,N_t/2]$; a relation that maps one channel onto another at
$-\Delta t$ leaves one independent complex function on the full $[0,N_t)$.

---

## 1. The three building-block identities -- and on which configs each holds

The whole analysis is built from three identities. The point of Sec. 3.3.5 is to track **which of
them holds per configuration and which only after gauge averaging**, separately for each mass setup.

### (I) Trace cyclicity $\Rightarrow$ $t_0\leftrightarrow t$ exchange -- **per configuration, every mass**

This is purely algebraic (cyclicity of the trace + commutativity of functions of a single operator),
so it holds for **every** $U$ and **every** mass.

- **Vector** (same operator $P$ at both ends, Eqs. (3.52)/(3.53), (3.63)-(3.66)):
$$
C_{V_{\pm\pm},c}(t_0,t;U)=\operatorname{tr}[K(t_0)P\,K(t)P]=\operatorname{tr}[K(t)P\,K(t_0)P]=C_{V_{\pm\pm},c}(t,t_0;U),
$$
and the disconnected piece is a product of two single-time traces, manifestly $t_0\leftrightarrow t$
symmetric. This is the paper's statement "$C^{nn'}_{V_{\pm\pm},c/d}$ are strictly symmetric in $t$ and $t_0$."
After averaging: $g_{\pm\pm}(\Delta t)=g_{\pm\pm}(-\Delta t)$ -- **even**.

- **Axial** (different structure at the two ends): cyclicity does *not* give $t_0\leftrightarrow t$
within one channel; instead it relates the two axial channels. Swapping $+\leftrightarrow-$ swaps
$K\leftrightarrow K^\dagger$ and the GW factors, and one finds
$$
C_{A_{+-}}(t,t_0;U)=\operatorname{tr}[K(t)(1-D^\dagger)\bar P\,K^\dagger(t_0)(1-D)P]
=\operatorname{tr}[K^\dagger(t_0)(1-D)P\,K(t)(1-D^\dagger)\bar P]=C_{A_{-+}}(t_0,t;U),
$$
i.e. Eq. (3.56), $\;C_{A_{-+}}(t_0,t;U)=C_{A_{+-}}(t,t_0;U)$. After averaging:
$$
g_{A_{-+}}(\Delta t)=g_{A_{+-}}(-\Delta t).
$$
This is the cyclicity-only relation; it holds for $m_P\neq0$ too (the $m_P$ axial uses
$\bar P=\tilde D_{m_P}^{\dagger-1}$, $P=D_{m_P}^{-1}$, but the manipulation is unchanged).

### (II) Per-configuration complex conjugation -- **vector massless / $m_F$ only**

Conjugating a trace gives $\operatorname{tr}[X]^*=\operatorname{tr}[X^\dagger]$. For the vector, with $P=D_m^{-1}$ and
$P^\dagger=D_m^{\dagger-1}$,
$$
C_{V_{++},c}(t_0,t;U)^*=\operatorname{tr}[K^\dagger(t_0)\,D_m^{\dagger-1}\,K^\dagger(t)\,D_m^{\dagger-1}]
=C_{V_{--},c}(t_0,t;U),
$$
because the $V_{--}$ contraction uses precisely the conjugate propagator $\eta\xi^\dagger=(D_m^{-1})^\dagger$
and $K^\dagger$. This is Eq. (3.54), $\;C_{V_{--},c/d}=C_{V_{++},c/d}^*$, **per configuration**.

- **Massless** ($D_m=D_{\rm ov}$) and **$m_F$** ($D_m=D_{m_F}$, $m_F\in\mathbb R$, so $D_{m_F}^\dagger=D_{\rm ov}^\dagger+m_F$):
  identity (II) holds -- the replacement $D_{\rm ov}\to D_{m_F}$ in (3.52)-(3.54) is exact.
- **$m_P$**: identity (II) **fails**. The $V_{--}$ correlator (3.65)/(3.66) is built from
  $\tilde D_{m_P}^{\dagger-1}$, not from $D_{m_P}^{\dagger-1}$. Concretely $m_P^*=-m_P$ gives
  $D_{m_P}^\dagger=D_{\rm ov}^\dagger-m_P$ in $C_{V_{++}}^*$, whereas $C_{V_{--}}$ carries
  $\tilde D_{m_P}^\dagger=D_{\rm ov}^\dagger-\tfrac{m_P}{1+m_P}$. Different shifts $\Rightarrow$
  $C_{V_{--}}\neq C_{V_{++}}^*$. This is the paper's "$C_{V_{++}}$ and $C_{V_{--}}$ are not simply related."

- **Axial, any mass**: there is **no** per-configuration reality relation. Daggering
  $C_{A_{+-}}$ produces
  $$
  C_{A_{+-}}(t_0,t;U)^*=\operatorname{tr}[K^\dagger(t_0)(1-D^\dagger)\bar G\,K(t)(1-D)G],\qquad G=D_{\rm ov}^{-1},\ \bar G=D_{\rm ov}^{\dagger-1},
  $$
  which matches **neither** $C_{A_{+-}}$ nor $C_{A_{-+}}$ (it is a distinct "dagger-on-both-kernels"
  object). The mixed $(1-D^\dagger)\bar G\,/\,(1-D)G$ structure spoils (II).

### (III) Time reversal (3.68) -- **operator identity always; usable after $\langle\cdot\rangle_U$ only if the measure is $T$-invariant**

Eq. (3.68) is a property of the overlap operator itself, valid for every $U$:
$$
\sigma_3\,D_{\rm ov}[U_{\cal T}]\,\sigma_3=D_{\rm ov}[U]^\dagger
\;\Longleftrightarrow\;
\begin{cases}
D_{\rm ov}[U]^{-1}=\sigma_3\,D_{\rm ov}[U_{\cal T}]^{\dagger-1}\,\sigma_3,\\[2pt]
(1-D_{\rm ov}[U])=\sigma_3\,(1-D_{\rm ov}[U_{\cal T}]^\dagger)\,\sigma_3,\\[2pt]
D_{\rm ov}[U]^{\dagger-1}=\sigma_3\,D_{\rm ov}[U_{\cal T}]^{-1}\,\sigma_3,
\end{cases}
$$
where $U_{\cal T}$ is the time-reflected gauge field. Differentiating w.r.t. the (real) link angle
gives the kernel transform
$$
\sigma_3\,K[U_{\cal T}](\bar t)\,\sigma_3=\varepsilon\,K[U](t)^\dagger,\qquad \varepsilon=\pm1,
$$
$\bar t$ being the reflected timeslice. The sign $\varepsilon$ ($-$ for the spatial current, $+$ for
the temporal/charge-density component) is **immaterial for every diagonal correlator**: the same
kernel sits at both ends, so $\varepsilon$ enters squared, $\varepsilon^2=1$. (It would matter only
for the off-diagonal $\langle j^0 j^i\rangle$, which is not measured.)

**(III) is usable on an observable only after the gauge average, and only if the gauge measure is
invariant under $U\to U_{\cal T}$.** This is the per-config-vs-after-average distinction that decides
the whole table:

| measure $T$-invariant? | why | (III) usable? |
|---|---|---|
| **massless** ($m_F=m_P=0$) | $\det D_{\rm ov}\det D_{\rm ov}^\dagger=\lvert\det D_{\rm ov}\rvert^2$ real-positive (even $N_f$, parity anomaly cancels), $R\to1$ | **yes** |
| **$m_F\neq0$** | $\sigma_{PS}=\bar\Psi\Psi$ is the **parity-even** mass; $\det(D_{\rm ov}+m_F)\det(D_{\rm ov}^\dagger+m_F)=\lvert\det(D_{\rm ov}+m_F)\rvert^2$ real-positive ($m_F\in\mathbb R$) | **yes** |
| **$m_P\neq0$** | $\sigma_{FS}$ is the **parity-odd** mass; the determinant / reweighting $R$ (Eq. (2.4)) carries a genuine phase (Chern-Simons / parity anomaly), $R$ is not real | **no** |

Time reflection in 3D Euclidean is a single-coordinate reflection, in the same class as parity; the
measure is $T$-invariant exactly when parity is unbroken. $m_F$ is parity-even, so it keeps $T$;
$m_P$ is parity-odd, so it breaks $T$ -- this is the same statement as "$R$ is a phase" below Eq. (2.5).

---

## 2. What (III) gives, channel by channel (massless / $m_F$, $T$-invariant measure)

Substitute the (III) rewrites of every factor (all in terms of $V\equiv U_{\cal T}$), let the
adjacent $\sigma_3$'s and the $\varepsilon^2=1$ collapse, then rename $V\to U$ under
$\langle f(U_{\cal T})\rangle_U=\langle f(U)\rangle_U$.

### Vector (two $P$'s, no daggers $\Rightarrow$ a genuine conjugation)
$$
C_{V_{++},c}(t_0,t;U)
=\operatorname{tr}\!\big[K[V](\bar t_0)^\dagger D[V]^{\dagger-1}K[V](\bar t)^\dagger D[V]^{\dagger-1}\big]
=\operatorname{tr}\!\big[D[V]^{-1}K[V](\bar t)D[V]^{-1}K[V](\bar t_0)\big]^{*},
$$
i.e. $\;C_{V_{++},c}(t_0,t;U)=C_{V_{++},c}(\bar t,\bar t_0;U_{\cal T})^{*}$. Averaging
($\bar t_0-\bar t=\Delta t$):
$$
\boxed{\,g_{++}(\Delta t)=g_{++}(\Delta t)^{*}\,}\quad\Rightarrow\quad \mathrm{Im}\,g_{++}=0\ \text{(exactly)}.
$$
Combined with cyclicity (I) ($g_{++}$ even) and (II) ($C_{V_{--}}=C_{V_{++}}^*\Rightarrow g_{--}=g_{++}^*=g_{++}$):
**the vector correlator is real, even, and $g_{++}=g_{--}$.** A single real even function.

### Axial (mixed $(1-D^\dagger)\bar G\,/\,(1-D)G$ $\Rightarrow$ NO conjugation)
The daggers in (III) are "used up" turning $D^\dagger$ into $D$; the result is a *direct* (un-starred)
trace:
$$
C_{A_{+-}}(t_0,t;U)=\operatorname{tr}\!\big[K[V](\bar t_0)^\dagger(1-D[V])D[V]^{-1}K[V](\bar t)(1-D[V]^\dagger)D[V]^{\dagger-1}\big]
=C_{A_{-+}}(\bar t_0,\bar t;U_{\cal T}).
$$
Averaging gives $g_{A_{+-}}(\Delta t)=g_{A_{-+}}(-\Delta t)$ -- **identical to cyclicity (3.56)**. So
for the axial **(III) yields nothing beyond (I)**: there is no reality. The gauge-averaged axial
correlator stays **complex**; the only relation is the channel reflection $g_{A_{+-}}(\Delta t)=g_{A_{-+}}(-\Delta t)$.

> Why the asymmetry between V and A: the vector trace has an *even* number of bare propagators and no
> $(1-D)$ factors, so (III) maps it to its own conjugate. The axial already mixes $D$ and $D^\dagger$
> (the Ginsparg-Wilson dressing $(1-D_{\rm ov})$), so (III) merely swaps $+\!-\leftrightarrow-\!+$
> with no leftover conjugation. The residual $\mathrm{Im}$ of the axial is a finite-$a$ (GW) effect
> that is **not** forced to zero by any exact lattice symmetry.

---

## 3. Case-by-case summary

### 3a. Massless ($m_F=m_P=0$) -- all currents conserved
- **Vector $V_{\pm\pm}$:** (I) even + (II) $C_{--}=C_{++}^*$ + (III) reality. Net: **real, even,
  $g_{++}=g_{--}$**. One real even function on $[0,N_t/2]$. The per-config $\mathrm{Im}\,C_{V_{++}}\neq0$
  averages to **exactly zero** (not $O(a^2)$): a free noise/consistency gauge.
- **Axial $A_{\pm\mp}$:** (I)/(3.56) only; (II) and (III) give no reality. Net: **complex**, with
  $g_{A_{+-}}(\Delta t)=g_{A_{-+}}(-\Delta t)$. One independent complex function over $[0,N_t)$.

### 3b. Flavor-breaking ($m_F\neq0$, $m_P=0$) -- $J_V,J_r$ conserved; $J_A$ NOT conserved
The $m_F$ term breaks the axial (chiral) symmetry, so there is **no exact conserved axial current**.
- **Vector $V_{\pm\pm}$:** identical structure to massless under $D_{\rm ov}\to D_{m_F}$. (II) holds
  ($m_F$ real) and the measure is $T$-invariant ($m_F$ parity-even), so (III) holds. Net: **real,
  even, $g_{++}=g_{--}$** -- same as massless.

#### Axial with the *massless kernel* (massless-limit estimator)
Even though $J_A$ is not conserved at $m_F\neq0$, we still measure an axial-axial correlator to
estimate the massless limit: we keep the **massless axial operator** -- the massless kernel $K$ and
the massless Ginsparg-Wilson dressing $(1-D_{\rm ov})$, exactly as in Eq. (3.55) -- and only replace
the propagators by the massive ones $D_{m_F}^{-1}$:
$$
C^{(m_F)}_{A_{+-}}(t_0,t;U)=\operatorname{tr}\!\big[K(t_0)\,(1-D_{\rm ov}^\dagger)\,D_{m_F}^{\dagger-1}\,K^\dagger(t)\,(1-D_{\rm ov})\,D_{m_F}^{-1}\big],
$$
$$
C^{(m_F)}_{A_{-+}}(t_0,t;U)=\operatorname{tr}\!\big[K^\dagger(t_0)\,(1-D_{\rm ov})\,D_{m_F}^{-1}\,K(t)\,(1-D_{\rm ov}^\dagger)\,D_{m_F}^{\dagger-1}\big].
$$
So the GW factors stay massless while the bare propagator slots carry $D_{m_F}^{\mp\dagger-1}$.
Its symmetry properties are **identical to the massless axial**:

- **(I) cyclicity / (3.56) -- per config.** Cyclicity is algebraic, independent of which operator
  sits where, so $C^{(m_F)}_{A_{-+}}(t_0,t;U)=C^{(m_F)}_{A_{+-}}(t,t_0;U)$ holds exactly.
- **(II) -- none.** The mixed $(1-D_{\rm ov}^\dagger)D_{m_F}^{\dagger-1}/(1-D_{\rm ov})D_{m_F}^{-1}$
  structure spoils per-config conjugation, just as in the massless axial.
- **(III) -- usable, and it reproduces (3.56).** This is the key point that makes the estimator a
  clean massless-limit probe: the measure is $T$-invariant ($m_F$ parity-even), **and both** the
  massless GW factor and the massive propagator obey the *same* relation (3.68),
  $$
  \sigma_3\,(1-D_{\rm ov}[U_{\cal T}])\,\sigma_3=1-D_{\rm ov}[U]^\dagger,\qquad
  \sigma_3\,D_{m_F}[U_{\cal T}]\,\sigma_3=D_{m_F}[U]^\dagger
  $$
  (the latter because $m_F\in\mathbb R$ is a $\sigma_3$-scalar shift). The substitution therefore
  goes through exactly as in the massless axial of Sec. 2, mapping $+\!-\to-\!+$ with **no leftover
  conjugation**:
  $$
  C^{(m_F)}_{A_{+-}}(t_0,t;U)=C^{(m_F)}_{A_{-+}}(\bar t_0,\bar t;U_{\cal T})
  \;\xrightarrow{\langle\cdot\rangle_U}\;
  g^{(m_F)}_{A_{+-}}(\Delta t)=g^{(m_F)}_{A_{-+}}(-\Delta t),
  $$
  i.e. (III) gives nothing beyond (3.56).

**Net:** the massless-kernel axial estimator at $m_F$ is **complex**, even-channel-paired by
$g^{(m_F)}_{A_{+-}}(\Delta t)=g^{(m_F)}_{A_{-+}}(-\Delta t)$ -- the **same** symmetry table row as the
genuine massless axial. Since parity is unbroken ($m_F$ parity-even, $T$-invariant measure), the
residual $\mathrm{Im}$ is again an unprotected finite-$a$ (GW) effect, **not** a parity signal; the
$m_F\to0$ extrapolation of this estimator inherits exactly the massless-axial symmetry structure.
(Had the propagator $D_{m_F}$ failed to satisfy (3.68) -- e.g. a parity-odd mass -- (III) would drop
and only (3.56) would survive; that is the $m_P$ case in 3c.)

### 3c. Parity-breaking ($m_P\neq0$, $m_F=0$) -- all currents conserved, but $T$ broken
- **Vector $V_{\pm\pm}$:** (I) even holds (cyclicity). **(II) fails** ($C_{--}$ uses $\tilde D_{m_P}$,
  not $D_{m_P}^\dagger$) and **(III) is unavailable** (measure not $T$-invariant: $R$ is a phase). Net:
  $C_{V_{++}}$ and $C_{V_{--}}$ are **two independent, even, complex** functions on $[0,N_t/2]$.
  The nonzero $\mathrm{Im}$ and the inequality $C_{V_{++}}\neq C_{V_{--}}^*$ are the **genuine
  parity-breaking signal** (the valence Chern-Simons effect, cf. Eqs. (2.4)/(2.5)).
- **Axial $A_{\pm\mp}$:** (I)/(3.56) still holds per config (cyclicity is mass-independent);
  (II)/(III) give no reality. Net: **complex**, $g_{A_{+-}}(\Delta t)=g_{A_{-+}}(-\Delta t)$ -- same
  form as the massless axial.

---

## 4. Master table

Per-configuration identities are marked "(cfg)"; identities needing the gauge average are marked
"($\langle\rangle_U$)". "even" $=g(\Delta t)=g(-\Delta t)$.

| Case | Channel | $t_0\!\leftrightarrow\! t$ (I, cfg) | conjugation (II) | $T$-reversal (III, $\langle\rangle_U$) | gauge-averaged result |
|---|---|---|---|---|---|
| **massless** | $V_{++},V_{--}$ | even (cfg) | $C_{--}=C_{++}^*$ (cfg) | usable $\Rightarrow$ real | **real, even, $g_{++}=g_{--}$** |
| | $A_{+-},A_{-+}$ | $A_{-+}(t_0,t)=A_{+-}(t,t_0)$ (cfg) | none | reproduces (3.56) only | **complex**, $g_{A_{+-}}(\Delta t)=g_{A_{-+}}(-\Delta t)$ |
| **$m_F\neq0$** | $V_{++},V_{--}$ | even (cfg) | $C_{--}=C_{++}^*$ (cfg) | usable $\Rightarrow$ real | **real, even, $g_{++}=g_{--}$** (as massless) |
| | $A_{+-},A_{-+}$ (massless kernel, $D_{m_F}^{-1}$ props) | $A_{-+}(t_0,t)=A_{+-}(t,t_0)$ (cfg) | none | reproduces (3.56) only | **complex**, $g_{A_{+-}}(\Delta t)=g_{A_{-+}}(-\Delta t)$ (as massless axial; $J_A$ not conserved -- massless-limit estimator) |
| **$m_P\neq0$** | $V_{++}$ | even (cfg) | fails | **not usable** ($T$ broken) | **complex, even**; independent |
| | $V_{--}$ | even (cfg) | fails | **not usable** | **complex, even**; independent of $V_{++}$ |
| | $A_{+-},A_{-+}$ | $A_{-+}(t_0,t)=A_{+-}(t,t_0)$ (cfg) | none | not usable | **complex**, $g_{A_{+-}}(\Delta t)=g_{A_{-+}}(-\Delta t)$ |

---

## 5. Independent time regimes

| Case | Channel(s) | independent data | fundamental domain |
|---|---|---|---|
| massless / $m_F$ | vector ($V_{++}=V_{--}$) | **1 real** function | $\Delta t\in[0,N_t/2]$ |
| massless | axial pair | **1 complex** function ($A_{-+}$ fixed by reflection of $A_{+-}$) | $\Delta t\in[0,N_t)$ |
| $m_F$ | axial pair (massless-kernel estimator) | **1 complex** function ($A_{-+}$ = reflected $A_{+-}$) -- as massless | $\Delta t\in[0,N_t)$ |
| $m_P$ | vector $V_{++}$ | **1 complex** function (Re + Im) | $\Delta t\in[0,N_t/2]$ |
| $m_P$ | vector $V_{--}$ | **1 complex** function, independent of $V_{++}$ | $\Delta t\in[0,N_t/2]$ |
| $m_P$ | axial pair | **1 complex** function ($A_{-+}$ = reflected $A_{+-}$) | $\Delta t\in[0,N_t)$ |

Folding rule: an **even** channel collapses to $[0,N_t/2]$; a **channel-reflection** relation
($A_{+-}\!\leftrightarrow\!A_{-+}$) leaves the full $[0,N_t)$ but ties the two channels, so the pair
carries one function's worth of independent data.

---

## 6. Data diagnostics implied by the table

- **Massless / $m_F$ vector:** $\langle\mathrm{Im}\,C_{V_{++}}\rangle_U$ is **exactly zero** in
  expectation (not an $O(a^2)$ residual) -- pure Monte-Carlo noise; a stable nonzero value flags a
  bug or a non-$T$-invariant measure. Check $g_{++}=g_{--}$ and reflection-evenness about $N_t/2$.
- **$m_P$ vector:** $\mathrm{Im}\,C_{V_{++}}$, $\mathrm{Im}\,C_{V_{--}}$, and the residual
  $C_{V_{++}}-C_{V_{--}}^*$ are **physical** -- the parity signal. Each channel is still even in
  $\Delta t$ (cyclicity), which remains a valid check.
- **Axial (massless, $m_F$ massless-kernel estimator, and $m_P$):** expect a nonzero $\mathrm{Im}$
  (no exact reality). The only exact check is the channel reflection
  $C_{A_{+-}}(\Delta t)=C_{A_{-+}}(-\Delta t)$ (Eq. (3.56)), which holds per configuration. For the
  $m_F$ massless-kernel estimator this reflection is additionally enforced by (III) (the measure is
  $T$-invariant), and the residual $\mathrm{Im}$ is a finite-$a$/noise quantity (parity unbroken) --
  it must extrapolate smoothly with $m_F\to0$, providing a check on the massless-limit fit.

---

### Provenance of each statement (separation by config)

- Holds **per configuration, all masses**: identity (I) -- vector $t_0\leftrightarrow t$ evenness and
  the axial relation (3.56).
- Holds **per configuration, vector massless/$m_F$ only**: identity (II) -- (3.54).
- Holds **only after $\langle\cdot\rangle_U$, and only for the $T$-invariant measures
  (massless, $m_F$)**: identity (III) -- vector reality; for the axial it degenerates to (3.56).
- For $m_P$: **(II) and (III) both drop**, leaving only (I); hence vector $++/--$ split into two
  independent complex channels and the imaginary parts become physical.
