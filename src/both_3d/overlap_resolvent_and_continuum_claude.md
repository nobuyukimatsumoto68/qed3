# The "resolvent" in the overlap conserved current, and the continuum limit

## 1. What is the resolvent

The overlap operator is built from the **sign function** of the Hermitian Wilson kernel. With
$X=D_W/\lambda_M$ (the two-component Wilson operator scaled by its largest eigenvalue), the overlap needs
the inverse square root $(X^\dagger X)^{-1/2}$ (this is what turns $X$ into $X/\sqrt{X^\dagger X}=\mathrm{sign}$).
It is evaluated by a **rational (Zolotarev) approximation**

$$ (X^\dagger X)^{-1/2} \;\approx\; p_0 \;+\; \sum_{\ell=1}^{\ell_\text{max}} p_\ell\,(X^\dagger X + q_\ell)^{-1}, $$

with positive Zolotarev coefficients $p_\ell, q_\ell$. The objects

$$ R_\ell \;\equiv\; (X^\dagger X + q_\ell)^{-1} $$

are the **resolvents** — the resolvent of the operator $X^\dagger X$ evaluated at the negative shift $-q_\ell$
(well outside the spectrum $X^\dagger X\ge 0$, so no singularity). These are exactly the shifted linear
systems the **multishift CG** solves; one $R_\ell\,\xi$ per pole $\ell$.

## 2. Where the resolvent enters the current

The conserved current is $K^{wz}=\partial D_\text{ov}/\partial\theta_{wz}$ (derivative w.r.t. a link phase).
Differentiating $D_\text{ov}\sim X(X^\dagger X)^{-1/2}$ gives two kinds of terms (Eq. 3.15 of the note):

$$ k^{wz} \;=\; \underbrace{W^{wz}\Big\{p_0+\sum_\ell p_\ell R_\ell\Big\}}_{\text{(I) bare hop, dressed}}
\;-\; \underbrace{X\sum_\ell p_\ell\,R_\ell\big(-W^{wz\dagger}X + X^\dagger W^{wz}\big)R_\ell}_{\text{(II) resolvent dressing}} . $$

- $W^{wz}=\partial X/\partial\theta_{wz}$ is the **bare hopping derivative** — the ultralocal nearest-neighbour
  piece. The **local current** keeps essentially only this ($W^{wz}$ with $\hat e\cdot\sigma$, no $\Omega$, no
  $-r\sigma_0$).
- Term (II) is the **response of the sign function itself** to the gauge field — it is built from the
  resolvents $R_\ell$ on both sides of $W^{wz}$. It is **non-local** (each $R_\ell$ connects distant sites,
  exponentially), and it is exactly the piece that makes $K$ an **exactly conserved** lattice current (the
  Ginsparg-Wilson correction). The local current **drops it entirely** and is therefore *not* conserved.

So "resolvent terms / resolvent dressing" = the parts of $K$ that carry the $R_\ell=(X^\dagger X+q_\ell)^{-1}$
factors — the overlap's non-locality, and the conserving correction.

## 3. Why this flips the sp sign at L=1 (recap)

A small correction cannot flip a sign. At $L=1$, $dt=1$: $|K^{sp}_\text{exact}|=5.3\times10^{-3}$ vs
$|K^{sp}_\text{local}|=8.9\times10^{-3}$ — **comparable** and **opposite**. So at this coarse lattice term (II)
is *not* subleading; it is as big as the bare term and opposite for sp, and it carries the sign. Conservation
ties $J_s$ to $J_t$; the local current breaks conservation (drops (II)), so its spatial sign is unconstrained
and comes out opposite. (For tp the bare term dominates, so there the sign survives.)

## 4. The continuum limit -- what happens as $a\to 0$ ($L\to\infty$)

The resolvent dressing (II) is a **lattice artifact** and vanishes in the continuum:

- The overlap's non-locality has a **fixed range in lattice units** (set by $1/\lambda_\text{min}$, a few
  sites). In **physical** units that range is $a\times(\text{lattice range})\to 0$ as $a\to0$ — the operator
  becomes ultralocal.
- The need for term (II) is the Ginsparg-Wilson defect: the bare point-split current fails to be conserved
  only by $O(a)$ (the GW relation $\{D_\text{ov},\gamma_5\}=a\,D_\text{ov}\gamma_5 D_\text{ov}$). So the
  conserving correction (II) is itself $O(a)$ and $\to 0$.
- Hence **both** currents converge to the same continuum conserved current $\bar\psi\gamma_\mu\psi$:
  $$ K_\text{local} \;\xrightarrow{a\to0}\; K_\text{exact} \;\xrightarrow{a\to0}\; \text{continuum } \bar\psi\gamma_\mu\psi . $$
  The exact one is conserved at every $a$ (so its sign is the physical, stable one); the local one is the
  bare current whose conserving defect shrinks like $O(a)$.

**Consequence (the sp sign discrepancy is a finite-$a$ cutoff effect) and what the DATA says.** As $a\to0$
both currents converge to the same continuum $\bar\psi\gamma_\mu\psi$, so the exact-vs-local sp sign
disagreement at $L=1$ must be a cutoff artifact. The remaining question is *which* sign is the continuum one.

**Measured ($L=1$ exact + local, $L=2$ local; tp is $+$ throughout):**
- $L=1$ exact sp: $+$  |  $L=1$ local sp: $-$  |  **$L=2$ local sp: $-$ (does NOT flip).**

So the **local sp sign is stable $-$** across $L=1,2$. The local current is the *bare* $\hat e\cdot\sigma$
current that $\to$ continuum as $a\to0$, so its stable $-$ is the **continuum-tracking** sign — and it agrees
with the CFT $G_s/G_t=-2$. It is therefore the **exact (conserved) current's $+$ at $L=1$ that is the cutoff
distortion**: at $L=1$ the whole $S^2$ is ~one lattice spacing across, term (II) (the $O(a)$ conserving
dressing) is $O(1)$ and flips the conserved-current correlator's sign away from the continuum; it should turn
$+\to-$ as $a\to0$ (not directly checkable — exact is $L=1$-only).

(EARLIER GUESS, now corrected by the $L=2$ data: I had predicted the *local* sign would flip $-\to+$ toward
the exact. The opposite is observed: local is stably $-$; the exact $+$ is the artifact.)

**Caveat -- magnitude is NOT yet converged.** The *sign* is the only clean readout at this coarseness; the
ratio $|K^{sp}/K^{tp}|$ is still large and noisy ($\approx 7$-$26$ over $dt$ at $L=1,2$, nowhere near $2$).
The magnitude $\to 2$ test needs $L=4$. Also still to rule out: an internal sign bug in
`jj_kbuild_exact`/`jj_contract_deter` for the exact sp (the $w_{sp}$ orientation / contraction) -- but the
stable local $L$-trend already points to $-$ (CFT-consistent) independently of the exact.
