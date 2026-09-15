# Continuum free propagator on $S^2\times\mathbb R$: the $k$-integral, per $(m,n,\iota_3)$

Companion derivation for `diagramA.nb`. Goal: eliminate the temporal (Matsubara) momentum $k$
**analytically** while **keeping $\iota_3$ un-summed**, so the eigenvalue-sign branch can be tracked through
the loop (it is not conserved at vertices). Summing $\iota_3$ afterwards reproduces Eq C.29.

## Modes and propagator

Eigenmode (labels: branch $\text{isc}=\pm1$, sign $\iota_3=\pm1$, $\iota_m=\mathrm{sgn}(m)$, temporal $k$,
spatial $(|m|,n)$), $z=\cos\theta$:
$$
\psi_{\text{isc},\iota_3,\iota_m,k,m,n}=e^{\,i(\text{isc}\,\iota_3\,k\,t+\iota_m|m|\phi)}\,\chi^{(\text{isc})}(z),\qquad
\chi^{(+)}=\begin{pmatrix}\cos\tfrac\gamma2\,A_0\\ \iota_3\sin\tfrac\gamma2\,B_0\end{pmatrix},\ 
\chi^{(-)}=\begin{pmatrix}\sin\tfrac\gamma2\,A_0\\ \iota_3\cos\tfrac\gamma2\,B_0\end{pmatrix},
$$
$$
A_0=\xi_{|m|,n}(\iota_m z),\quad B_0=\iota_m i(-1)^n\xi_{|m|,n}(-\iota_m z),\quad
\lambda=n+|m|+\tfrac12,\ \Lambda=\sqrt{k^2+\lambda^2},\ \cos\gamma=\tfrac k\Lambda,
$$
$$
\cos^2\tfrac\gamma2=\tfrac12(1+\tfrac k\Lambda),\ \sin^2\tfrac\gamma2=\tfrac12(1-\tfrac k\Lambda),\ \cos\tfrac\gamma2\sin\tfrac\gamma2=\tfrac{\lambda}{2\Lambda}.
$$
Verified: $D\psi=i\,\iota_3\Lambda\,\psi$, so $G=\sum_K\psi_K\psi_K^\dagger/(i\,\iota_3\Lambda\,\mathcal N)$,
$\mathcal N=4\pi c^2_{|m|,n}$. Absorb the phase into
$A=A_0e^{im\phi}$, $B=B_0e^{im\phi}$ ($A$ real$\times$phase, $B$ imaginary$\times$phase, $B_0^*=-B_0$).

## Outer products $\chi(z)\chi(z')^\dagger$ (temporal phase $e^{\,i\,\text{isc}\,\iota_3 k\tau}$, $\tau=t-t'$)

| element | $\text{isc}=+$ coeff | $\text{isc}=-$ coeff | spatial |
|---|---|---|---|
| $(1,1)$ | $\cos^2\tfrac\gamma2$ | $\sin^2\tfrac\gamma2$ | $AA'^*$ |
| $(2,2)$ | $-\sin^2\tfrac\gamma2$ | $-\cos^2\tfrac\gamma2$ | $BB'^*$ |
| $(1,2)$ | $-\iota_3\cos\tfrac\gamma2\sin\tfrac\gamma2$ | $-\iota_3\cos\tfrac\gamma2\sin\tfrac\gamma2$ | $AB'^*$ |
| $(2,1)$ | $+\iota_3\cos\tfrac\gamma2\sin\tfrac\gamma2$ | $+\iota_3\cos\tfrac\gamma2\sin\tfrac\gamma2$ | $BA'^*$ |

## The $\text{isc}$ sum only ($\iota_3$ kept fixed)

Weight $1/(i\,\iota_3\Lambda)$. Sum $\text{isc}=\pm$, $\iota_3$ **fixed**.

**Diagonal $(1,1)$:** $\displaystyle\sum_\text{isc}\text{coeff}\,e^{i\,\text{isc}\,\iota_3 k\tau}
=\cos^2\tfrac\gamma2 e^{i\iota_3k\tau}+\sin^2\tfrac\gamma2 e^{-i\iota_3k\tau}
=\cos(k\tau)+i\,\tfrac k\Lambda\,\iota_3\sin(k\tau)$. Dividing by $i\iota_3\Lambda$:
$$
(1,1)=AA'^*\!\int\!\frac{dk}{2\pi}\Big[\underbrace{\frac{\cos(k\tau)}{i\iota_3\Lambda}}_{\text{Bessel},\ \propto\iota_3}
+\underbrace{\frac{k\sin(k\tau)}{\Lambda^2}}_{\text{exp}}\Big].
$$
$(2,2)$ is identical with $\cos\to-\cos$ on the Bessel piece (opposite $K_0$ sign), spatial $BB'^*$.

**Off-diagonal $(1,2)$:** the coefficient's $\iota_3$ cancels the $1/\iota_3$, and $\sum_\text{isc}$ gives
$2\cos(k\tau)$ ($\iota_3$ drops from $\cos$). So **the off-diagonal is $\iota_3$-independent**:
$$
(1,2)=AB'^*\!\int\!\frac{dk}{2\pi}\frac{-2\cos\tfrac\gamma2\sin\tfrac\gamma2\cos(k\tau)}{i\Lambda}
=AB'^*\!\int\!\frac{dk}{2\pi}\frac{-\lambda\cos(k\tau)}{i\Lambda^2}.
$$

## $k$-integrals ($N_t\to\infty$)

$$
\int\!\frac{dk}{2\pi}\frac{\cos(k\tau)}{\Lambda}=\frac{K_0(\lambda|\tau|)}{\pi}\ \ (\text{Bessel}),\qquad
\int\!\frac{dk}{2\pi}\frac{k\sin(k\tau)}{\Lambda^2}=\tfrac12\operatorname{sgn}(\tau)e^{-\lambda|\tau|},\qquad
\int\!\frac{dk}{2\pi}\frac{\cos(k\tau)}{\Lambda^2}=\frac{e^{-\lambda|\tau|}}{2\lambda}.
$$

## Result — per-$(m,n,\iota_3)$, NOT summed over $\iota_3$

$$
\boxed{\;G^{(m,n,\iota_3)}(x,x';\tau)=\frac{1}{4\pi c^2_{|m|,n}}
\begin{pmatrix}
\tfrac12\operatorname{sgn}(\tau)e^{-\lambda|\tau|}AA'^*-\dfrac{i\,\iota_3}{\pi}K_0(\lambda|\tau|)\,AA'^* &
\ \dfrac{i}{2}\,e^{-\lambda|\tau|}AB'^*\\[10pt]
\dfrac{i}{2}\,e^{-\lambda|\tau|}BA'^* &
\ \tfrac12\operatorname{sgn}(\tau)e^{-\lambda|\tau|}BB'^*+\dfrac{i\,\iota_3}{\pi}K_0(\lambda|\tau|)\,BB'^*
\end{pmatrix}\;}
$$
($1/(i\iota_3)=-i\iota_3$ and $-1/i=i$ used; $-\lambda/i\cdot\tfrac{1}{2\lambda}=\tfrac i2$.)

- **Diagonal** carries a **Bessel $K_0(\lambda|\tau|)$** term $\propto\iota_3$ (imaginary, $\iota_3$-odd), on
  top of the $\tfrac12\operatorname{sgn}(\tau)e^{-\lambda|\tau|}$ propagating part.
- **Off-diagonal** is $\iota_3$-**independent**, Bessel-free, even in $\tau$: $\tfrac i2 e^{-\lambda|\tau|}$.
- **Summing $\iota_3=\pm$**: the $K_0$ pieces cancel, the rest doubles → Eq C.29
  $\big(\tfrac{e^{-\lambda|\tau|}}{4\pi c^2}[\operatorname{sgn}(\tau)AA'^*,-iAB'^*;-iBA'^*,-\operatorname{sgn}(\tau)BB'^*]\big)$.
  Keeping $\iota_3$ (as here) retains the $K_0$.

## Consequences for the loop

- $\iota_3$ is a **fixed label per propagator segment** here; it changes only through vertex overlaps
  (which are $\iota_3$-mixing — the off-diagonal $A\!\leftrightarrow\!B$ and $\tilde\tau$), consistent with
  "$\iota_3$ is not conserved around the loop".
- The **$m$-conservation** at each $\ell=0$ vertex is exact ($\int d\phi\,e^{i(m'-m)\phi}$); the $z$-integrals
  give the $\xi$ overlaps.
- **Equal-time $\tilde\tau=G_\text{eq}-\tfrac12$**: at $\tau\to0$, $\operatorname{sgn}(0)$ on the diagonal is the
  contact ($\tfrac12$) that $\tilde\tau$ removes; $K_0(\lambda|\tau|)$ **diverges** as $\tau\to0$
  ($\sim-\ln|\tau|$), so the per-$\iota_3$ equal-time diagonal needs its short-distance piece handled with
  care (this divergence is exactly what cancels between $\iota_3=\pm$; keep it symbolic and let it cancel
  when the loop closes). The off-diagonal is finite at $\tau=0$.

## Finite $N_t$ (thermal) version

Antiperiodic $k=\tfrac{2\pi}{N_t a_t}(\mathbb Z+\tfrac12)$: replace the continuum $k$-integrals by the
fermionic Matsubara sums (closed form $\cosh/\sinh$ ratios about $\lambda\tfrac{N_t a_t}{2}$), reducing to the
above as $N_t\to\infty$. Use finite $N_t=128$ to compare with the perambulator; use $N_t\to\infty$ for the
clean diagram-A cut.

Refs: `qed3_v2-6.pdf` App C (C.10, C.17, C.28, C.29); `diagramA.nb`.
