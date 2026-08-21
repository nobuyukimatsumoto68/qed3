# Leading masses for the SP channel, $\ell=0,1$ (from qed3int_v2-15 Eqs. 4.26, 4.27, 4.33)

Goal: by analogy with the paper's TP derivation (Eq. 4.33 $\to$ 4.34), get the leading exponents
(masses) of the **spatially projected** correlator $G_s$ in the $\ell=0$ and $\ell=1$ channels.

## Setup (paper conventions)

Cylinder ($S^2\times\mathbb{R}$) correlator of the conserved current, Eq. (4.26):
$$
f_{\rm cyl}^{ab}(t;\hat n_1,\hat n_2)=C_j\,\Lambda^a{}_\alpha(\hat n_1)\Lambda^b{}_\beta(\hat n_2)\,
\Big[\delta_{\alpha\beta}-\frac{2\,v_\alpha v_\beta}{w}\Big]\,w^{-\Delta}\,e^{-\Delta t},
$$
with $v=\hat n_1-e^{-t}\hat n_2$, $w=|v|^2=1-2\,\hat n_{12}e^{-t}+e^{-2t}$, $\hat n_{12}=\hat n_1\!\cdot\!\hat n_2\equiv x$.
Frame: $\hat e^{1}=\hat\theta,\ \hat e^{2}=\hat\phi,\ \hat e^{3}=\hat n$ (the "$3$" axis is radial $=$ temporal).
$\Delta$ = current dimension ($\Delta=D-1=2$ in $D=3$ for the conserved vector current).

- TP (Eq. 4.30): $G_t=e^a_3 e^b_3 f^{ab}=f^{33}$.
- SP (Eq. 4.27): $G_s=\dfrac{\delta^{ab}-e^a_3 e^b_3}{D-1}f^{ab}=\dfrac{1}{2}\big(f^{11}+f^{22}\big)$.

So
$$
G_s=\tfrac12\,S_{\alpha\beta}\,T_{\alpha\beta}\,w^{-\Delta}e^{-\Delta t}\,C_j ,\qquad
S_{\alpha\beta}\equiv\hat\theta_{1\alpha}\hat\theta_{2\beta}+\hat\phi_{1\alpha}\hat\phi_{2\beta},\quad
T_{\alpha\beta}=\delta_{\alpha\beta}-\frac{2v_\alpha v_\beta}{w}.
$$
$S$ is the **tangential frame overlap** (sum over the two spatial frame legs).

## Key contractions

Because $\hat e^{i}(\hat n)\!\cdot\!\hat n=0$ for $i=\theta,\phi$:
$$
S_{\alpha\beta}\,\hat n_{1\alpha}\hat n_{1\beta}=0,\quad
S_{\alpha\beta}\,\hat n_{2\alpha}\hat n_{2\beta}=0,\quad
S_{\alpha\beta}\,\hat n_{1\alpha}\hat n_{2\beta}=0,
$$
$$
S_{\alpha\alpha}=\hat\theta_1\!\cdot\!\hat\theta_2+\hat\phi_1\!\cdot\!\hat\phi_2\equiv Q,\qquad
S_{\alpha\beta}\,\hat n_{2\alpha}\hat n_{1\beta}\equiv W .
$$
Only the $\hat n_2\hat n_1$ structure survives the contraction with $S$.

Let $\varepsilon\equiv e^{-t}$. Expanding $T$ and $w^{-\Delta}$:
$$
S_{\alpha\beta}T_{\alpha\beta}=Q+2W\varepsilon+4Wx\,\varepsilon^2+O(\varepsilon^3),\qquad
w^{-\Delta}=1+2\Delta x\,\varepsilon+\big(2\Delta(\Delta{+}1)x^2-\Delta\big)\varepsilon^2+O(\varepsilon^3).
$$
Hence
$$
G_s=\tfrac{C_j}{2}\,\varepsilon^{\Delta}\Big[\,Q+\big(2W+2\Delta x\,Q\big)\varepsilon
+O(\varepsilon^2)\Big].
$$

## $\ell$-projection (Eq. 4.33 with scalar $Y_{\ell m}$)

The diagonal projection $G^{\ell}\equiv\int\!\frac{d^2\hat n_1 d^2\hat n_2}{4\pi}\,Y_{\ell m}(\hat n_1)Y_{\ell m}(\hat n_2)\,G_s$.

**Leading term $\propto Q\,\varepsilon^\Delta$** (with $Q=\hat\theta_1\!\cdot\!\hat\theta_2+\hat\phi_1\!\cdot\!\hat\phi_2$):
- $\ell=1$: the tangential (VSH) overlap $Q$ has a nonzero scalar-$Y$ $\ell=1$ coefficient.
  **SP $\ell=1$ leads at $e^{-\Delta t}$, mass $\Delta$.**
- $\ell=0$: $G^{0}\propto\int\!\!\int Q=\big(\int\hat\theta_1\big)\!\cdot\!\big(\int\hat\theta_2\big)+\big(\int\hat\phi_1\big)\!\cdot\!\big(\int\hat\phi_2\big)$.
  **CORRECTION (was wrong): this is NOT zero.** A frame vector does *not* integrate to zero —
  $\int d^2\hat n\,\hat\theta_z=2\pi\int_0^\pi(-\sin\theta)\sin\theta\,d\theta=-\pi^2\neq0$ (the $\hat\theta,\hat\phi$
  frame is singular at the poles, so $\hat\theta_z=-\sin\theta$ is even under $\theta\to\pi-\theta$ and
  survives). Hence the $\ell=0$ coefficient of $Q$ is **nonzero** and **SP $\ell=0$ ALSO leads at
  $e^{-\Delta t}$, mass $\Delta$.**

## RESULT (leading masses)

$$
\boxed{\;m_{\rm sp,\ \ell=0}=m_{\rm sp,\ \ell=1}=\Delta\;}
$$
i.e. for the current ($\Delta=2$): **SP $\ell=0\to2$ and SP $\ell=1\to2$**, both equal to TP $\ell=1$.
This is consistent with the measured $\mathrm{sp}\,\ell0/\mathrm{tp}\,\ell1\approx1$ and
$\mathrm{sp}\,\ell1/\mathrm{tp}\,\ell1\approx1$.

Contrast with TP (Eq. 4.34): $m_{\rm tp,\ell=1}=\Delta$, $m_{\rm tp,\ell=2}=\Delta+1$, and
$m_{\rm tp,\ell=0}$ *absent* (coefficient $\propto(2-\Delta)=0$; charge conservation $Q|0\rangle=0$).
The SP $\ell=0$ is not protected by charge conservation and appears at the *same* leading order $\Delta$.

(EARLIER ERROR, now fixed: I claimed $\int\hat e^i=0\Rightarrow$ SP $\ell=0$ at $\Delta+1$. Wrong — the
$z$-component survives, so SP $\ell=0$ is at $\Delta$. The nonzero $\int\hat\theta$ is the coordinate-frame
artifact behind the "VSH mix", but it lands $\ell=0$ at $\Delta$, not $\Delta+1$.)

## Comparison to the measured ratios & caveats

Measured (this session, `ratio_axial_sp_over_tp_vs_a2`):
$\ \mathrm{sp}\,\ell1/\mathrm{tp}\,\ell1\approx1.0$ and $\ \mathrm{sp}\,\ell0/\mathrm{tp}\,\ell1\approx1.0$.
Both AGREE with the (corrected) prediction $\Delta/\Delta=1$: sp $\ell=0$ and $\ell=1$ are both $\Delta=2$.

Two things to reconcile before trusting either number:
1. **Vector vs axial.** Eqs. (4.26)-(4.36) are for the *conserved* current ($\Delta=2$, $\ell=0$ TP
   protected). The lattice `ylm_axial` is the **axial** current, not conserved; its $\Delta_A$ and the
   $\ell=0$ (non-)vanishing need the axial version of the block.
2. **Sign / naming.** Paper $G_s>0$ (4.28), $G_t<0$ (4.31); lattice `tp=s3` is **positive**, `sp=s1+s2`
   is **negative** (verified in the fold diagnostic). The signs line up as $\text{lattice tp}\leftrightarrow G_s$,
   $\text{lattice sp}\leftrightarrow G_t$ — i.e. the lattice `tp`/`sp` labels may be *opposite* to the
   paper's $G_t$/$G_s$. If so the roles swap and the "sp $\ell=0$" we measured is really the TP-type
   $\ell=0$, whose leading $e^{-\Delta t}$ analysis must be redone with the axial block.
3. The scalar-$Y$ projection of the *tangential* (VSH) sp correlator is a **VSH mix**, so the measured
   scalar-$\ell=0$ can be contaminated by $\ell=1$-like pieces (mass $\Delta$), which would pull the
   ratio toward 1 regardless of the clean channel's $\Delta+1$.

NEXT to settle it: (a) redo the block for the axial current dimension; (b) confirm the s1,s2,s3 $\to$
temporal/spatial frame mapping against the paper's $\Lambda$ convention; (c) if a clean $\ell=0$ is
wanted, use the VSH (not scalar-$Y$) projection.
