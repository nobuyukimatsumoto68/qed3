# Reality and $t$-symmetry of $\langle J_V J_V\rangle$ from Eq. (3.58) -- working note

Reference: `qed3int_v2-11.pdf`, massless limit (Sec. 3.3.2).  The text around Eq. (3.58) is cut off
("...after the gauge averaging $\langle\cdot\rangle_U$"), so this note reconstructs what (3.58) can and
cannot give for the reality of the vector correlator.  Please correct the current $T$-transform if my
kernel assumption below is off.

## CONVERGED CONCLUSION

The vector correlator $\langle J_V J_V\rangle$ is **T-even** (exactly) and **real** -- the gauge-averaged
imaginary part is **strictly zero** (exact), NOT an $O(a^2)$ residual.  Precisely:

- **T-even -- exact.**  By trace cyclicity ($C(t_0\to t)=C(t\to t_0)$) plus translation invariance from gauge
  averaging, $g(\Delta t)\equiv\langle C_V(\Delta t)\rangle_U$ satisfies $g(\Delta t)=g(-\Delta t)$.  The
  WHOLE complex correlator is reflection-even, so **both Re and Im are T-even** (the Im is not T-odd).  This
  is kinematic (algebra + translation), independent of (3.58) and of the cutoff.

- **Real -- the gauge-averaged $\mathrm{Im}$ is STRICTLY ZERO (exact), not $O(a^2)$.**  For the massless,
  even-$N_f$ setup the gauge measure is $T$-invariant (the parity anomaly cancels among the even number of
  2-component fermions), so exact lattice $T$ (3.58) $+$ exact cyclicity give $\langle C\rangle=\langle
  C\rangle^{*}$ EXACTLY: $\mathrm{Im}\,g(\Delta t)=0$ with NO $O(a)/O(a^2)$ residual.  Per configuration
  $\mathrm{Im}\,C\neq0$, but it averages to exactly zero, so in the **data the Im is pure noise**.
  (Reality here is kinematic for even $N_f$.  In the parity sector $m_P$ the measure is $T$-broken by the $R$
  phase, so there $\mathrm{Im}\neq0$ -- the genuine parity signal.)

- **What this gives.**  $\mathrm{Im}\,g$ is an exactly-zero, signal-free channel: $\langle\mathrm{Im}\rangle$
  beyond noise = a consistency/bug failure (strong check) for massless/$m_F$, and a free noise gauge for the
  $|G|$ estimate (see `jj_abs_estimator_claude.md`).  A T-odd $\mathrm{Im}$ would also violate cyclicity, so
  any T-odd appearance is statistical noise ($\sim1/\sqrt N$).  $\mathrm{Re}\,g$ is the clean T-even physical
  signal; the CFT content lives in its functional form (and the $m_P$ $R$-phase), not in any massless Im.

Details and the derivations behind each bullet are below.

## Ingredients (massless)

Per gauge configuration $U$ (Eqs. 3.52, 3.54):
$$
A(t_0,t;U)\;\equiv\;C^{nn'}_{V++,c}(t_0\to t;U)\;=\;\mathrm{tr}\!\big[K^{nn'}(t_0)\,D_\text{ov}^{-1}\,K^{nn'}(t)\,D_\text{ov}^{-1}\big],
$$
$$
C^{nn'}_{V--,c}(t_0\to t;U)\;=\;A(t_0,t;U)^{*}. \tag{3.54}
$$
Eq. (3.54) is an **exact, per-configuration** identity (complex-conjugating the Wick contraction sends
$D_\text{ov}^{-1}\!\to\!D_\text{ov}^{\dagger-1}$ and $K\!\to\!K^\dagger$).  The disconnected pieces obey the
same conjugation (3.53/3.54).

Time reversal (Eq. 3.58), valid **only after** $\langle\cdot\rangle_U$ (the gauge measure is $T$-invariant):
$$
\sigma_3\,D_\text{ov}[U_T]\,\sigma_3 = D_\text{ov}[U]^\dagger
\;\Longrightarrow\;
D_\text{ov}[U]^{-1} = \sigma_3\,\big(D_\text{ov}[U_T]^{-1}\big)^{\dagger}\,\sigma_3 .
$$

## What (3.58) yields

Insert the $T$-rewrite of $D_\text{ov}^{-1}$ into $A$ and use $\langle\,f(U_T)\,\rangle_U=\langle\,f(U)\,\rangle_U$.
The kernel carries a definite $T$-parity (this is the step the cut-off text should state precisely).  With
$T$ the antiunitary time reversal, the **spatial** current is $T$-odd and the **temporal** ($j^0$, charge
density) is $T$-even:
$$
\sigma_3\,K^{nn'}|_T\,\sigma_3 = -\,K^{nn'\dagger}\ (\text{spatial, }T\text{-odd}),\qquad
\sigma_3\,K^{t}|_T\,\sigma_3 = +\,K^{t\dagger}\ (\text{temporal, }T\text{-even}).
$$
(My earlier "$K\to K^\dagger$" dropped this $\pm$ sign.)  In a DIAGONAL correlator the same kernel sits at
both ends, so the $T$-parity sign enters **squared** ($(\pm1)^2=+1$) -- it cancels and does not affect the
box below; it would only matter for the off-diagonal $\langle j^0 j^i\rangle$ (one of each $\Rightarrow$ net
$T$-odd, not measured).  The $\dagger$'s then turn the trace into its complex conjugate and swap
$t_0\leftrightarrow t$:
$$
\boxed{\;\langle A(t_0,t)\rangle_U \;=\; \langle A(t,t_0)\rangle_U^{*}\;}
\qquad\text{i.e.}\qquad
\langle C_V(\Delta t)\rangle_U \;=\; \langle C_V(-\Delta t)\rangle_U^{*},\quad \Delta t=t-t_0 .
$$
So **(3.58) is a reflection$\,+\,$conjugation relation** (a $t$-symmetry statement), tying the correlator at
$+\Delta t$ to the *conjugate* at $-\Delta t$ after gauge averaging.  It is the vector analogue of the axial
$C_{A-+}(t_0\to t)=C_{A+-}(t\to t_0)$ (3.57), but with an extra complex conjugation.  (If $T$ is instead a
pure time-axis reflection, temporal$\leftrightarrow$spatial parities swap; immaterial for the diagonal box
since the sign squares out.)

## Reflection-evenness is automatic (trace cyclicity), NOT a Lorentz-parity statement

[CORRECTION to an earlier draft that invoked "temporal even / mixed odd" Lorentz parity -- that was wrong.]
The connected estimator (3.52) has the SAME kernel $K$ at both ends, so by **trace cyclicity**, per config,
$$
C_{V++,c}(t_0\to t;U)=\mathrm{tr}[K(t_0)D_\text{ov}^{-1}K(t)D_\text{ov}^{-1}]
   =\mathrm{tr}[K(t)D_\text{ov}^{-1}K(t_0)D_\text{ov}^{-1}]=C_{V++,c}(t\to t_0;U),
$$
and the disc piece (3.53) is a product of two single-time traces, also $t_0\leftrightarrow t$ symmetric.
After gauge averaging restores translation invariance ($C\to g(\Delta t)$), the midpoint exchange becomes
$$
g(\Delta t)=g(-\Delta t).
$$
This holds for **any diagonal projection** (tp, sp, ylm; same kernel both ends) -- it is NOT tp-specific and
NOT a continuum reflection-parity argument.  "Always reflection-even" is justified ONLY at this level
(same kernel both ends $+$ gauge-averaged translation invariance), with the usual lattice link-placement /
contact caveats affecting magnitudes near $\Delta t=0$, not the symmetry.

## Does this give reality?  -> it reduces to the $\sigma_3$/channel action of $T$

Reflection-evenness (cyclicity) + the box (3.58) would give $g=g^{*}$ (real) **iff $T$ maps the measured
channel to ITSELF**:
$$
g(\Delta t)\overset{(3.58)}{=}\langle C(-\Delta t)\rangle^{*}\overset{\text{cyclicity}}{=}g(\Delta t)^{*}
\;\Rightarrow\; g\in\mathbb{R}.
$$
But the $\sigma_3$ in (3.58) acts on exactly the components that distinguish the two channels, so if $T$
**swaps** $+\leftrightarrow-$, then (3.58) only reproduces $C_{V--}=C_{V++}^{*}$ (Eq. 3.54) -- i.e. it
re-derives reflection-evenness and says **nothing new about reality**.  Which case holds is precisely the
content of the sentence the v2-11 text cuts off after "...after the gauge averaging $\langle\cdot\rangle_U$".

## Bottom line

- **Reflection-evenness** $g(\Delta t)=g(-\Delta t)$ is automatic from **trace cyclicity** (+ translation),
  generic to the diagonal correlators -- not from (3.58), not from Lorentz parity, not tp-special.
- Whether **(3.58) yields reality** depends entirely on the $\sigma_3$/channel action of $T$:
  channel-preserving $\Rightarrow$ real; channel-swapping $\Rightarrow$ only re-derives (3.54), no reality.
- The empirical massless tp data ($\mathrm{Im}\approx0$) is consistent with the channel-preserving case,
  but that identification is the assumption to nail down (the cut-off step).
- A **single channel** alone (`Vpp`) is real only in the channel-preserving case; otherwise its imaginary
  part survives and is tied to $C_{V--}=C_{V++}^{*}$.

Consistent with the paper's remark (below 3.70) that there is "no simple relation among the vector currents"
(only the axials get the clean (3.57)); for the vectors the strongest clean statement is the reflection
(cyclicity) + the conjugation (3.58/3.54), whose collapse to reality needs the channel-action input.

## Reality as a CFT / continuum test (observable)

[CORRECTION: an earlier draft claimed reality is "guaranteed in the continuum by Hermiticity" -- WRONG.]
Hermiticity of the current does NOT by itself make the Euclidean two-point real; that needs the Euclidean
MEASURE to be real (reflection-positive, conjugation-compatible reflection).  In QED3 this is exactly where
it can fail: the **parity anomaly** -- in odd dimensions the fermion determinant carries a phase (induced
Chern-Simons / $\eta$-invariant), so $\det D$ is generically COMPLEX.  That phase is the parity-breaking
effect the paper tracks via the reweighting factor $R$ (Eq. 2.5; cf. the remark below it that the $R$
fluctuation becomes negligible once parity is restored at large $N_f$).  With a complex weight, there is no
kinematic reason for $\langle C_{V++}\rangle_U$ to be real.

So reality of $\langle JJ\rangle$ is NOT guaranteed -- not by the lattice symmetry (needs the $T$
channel-action) and NOT kinematically in the continuum (needs the parity-breaking phase to be absent/cancel,
i.e. parity $T$ effectively unbroken).  The Sec. 4 CFT forms
$\langle j^\alpha j^\beta\rangle = C_J I^{\alpha\beta}/|x-y|^{2\Delta}$ (real $C_J$, real tensor) are real
because the presumed IR fixed point is **parity-invariant**.  Hence the **imaginary part of the single
channel is a quantitative test of parity restoration / conformality of the IR**:
$$
\mathrm{Im}\,\langle C_{V++}(\Delta t)\rangle_U
   = \tfrac{1}{2i}\big\langle C_{V++}-C_{V--}\big\rangle_U
   \;\xrightarrow{\ \text{conformal IR}\ }\ 0 \quad(\text{away from contact}).
$$
## Continuum: per-config reality from overlap anti-Hermiticity

[SCOPE: this section is a PER-CONFIGURATION statement.  The GAUGE-AVERAGED $\mathrm{Im}$ is a stronger
matter and is **exactly zero** (exact lattice $T$ (3.58) $+$ cyclicity $+$ $T$-invariant even-$N_f$ massless
measure) -- see the converged conclusion.  So the "cutoff effect" language below describes only the
per-config Im (finite-$a$ noise that averages to exactly zero); it is NOT an $O(a^2)$ residual in the
observable $\langle\mathrm{Im}\rangle_U$.]

A per-config statement.  The kernel is
$K^{wz}=\partial D_\text{ov}/\partial\theta_{wz}$ (derivative w.r.t. the REAL link angle), so
$K^\dagger=\partial D_\text{ov}^\dagger/\partial\theta$.  The massless continuum Dirac operator is
anti-Hermitian, so FORMALLY $D_\text{ov}^\dagger\to -D_\text{cont}$, and hence $K^\dagger\to -K$.  Then
$$
C_{V++}^{*}=C_{V--}=\mathrm{tr}[K^\dagger(t_0)D_\text{ov}^{\dagger-1}K^\dagger(t)D_\text{ov}^{\dagger-1}]
   \xrightarrow{\ \text{cont}\ }\mathrm{tr}[(-K)(-D^{-1})(-K)(-D^{-1})]=C_{V++},
$$
four sign flips cancelling pairwise (disc: each trace $\mathrm{tr}(D^{\dagger-1}K^\dagger)\to\mathrm{tr}(D^{-1}K)$).
So $C_{V++}^{*}\to C_{V++}$: **real PER CONFIGURATION** in the continuum.  The imaginary part is therefore
controlled by the finite-$a$ **failure of overlap anti-Hermiticity** $D_\text{ov}^\dagger+D_\text{cont}\neq0$
(the GW-type term) -- a **cutoff effect**, vanishing as $a\to0$ without gauge averaging.

IMPORTANT: this per-config non-anti-Hermiticity does NOT survive in the OBSERVABLE.  After gauge averaging,
the exact lattice $T$ (3.58) $+$ cyclicity force $\langle\mathrm{Im}\,C_{V++}\rangle_U=0$ EXACTLY (for the
$T$-invariant even-$N_f$ massless measure) -- there is no $O(a)/O(a^2)$ residual in $\langle\mathrm{Im}\rangle_U$.
The per-config $\mathrm{Im}\,C_{V++}\neq0$ (by 3.54, $C_{V--}=C_{V++}^{*}$) is finite-$a$ noise that averages
to exactly zero.  (The parity anomaly / $R$ phase shows up instead in the MASS sector $m_P$, where the
measure is $T$-broken; there $\langle\mathrm{Im}\rangle\neq0$ is the genuine signal.)

So in the DATA, $\langle\mathrm{Im}\rangle_U$ is **pure noise** (consistent with exactly zero); the notebooks
carry $\mathrm{Im}\,$Vpp and the Vmm$-\mathrm{conj}($Vpp$)$ residual purely as a consistency / noise gauge.

## Parity (in $\Delta t$) of the finite-$a$ Im: EVEN, by cyclicity (not odd)

At finite cutoff two relations act on $g(\Delta t)\equiv\langle C_{V++}(\Delta t)\rangle_U$:
- **Box (3.58):** $g(\Delta t)=g(-\Delta t)^{*}$  $\Rightarrow$  $\mathrm{Re}\,g$ even, $\mathrm{Im}\,g$ **odd**.
- **Trace cyclicity (EXACT, per config, matched kernels):**
  $C(t_0\to t)=\mathrm{tr}[K(t_0)D^{-1}K(t)D^{-1}]=\mathrm{tr}[K(t)D^{-1}K(t_0)D^{-1}]=C(t\to t_0)$, so after
  translation $g(\Delta t)=g(-\Delta t)$  $\Rightarrow$  $\mathrm{Re}\,g$ even, $\mathrm{Im}\,g$ **even**.

Cyclicity is algebraic (independent of 3.58 and of $a$) and the estimator is matched
(`psi`$=D^{-\dagger}K^\dagger(t_0)\eta$, sink $K(t)D^{-1}\eta$, and $(K^\dagger)^\dagger(D^{-\dagger})^\dagger
=K D^{-1}$ reproduces $\mathrm{tr}[K(t_0)D^{-1}K(t)D^{-1}]$).  So **cyclicity forces $\mathrm{Im}\,g$ to be
reflection-EVEN** -- the same symmetry as $\mathrm{Re}\,g$.

Consequences (3.58 is EXACT for the $T$-invariant even-$N_f$ massless measure):
- Im even (cyclicity) $+$ Im odd (box) $\Rightarrow$ $\mathrm{Im}\,g=0$ **EXACTLY** (reality); NO $O(a^2)$
  residual in the gauge average.  Per config $\mathrm{Im}\,C\neq0$ but it averages to exactly zero.
- A genuinely **T-ODD** Im would *violate cyclicity* -> it cannot be physical; it is the **statistical
  noise** (cyclicity/translation symmetry is only restored after $\langle\cdot\rangle_U$, so finite-$N$ noise
  carries an odd component that must average away $\sim 1/\sqrt N$).

(Channel-swap caveat: if $T$ in (3.58) swaps $+\leftrightarrow-$, the box reduces via (3.54) to
$g(\Delta t)=g(-\Delta t)$ = cyclicity again -> no reality, but Im still EVEN.  Either channel action gives
Im EVEN.)

**Diagnostic (to add to the notebook):** split $\mathrm{Im}\,g(\Delta t)$ into even/odd parts about
$\Delta t=N_t/2$.  Prediction: BOTH parts are consistent with zero -- $\langle\mathrm{Im}\rangle$ is exactly
zero in expectation, so the **even** part is statistical (no physical residual for even-$N_f$ massless) and
the **odd** part is noise ($\sim1/\sqrt N$).  A stable nonzero EVEN part beyond noise would flag a bug or a
non-$T$-invariant measure (e.g. the $m_P$ sector); a stable ODD part flags a source/sink ($K^\dagger$ vs $K$)
asymmetry.  ($\mathrm{Re}\,g$ is the clean T-even physical signal, from cyclicity and the box.)
